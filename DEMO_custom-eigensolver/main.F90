program main
  use cudafor
  use cublas
  use cusolverDn
  use nvtx
  use zhegst_jdr
  use funcs
  use compare_utils
  implicit none
  
  integer :: N,i,j, info, lda, istat, Lwork
  character(len=20) :: arg
  real(8) :: wallclock
  real(8) :: rv, iv, ts, te
  real(8), allocatable, dimension(:) :: W, W1, W2, W3, workspace, rworkspace
  complex(8), dimension(:,:), allocatable :: A, B
  complex(8), dimension(:,:), allocatable, device :: A_d, B_d, A1_d, A2_d, A3_d, B1_d, B2_d, B3_d
  complex(8), device, dimension(:), allocatable :: workspace_d
  integer, device :: devInfo_d
  type(cusolverDnHandle) :: h


  ! Parse command line arguments
  i = command_argument_count()

  if (i == 1) then
    ! If N is provided, generate random hermetian matrices for A and B
    print*, "Using randomly-generated matrices..."
    call get_command_argument(1, arg)
    read(arg, *)  N
    print*, "Running with N = ", N
    lda = N

    ! Create random positive-definite hermetian matrices on host
    call create_random_hermetian_pd(A, N)
    call create_random_hermetian_pd(B, N)

  else if (i == 0) then
    ! If no argument provided, read QEspresso matrices A_input and B_input
    print*, "Using QEspresso matrices..."
    N = 1600
    print*, "Running with N = ", N
    lda = 3200

    call read_QEspresso_matrices(A, B, N, lda)

  else
    print*, "Usage:\n\t ./main [N]"
    call exit
  endif


  ! Copy matrices to device
  allocate(A_d, source = A)
  allocate(A1_d, source = A)
  allocate(A2_d, source = A)
  allocate(A3_d, source = A)

  allocate(B_d, source = B)
  allocate(B1_d, source = B)
  allocate(B2_d, source = B)

  ! Initialize magma, cublas, etc
  call magmaf_init
  istat = cublasInit
  if (istat /= CUBLAS_STATUS_SUCCESS) write(*,*) 'cublas intialization failed'

  ! TASK 1:  Perform cholesky factorization of B, B := L ---------------------
  print*, "Performing cholesky factorizations...."
  ! CASE 1: Factor using Intel MKL 
  ts = wallclock()
  call nvtxStartRange("MKL zpotrf")
  call zpotrf('L', N, B, lda, istat)
  call nvtxEndRange
  te = wallclock()
  if (istat /= 0) write(*,*) 'MKL zpotrf failed'
  print*, "\tTime for MKL zpotrf = ", te - ts
  print*

  ! CASE 2: Factor using MAGMA
  ts = wallclock()
  call nvtxStartRange("magmaf_zpotrf_gpu")
  call magmaf_zpotrf_gpu('L', N, B1_d, lda, istat)
  call nvtxEndRange
  te = wallclock()
  if (istat /= 0) write(*,*) 'magmaf_zpotrf failed'
  print*, "\tTime for magmaf_zpotrf_gpu = ", te - ts
  print*

  ! CASE 3: Factor using CUSOLVER
  istat = cusolverDnCreate(h)
  if (istat /= CUSOLVER_STATUS_SUCCESS) write(*,*) 'handle creation failed'

  istat = cusolverDnZpotrf_bufferSize(h, CUBLAS_FILL_MODE_LOWER, N, B2_d, lda, Lwork)
  if (istat /= CUSOLVER_STATUS_SUCCESS) write(*,*) 'cusolverDnZpotrf_buffersize failed'

  allocate(workspace_d(Lwork))

  ts = wallclock()
  call nvtxStartRange("cusolverDnZpotrf")
  istat = cusolverDnZpotrf(h, CUBLAS_FILL_MODE_LOWER, N, B2_d, lda, workspace_d, Lwork, devInfo_d)
  call nvtxEndRange
  te = wallclock()
  print*, "\tTime for cusolverDnZpotrf = ", te - ts
  print*

  if (istat /= CUSOLVER_STATUS_SUCCESS) write(*,*) 'cusolverDnZpotrf failed'

  istat = devInfo_d
  if (istat /= 0) write(*,*) 'Cholesky factorization failed'

  istat = cusolverDnDestroy(h)
  if (istat /= CUSOLVER_STATUS_SUCCESS) write(*,*) 'handle destruction failed'


  ! TASK 2:  Reduce generalized Eigenproblem to standard form ----------------
  print*, "Performing reduction to generalized Eigenproblem...."
  ! CASE 1: Perform computation on CPU using MKL
  ts = wallclock()
  call nvtxStartRange("MKL")
  call zhegst(1, 'L', N, A, lda, B, lda, istat)
  call nvtxEndRange
  te = wallclock()
  if (istat /= 0) write(*,*) 'MKL zhegst failed'
  print*, "\tTime for MKL zhegst = ", te - ts
  print*

  A_d = A

  ! CASE 2: Perform computation using MAGMA routine
  ts = wallclock()
  call nvtxStartRange("MAGMA", 0)
  call magmaf_zhegst_gpu(1, 'L', N, A1_d, lda, B1_d, lda, istat)
  call nvtxEndRange
  te = wallclock()
  if (istat /= 0) write(*,*) 'magma zhegst failed'
  print*, "\tTime for magmaf_zhegst_gpu = ", te - ts
#ifdef ACC
  call compare_ltri(A_d, A1_d, N, N)
#endif
  print*

  ! CASE 3: Perform computation using 2 cublas ZTRSM calls
  ts = wallclock()
  call nvtxStartRange("GPU V1", 1)
  call zhegst_gpu_v1(1, 'L', N, A2_d, lda, B2_d, lda)
  call nvtxEndRange
  te = wallclock()

  print*, "\tTime for zhegst_gpu_v1 = ", te - ts
#ifdef ACC
  call compare_ltri(A_d, A2_d, N, N)
#endif
  print*

  ! CASE 4: Perform computation using 2 cublas ZTRSM calls for subblock
  ts = wallclock()
  ts = wallclock()
  call nvtxStartRange("GPU V2", 0)
  call zhegst_gpu_v2(1, 'L', N, A3_d, lda, B2_d, lda, 448)
  call nvtxEndRange
  te = wallclock()
  print*, "\tTime for zhegst_gpu_v2 = ", te - ts
#ifdef ACC
  call compare_ltri(A_d, A3_d, N, N)
#endif
  print*


#ifdef ACC
  !! TASK 3: Compute eigenvalues (all done on CPU using MKL) and compare ------
  allocate(workspace(N*N))
  allocate(rworkspace(3*N - 2))
  allocate(W(N))
  allocate(W1(N))
  allocate(W2(N))
  allocate(W3(N))
  
  ! Compute eigenvalues using MKL for baseline
  call zheev('N', 'L', N, A, lda, W, workspace, N*N, rworkspace, istat)
  if (istat /= 0) write(*,*) 'MKL zheev failed'

  A = A1_d
  call zheev('N', 'L', N, A, lda, W1, workspace, N*N, rworkspace, istat)
  if (istat /= 0) write(*,*) 'MKL zheev failed'

  A = A2_d
  call zheev('N', 'L', N, A, lda, W2, workspace, N*N, rworkspace, istat)
  if (istat /= 0) write(*,*) 'MKL zheev failed'

  A = A3_d
  call zheev('N', 'L', N, A, lda, W3, workspace, N*N, rworkspace, istat)
  if (istat /= 0) write(*,*) 'MKL zheev failed'

  call compare(W, W1, N)
  call compare(W, W2, N)
  call compare(W, W3, N)

  call write_array(W, N, "W_MKL.dat")
  call write_array(W1, N, "W_MAGMA.dat")
  call write_array(W2, N, "W_GPU_V1.dat")
  call write_array(W3, N, "W_GPU_V2.dat")

#endif


end program main
