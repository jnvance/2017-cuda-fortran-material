module funcs
  contains
  
  ! Subroutine to create a random positive-definite Hermitian matrix
  subroutine create_random_hermetian_pd(A, N)
    use cudafor
    use cublas
    complex(8), allocatable, dimension(:,:) :: A, temp
    complex(8), allocatable, dimension(:,:), device :: A_d, temp_d
    complex(8) :: val
    real(8) :: rv, iv
    integer :: i, j, N

    allocate(A(N,N))
    allocate(temp(N,N))

    ! Create general hermetian temp
    do j = 1, N
      do i = 1, N
        if (i > j) then
          call random_number(rv)
          call random_number(iv)
          temp(i,j) = cmplx(rv, iv, 8)
          temp(j,i) = conjg(temp(i,j))
        else if (i == j) then
          call random_number(rv)
          temp(i,j) = rv
        end if
      end do
    end do

    allocate(A_d, source = A)
    allocate(temp_d, source = temp)

    ! Multiply temp by conjugate transpose of temp to get positive definite hermetian A
    call cublaszgemm('N', 'C', N, N, N, cmplx(1.0, 0.0, 8), temp_d, N, temp_d, N, cmplx(0.0, 0.0, 8), A_d, N)

    A = A_d
    deallocate(temp)
    deallocate(A_d)
    deallocate(temp_d)
        
  end subroutine


  subroutine read_QEspresso_matrices(A, B, N, lda)
    complex(8), allocatable, dimension(:,:) :: A, B
    integer :: i, j

    allocate(A(lda, N))
    allocate(B(lda, N))

    open(unit=10, status='old', file="A_input", form="unformatted")
    read(10) A
    close(10)

    open(unit=10, status='old', file="B_input", form="unformatted")
    read(10) B
    close(10)

    ! Populate missing entries (QEspresso matrices contain upper triangular portion only)
    do j = 1, N
      do i = 1, N
        if (i > j) then
          A(i,j) = conjg(A(j,i))
          B(i,j) = conjg(B(j,i))
        end if
      end do
    end do

  end subroutine
end module funcs

module compare_utils

  interface compare
    module procedure compare_real_1d_cpu
  end interface compare

  interface compare_ltri
    module procedure compare_ltri_complex_2d_cpu
    module procedure compare_ltri_complex_2d_gpu
  end interface compare_ltri

  interface write_array
    module procedure write_real_1d_cpu
  end interface write_array

  contains

  subroutine compare_real_1d_cpu(A1_h,A2_h,N)
    implicit none
    real(8), dimension(:) :: A1_h, A2_h
    real(8), dimension(:), allocatable :: A1, A2
    real(8) :: maxerr,perr,l2normerr,norm,buf
    integer :: i,j,k,N,imax
    character (len=4) :: itcount
    character (len=1) :: proc

    allocate(A1, source = A1_h)
    allocate(A2, source = A2_h)

    l2normerr = 0.d0
    norm = 0.d0
    maxerr = 0.d0
    imax=1
    do i=1, N
        if(abs(A1(i)) >= 1e-10) then
          perr = abs(A1(i) - A2(i))/abs(A1(i))*100.d0
          norm = norm + abs(A1(i)*A1(i));
          l2normerr = l2normerr + abs((A1(i) - A2(i))*(A1(i) - A2(i)))
        else
          perr = 0.d0
        endif
        if(perr>maxerr .and. A1(i)/=0.0d0 .and. A2(i)/=0.0d0) then
          maxerr = perr
          imax = i
        endif
    enddo

  norm = sqrt(norm)
  l2normerr = sqrt(l2normerr)
  if(l2normerr /= 0.d0) then
    l2normerr = l2normerr/norm
    write(*,"(A16,2X,ES10.3,A12,ES10.3,A6,I5,A6,2X,E20.14,2X,A6,2X,E20.14)") &
    "l2norm error",l2normerr,"max error",maxerr,"% at",imax,"A1=",A1(imax),"A2=",A2(imax)
  else
    write(*,"(A16)") "EXACT MATCH"
  endif

  deallocate(A1, A2)

  end subroutine compare_real_1d_cpu

  subroutine compare_ltri_complex_2d_cpu(A1_h,A2_h, N, M)
    implicit none
    complex(8), dimension(:,:) :: A1_h, A2_h
    complex(8), dimension(:,:), allocatable :: A1, A2
    real(8) :: maxerr,perr,l2normerr,norm,buf
    integer :: i,j,k,imax,jmax,kmax, N, M
    character (len=4) :: itcount
    character (len=1) :: proc

    allocate(A1, source = A1_h)
    allocate(A2, source = A2_h)

    l2normerr = 0.d0
    norm = 0.d0
    maxerr = 0.d0
    imax=1
    jmax=1
    kmax=1
    do j= 1, M
      do i= 1, N
        if (j .le. i) then
          !print*, "i, j: " , i, " ", j, " A1, A2: ", A1(i,j), " ", A2(i,j)
          if(abs(A1(i,j)) >= 1e-10) then
            perr = abs(A1(i,j) - A2(i,j))/abs(A1(i,j))*100.d0
            norm = norm + abs(A1(i,j)*A1(i,j));
            l2normerr = l2normerr + abs((A1(i,j) - A2(i,j))*(A1(i,j) - A2(i,j)))
          else
            perr = 0.d0
          endif
          if(perr>maxerr .and. A1(i,j)/=0.0d0 .and. A2(i,j)/=0.0d0) then
            maxerr = perr
            imax = i
            jmax = j
          endif
        endif
      enddo
   enddo

  norm = sqrt(norm)
  l2normerr = sqrt(l2normerr)
  if(l2normerr /= 0.d0) then
    l2normerr = l2normerr/norm
    write(*,"(A16,2X,ES10.3,A12,ES10.3,A6,I5,I5,A6,2X,E20.14,1X,E20.14,2X,A6,2X,E20.14,1X,E20.14)") &
    "l2norm error",l2normerr,"max error",maxerr,"% at",imax,jmax,"A1=",REAL(A1(imax,jmax)),AIMAG(A1(imax,jmax)),"A2=",REAL(A2(imax,jmax)),AIMAG(A2(imax,jmax))
  else
    write(*,"(A16)") "EXACT MATCH"
  endif

  deallocate(A1, A2)

  end subroutine compare_ltri_complex_2d_cpu


  subroutine compare_ltri_complex_2d_gpu(A1_d,A2_d, N, M)
    implicit none
    complex(8), device, dimension(:,:) :: A1_d, A2_d
    complex(8), dimension(:,:), allocatable :: A1, A2
    real(8) :: maxerr,perr,l2normerr,norm,buf
    integer :: i,j,k,imax,jmax,kmax, N, M
    character (len=4) :: itcount
    character (len=1) :: proc

    allocate(A1, source = A1_d)
    allocate(A2, source = A2_d)

    l2normerr = 0.d0
    norm = 0.d0
    maxerr = 0.d0
    imax=1
    jmax=1
    kmax=1
    do j= 1, M
      do i= 1, N
        if (j .le. i) then
          !print*, "i, j: " , i, " ", j, " A1, A2: ", A1(i,j), " ", A2(i,j)
          if(abs(A1(i,j)) >= 1e-10) then
            perr = abs(A1(i,j) - A2(i,j))/abs(A1(i,j))*100.d0
            norm = norm + abs(A1(i,j)*A1(i,j));
            l2normerr = l2normerr + abs((A1(i,j) - A2(i,j))*(A1(i,j) - A2(i,j)))
          else
            perr = 0.d0
          endif
          if(perr>maxerr .and. A1(i,j)/=0.0d0 .and. A2(i,j)/=0.0d0) then
            maxerr = perr
            imax = i
            jmax = j
          endif
        endif
      enddo
   enddo

  norm = sqrt(norm)
  l2normerr = sqrt(l2normerr)
  if(l2normerr /= 0.d0) then
    l2normerr = l2normerr/norm
    write(*,"(A16,2X,ES10.3,A12,ES10.3,A6,I5,I5,A6,2X,E20.14,1X,E20.14,2X,A6,2X,E20.14,1X,E20.14)") &
    "l2norm error",l2normerr,"max error",maxerr,"% at",imax,jmax,"A1=",REAL(A1(imax,jmax)),AIMAG(A1(imax,jmax)),"A2=",REAL(A2(imax,jmax)),AIMAG(A2(imax,jmax))
  else
    write(*,"(A16)") "EXACT MATCH"
  endif

  deallocate(A1, A2)

  end subroutine compare_ltri_complex_2d_gpu

  subroutine write_real_1d_cpu(A, N, filename)
    real(8), dimension(:) :: A
    integer :: N
    character(len=*) :: filename

    open(unit=10, status='replace',file=filename)

    do i = 1,N
      write(10, "(ES22.14)") A(i)
    end do

    close(10)
  end subroutine

end module compare_utils
! ----
! nvtx
! ----

module nvtx
  use iso_c_binding
  use cudafor
  implicit none

  integer,private :: col(7) = [ Z'0000ff00', Z'000000ff', Z'00ffff00',Z'00ff00ff',Z'0000ffff', &
       Z'00ff0000', Z'00ffffff']
  character(len=256),private :: tempName

  type, bind(C):: nvtxEventAttributes
     integer(C_INT16_T):: version=1
     integer(C_INT16_T):: size=48 !
     integer(C_INT):: category=0
     integer(C_INT):: colorType=1 ! NVTX_COLOR_ARGB = 1
     integer(C_INT):: color
     integer(C_INT):: payloadType=0 ! NVTX_PAYLOAD_UNKNOWN = 0
     integer(C_INT):: reserved0
     integer(C_INT64_T):: payload   ! union uint,int,double
     integer(C_INT):: messageType=1  ! NVTX_MESSAGE_TYPE_ASCII     = 1 
     type(C_PTR):: message  ! ascii char
  end type nvtxEventAttributes
  
  interface nvtxRangePush
     ! push range with custom label and standard color
     subroutine nvtxRangePushA(name) bind(C, name='nvtxRangePushA')
       use iso_c_binding
       character(kind=C_CHAR,len=*) :: name
     end subroutine nvtxRangePushA
     
     ! push range with custom label and custom color
     subroutine nvtxRangePushEx(event) bind(C, name='nvtxRangePushEx')
       use iso_c_binding
       import:: nvtxEventAttributes
       type(nvtxEventAttributes):: event
     end subroutine nvtxRangePushEx
  end interface nvtxRangePush
  
  interface nvtxRangePop
     subroutine nvtxRangePop() bind(C, name='nvtxRangePop')
     end subroutine nvtxRangePop
  end interface nvtxRangePop
  
contains
  
  subroutine nvtxStartRange(name,id)
    character(kind=c_char,len=*) :: name
    integer, optional:: id
    type(nvtxEventAttributes):: event
    integer :: istat
    istat = cudaDeviceSynchronize()

    tempName=trim(name)//c_null_char

    if ( .not. present(id)) then
       call nvtxRangePush(tempName)
    else
       event%color=col(mod(id,7)+1)
       event%message=c_loc(tempName)
       call nvtxRangePushEx(event)
    end if
  end subroutine nvtxStartRange

  subroutine nvtxStartRangeAsync(name,id)
    character(kind=c_char,len=*) :: name
    integer, optional:: id
    type(nvtxEventAttributes):: event

    tempName=trim(name)//c_null_char

    if ( .not. present(id)) then
       call nvtxRangePush(tempName)
    else
       event%color=col(mod(id,7)+1)
       event%message=c_loc(tempName)
       call nvtxRangePushEx(event)
    end if
  end subroutine nvtxStartRangeAsync

  
  subroutine nvtxEndRange
    integer :: istat
    istat = cudaDeviceSynchronize()
    call nvtxRangePop
  end subroutine nvtxEndRange

  subroutine nvtxEndRangeAsync
    call nvtxRangePop
  end subroutine nvtxEndRangeAsync
  
end module nvtx

