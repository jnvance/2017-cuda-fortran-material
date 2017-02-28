subroutine vector_add(a,b,N)
  real:: a(*)
  real:: b
  integer:: i,N
  real :: a_d(N)

#ifdef __CUDA__
  attributes(device) :: a_d
#endif
 
  a_d(1:N) = a(1:N)
  print *,"N=",N
  !$cuf kernel do <<<*,*>>>
   do i=1,N
     a_d(i)=a_d(i)+b
   end do
  
  a(1:N) = a_d(1:N)

end subroutine vector_add


