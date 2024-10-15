program hw2_q5
   use spectral

   ! Utility
   integer :: i
   
   ! Size
   integer, parameter :: N=16

   ! Wavenumbers
   integer, dimension(N) :: k1v
   
   ! Function and fourier coefficients
   double complex, dimension(N) :: f1, f1hat 

   ! Space
   real(8) :: xmin, xmax, dx, x

   xmin = 0.
   xmax = 2. * PI
   dx   = (xmax - xmin)/N

   ! Get function and wavenumber
   k1v = get_kv(N)
   do i=1, N
      x = (i - 1)*dx
      f1(i) = 0.5*COS(x) + 0.25*COS(6*x) + 0.125*COS(12*x) 
   end do

   call FFT_1D(f1, f1hat, N)

   call write_complex_array(f1hat, N, "./outs/q5_16.txt")

end program hw2_q5