program hw2_q4
   use spectral

   ! Utility
   integer :: i
   
   ! Size
   integer, parameter :: N=256

   ! Wavenumbers
   integer, dimension(N) :: k1v
   
   ! Function and fourier coefficients
   double complex, dimension(N) :: f1, f1hat 
   double complex, dimension(N) :: ac, achat 

   ! Space
   real(8) :: xmin, xmax, dx, x

   ! Problem-specific parameters
   integer :: m=0
   real(8) :: phi=0

   xmin = 0.
   xmax = 2. * PI
   dx   = (xmax - xmin)/N

   ! Get function and wavenumber
   k1v = get_kv(N)
   do i=1, N
      x = (i - 1)*dx
      f1(i) = COS(2*PI*m/xmax*(x + phi))
   end do

   call FFT_1D(f1, f1hat, N)
   achat = compute_autocorr(f1hat, N)
   call iFFT_1D(ac, achat, N)

   call write_real_array(ac, N, "./outs/q4a_m0_p0.txt")

end program hw2_q4