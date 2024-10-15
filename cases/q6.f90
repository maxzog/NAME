program hw2_q6
   use spectral

   ! Utility
   integer :: i
   
   ! Size
   integer, parameter :: N=256

   ! Wavenumbers
   integer, dimension(N) :: k1v
   
   ! Function and fourier coefficients
   double complex, dimension(N) :: f1, f1hat 
   double complex, dimension(N) :: Ek

   ! Space
   real(8) :: xmin, xmax, dx, x

   xmin = 0.
   xmax = 2. * PI
   dx   = (xmax - xmin)/N

   ! Get function and wavenumber
   k1v = get_kv(N)
   f1=0.
   do i=1, N
      x = (i - 1)*dx
      ! f1(i) = EXP(-x**2) 
      if (x<=PI) f1(i) = SIN(x)
   end do

   call FFT_1D(f1, f1hat, N)
   Ek = compute_spectrum(f1hat, N)
   call write_complex_array(Ek, N, "./outs/q6b_256_Ek.txt")

end program hw2_q6