program hw2_q3
   use spectral
   use finite_difference

   ! Utility
   integer :: i
   
   ! Size
   integer, parameter :: N=32

   ! Wavenumbers
   integer, dimension(N) :: k1v
   
   ! Function and fourier coefficients
   double complex, dimension(N) :: f1, f1hat_k 
   double complex, dimension(N) :: df1

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
      f1(i) = SIN(2*x) + COS(6*x)
      ! if (x<=PI) f1(i) = SIN(x)
   end do

   call cdf2(f1, df1, dx, N)
   call write_real_array(df1, N, "./outs/q3a_cdf2.txt")

   call FFT_1D(f1, f1hat_k, N)
   call ddx_1D(f1hat_k, k1v, N)
   call iFFT_1D(f1, f1hat_k, N)
   call write_real_array(f1, N, "./outs/q3a.txt")
   call write_complex_array(f1hat_k, N, "./outs/q3a_hat.txt")

end program hw2_q3