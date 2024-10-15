program hw2_q1
   use spectral

   ! Utility
   integer :: i
   
   ! Size
   integer, parameter :: N=8
   integer, parameter :: M=8
   integer, parameter :: L=8

   ! Wavenumbers
   integer, dimension(N) :: k1v
   integer, dimension(M) :: k2v
   integer, dimension(L) :: k3v
   
   ! Function and fourier coefficients
   double complex, dimension(N, M, L) :: f3, f3hat_k 
   double complex, dimension(N, M)    :: f2, f2hat_k 
   double complex, dimension(N)       :: f1, f1hat_k 

   ! Space
   real(8) :: xmin, xmax, dx, x
   real(8) :: ymin, ymax, dy, y
   real(8) :: zmin, zmax, dz, z

   xmin = 0.
   xmax = 2. * PI
   dx   = (xmax - xmin)/N

   ymin = 0.
   ymax = 2. * PI
   dy   = (ymax - ymin)/M

   zmin = 0.
   zmax = 2. * PI
   dz   = (zmax - zmin)/L

   ! Get function and wavenumber
   do k=1, L
      k3v(k) = k - L/2 - 1
      z = (k - 1)*dz
      do j=1, M
         k2v(j) = j - M/2 - 1
         y = (j - 1)*dy
         do i=1, N
            k1v(i) = i - N/2 - 1
            x = (i - 1)*dx
            ! Function evals
            f2(i,j)   = SIN(2*x)*COS(4*y) + COS(6*x)*SIN(5*y)
            f3(i,j,k) = SIN(2*x)*COS(4*y)*COS(3*z) + COS(4*x)*SIN(y)
         end do
      end do
   end do

   call FFT_2D(f2, f2hat_k, N, M)

   call FFT_3D(f3, f3hat_k, N, M, L)

   call write_complex_array_2D(f2hat_k, N, M, "./outs/q1_2d_8.txt")
   call write_complex_array_3D(f3hat_k, N, M, L, "./outs/q1_3d_8.txt")
end program hw2_q1