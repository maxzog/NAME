program hw3_q2
   use spectral
   use grid_class
   implicit none

   ! Grid type
   type(grid) :: mesh

   ! Size
   integer, parameter :: N=32
   integer, parameter :: M=32

   ! Wavenumbers
   integer, dimension(N) :: k1v
   integer, dimension(M) :: k2v
   
   ! Function and fourier coefficients
   double complex, dimension(N,M) :: P, Q

   ! Initialize domain
   mesh=grid(Nx=N,Ny=M,Lx=2*PI,Ly=2*PI)
   call mesh%print

   init_problem: block
      integer :: i, j
      do j=1,mesh%Ny
         do i=1,mesh%Nx
            Q(i,j)=-5.0*SIN(mesh%xv(i))*COS(2.0*mesh%yv(j))
         end do
      end do
      P=0.0
   end block init_problem

   ! Transform RHS
   call FFT_2D(Q, Q, mesh%nx, mesh%ny)
   ! Solve
   call mesh%solve_poisson_2D(P, Q)
   ! Output the coeffs
   call write_complex_array_2D(P, mesh%nx, mesh%ny, "./outs/q2_phat.txt")
   ! Transform
   call iFFT_2D(P, P, mesh%nx, mesh%ny)
   ! Output the function
   call write_complex_array_2D(P, mesh%nx, mesh%ny, "./outs/q2_p.txt")

end program hw3_q2