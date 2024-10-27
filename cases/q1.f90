program hw3_q1
   use spectral
   use grid_class
   implicit none

   ! Grid type
   type(grid) :: mesh

   ! Size
   integer, parameter :: N=64

   ! Wavenumbers
   integer, dimension(N) :: k1v
   
   ! Function and fourier coefficients
   double complex, dimension(N) :: f
   double complex, dimension(N) :: E

   ! Initialize domain
   mesh=grid(Nx=N,Lx=4*PI)
   call mesh%print

   !> Compute un-windowed spectrum   init_problem: block
   init_problem: block
      integer :: i
      do i=1,mesh%Nx
         f(i)=cos(2.0*PI*9.5/mesh%Lx*mesh%xv(i))
      end do
   end block init_problem

   ! Transform RHS
   call FFT_1D(f, f, mesh%nx)
   ! Compute
   E=compute_spectrum(f,mesh%nx)
   ! Output the spectrum
   call write_complex_array_1D(E, mesh%nx, "./outs/q1_Eff.txt")

   !> Compute windowed spectrum
   init_problem_window: block
      integer :: i
      do i=1,mesh%Nx
         f(i)=cos(2.0*PI*9.5/mesh%Lx*mesh%xv(i))*SIN(PI*mesh%xv(i)/mesh%Lx)
      end do
   end block init_problem_window

   ! Transform RHS
   call FFT_1D(f, f, mesh%nx)
   ! Compute
   E=compute_spectrum(f,mesh%nx)
   ! Output the spectrum
   call write_complex_array_1D(E, mesh%nx, "./outs/q1_Egg.txt")
end program hw3_q1