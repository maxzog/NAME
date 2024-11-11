program hw4_q1
   use spectral
   use grid_class
   use physics_class
   implicit none

   ! Grid type
   type(grid) :: mesh

   ! State type
   type(problem) :: state

   ! Size
   integer, parameter :: N=64
   
   ! State 
   double complex, dimension(N) :: u

   ! Initialize domain
   mesh=grid(Nx=N,Lx=2*PI)
   call mesh%print

   state=problem_constructor(mesh=mesh, eqn=Burgers, scheme=RK4, alias=.true., visc=0.05d0)

   init_timer: block
      real(8) :: dt_max
      ! Calculate the maximum allowable dt for stability
      dt_max = 0.0005d0 * (mesh%dx**2) / (4.0d0 * state%visc)
      state%tf=2.0d0
      state%t=0.0
      state%dt=dt_max
      state%step=0
      state%stepf=9999999
   end block init_timer

   !> Compute un-windowed spectrum   init_problem: block
   init_problem: block
      integer :: i
      do i=1,mesh%Nx
         state%U(i)=sin(mesh%xv(i))
      end do
   end block init_problem

   ! Transform the initial condition
   call FFT_1D(state%U, state%U, mesh%nx)

   call write_complex_array_1D(state%U, mesh%nx, "./outs/state_ic.txt")

   ! Advance Fourier coeffs
   do while(.not.(state%sdone.or.state%tdone))
      call state%advance(mesh)
      call state%adjust_time()
      if (mod(state%step, 250) .eq. 0) call output(state, mesh)
   end do

   call output(state, mesh)

   contains

   subroutine output(state, mesh)
      implicit none
      type(problem), intent(in) :: state
      type(grid), intent(in) :: mesh

      ! Define the filename based on the current step
      character(len=100) :: filename
      if (state%dealias) write(filename, '(A,I5.5,A)') "./outs/state_rk4_", state%step, ".txt"
      if (.not.state%dealias) write(filename, '(A,I5.5,A)') "./outs_no/state_rk4_", state%step, ".txt"
   
      ! Write the complex array to the generated filename
      call write_complex_array_1D(state%U, mesh%nx, filename)
   
      ! Print the current time and step
      print *, "Time :: ", state%t, " Step :: ", state%step
   end subroutine output

end program hw4_q1
