program hw5
   use spectral
   use grid_class
   use navierstokes_class
   implicit none

   ! Variables
   integer :: N
   logical :: alias

   ! Command-line arguments
   integer :: iarg
   character(len=20) :: arg
   real(8) :: visc, k0, u0

   ! Type inits for grid and state 
   type(grid) :: mesh
   type(problem) :: state

   ! Energy spectrum
   real(8), dimension(:), allocatable :: Ek,Ekv

   ! Defaults
   N = 32                  ! Default grid size
   alias = .true.          ! Default to dealiasing enabled
   visc=1.0d0

   ! Read command-line arguments
   do iarg = 1, command_argument_count()
      call get_command_argument(iarg, arg)
      select case (iarg)
         case (1)
            read(arg, *) N
         case (2)
            read(arg, *) alias
         case (3)
            read(arg, *) visc
      end select
   end do

   ! Prep spectrum and axis
   allocate(Ek (1:N)); Ek=0.0d0
   allocate(Ekv(1:N)); Ekv=0.0d0

   ! Initialize domain
   mesh=grid(Nx=N,Ny=N,Nz=N,Lx=2*PI,Ly=2*PI,Lz=2*PI)
   mesh%k1v=get_kv(mesh%Nx)
   mesh%k2v=get_kv(mesh%Ny)
   mesh%k3v=get_kv(mesh%Nz)

   state=problem_constructor(mesh=mesh, scheme=RK4, alias=alias, visc=visc)

   init_timer: block
      real(8) :: dt_max
      ! Calculate the maximum allowable dt for stability
      state%tf=5.0
      state%t=0.0
      state%dt=0.001
      state%step=0
      state%stepf= 50 !99999999
   end block init_timer

   k0=5.0d0
   u0=1.0d0

   call state%initialize_hit(mesh,k0,u0)

   call write_complex_array_3D(state%U(1,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/Uh_ic.txt")
   call write_complex_array_3D(state%U(2,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/Vh_ic.txt")
   call write_complex_array_3D(state%U(3,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/Wh_ic.txt")

   call compute_radial_spectrum(state%U, mesh, Ek, Ekv)

   call write_double_array(Ek, N, "./outs/Ek_ic.txt")
   call write_double_array(Ekv, N, "./outs/Ekax_ic.txt")

   ! Transform the initial condition
   call state%transform_vel(mesh, BACKWARD)

   call write_complex_array_3D(state%U(1,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/U_ic.txt")
   call write_complex_array_3D(state%U(2,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/V_ic.txt")
   call write_complex_array_3D(state%U(3,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/W_ic.txt")

   ! Transform the initial condition
   call state%transform_vel(mesh, FORWARD)

   ! Advance Fourier coeffs
   do while(.not.(state%sdone.or.state%tdone))
      call state%advance(mesh)
      call state%adjust_time()
      if (mod(state%step,10).eq.1) then
         call compute_radial_spectrum(state%U, mesh, Ek, Ekv)
         call state%transform_vel(mesh, BACKWARD)
         call output(state, mesh)
         call state%transform_vel(mesh, FORWARD)
      end if
   end do

   ! Transform the initial condition
   call state%transform_vel(mesh, BACKWARD)

   call write_complex_array_3D(state%U(1,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/U_later.txt")
   call write_complex_array_3D(state%U(2,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/V_later.txt")
   call write_complex_array_3D(state%U(3,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/W_later.txt")

   deallocate(Ek)

   contains

   subroutine output(state, mesh)
      implicit none
      type(problem), intent(in) :: state
      type(grid), intent(in) :: mesh
      integer :: i

      character(len=100) :: filename
      character(len=10) :: n_dir, visc_dir
      character(len=5) :: n_str, visc_str
   
      ! Convert N and visc to character strings
      write(n_str, '(I0)') N
      write(visc_str, '(F3.1)') visc
   
      ! Create directory names based on N and viscosity
      write(n_dir, '(A)') "n" // trim(adjustl(n_str))
      write(visc_dir, '(A)') "nu" // trim(adjustl(visc_str)) 
   
      write(filename, '(A,I5.5,A)') "./outs/U_", state%step, ".txt"
      ! Write the complex array to the generated filename
      call write_complex_array_3D(state%U(1,:,:,:), mesh%nx, mesh%ny, mesh%nz, filename)
      write(filename, '(A,I5.5,A)') "./outs/V_", state%step, ".txt"
      ! Write the complex array to the generated filename
      call write_complex_array_3D(state%U(2,:,:,:), mesh%nx, mesh%ny, mesh%nz, filename)
      write(filename, '(A,I5.5,A)') "./outs/W_", state%step, ".txt"
      ! Write the complex array to the generated filename
      call write_complex_array_3D(state%U(3,:,:,:), mesh%nx, mesh%ny, mesh%nz, filename)

      write(filename, '(A,I5.5,A)') "./outs/Ek_", state%step, ".txt"
      call write_double_array(Ek, N, filename)
      write(filename, '(A,I5.5,A)') "./outs/Ekax_", state%step, ".txt"
      call write_double_array(Ekv, N, filename)
   
      ! Print the current time and step
      print *, "Time :: ", state%t, " Step :: ", state%step
   end subroutine output
end program hw5
