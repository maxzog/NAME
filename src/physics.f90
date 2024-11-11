module physics_class

   public :: problem

   !> Equations to be solved
   integer, parameter :: diffusion=1
   integer, parameter :: advection=2
   integer, parameter :: advection_diffusion=3
   integer, parameter :: Burgers=4
   integer, parameter :: NavierStokes=5 ! Not available

   !> Integration scheme
   integer, parameter :: euler=1
   integer, parameter :: RK2=2
   integer, parameter :: RK4=3

   type problem
      double complex, dimension(:), allocatable :: U
      double complex, dimension(:), allocatable :: RHS
      real(8) :: visc = 0.0
      real(8) :: c = 0.0
      integer :: equation = diffusion
      integer :: scheme = euler
      real(8) :: CFL=0.0
      real(8) :: dt=0.0
      real(8) :: tf=0.0
      real(8) :: t =0.0
      integer :: step=0
      integer :: stepf=0
      logical :: dealias=.true.
      logical :: sdone=.false.
      logical :: tdone=.false.
   contains
      procedure, public :: advance
      procedure, public :: get_rhs
      procedure, public :: adjust_time
   end type problem

   interface problem
      procedure problem_constructor
   end interface problem

   contains

   function problem_constructor(mesh, eqn, scheme, alias, visc, c, CFL, dt, tf, t, step, stepf) result(self)
      use grid_class, only: grid
      implicit none
      type(problem) :: self
      type(grid) :: mesh
      integer :: eqn, scheme
      real(8), optional :: visc,c,CFL,dt,tf,t
      integer, optional :: step,stepf
      logical, optional :: alias
      
      allocate(self%U  (1:mesh%nx)); self%U  =0.0
      allocate(self%RHS(1:mesh%nx)); self%RHS=0.0

      ! Maker sure we've got everything
      if (.not.present(visc) .and. eqn == diffusion) call abort
      if (.not.present(visc) .and. eqn == advection_diffusion) call abort
      if (.not.present(visc) .and. eqn == Burgers) call abort
      if (.not.present(c) .and. eqn == advection) call abort
      if (.not.present(c) .and. eqn == advection_diffusion) call abort

      if (present(alias)) self%dealias=alias

      self%scheme=scheme
      self%equation=eqn

      if (present(c)) then
         self%c=c
      else
         self%c=0.0
      end if
      if (present(visc)) then
         self%visc=visc
      else
         self%visc=0.0
      end if
   end function problem_constructor

   subroutine advance(this, mesh)
      use grid_class
      implicit none
      class(problem), intent(inout) :: this
      class(grid), intent(in) :: mesh
      
      select case(this%scheme)
      case(euler)
         call advance_euler(this, mesh)
      case(RK2)
         call advance_rk2(this, mesh)
      case(RK4)
         call advance_rk4(this, mesh)
      end select
      
      contains
      
      subroutine advance_euler(this, mesh)
         implicit none 
         class(problem), intent(inout) :: this
         class(grid), intent(in) :: mesh

         call this%get_rhs(mesh)
         this%U = this%U + this%RHS*this%dt
      end subroutine advance_euler

      subroutine advance_rk2(this, mesh)
         implicit none
         class(problem), intent(inout) :: this
         class(grid), intent(in) :: mesh
         double complex, dimension(:), allocatable :: Utemp

         allocate(Utemp(1:size(this%U))); Utemp=this%U

         call this%get_rhs(mesh)
         this%U = this%U + 0.5d0*this%RHS*this%dt
         call this%get_rhs(mesh)
         this%U = Utemp + this%RHS*this%dt

         deallocate(Utemp)
      end subroutine advance_rk2

      subroutine advance_rk4(this, mesh)
         implicit none
         class(problem), intent(inout) :: this
         class(grid), intent(in) :: mesh
         double complex, dimension(:), allocatable :: Utemp
         double complex, dimension(:), allocatable :: k1,k2,k3,k4

         ! Store state at step n
         allocate(Utemp(1:size(this%U))); Utemp=this%U

         ! Create storage for intermediate steps
         allocate(k1(1:size(this%U))); k1=0.0d0
         allocate(k2(1:size(this%U))); k2=0.0d0
         allocate(k3(1:size(this%U))); k3=0.0d0

         ! Call RHS for first stage
         call this%get_rhs(mesh)
         k1=this%RHS
         ! Advance using first stage
         this%U = Utemp + 0.5d0*this%dt*k1
         call this%get_rhs(mesh)
         k2=this%RHS
         ! Advance using second stage
         this%U = Utemp + 0.5d0*this%dt*k2
         call this%get_rhs(mesh)
         k3=this%RHS
         this%U = Utemp + this%dt*k3
         ! Call RHS for fourth stage
         call this%get_rhs(mesh)

         this%U = Utemp + this%dt/6.0d0 * (k1 + 2.0d0*k2 + 2.0d0*k3 + this%RHS)

         deallocate(Utemp,k1,k2,k3)
      end subroutine advance_rk4
   end subroutine advance


   subroutine get_rhs(this, mesh)
      use grid_class
      implicit none
      class(problem), intent(inout) :: this
      class(grid), intent(in) :: mesh

      !> Call the subroutine based on eqn specified
      select case(this%equation)
      case(diffusion)
         call get_rhs_diffusion(this)
      case(advection)
         call get_rhs_advection(this)
      case(advection_diffusion)
         call get_rhs_advection_diffusion(this)
      case(Burgers)
         if (this%dealias) then
            call get_rhs_burgers_a(this)
         else
            call get_rhs_burgers(this)
         end if
      end select

      contains
      
      !> Compute RHS of 1D diffusion equation
      subroutine get_rhs_diffusion(this)
         implicit none
         class(problem), intent(inout) :: this
         integer :: i
         do i=1,mesh%Nx
            this%RHS(i) = -this%visc*mesh%k1v(i)**2*this%U(i)
         end do
      end subroutine get_rhs_diffusion

      !> Compute RHS of 1D advection
      subroutine get_rhs_advection(this)
         use spectral, only: i_unit
         implicit none
         class(problem), intent(inout) :: this
         integer :: i
         do i=1,mesh%Nx
            this%RHS(i) = -i_unit*this%c*mesh%k1v(i)*this%U(i)
         end do
      end subroutine get_rhs_advection

      !> Compute RHS of 1D advection-diffusion
      subroutine get_rhs_advection_diffusion(this)
         use spectral, only: i_unit
         implicit none
         class(problem), intent(inout) :: this
         integer :: i
         do i=1,mesh%Nx
            this%RHS(i) = -i_unit*this%c*mesh%k1v(i)*this%U(i) & 
            &             - this%visc*mesh%k1v(i)**2*this%U(i)
            print *, this%U(i)
         end do
      end subroutine get_rhs_advection_diffusion

      !> Compute RHS of 1D Burgers
      subroutine get_rhs_burgers(this)
         use spectral
         implicit none
         class(problem), intent(inout) :: this
         double complex, dimension(:), allocatable :: U2
         integer :: i

         allocate(U2(1:size(this%U))); U2=0.0d0

         ! Compute non-linear term in physical space
         call iFFT_1D(U2, this%U, mesh%nx)
         do i=1,mesh%Nx
            U2(i)=0.5d0*U2(i)**2
         end do
         ! Transform back
         call FFT_1D(U2, U2, mesh%nx)
        
         do i=1,mesh%Nx
            this%RHS(i) = -i_unit*mesh%k1v(i)*U2(i) & 
            &             - this%visc*mesh%k1v(i)**2*this%U(i)
         end do

         deallocate(U2)
      end subroutine get_rhs_burgers

      !> Compute RHS of 1D Burgers
      subroutine get_rhs_burgers_a(this)
         use spectral
         implicit none
         class(problem), intent(inout) :: this
         double complex, dimension(:), allocatable :: U2
         double complex, dimension(:), allocatable :: Upad
         integer :: i,j,Npad

         Npad=size(this%U)+mesh%nx
         allocate(Upad(1:Npad)); Upad=0.0d0

         j=0
         do i=1,Npad
            if (i < mesh%nx/2) then
               j=j+1
               Upad(i)=this%U(j)
            elseif (i < 1.5*mesh%nx) then
               Upad(i)=0.0d0
            else
               j=j+1
               Upad(i)=this%U(j)
            end if
         end do

         allocate(U2(1:Npad)); U2=0.0d0

         ! Compute non-linear term in physical space
         call iFFT_1D(U2, Upad, Npad)
         do i=1,Npad
            U2(i)=0.5d0*U2(i)**2
         end do
         ! Transform back
         call FFT_1D(U2, U2, Npad)

         do i=1,mesh%Nx
            if (i < mesh%Nx/2) then
               Upad(i) = U2(i)
            else 
               Upad(i) = U2(i+mesh%Nx)
            end if
         end do
        
         do i=1,mesh%Nx
            this%RHS(i) = -i_unit*mesh%k1v(i)*Upad(i) & 
            &             - this%visc*mesh%k1v(i)**2*this%U(i)
         end do

         deallocate(U2, Upad)
      end subroutine get_rhs_burgers_a
   end subroutine get_rhs

   subroutine adjust_time(this)
      implicit none
      class(problem), intent(inout) :: this

      this%t = this%t + this%dt
      this%step = this%step + 1

      this%tdone = this%t > this%tf
      this%sdone = this%step > this%stepf
   end subroutine adjust_time


end module physics_class