module navierstokes_class
   use grid_class, only: grid

   public :: problem

   !> Integration scheme
   integer, parameter :: euler=1
   integer, parameter :: RK2=2
   integer, parameter :: RK4=3

   !> Transform flags
   integer, parameter :: FORWARD=1
   integer, parameter :: BACKWARD=2

   type problem
      double complex, dimension(:,:,:,:), allocatable :: U
      double complex, dimension(:,:,:,:), allocatable :: U2
      !!      (uu uv uw)   (1 4 6)
      !! U2 = (uv vv vw) = (4 2 5)  
      !!      (uw vw ww)   (6 5 3)
      double complex, dimension(:,:,:,:), allocatable :: RHS
      double complex, dimension(:,:,:), allocatable :: P
      real(8) :: visc = 0.0
      real(8) :: c = 0.0
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
      procedure, public :: transform_vel
      procedure, public :: transform_U2
      procedure, public :: compute_U2
      procedure, public :: compute_phat
      procedure, public :: adjust_time
      procedure, public :: initialize_hit
   end type problem

   interface problem
      procedure problem_constructor
   end interface problem

   contains

   function problem_constructor(mesh, scheme, alias, visc, c, CFL, dt, tf, t, step, stepf) result(self)
      implicit none
      type(problem) :: self
      type(grid) :: mesh
      integer :: scheme
      real(8), optional :: visc,c,CFL,dt,tf,t
      integer, optional :: step,stepf
      logical, optional :: alias
      
      allocate(self%P  (  1:mesh%nx,1:mesh%ny,1:mesh%nz)); self%P  =0.0
      allocate(self%U  (3,1:mesh%nx,1:mesh%ny,1:mesh%nz)); self%U  =0.0
      allocate(self%U2 (6,1:mesh%nx,1:mesh%ny,1:mesh%nz)); self%U2 =0.0
      allocate(self%RHS(3,1:mesh%nx,1:mesh%ny,1:mesh%nz)); self%RHS=0.0

      if (present(alias)) self%dealias=alias

      self%scheme=scheme

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
         double complex, dimension(:,:,:,:), allocatable :: Utemp

         allocate(Utemp, MOLD=this%U); Utemp=this%U

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
         double complex, dimension(:,:,:,:), allocatable :: Utemp
         double complex, dimension(:,:,:,:), allocatable :: k1,k2,k3,k4

         ! Store state at step n
         allocate(Utemp, MOLD=this%U); Utemp=this%U

         ! Create storage for intermediate steps
         allocate(k1, MOLD=this%U); k1=0.0d0
         allocate(k2, MOLD=this%U); k2=0.0d0
         allocate(k3, MOLD=this%U); k3=0.0d0

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
      implicit none
      class(problem), intent(inout) :: this
      class(grid), intent(in) :: mesh

      !> Call the subroutine based on eqn specified
      if (this%dealias) then
         ! call get_rhs_ns_a(this)
      else
         call get_rhs_ns(this)
      end if

      contains
      !> 3D Navier-Stokes RHS without dealiasing
      subroutine get_rhs_ns(this)
         use spectral
         implicit none
         class(problem), intent(inout) :: this
         integer :: i,j,k
         real(8) :: kmag

         ! Bring velocity to physical space
         call this%transform_vel(mesh, BACKWARD)
         ! Compute non-linear term
         call this%compute_U2(mesh)
         ! Bring velocity back to spectral space
         call this%transform_vel(mesh, FORWARD)
         ! Bring non-linear term to spectral space
         call this%transform_U2(mesh, FORWARD)

         ! Compute pressure in spectral space
         call this%compute_phat(mesh)
         
         do k=1,mesh%Nz
            do j=1,mesh%Ny
               do i=1,mesh%Nx
                  kmag=sqrt(real(mesh%k1v(i)**2 + mesh%k2v(j)**2 + mesh%k3v(k)**2, 8))
                  ! RHSi = - d/dxj (uiuj) + d/dxi P + nu d^2/dxj^2 u_i
                  ! RHS - x
                  this%RHS(1,i,j,k) = -i_unit*(mesh%k1v(i)*this%U2(1,i,j,k) + mesh%k2v(j)*this%U2(4,i,j,k) + mesh%k3v(k)*this%U2(6,i,j,k)) &
                  &                   -i_unit*mesh%k1v(i)*this%P(i,j,k)                                                                    &
                  &                   -this%visc*kmag**2*this%U(1,i,j,k)
                  ! RHS - y
                  this%RHS(2,i,j,k) = -i_unit*(mesh%k1v(i)*this%U2(4,i,j,k) + mesh%k2v(j)*this%U2(2,i,j,k) + mesh%k3v(k)*this%U2(5,i,j,k)) &
                  &                   -i_unit*mesh%k2v(j)*this%P(i,j,k)                                                                    &
                  &                   -this%visc*kmag**2*this%U(2,i,j,k)
                  ! RHS - z 
                  this%RHS(3,i,j,k) = -i_unit*(mesh%k1v(i)*this%U2(1,i,j,k) + mesh%k2v(j)*this%U2(4,i,j,k) + mesh%k3v(k)*this%U2(6,i,j,k)) &
                  &                   -i_unit*mesh%k3v(k)*this%P(i,j,k)                                                                    &
                  &                   -this%visc*kmag**2*this%U(3,i,j,k)
                  this%RHS(3,i,j,k) = 0.0d0
               end do
            end do
         end do

      end subroutine get_rhs_ns
   end subroutine get_rhs

   !> Compute the pressure in spectral space
   subroutine compute_phat(this, mesh)
      implicit none
      class (problem), intent(inout) :: this
      class (grid), intent(in) :: mesh
      integer :: i,j,k
      real(8) :: kmag

      this%P=0.0
      do k=1,mesh%nz
         do j=1,mesh%ny
            do i=1,mesh%nx
               kmag=sqrt(real(mesh%k1v(i)**2 + mesh%k2v(j)**2 + mesh%k3v(k)**2, 8))
               if (kmag.lt.epsilon(1.0d0)) cycle
               this%P(i,j,k) = this%P(i,j,k) - mesh%k1v(i)**2 * this%U2(1,i,j,k)
               this%P(i,j,k) = this%P(i,j,k) - mesh%k2v(j)**2 * this%U2(2,i,j,k)
               this%P(i,j,k) = this%P(i,j,k) - mesh%k3v(k)**2 * this%U2(3,i,j,k)
               this%P(i,j,k) = this%P(i,j,k) - 2.0*mesh%k1v(i)*mesh%k2v(j) * this%U2(4,i,j,k)
               this%P(i,j,k) = this%P(i,j,k) - 2.0*mesh%k2v(j)*mesh%k3v(k) * this%U2(5,i,j,k)
               this%P(i,j,k) = this%P(i,j,k) - 2.0*mesh%k1v(i)*mesh%k3v(k) * this%U2(6,i,j,k)
               this%P(i,j,k) = this%P(i,j,k) / kmag**2
            end do
         end do
      end do
   end subroutine compute_phat

   !> Compute stress tensor
   subroutine compute_U2(this, mesh)
      implicit none
      class (problem), intent(inout) :: this
      class (grid), intent(in) :: mesh
      integer :: i,j,k

      ! Compute non-linear term in physical space
      do k=1,mesh%Nz
         do j=1,mesh%Ny
            do i=1,mesh%Nx
               this%U2(1,i,j,k) = this%U(1,i,j,k)*this%U(1,i,j,k)
               this%U2(2,i,j,k) = this%U(2,i,j,k)*this%U(2,i,j,k)
               this%U2(3,i,j,k) = this%U(3,i,j,k)*this%U(3,i,j,k)
               
               this%U2(4,i,j,k) = this%U(1,i,j,k)*this%U(2,i,j,k)
               this%U2(5,i,j,k) = this%U(2,i,j,k)*this%U(3,i,j,k)
               this%U2(6,i,j,k) = this%U(1,i,j,k)*this%U(3,i,j,k)
            end do   
         end do   
      end do   
   end subroutine compute_U2

   !> For convenience
   !> Performs three 3D transforms on the velocity field
   subroutine transform_vel(this, mesh, direction)
      use spectral
      implicit none
      class (problem), intent(inout) :: this
      class (grid), intent(in) :: mesh
      double complex, dimension(:,:,:), allocatable :: temp_arr
      integer :: direction

      allocate(temp_arr, MOLD=this%U(1,:,:,:))

      select case(direction)
      case(FORWARD)
         temp_arr=this%U(1,:,:,:)
         call FFT_3D(temp_arr, temp_arr, mesh%nx, mesh%ny, mesh%nz)
         this%U(1,:,:,:)=temp_arr
         temp_arr=this%U(2,:,:,:)
         call FFT_3D(temp_arr, temp_arr, mesh%nx, mesh%ny, mesh%nz)
         this%U(2,:,:,:)=temp_arr
         temp_arr=this%U(3,:,:,:)
         call FFT_3D(temp_arr, temp_arr, mesh%nx, mesh%ny, mesh%nz)
         this%U(3,:,:,:)=temp_arr
      case(BACKWARD)
         temp_arr=this%U(1,:,:,:)
         call iFFT_3D(temp_arr, temp_arr, mesh%nx, mesh%ny, mesh%nz)
         this%U(1,:,:,:)=temp_arr
         temp_arr=this%U(2,:,:,:)
         call iFFT_3D(temp_arr, temp_arr, mesh%nx, mesh%ny, mesh%nz)
         this%U(2,:,:,:)=temp_arr
         temp_arr=this%U(3,:,:,:)
         call iFFT_3D(temp_arr, temp_arr, mesh%nx, mesh%ny, mesh%nz)
         this%U(3,:,:,:)=temp_arr
      end select

      deallocate(temp_arr)
   end subroutine transform_vel

   !> For convenience
   !> Performs six 3D transforms on the stress tensor
   subroutine transform_U2(this, mesh, direction)
      use spectral
      implicit none
      class (problem), intent(inout) :: this
      class (grid), intent(in) :: mesh
      double complex, dimension(:,:,:), allocatable :: temp_arr
      integer :: direction

      allocate(temp_arr, MOLD=this%U2(1,:,:,:))
      select case(direction)
      case(FORWARD)
         temp_arr=this%U2(1,:,:,:)
         call FFT_3D(temp_arr, temp_arr, mesh%nx, mesh%ny, mesh%nz)
         this%U2(1,:,:,:)=temp_arr
         temp_arr=this%U2(2,:,:,:)
         call FFT_3D(temp_arr, temp_arr, mesh%nx, mesh%ny, mesh%nz)
         this%U2(2,:,:,:)=temp_arr
         temp_arr=this%U2(3,:,:,:)
         call FFT_3D(temp_arr, temp_arr, mesh%nx, mesh%ny, mesh%nz)
         this%U2(3,:,:,:)=temp_arr
         temp_arr=this%U2(4,:,:,:)
         call FFT_3D(temp_arr, temp_arr, mesh%nx, mesh%ny, mesh%nz)
         this%U2(4,:,:,:)=temp_arr
         temp_arr=this%U2(5,:,:,:)
         call FFT_3D(temp_arr, temp_arr, mesh%nx, mesh%ny, mesh%nz)
         this%U2(5,:,:,:)=temp_arr
         temp_arr=this%U2(6,:,:,:)
         call FFT_3D(temp_arr, temp_arr, mesh%nx, mesh%ny, mesh%nz)
         this%U2(6,:,:,:)=temp_arr
      case(BACKWARD)
         temp_arr=this%U2(1,:,:,:)
         call iFFT_3D(temp_arr, temp_arr, mesh%nx, mesh%ny, mesh%nz)
         this%U2(1,:,:,:)=temp_arr
         temp_arr=this%U2(2,:,:,:)
         call iFFT_3D(temp_arr, temp_arr, mesh%nx, mesh%ny, mesh%nz)
         this%U2(2,:,:,:)=temp_arr
         temp_arr=this%U2(3,:,:,:)
         call iFFT_3D(temp_arr, temp_arr, mesh%nx, mesh%ny, mesh%nz)
         this%U2(3,:,:,:)=temp_arr
         temp_arr=this%U2(4,:,:,:)
         call iFFT_3D(temp_arr, temp_arr, mesh%nx, mesh%ny, mesh%nz)
         this%U2(4,:,:,:)=temp_arr
         temp_arr=this%U2(5,:,:,:)
         call iFFT_3D(temp_arr, temp_arr, mesh%nx, mesh%ny, mesh%nz)
         this%U2(5,:,:,:)=temp_arr
         temp_arr=this%U2(6,:,:,:)
         call iFFT_3D(temp_arr, temp_arr, mesh%nx, mesh%ny, mesh%nz)
         this%U2(6,:,:,:)=temp_arr
      end select
      deallocate(temp_arr)
   end subroutine transform_U2

   subroutine initialize_hit(this, mesh)
      implicit none
      class (problem), intent(inout) :: this
      class (grid), intent(in) :: mesh
      double complex :: top, bot, alpha, beta
      integer :: i,j,k
      real(8) :: kmag
      
      this%U=0.0d0
      do k=1,mesh%Nz/2
         do j=1,mesh%Ny/2
            do i=1,mesh%Nx/2
               kmag=sqrt(real(mesh%k1v(i)**2 + mesh%k2v(j)**2 + mesh%k3v(k)**2, 8))
               if (kmag.lt.epsilon(1.0d0)) cycle
               ! u (x) component
               top = alpha*kmag*mesh%k2v(j) + beta*mesh%k1v(i)*mesh%k3v(k)
               bot = kmag*sqrt(real(mesh%k1v(i)**2 + mesh%k2v(j)**2, 8))
               this%U(1,mesh%nx-(i-2),mesh%ny-(j-2),mesh%nz-(k-2)) = top/bot
               this%U(1,i,j,k) = conjg(top/bot)
               ! v (y) component
               top = beta*mesh%k2v(j)*mesh%k3v(k) - alpha*kmag*mesh%k1v(i) 
               bot = kmag*sqrt(real(mesh%k1v(i)**2 + mesh%k2v(j)**2, 8))
               this%U(2,mesh%nx-(i-2),mesh%ny-(j-2),mesh%nz-(k-2)) = top/bot
               this%U(2,i,j,k) = conjg(top/bot)
               ! w (z) component
               top = beta*sqrt(real(mesh%k1v(i)**2 + mesh%k2v(j)**2, 8)) 
               bot = kmag 
               this%U(3,mesh%nx-(i-2),mesh%ny-(j-2),mesh%nz-(k-2)) = top/bot
               this%U(3,i,j,k) = conjg(top/bot)
            end do
         end do
      end do
   end subroutine initialize_hit

   function get_alpha_beta(kmag) result(alpha, beta)
      real(8) :: kmag
      double complex :: alpha, beta
      alpha=sqrt(spectrum(kmag) / (4*PI*kmag**2))*exp(i_unit(theta1))*cos(phi)
      alpha=sqrt(spectrum(kmag) / (4*PI*kmag**2))*exp(i_unit(theta2))*sin(phi)
   end function get_alpha

   subroutine adjust_time(this)
      implicit none
      class(problem), intent(inout) :: this

      this%t = this%t + this%dt
      this%step = this%step + 1

      this%tdone = this%t > this%tf
      this%sdone = this%step > this%stepf
   end subroutine adjust_time

end module navierstokes_class