module finite_difference

   contains

      ! One dimensional second-order central difference
      ! assumes periodic domain
      subroutine cdf2(f, df, dx, nx)
         implicit none
         double complex, intent(inout), dimension(nx) :: f, df
         real(8), intent(in) :: dx
         integer, intent(in) :: nx
         integer :: i

         do i=2,nx-1
            df(i) = 0.5*f(i+1) - 0.5*f(i-1)
         end do
         df(1)  = 0.5*f(2) - 0.5*f(nx)
         df(nx) = 0.5*f(1) - 0.5*f(nx-1)
         df = df/dx
      end subroutine cdf2

end module finite_difference