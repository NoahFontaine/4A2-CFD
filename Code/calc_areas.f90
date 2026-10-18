      
      subroutine calc_areas(g)

!     Calculate the area of the quadrilateral cells and the lengths of the side
!     facets

!     Explicitly declare the required variables
      use types
      use routines
      implicit none
      type(t_grid), intent(inout) :: g
      integer :: ni, nj

!     Declare integers or any extra variables you need here
!     INSERT
      integer :: i, j
      real :: l_i_min, l_j_min
!     END INSERT

!     Get the size of the mesh and store locally for convenience
      ni = g%ni; nj = g%nj;

!     Calculate the areas of the cells and store in g%area. The area of any
!     quadrilateral is half of the magnitude of the cross product of the two
!     vectors that form the diagonals. Check the order of your product so that
!     the values come out positive! You can do this using two nested loops in
!     the i and j-directions or in a vectorised way by indexing the coordinate
!     arrays with lists of indices
!     INSERT
      do i = 1, ni-1
            do j = 1, nj-1
                  g%area(i,j) = 0.5*abs((g%x(i+1,j+1)-g%x(i,j))*(g%y(i+1,j)-g%y(i,j+1)) - &
                                        (g%y(i+1,j+1)-g%y(i,j))*(g%x(i+1,j)-g%x(i,j+1)))
                  end do
      end do
!     END INSERT

!     Calculate the projected lengths in the x and y-directions on all of the
!     "i = const" facets and store them in g%lx_i and g%ly_i. When combined
!     together these two components define a vector that is normal to the facet,
!     pointing inwards towards the centre of the cell. This is only the case for
!     the left hand side of the cell, the vector stored in position i,j points
!     towards the centre of the i,j cell
!     INSERT
      do i = 1, ni
            do j = 1, nj-1
!                 Looking at i=const facets (i.e. the left hand side of the cell).
!                 The normal vector is given by n = (y(i,j+1) - y(i,j), -(x(i,j+1) - x(i,j))).
!                 This is the [lx_i, ly_i] vector.
                  g%lx_i(i,j) = g%y(i,j+1) - g%y(i,j)
                  g%ly_i(i,j) = -(g%x(i,j+1) - g%x(i,j))
            end do
      end do
!     END INSERT

!     Now repeat the calculation for the project lengths on the "j=const"
!     facets. 
!     INSERT
      do i = 1, ni-1
            do j = 1, nj
!                 Looking at j=const facets (i.e. the bottom side of the cell).
!                 The normal vector is given by n = (-(y(i+1,j) - y(i,j)), x(i+1,j) - x(i,j)).
!                 This is the [lx_j, ly_j] vector.
                  g%lx_j(i,j) = -(g%y(i+1,j) - g%y(i,j))
                  g%ly_j(i,j) = g%x(i+1,j) - g%x(i,j)
            end do
      end do
!     END INSERT

!     Find the minimum length scale in the mesh, this is defined as the length
!     of the shortest side of all the cells. Call this length "l_min", it is used
!     to set the timestep from the CFL number. Start by calculating the lengths
!     of the i and j facets by using the intrinsic function "hypot", this avoids
!     underflow and overflow errors. Then find the overal minimum value using
!     both the "min" and "minval" functions.
!     INSERT
      l_i_min = minval(hypot(g%lx_i, g%ly_i))
      l_j_min = minval(hypot(g%lx_j, g%ly_j))
      g%l_min = min(l_i_min, l_j_min)
!     END INSERT
!
!     Print the overall minimum length size that has been calculated
      write(6,*) 'Calculated cell areas and facet lengths'
      write(6,*) '  Overall minimum element size = ', g%l_min
      write(6,*)

      end subroutine calc_areas
