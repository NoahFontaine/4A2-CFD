
      module smooth_stencil

!     Packaging a subroutine in a module allows it to recieve the data
!     conveniently as assumed shape arrays
      
      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine smooth_array(av,prop)

!     This subroutine smooths "prop" to stabilise the calculation, the basic 
!     solver uses second order smoothing, many improvements are possible.

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      real, intent(inout) :: prop(:,:)
      real, dimension(size(prop,1),size(prop,2)) :: prop_avg
      integer :: ni, nj

!     Get the block size and store locally for convenience
      ni = size(prop,1); nj = size(prop,2)

!     Calculate the average values at the nodes in the interior region of the
!     mesh, use the four neighbouring nodes in the plus and minus i and 
!     j-directions.
!     INSERT
      prop_avg(2:ni-1,2:nj-1) = 0.25 * (prop(1:ni-2,2:nj-1) + prop(3:ni,2:nj-1) + &
          prop(2:ni-1,1:nj-2) + prop(2:ni-1,3:nj))
!     END INSERT

!     Edge values are also averaged in both the i and j-directions. Parallel to
!     the boundary the averaging is centred, the averages of two nodes are taken
!     either side of the current point. Perpendicular to the boundary the
!     algorithm is one-sided, the value at the current point is extrapolated
!     from the values at two nodes away from the boundary point.
!     INSERT
      prop_avg(1,2:nj-1) = (prop(1,1:nj-2) + prop(1,3:nj) + &
          2.0 * prop(2,2:nj-1) - prop(3,2:nj-1)) / 3.0
      prop_avg(ni,2:nj-1) = (prop(ni,1:nj-2) + prop(ni,3:nj) + &
          2.0 * prop(ni-1,2:nj-1) - prop(ni-2,2:nj-1)) / 3.0
      prop_avg(2:ni-1,1) = (prop(1:ni-2,1) + prop(3:ni,1) + &
          2.0 * prop(2:ni-1,2) - prop(2:ni-1,3)) / 3.0
      prop_avg(2:ni-1,nj) = (prop(1:ni-2,nj) + prop(3:ni,nj) + &
          2.0 * prop(2:ni-1,nj-1) - prop(2:ni-1,nj-2)) / 3.0
!     END INSERT

!     The corner values are not currently smoothed
      prop_avg([1,ni],[1,nj]) = prop([1,ni],[1,nj])

!     Now apply the artificial viscosity by smoothing "prop" towards "prop_avg",
!     take (1-sfac) * the calculated value of the property + sfac * the average 
!     of the surrounding values. 
!     INSERT
      prop = (1.0 - av%sfac) * prop + av%sfac * prop_avg
!     END INSERT

      end subroutine smooth_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module smooth_stencil

