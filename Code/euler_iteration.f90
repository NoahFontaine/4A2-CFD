
      subroutine euler_iteration(av,g)

!     This subroutine calculates the fluxes into each cell and then sums them to
!     update the primary flow properties

!     Explicitly declare the required variables
      use types
      use flux_stencil
      use smooth_stencil
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), intent(inout) :: g
      real, dimension(g%ni,g%nj-1) :: mass_i, flux_i
      real, dimension(g%ni-1,g%nj) :: mass_j, flux_j
      integer :: i, j, ni, nj

!     Get the block size and store locally for convenience
      ni = g%ni; nj = g%nj

!     Setup the continuity equation by calculating the mass flow through
!     the facets in both the i and j-directions. Store these values in
!     "mass_i" and "mass_j"
!     INSERT
      mass_i = 0.5 * (g%rovx(1:ni,1:nj-1) + g%rovx(1:ni,2:nj)) * g%lx_i(1:ni,1:nj-1) + &
      0.5 * (g%rovy(1:ni,1:nj-1) + g%rovy(1:ni,2:nj)) * g%ly_i(1:ni,1:nj-1)
      mass_j = 0.5 * (g%rovx(1:ni-1,1:nj) + g%rovx(2:ni,1:nj)) * g%lx_j(1:ni-1,1:nj) + &
      0.5 * (g%rovy(1:ni-1,1:nj) + g%rovy(2:ni,1:nj)) * g%ly_j(1:ni-1,1:nj)
!     END INSERT

!     Apply the wall boundary condition by checking that two nodes at the
!     end of a facet are both on a wall, if so then the appropriate mass
!     flow array is set to have zero flow through that facet
      where(g%wall(1:ni-1,:) .and. g%wall(2:ni,:)) mass_j = 0 
      where(g%wall(:,1:nj-1) .and. g%wall(:,2:nj)) mass_i = 0 

!     Update the density with mass fluxes by calling "sum_fluxes"
!     INSERT
      call sum_fluxes(av,mass_i,mass_j,g%area,g%ro,g%dro)
!     END INSERT

!     Setup the conservation of energy equation by calculated the enthalpy flux
!     and storing the values in "flux_i" and "flux_j", you will need "mass_i"
!     and "mass_j" from before
!     INSERT

!     Update the internal energy with enthalpy fluxes
!     INSERT

!     Setup the x-momentum equation including momentum flux and pressure forces
!     INSERT

!     Update the x-momentum with momentum flux
!     INSERT

!     Setup the y-momentum equation including momentum flux and pressure forces
!     INSERT

!     Update the y-momentum with momentum flux
!     INSERT

!     Add artificial viscosity by smoothing all of the primary flow variables
      call smooth_array(av,g%ro)
      call smooth_array(av,g%roe)
      call smooth_array(av,g%rovx)
      call smooth_array(av,g%rovy)
      

      end subroutine euler_iteration

