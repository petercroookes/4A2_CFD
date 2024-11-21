
      subroutine euler_iteration(av,g)

!     This subroutine calculates the fluxes into each cell and then sums them to
!     update the primary flow properties

!     Explicitly declare the required variables
      use ieee_arithmetic  
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

      !print *, 'Check first and second values'
      !print *, 'rovx1', g%rovx(1,1), 'rovx2', g%rovx(1,2), g%rovx(2,1)
      !print *, 'rovy2', g%rovy(1,1), 'rovy2', g%rovy(1,2), g%rovy(2,1)
      !print *, 'First Density', g%ro(1,:)

      mass_i(:,:) = 0.5 * (((g%rovx(:,1:nj-1) + g%rovx(:, 2:nj)) * &
           g%lx_i(:, 1:nj-1)) + (g%rovy(:, 1:nj-1) + g%rovy(:, 2:nj)) * g%ly_i(:, 1:nj-1))
      mass_j(:,:) = 0.5 * (((g%rovx(1:ni-1, :) + g%rovx(2:ni, :)) * &
           g%lx_j(1:ni-1,:)) + (g%rovy(1:ni-1, :) + g%rovy(2:ni, :)) * g%ly_j(:ni-1,:))  
      
      ! Check if any NaN values are in mass_i or mass_j
      !if (any(ieee_is_nan(mass_i)) .or. any(ieee_is_nan(mass_j))) then
      !   print *, "NaN detected in mass_i / mass_j"
      !   print *,  'mass_i', mass_i, 'mass_j', mass_j
      !end if

      !print *, 'Find sum of mass fluxes at i = 1 cell:'
      !do j = 1, nj - 1
      !   print *, 'Cell : i = 1, j = ', j
      !   print *, mass_i(1,j) + mass_j(1,j) - mass_i(2,j) - mass_j(1,j+1)
      !end do
      
!     Apply the wall boundary condition by checking that two nodes at the
!     end of a facet are both on a wall, if so then the appropriate mass
!     flow array is set to have zero flow through that facet
      where(g%wall(1:ni-1,:) .and. g%wall(2:ni,:)) mass_j = 0 
      where(g%wall(:,1:nj-1) .and. g%wall(:,2:nj)) mass_i = 0 

!     Update the density with mass fluxes by calling "sum_fluxes"
!     INSERT

      call sum_fluxes(av, mass_i, mass_j, g%area, g%ro, g%dro)

      ! Check if any NaN values in dro
      if (any(ieee_is_nan(g%dro))) then
         print *, "NaN detected in dro"
      end if
      
!     Setup the conservation of energy equation by calculated the enthalpy flux
!     and storing the values in "flux_i" and "flux_j", you will need "mass_i"
!     and "mass_j" from before
!     INSERT

      flux_i(:,:) = 0.5 * mass_i(:,:) * (g%hstag(:, 1:nj-1) + g%hstag(:, 2:nj))
      flux_j(:,:) = 0.5 * mass_j(:,:) * (g%hstag(1:ni-1, :) + g%hstag(2:ni, :))
      
!     Update the internal energy with enthalpy fluxes
!     INSERT

      call sum_fluxes(av, flux_i, flux_j, g%area, g%roe, g%droe)

      ! Check if any NaN values in droe
      if (any(ieee_is_nan(g%droe))) then
         print *, "NaN detected in droe"
      end if
      
!     Setup the x-momentum equation including momentum flux and pressure forces
!     INSERT

      flux_i(:,:) = 0.5 * (mass_i(:,:) * (g%vx(:,1:nj-1) + g%vx(:,2:nj)) &
           + (g%p(:, 1:nj-1) + g%p(:, 2:nj)) * g%lx_i(:,:)) 
      flux_j(:,:) = 0.5 * (mass_j(:,:) * (g%vx(1:ni-1,:) + g%vx(2:ni,:)) &
           + (g%p(1:ni-1,:) + g%p(2:ni,:)) * g%lx_j(:,:))

!     Update the x-momentum with momentum flux
!     INSERT

      call sum_fluxes(av, flux_i, flux_j, g%area, g%rovx, g%drovx)

      ! Check if any NaN values in drovx or drovy
      if (any(ieee_is_nan(g%drovx))) then
         print *, "NaN detected in dro"
      end if

      if (any(ieee_is_nan(g%drovy))) then
         print *, "NaN detected in drovy"
      end if
      
!     Setup the y-momentum equation including momentum flux and pressure forces
!     INSERT

      flux_i(:,:) = 0.5 * (mass_i(:,:) * (g%vy(:,1:nj-1) + g%vy(:,2:nj)) &
           + (g%p(:, 1:nj-1) + g%p(:, 2:nj)) * g%ly_i(:,:)) 
      flux_j(:,:) = 0.5 * (mass_j(:,:) * (g%vy(1:ni-1,:) + g%vy(2:ni,:)) &
           + (g%p(1:ni-1,:) + g%p(2:ni,:)) * g%ly_j(:,:))
      
!     Update the y-momentum with momentum flux
!     INSERT

      call sum_fluxes(av, flux_i, flux_j, g%area, g%rovy, g%drovy)
      
!     Add artificial viscosity by smoothing all of the primary flow variables
      call smooth_array(av,g%ro)
      call smooth_array(av,g%roe)
      call smooth_array(av,g%rovx)
      call smooth_array(av,g%rovy)

      end subroutine euler_iteration


