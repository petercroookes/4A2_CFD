
      module smooth_stencil

!     Packaging a subroutine in a module allows it to recieve the data
!     conveniently as assumed shape arrays
      
      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine smooth_array(av,prop,corr)

!     This subroutine smooths "prop" to stabilise the calculation, the basic 
!     solver uses second order smoothing, many improvements are possible.

!     Explicitly declare the required variables
      use types
      use ieee_arithmetic 
      implicit none
      type(t_appvars), intent(in) :: av
      real, intent(inout) :: prop(:,:), corr(:,:)
      real, dimension(size(prop,1),size(prop,2)) :: prop_avg, corr_total
      real :: fcorr = 0.9
      integer :: ni, nj, i, j

!     Get the block size and store locally for convenience
      ni = size(prop,1); nj = size(prop,2)
      
!     Calculate the average values at the nodes in the interior region of the
!     mesh, use the four neighbouring nodes in the plus and minus i and 
!     j-directions.
!     INSERT

      do i = 2, ni - 1
         do j = 2, nj - 1
            prop_avg(i,j) = 0.25 * (prop(i-1,j) + prop(i,j-1) &
                 + prop(i+1,j) + prop(i,j+1))
         end do
      end do
      
!     Edge values are also averaged in both the i and j-directions. Parallel to
!     the boundary the averaging is centred, the averages of two nodes are taken
!     either side of the current point. Perpendicular to the boundary the
!     algorithm is one-sided, the value at the current point is extrapolated
!     from the values at two nodes away from the boundary point.
!     INSERT

      do i = 2, ni - 1
         prop_avg(i,1) = (prop(i-1,1) + prop(i+1,1) + 2*prop(i,2) &
              - prop(i,3)) / 3.0
         prop_avg(i,nj) = (prop(i-1,nj) + prop(i+1,nj) + 2*prop(i,nj-1) &
              - prop(i,nj-2)) / 3.0
      end do

      do j = 2, nj - 1
         prop_avg(1,j) = (prop(1,j-1) + prop(1,j+1) + 2*prop(2,j) &
              - prop(3,j)) / 3.0
         prop_avg(ni,j) = (prop(ni,j-1) + prop(ni,j+1) + 2*prop(ni-1,j) &
              - prop(ni-2,j)) / 3.0
      end do
      
!     The corner values are not currently smoothed
      prop_avg([1,ni],[1,nj]) = prop([1,ni],[1,nj])

      ! Check if any NaN values are in mass_i or mass_j
      if (any(ieee_is_nan(corr)) .or. any(ieee_is_nan(prop_avg)) .or. any(ieee_is_nan(prop_avg))) then
         print *, "NaN detected in  prop, prop_avg or corr"
         print *,  'corr', corr, 'prop', prop, 'prop_avg', prop_avg
         stop
      end if
      
!     Calculated deferred corrections. corr is typically updated by only
!     small proportion of every timestep

      corr_total = fcorr * (prop - prop_avg)
      corr = 0.99 * corr + 0.01 * corr_total
      
!     Now apply the artificial viscosity by smoothing "prop" towards "prop_avg",
!     take (1-sfac) * the calculated value of the property + sfac * the average 
!     of the surrounding values. 
!     INSERT

!     Deferred correction applied, prop smoothed towards prop_avg + corr
!     instead of prop_avg

      prop = av%sfac * (prop_avg + corr) + (1 - av%sfac) * prop
      
      end subroutine smooth_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module smooth_stencil


