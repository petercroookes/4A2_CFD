      
      subroutine set_timestep(av,g,bcs)

!     This subroutine sets a single value for the time step based on the 
!     stagnation speed of sound and the minimum length scale of any element

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(inout) :: av
      type(t_grid), intent(in) :: g
      type(t_bconds), intent(in) :: bcs
!     Add 2D-arrays for local speed of sound and local velocity
      real, dimension(:,:), allocatable :: a_local, v_local, p_avg, ro_avg
      real, dimension(:,:), allocatable :: vx_avg, vy_avg
      integer :: ni, nj, i, j

      ni = g%ni; nj = g%nj
      
!     Calculate the stagnation speed of sound from the inlet stagnation
!     temperature and gas constants
!     INSERT

!     astag = sqrt(av%gam * av%rgas * bcs%tstag)
      
!     Assume that the maximum flow speed is also equal to "astag". This will 
!     be pessimistic for subsonic flows but may be optimistic for supersonic 
!     flows. In the latter case the length of the time step as determined by 
!     may need to be reduced by improving this routine or varying the CFL number
!     INSERT

!      v_max = astag * 2
      
!     Calculate local speed of sound and local velocity across grid
      allocate(p_avg(ni-1,nj-1),ro_avg(ni-1,nj-1),a_local(ni-1,nj-1))
      allocate(vx_avg(ni-1,nj-1),vy_avg(ni-1,nj-1),v_local(ni-1,nj-1))
!      p_avg(1:ni-1,1:nj-1) = 0.25 * (g%p(1:ni-1,1:nj-1) + g%p(2:ni,1:nj-1) &
!           g%p(1:ni-1,2:nj) + g%p(2:ni,2:nj))
!      ro_avg(1:ni-1,1:nj-1) = 0.25 * (g%ro(1:ni-1,1:nj-1) + g%ro(2:ni,1:nj-1) &
!           + g%ro(1:ni-1,2:nj) + g%ro(2:ni,2:nj))
      
      do i = 1, ni-1
         do j = 1, nj-1
            p_avg(i,j) = 0.25 * (g%p(i,j) + g%p(i+1,j) + g%p(i,j+1) + g%p(i+1,j+1))
            ro_avg(i,j) = 0.25 * (g%ro(i,j) + g%ro(i+1,j) + g%ro(i,j+1) + g%ro(i+1,j+1))
         end do
      end do
      a_local(:,:) = sqrt(av%gam * p_avg(:,:) / ro_avg(:,:))
      
!      vx_avg(1:ni-1,1:nj-1) = 0.25 * (g%vx(1:ni-1,1:nj-1) + g%vx(2:ni,1:nj-1) &
!           + g%vx(1:ni-1,2:nj) + g%vx(2:ni,2:nj))
!      vy_avg(1:ni-1,1:nj-1) = 0.25 * (g%vy(1:ni-1,1:nj-1) + g%vy(2:ni,1:nj-1) &
!           + g%vy(1:ni-1,2:nj) + g%vy(2:ni,2:nj))
!      v_local(:,:) = hypot(vx_avg(:,:), vy_avg(:,:))

      do i = 1, ni-1
         do j = 1, nj-1
            vx_avg(i,j) = 0.25 * (g%vx(i,j) + g%vx(i+1,j) + g%vx(i,j+1) + g%vx(i+1,j+1))
            vy_avg(i,j) = 0.25 * (g%vy(i,j) + g%vy(i+1,j) + g%vy(i,j+1) + g%vy(i+1,j+1))
         end do
      end do
      v_local(:,:) = hypot(abs(vx_avg(:,:)),abs(vy_avg(:,:)))
      
!     Calculate the timestep using the CFL number and store it in "av%dt"
!     dt swapped for dt_total here to implement Runge-Kutta
!     dt_total formula updated for spatially varying timestep

      av%dt_total(:,:) = av%cfl * g%l_min(:,:) / (a_local(:,:) + v_local(:,:))
      
!     Print the calculated timestep and some intermediate values
!     INSERT

!      write(6,*) 'Smallest timestep =', minval(av%dt_total)
!      write(6,*) 'Largest timestep =', maxval(av%dt_total)
      end subroutine set_timestep


