      
      subroutine apply_bconds(av,g,bcs)

!     This subroutine applies both the inlet and outlet boundary conditions, as
!     it modifies both the primary and secondary flow variables they must be
!     calculated first

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), intent(inout) :: g
      type(t_bconds), intent(inout) :: bcs

!     Declare the other variables you need here
!     INSERT

      real :: tstat(g%ni,g%nj), v(g%ni,g%nj)
      
!     At the inlet boundary the change in density is driven towards "rostag",
!     which is then used to obtain the other flow properties to match the
!     specified stagnation pressure, temperature and flow angle. 

!     To help prevent instabilities forming at the inlet boundary condition the 
!     changes in inlet density are relaxed by a factor "rfin" normally set to 
!     0.25 but it can be reduced further.

!     It is also worth checking if "ro" is greater than "rostag" and limiting 
!     the values to be slightly less than "rostag". This can prevent the solver 
!     crashing during severe transients.
      if(av%nstep == 1) then
          bcs%ro = g%ro(1,:)
      else
          bcs%ro = bcs%rfin * g%ro(1,:) + (1 - bcs%rfin) * bcs%ro
      endif
      bcs%ro = min(bcs%ro,0.9999 * bcs%rostag)

!     Calculate "p(1,:)", "rovx(1,:)", "rovy(1,:)" and "roe(1,:)" from the inlet 
!     "ro(:)", "pstag", "tstag" and "alpha". Also set "vx(1,:)", "vy(1,:)" and 
!     "hstag(1,:)"
!     INSERT

      tstat(1,:) = bcs%tstag * (bcs%ro(:) / bcs%rostag)**(av%gam - 1.0)
      g%p(1,:) = bcs%pstag * (bcs%ro(:) / bcs%rostag ) ** av%gam
      v(1,:) = sqrt(2.0 * av%cp * (bcs%tstag - tstat(1,:)))
      g%vx(1,:) = v(1,:) * cos(bcs%alpha); g%vy(1,:) = v(1,:) * sin(bcs%alpha)
      g%rovx(1,:) = bcs%ro(:) * g%vx(1,:); g%rovy(1,:) = bcs%ro * g%vy(1,:)
      g%roe(1,:) = bcs%ro(:) * (av%cv * tstat(1,:) + 0.5 * (v(1,:))**2.0)
      g%hstag(1,:) = av%cp * bcs%tstag

!     For the outlet boundary condition set the value of "p(ni,:)" to the
!     specified value of static pressure "p_out" in "bcs"
!     INSERT

      g%p(g%ni,:) = bcs%p_out
      
      end subroutine apply_bconds


