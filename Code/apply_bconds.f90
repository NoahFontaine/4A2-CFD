      
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
      real :: t_stat(g%nj)
!     INSERT

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
      g%p(1,:) = bcs%pstag * (bcs%ro / bcs%rostag)**av%fgam
      t_stat = g%p(1,:) / (av%rgas * bcs%ro)
      g%rovx(1,:) = bcs%ro * sqrt(2 * av%cp * (bcs%tstag - t_stat)) * cos(bcs%alpha)
      g%rovy(1,:) = bcs%ro * sqrt(2 * av%cp * (bcs%tstag - t_stat)) * sin(bcs%alpha)
      g%roe(1,:) = bcs%ro * (av%cv * bcs%tstag + 0.5 * (g%rovx(1,:)**2 + g%rovy(1,:)**2) / bcs%ro)
      g%vx(1,:) = g%rovx(1,:) / bcs%ro
      g%vy(1,:) = g%rovy(1,:) / bcs%ro
      g%hstag(1,:) = av%cp * (g%p(1,:) / (av%rgas * bcs%ro)) + 0.5 * (g%vx(1,:)**2 + g%vy(1,:)**2)
!     END INSERT

!     For the outlet boundary condition set the value of "p(ni,:)" to the
!     specified value of static pressure "p_out" in "bcs"
!     INSERT
      g%p(g%ni,:) = bcs%p_out
!     END INSERT

      end subroutine apply_bconds


