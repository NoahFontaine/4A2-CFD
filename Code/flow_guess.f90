      
      subroutine flow_guess(av,g,bcs,guesstype)

!     This calculates an initial guess of the primary flowfield variables and
!     stores them at the nodes within the mesh datatype

!     Explicitly declare the required variables
      use types
      use routines
      implicit none
      type(t_appvars), intent(inout) :: av
      type(t_grid), intent(inout) :: g
      type(t_bconds), intent(in) :: bcs
      integer, intent(in) :: guesstype
      integer :: i, j, ni, nj, j_mid
      
!     Variables required for the crude guess
      real :: t_out, v_out, ro_out, lx, ly, l

!     Variables required for the improved guess, you will need to add to these
      real :: l_i(g%ni), m_out, mach_lim, t_lim, v_guess(g%ni), ro_guess(g%ni), t_guess(g%ni), p_guess(g%ni)
!     INSERT

!     Get the size of the mesh and store locally for convenience
      ni = g%ni; nj = g%nj;

!     Assuming isentropic flow to the exit plane calculate the static
!     temperature and the exit velocity
      t_out = bcs%tstag * (bcs%p_out / bcs%pstag)**av%fgam
      v_out = (2 * av%cp * (bcs%tstag - t_out))**0.5
      ro_out = bcs%p_out / (av%rgas * t_out)

!     Determine which guess calculation method to use by the value of "guesstype"
      if(guesstype == 1) then

!         Store the exit density and internal energy as if they were uniform 
          g%ro = ro_out 
          g%roe  = g%ro * (av%cv * t_out + 0.5 * v_out**2)

!         Calculate the gradient of the mesh lines in the centre of the domain
!         to determine the assumed direction of the flow
          j_mid = nj / 2
          do i = 1,ni-1
              lx = g%lx_j(i,j_mid); ly = g%ly_j(i,j_mid); 
              l = hypot(lx,ly)
              g%rovx(i,:) = g%ro(i,:) * v_out * ly / l
              g%rovy(i,:) = -g%ro(i,:) * v_out * lx / l
          end do

!         Copy the values to the "i = ni" nodes as an approximation
          g%rovx(ni,:) = g%rovx(ni-1,:)
          g%rovy(ni,:) = g%rovy(ni-1,:)

!         Print the guess that has been calculated
          write(6,*) 'Crude flow guess calculated'
          write(6,*) '  At first point ro =', g%ro(1,1), 'roe =', &
              g%roe(1,1), 'rovx =', g%rovx(1,1), 'rovy =', g%rovy(1,1)
          write(6,*)

      else if(guesstype == 2) then 

!         Calculate the length of each "i = const" line between the "j = 1" and 
!         "j = nj" boundaries of the domain and store it in the local variable
!         "l_i". You could calculate the length along each i-facet from the x 
!         and y projected lengths with "hypot" and then sum them up in the
!         second dimension with "sum". 
!         INSERT
          l_i = sum(hypot(g%lx_i(1:ni,1:nj-1), g%ly_i(1:ni,1:nj-1)), dim=2)
!         END INSERT

!         Use the exit temperature, density and velocity calculated for the 
!         crude guess with "l_i" to estimate the mass flow rate at the exit
!         INSERT
          m_out = ro_out * v_out * l_i(ni)
!         END INSERT

!         Set a limit to the maximum allowable mach number in the initial
!         guess, call this "mach_lim", calculate the corresponding temperature,
!         called "t_lim"
!         INSERT
          mach_lim = 0.1
          t_lim = bcs%tstag / (1 + 0.5 * (mach_lim**2) * (av%gam - 1))
!         END INSERT
!         Now estimate the velocity and density at every "i = const" line, call 
!         the velocity "v_guess(i)" and the density "ro_guess(i)":
!             1. Assume density is constant at the exit value
!             2. Use continuity and "l_i(i)" to estimate the velocity
!             3. Assume stagnation temperature is constant for static temp
!             4. Limit the static temperature, lookup intrinsic "max"
!             5. Calculate the density throughout "ro_guess(i)"
!             6. Update the estimate of the velocity "v_guess(i)" 
!         INSERT
          do i = 1, ni
              v_guess(i) = m_out / (ro_out * l_i(i))
              t_guess(i) = bcs%tstag - 0.5 * v_guess(i)**2 / av%cp
              t_guess(i) = max(t_guess(i), t_lim)
              p_guess(i) = bcs%pstag * (t_guess(i) / bcs%tstag)**(av%gam/(av%gam-1))
              ro_guess(i) = p_guess(i) / (av%rgas * t_guess(i))
              v_guess(i) = m_out / (ro_guess(i) * l_i(i))
          end do
!         END INSERT

!         Direct the calculated velocity to be parallel to the "j = const"
!         gridlines for all values of i and j. This can be achieved with a 
!         similar calculation to the "j = nj/2" one that was performed in the 
!         crude guess. Then set all of ro, roe, rovx and rovy, note that roe 
!         includes the kinetic energy component of the internal energy.
!         INSERT
          do i = 1, ni-1
              g%ro(i,:) = ro_guess(i)
              g%roe(i,:) = g%ro(i,:) * (av%cv * t_guess(i) + 0.5 * v_guess(i)**2)
          end do

          j_mid = nj / 2
          do i = 1,ni-1
              lx = g%lx_j(i,j_mid); ly = g%ly_j(i,j_mid); 
              l = hypot(lx,ly)
              g%rovx(i,:) = g%ro(i,:) * v_guess(i) * ly / l
              g%rovy(i,:) = -g%ro(i,:) * v_guess(i) * lx / l
          end do
!         END INSERT

!         Make sure the guess has been copied for the "i = ni" values too
!         INSERT
          g%ro(ni,:) = ro_guess(ni)
          g%roe(ni,:) = g%ro(ni,:) * (av%cv * t_guess(ni) + 0.5 * v_guess(ni)**2)
          g%rovx(ni,:) = g%rovx(ni-1,:)
          g%rovy(ni,:) = g%rovy(ni-1,:)
!         END INSERT

!         Print the first elements of the guess like for the crude guess
!         INSERT
          write(6,*) 'Improved flow guess calculated'
          write(6,*) '  At first point ro =', g%ro(1,1), 'roe =', &
                g%roe(1,1), 'rovx =', g%rovx(1,1), 'rovy =', g%rovy(1,1)
          write(6,*)
!         END INSERT

      end if

!     The initial guess values derived from the boundary conditions are also
!     useful as a reference to non-dimensionalise the convergence
      av%ro_ref = sum(g%ro(1,:)) / nj
      av%roe_ref = sum(g%roe(1,:)) / nj
      av%rov_ref = max(sum(g%rovx(1,:)),sum(g%rovy(1,:))) / nj

      end subroutine flow_guess


