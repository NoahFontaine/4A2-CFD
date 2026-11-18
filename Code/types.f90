
      module types

!     Define the derived data types used in the main program and subroutines,
!     you many need to create further variables for specific extensions

!     Appvars type contains all application variables and gas constants
      type t_appvars

!         Case name
          character(len=:), allocatable :: casename

!         Gas Properties
          real :: rgas, gam, cp, cv, fgam

!         Timestepping, smoothing and other run options
          real ::  cfl, sfac, dt, dt_total, d_max, d_avg
          integer :: nsteps, nstep

!         Reference values of the primary flow variables
          real :: ro_ref, roe_ref, rov_ref

!         Size of the mesh for single block test cases
          integer :: ni, nj

!         Number of blocks and matching patches for multi-block extension
          integer :: nn, nm

      end type t_appvars

!     Boundary condition type contains inlet and outlet data
      type t_bconds

!         Single value floats at the inlet
          real :: pstag, tstag, alpha, rfin, rostag

!         Vectors of floats at the inlet
          real, dimension(:), allocatable :: p, ro

!         Single value float at the outlet
          real :: p_out

!         Block numbers of the inlet and outlet in multi-block extension
          integer :: n_in, n_out

      end type t_bconds

!     Matching patch type contains block to block interface data
      type t_match

!         Length of the matching patch and adjoining block numbers
          integer :: nk, n_1, n_2

!         Lists of indices that are coincident on both sides
          integer, dimension(:), allocatable :: i_1, j_1, i_2, j_2

      end type t_match

!     Geometry type contains the boundaries of the case geometry
      type t_geometry

!         Curve length
          integer :: ni_a, ni_b

!         Coordinate data for upper and lower domain boundaries
          real, dimension(:), allocatable :: x_a, y_a, x_b, y_b

      end type t_geometry

!     Grid type contains coordinate data and all flow variables
      type t_grid

!         Mesh size for all cases
          integer :: ni, nj

!         Mesh coordinate data in 2D matrices
          real, dimension(:,:), allocatable :: x, y, area, lx_i, ly_i, &
              lx_j, ly_j
          real  ::  l_min

!         Primary variables at nodes
          real, dimension(:,:), allocatable :: ro, roe, rovx, rovy

!         Variables to hold cell increments
          real, dimension(:,:), allocatable :: dro, droe, drovx, drovy

!         Secondary variables at nodes
          real, dimension(:,:), allocatable :: p, hstag, vx, vy 

!         Logical array to store wall locations for the nodes
          logical, dimension(:,:), allocatable :: wall

      end type t_grid

      end module types


