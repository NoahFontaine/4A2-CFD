      
      subroutine read_geom(av,geom)

!     Read in the curves specifying the case geometry

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_geometry), intent(out) :: geom
      
!     Declare an integer to count the number of points in the curve as you loop
!     INSERT
      integer :: i
!     END INSERT

!     Open the file and assign to unit 1
      open(1,file='geom_' // av%casename // '.txt')

!     Read the length of the first curve and allocate the length of the arrays
!     required to store its geometry within memory in the "geom" variable
      read(1,*) geom%ni_a
      allocate(geom%x_a(geom%ni_a),geom%y_a(geom%ni_a))

!     Read the x and y coordinates of the curve from the file, create a do loop
!     that iterates from 1 to "geom%ni_a", the length of the first domain curve
!     INSERT 
      do i = 1, geom%ni_a
         read(1,*) geom%x_a(i), geom%y_a(i)
      end do
!     END INSERT

!     Repeat the process for the second curve, read its length, allocate the
!     memory, then read the coordinates line by line from the file
!     INSERT
      read(1,*) geom%ni_b
      allocate(geom%x_b(geom%ni_b),geom%y_b(geom%ni_b))
      do i = 1, geom%ni_b
         read(1,*) geom%x_b(i), geom%y_b(i)
      end do
!     END INSERT

!     Print the lengths of the curves that have been successfully read
      write(6,*) 'Read domain curves from file'
      write(6,*) '  Curve lengths ni_a =', geom%ni_a, 'ni_b = ', geom%ni_b
      write(6,*)

!     Close the unit now everything has been read
!     INSERT
      close(1)
!     END INSERT

      end subroutine read_geom


