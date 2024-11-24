      
      subroutine calc_areas(g)

!     Calculate the area of the quadrilateral cells and the lengths of the side
!     facets

!     Explicitly declare the required variables
      use types
      use routines
      implicit none
      type(t_grid), intent(inout) :: g
      integer :: ni, nj

!     Declare integers or any extra variables you need here
!     INSERT

      integer :: i, j
      real :: x1, x2, x3, x4, y1, y2, y3, y4, d1_x, d1_y, d2_x, d2_y

      ! Declare lengths_i and lengths_j as allocatable arrays
      real, allocatable :: lengths_i(:,:), lengths_j(:,:)
      
!     Get the size of the mesh and store locally for convenience
      ni = g%ni; nj = g%nj;

!     Calculate the areas of the cells and store in g%area. The area of any
!     quadrilateral is half of the magnitude of the cross product of the two
!     vectors that form the diagonals. Check the order of your product so that
!     the values come out positive! You can do this using two nested loops in
!     the i and j-directions or in a vectorised way by indexing the coordinate
!     arrays with lists of indices
!     INSERT

      do j = 1, nj-1
         do i = 1, ni-1

            ! Extract corner points of the current cell
            x1 = g%x(i, j)
            y1 = g%y(i, j)
            x2 = g%x(i+1, j)
            y2 = g%y(i+1, j)
            x3 = g%x(i+1, j+1)
            y3 = g%y(i+1, j+1)
            x4 = g%x(i, j+1)
            y4 = g%y(i, j+1)

            ! Compute diagonal vectors
            d1_x = x3 - x1
            d1_y = y3 - y1
            d2_x = x4 - x2
            d2_y = y4 - y2

            ! Calculate the area of the quadrilateral
            g%area(i, j) = 0.5 * abs(d1_x * d2_y - d1_y * d2_x)

         end do
      end do
      
!     Calculate the projected lengths in the x and y-directions on all of the
!     "i = const" facets and store them in g%lx_i and g%ly_i. When combined
!     together these two components define a vector that is normal to the facet,
!     pointing inwards towards the centre of the cell. This is only the case for
!     the left hand side of the cell, the vector stored in position i,j points
!     towards the centre of the i,j cell
!     INSERT

      g%lx_i(1:ni,1:nj-1) = g%y(1:ni,2:nj) - g%y(1:ni,1:nj-1)
      g%ly_i(1:ni,1:nj-1) = g%x(1:ni,1:nj-1) - g%x(1:ni,2:nj)

      
!     Now repeat the calculation for the project lengths on the "j=const"
!     facets. 
!     INSERT

      g%lx_j(1:ni-1,1:nj) = g%y(1:ni-1,1:nj) - g%y(2:ni,1:nj)
      g%ly_j(1:ni-1,1:nj) = g%x(2:ni,1:nj) - g%x(1:ni-1,1:nj)

      
!     Find the minimum length scale in the mesh, this is defined as the length
!     of the shortest side of all the cells. Call this length "l_min", it is used
!     to set the timestep from the CFL number. Start by calculating the lengths
!     of the i and j facets by using the intrinsic function "hypot", this avoids
!     underflow and overflow errors. Then find the overal minimum value using
!     both the "min" and "minval" functions.
!     INSERT

      ! Allocate temporarily arrays holding mesh cell side lengths in i and j directions
      allocate(lengths_i(ni-1,nj-1), lengths_j(ni-1,nj-1))
      
      do j = 1, nj-1
         do i = 1, ni-1

            ! Extract corner points of the current cell
            x1 = g%x(i, j)
            y1 = g%y(i, j)
            x2 = g%x(i+1, j)
            y2 = g%y(i+1, j)
            x3 = g%x(i+1, j+1)
            y3 = g%y(i+1, j+1)
            x4 = g%x(i, j+1)
            y4 = g%y(i, j+1)

            ! Calculate the side lenghts
            lengths_i(i,j) = hypot((x2-x1),(y2-y1))
            lengths_j(i,j) = hypot((x4-x1),(y4-y1)) 
            
         end do
      end do

      ! Calculate minimum length and deallocate the lengths arrays
      ! Spatially varying timestep implemented:
      
      do j = 1, nj-2
         do i = 1, ni-2
            g%l_min(i,j) = min(lengths_i(i,j),lengths_i(i,j+1),&
                 lengths_j(i,j), lengths_j(i+1,j))
         end do
      end do

      do i = 1, ni-2
         do j = nj-1, nj-1
            g%l_min(i,j) = min(lengths_i(i,j),lengths_j(i,j),lengths_j(i+1,j))
         end do
      end do

      do j = 1, nj-2
         do i = ni-1, ni-1
            g%l_min(i,j) = min(lengths_i(i,j),lengths_j(i,j),lengths_i(i,j+1))
         end do
      end do

      g%l_min(ni-1,nj-1) = min(lengths_i(ni-1,nj-1),lengths_j(ni-1,nj-1))
         
      deallocate(lengths_i, lengths_j)
      
!     Print the overall minimum length size that has been calculated
      write(6,*) 'Calculated cell areas and facet lengths'
      write(6,*) '  Overall minimum element size = ', minval(g%l_min)
      write(6,*) 'First values of projected lengths:', g%lx_i(1,1), g%ly_i(1,1), &
           g%lx_j(1,1), g%ly_j(1,1)
      write(6,*) 'Shape of projected length arrays:', shape(g%lx_i), shape(g%ly_i), &
           shape(g%lx_j), shape(g%ly_j)
      
      end subroutine calc_areas
