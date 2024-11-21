      
      subroutine check_mesh(g)

!     Check the cell area and facet length calculations before attempting to
!     solve the flow, make sure you do this for both the "bump" and "bend" test
!     cases

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_grid), intent(inout) :: g
      real :: area_min, dx_error, dy_error, tol
      integer :: ni, nj, i, j
      real :: x1, y1, x2, y2, x3, y3, x4, y4
      real :: lx_i2, ly_i2, lx_j2, ly_j2
      real, allocatable :: proj_error(:,:)

!     Get the size of the mesh and store locally for convenience
      ni = g%ni; nj = g%nj;

!     Exact checking of floating point numbers never goes well, define a
!     small tolerance to use for a comparative operation instead
      tol = 1e-3 * g%l_min

!     Check that all of the cell areas are positive, either with the intrinsic
!     "minval" function or with nested do loops. Print the output to the screen
!     and flag negative numbers as an error with an if statement to "stop" the
!     program
      !     INSERT

      do i = 1, ni-1
         do j = 1, nj-1
            if (g%area(i,j) < 0) then
               print *, 'Area value is negative. Cell Area is', g%area(i,j), &
                    'at i, j =', i, j
               stop
            end if
         end do
      end do

!     Next check that the sum of the edge vectors around every quadrilateral is 
!     very nearly zero in both the x and y-coordinate directions. You can
!     complete this with some elementwise addition of the arrays and use of the
!     "maxval" and "abs" intrinsic functions.
!     INSERT

      allocate(proj_error(ni-1, nj-1))
      
      do i = 1, ni-2
         do j = 1, nj-2
            dx_error = g%lx_i(i,j+1) - g%lx_i(i,j) + g%lx_j(i+1,j) - g%lx_j(i,j)
            dy_error = g%ly_i(i,j+1) - g%ly_i(i,j) + g%ly_j(i+1,j) - g%ly_j(i,j)
            proj_error(i,j) = max(abs(dx_error),abs(dy_error))
            
            if (proj_error(i,j) > tol) then
               !print *, 'Warning: Sum of projected facet vectors too large. &
               !Length =', proj_error(i,j), 'i,j =', i, j
            end if
         end do
      end do

!     It may be worthwhile to complete some other checks, the prevous call to
!     the "write_output" subroutine has written a file that you can read and
!     postprocess using the Python script plot_mesh.py. This program also has
!     access to all of the mesh parameters used within the solver that you could
!     inspect graphically.

!     Print a blank line
      write(6,*)

      deallocate(proj_error)
      end subroutine check_mesh
