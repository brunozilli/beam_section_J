c=======================================================================
c mesh_checker.f
c
c Check and fix triangular mesh: orientation, degeneracy, validity.
c Industrial-grade validation for FEM pipelines.
c
c Authors: Bruno Zilli & DeepSeek
c Licence: MIT
c=======================================================================

      subroutine check_and_fix_mesh(
     &     nnode, ntri,
     &     conn, y, z,
     &     A_total,
     &     n_flipped, n_degenerate, ierr)

      implicit none

c     INPUT
      integer, intent(in) :: nnode, ntri
      integer, intent(inout) :: conn(3, ntri)
      double precision, intent(in) :: y(nnode), z(nnode)

c     OUTPUT
      double precision, intent(out) :: A_total
      integer, intent(out) :: n_flipped, n_degenerate, ierr

c     LOCALS
      integer :: i, i1, i2, i3, tmp
      double precision :: y1, y2, y3, z1, z2, z3
      double precision :: det, area
      double precision, parameter :: tol = 1.0d-14

c     Initialise
      A_total = 0.0d0
      n_flipped = 0
      n_degenerate = 0
      ierr = 0

      if (ntri .le. 0) then
         write(*,*) 'ERROR: no triangles in mesh'
         ierr = 1
         return
      endif

c     Loop over all triangles
      do i = 1, ntri

c       Connectivity check
        i1 = conn(1,i)
        i2 = conn(2,i)
        i3 = conn(3,i)

        if (i1 .lt. 1 .or. i1 .gt. nnode .or.
     &      i2 .lt. 1 .or. i2 .gt. nnode .or.
     &      i3 .lt. 1 .or. i3 .gt. nnode) then
            write(*,*) 'ERROR: invalid node index in element', i
            ierr = 2
            return
        endif

c       Coordinates
        y1 = y(i1)
        z1 = z(i1)
        y2 = y(i2)
        z2 = z(i2)
        y3 = y(i3)
        z3 = z(i3)

c       Determinant (2 x area)
        det = (y2 - y1)*(z3 - z1) - (y3 - y1)*(z2 - z1)

c       Check for degenerate triangle (zero area)
        if (dabs(det) .lt. tol) then
            n_degenerate = n_degenerate + 1
            cycle
        endif

c       Fix orientation: ensure all triangles have positive area
        if (det .lt. 0.0d0) then
c         Swap node 2 and node 3 to flip orientation
            tmp = conn(2,i)
            conn(2,i) = conn(3,i)
            conn(3,i) = tmp

            n_flipped = n_flipped + 1
            det = -det
        endif

c       Accumulate total area
        area = 0.5d0 * det
        A_total = A_total + area

      end do

c     Final sanity check
      if (A_total .lt. tol) then
         write(*,*) 'WARNING: total area near zero'
         ierr = 3
      endif

c     Diagnostic output
      if (n_flipped .gt. 0) then
         write(*,*) 'MESH CHECKER: flipped', n_flipped, 'triangles'
      endif
      if (n_degenerate .gt. 0) then
         write(*,*) 'MESH CHECKER: skipped', n_degenerate, 
     &              'degenerate triangles'
      endif
      write(*,*) 'MESH CHECKER: total area =', A_total

      return
      end subroutine check_and_fix_mesh
