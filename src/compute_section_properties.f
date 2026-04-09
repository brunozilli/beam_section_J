c=======================================================================
c compute_section_properties.f
c
c Compute cross-section properties (area, centroid, moments of inertia)
c using a triangular mesh. Uses exact analytical formulae for triangles.
c
c Authors: Bruno Zilli & DeepSeek
c Licence: MIT
c=======================================================================

      subroutine compute_section_properties(nnode_sec, ntri_sec,
     &                                      conn_sec, y_sec, z_sec,
     &                                      A, y_c, z_c,
     &                                      I_y, I_z, I_yz, J_polar)

      implicit none

c     INPUT
      integer, intent(in) :: nnode_sec, ntri_sec
      integer, intent(in) :: conn_sec(3, ntri_sec)
      double precision, intent(in) :: y_sec(nnode_sec)
      double precision, intent(in) :: z_sec(nnode_sec)

c     OUTPUT
      double precision, intent(out) :: A, y_c, z_c
      double precision, intent(out) :: I_y, I_z, I_yz
      double precision, intent(out) :: J_polar

c     LOCALS
      integer :: i, i1, i2, i3
      double precision :: y1, y2, y3, z1, z2, z3
      double precision :: area_tri
      double precision :: S_y, S_z
      double precision :: Iyy0, Izz0, Iyz0
      double precision :: det

c     Initialise output variables
      A   = 0.0d0
      S_y = 0.0d0
      S_z = 0.0d0
      I_y = 0.0d0
      I_z = 0.0d0
      I_yz = 0.0d0

      if (ntri_sec .le. 0) return

c     Loop over all triangular elements
      do i = 1, ntri_sec

        i1 = conn_sec(1,i)
        i2 = conn_sec(2,i)
        i3 = conn_sec(3,i)

        y1 = y_sec(i1)
        z1 = z_sec(i1)
        y2 = y_sec(i2)
        z2 = z_sec(i2)
        y3 = y_sec(i3)
        z3 = z_sec(i3)

c       Oriented area (2 x area) - sign indicates orientation
        det = (y2 - y1)*(z3 - z1) - (y3 - y1)*(z2 - z1)
        area_tri = 0.5d0 * det

c       Accumulate total area and first moments
        A = A + area_tri
        S_y = S_y + area_tri * (y1 + y2 + y3) / 3.0d0
        S_z = S_z + area_tri * (z1 + z2 + z3) / 3.0d0

c       Second moments about origin (0,0) using exact triangle formulae
c       Iyy = ∫ z² dA
        Iyy0 = area_tri *
     &        (z1*z1 + z2*z2 + z3*z3 +
     &         z1*z2 + z2*z3 + z3*z1) / 6.0d0

c       Izz = ∫ y² dA
        Izz0 = area_tri *
     &        (y1*y1 + y2*y2 + y3*y3 +
     &         y1*y2 + y2*y3 + y3*y1) / 6.0d0

c       Iyz = ∫ y*z dA
        Iyz0 = area_tri *
     &        (2.0d0*(y1*z1 + y2*z2 + y3*z3) +
     &         y1*z2 + y2*z1 +
     &         y2*z3 + y3*z2 +
     &         y3*z1 + y1*z3) / 12.0d0

        I_y  = I_y  + Iyy0
        I_z  = I_z  + Izz0
        I_yz = I_yz + Iyz0

      end do

c     Check for valid area
      if (dabs(A) .lt. 1.0d-16) then
         write(*,*) 'ERROR: zero or invalid section area'
         return
      endif

c     Section centroid
      y_c = S_y / A
      z_c = S_z / A

c     Parallel axis theorem (Steiner) - shift to centroid
      I_y  = I_y  - A * z_c*z_c
      I_z  = I_z  - A * y_c*y_c
      I_yz = I_yz - A * y_c*z_c

c     Polar moment of inertia (not the torsion constant!)
      J_polar = I_y + I_z

      return
      end subroutine compute_section_properties
