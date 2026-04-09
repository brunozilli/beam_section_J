c src/torsion_j.f
c Compute torsional constant J using FEM (Prandtl's stress function)
c Solves: ∇²φ = -2 on Ω, φ = 0 on boundary, then J = 2∫φ dA
c
c Authors: Bruno Zilli & DeepSeek
c Licence: MIT
c=======================================================================

      module torsion_j
      
      implicit none
      
      contains
      
c=======================================================================
      subroutine compute_torsion_J(nnode, ntri, conn, y, z, J_torsion,
     &                              verbose)
c=======================================================================
      implicit none
      
      integer, intent(in) :: nnode, ntri
      integer, intent(in) :: conn(3, ntri)
      double precision, intent(in) :: y(nnode), z(nnode)
      double precision, intent(out) :: J_torsion
      logical, intent(in), optional :: verbose
      
c     Local variables
      integer, allocatable :: bound(:), map(:)
      double precision, allocatable :: Kmat(:,:), rhs(:), phi_free(:)
      double precision, allocatable :: phi(:)
      integer :: nfree, i, j, k, m, n, i1, i2, i3, info
      double precision :: area, Ke(3,3), fe(3)
      double precision :: b1, c1, b2, c2, b3, c3
      double precision :: y1, y2, y3, z1, z2, z3
      logical :: verb
      double precision :: t1, t2
      
      verb = .false.
      if (present(verbose)) verb = verbose
      
      if (verb) then
        call cpu_time(t1)
        write(*,*) ''
        write(*,*) '=== Torsional Constant J (FEM) ==='
        write(*,*) '  Nodes:      ', nnode
        write(*,*) '  Triangles:  ', ntri
      end if
      
c     Step 1: find boundary nodes
      allocate(bound(nnode))
      call find_boundary_nodes(nnode, ntri, conn, bound, verb)
      
c     Step 2: count free nodes (interior)
      nfree = 0
      do i = 1, nnode
        if (bound(i) .eq. 0) nfree = nfree + 1
      end do
      
c     Step 3: build mapping from global node to free node index
      allocate(map(nnode))
      j = 0
      do i = 1, nnode
        if (bound(i) .eq. 0) then
          j = j + 1
          map(i) = j
        else
          map(i) = 0
        end if
      end do
      
      if (verb) then
        write(*,*) '  Free nodes:   ', nfree
        write(*,*) '  Boundary nodes:', nnode - nfree
      end if
      
c     Step 4: if no interior nodes, use analytical approximation
      if (nfree .eq. 0) then
        if (verb) write(*,*) '  WARNING: No free nodes, using approximation'
        call analytical_J_approx(nnode, ntri, conn, y, z, J_torsion, verb)
        deallocate(bound, map)
        return
      end if
      
c     Step 5: allocate global system
      allocate(Kmat(nfree, nfree), rhs(nfree))
      allocate(phi_free(nfree), phi(nnode))
      
      do i = 1, nfree
        rhs(i) = 0.0d0
        do j = 1, nfree
          Kmat(i,j) = 0.0d0
        end do
      end do
      
c     Step 6: assemble element matrices
      do k = 1, ntri
        i1 = conn(1, k)
        i2 = conn(2, k)
        i3 = conn(3, k)
        
        y1 = y(i1)
        y2 = y(i2)
        y3 = y(i3)
        z1 = z(i1)
        z2 = z(i2)
        z3 = z(i3)
        
c       Element area
        area = 0.5d0 * abs((y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))
        
        if (area .lt. 1.0d-12) then
          if (verb) write(*,*) '  Warning: degenerate triangle', k, ' skipped'
          cycle
        end if
        
c       Coefficients for stiffness matrix
        b1 = y2 - y3
        c1 = z3 - z2
        b2 = y3 - y1
        c2 = z1 - z3
        b3 = y1 - y2
        c3 = z2 - z1
        
c       Element stiffness matrix (3x3)
        Ke(1,1) = (b1*b1 + c1*c1) / (4.0d0*area)
        Ke(1,2) = (b1*b2 + c1*c2) / (4.0d0*area)
        Ke(1,3) = (b1*b3 + c1*c3) / (4.0d0*area)
        Ke(2,1) = Ke(1,2)
        Ke(2,2) = (b2*b2 + c2*c2) / (4.0d0*area)
        Ke(2,3) = (b2*b3 + c2*c3) / (4.0d0*area)
        Ke(3,1) = Ke(1,3)
        Ke(3,2) = Ke(2,3)
        Ke(3,3) = (b3*b3 + c3*c3) / (4.0d0*area)
        
c       Element load vector (RHS) for f = 2
        fe(1) = 2.0d0 * area / 3.0d0
        fe(2) = 2.0d0 * area / 3.0d0
        fe(3) = 2.0d0 * area / 3.0d0
        
c       Assemble into global matrix
        do i = 1, 3
          m = map(conn(i, k))
          if (m .gt. 0) then
            rhs(m) = rhs(m) + fe(i)
            do j = 1, 3
              n = map(conn(j, k))
              if (n .gt. 0) then
                Kmat(m, n) = Kmat(m, n) + Ke(i, j)
              end if
            end do
          end if
        end do
      end do
      
      if (verb) write(*,*) '  Solving with LAPACK DPOSV...'
      
c     Add small epsilon to diagonal to ensure positive definiteness
      do i = 1, nfree
        Kmat(i,i) = Kmat(i,i) + 1.0d-12
      end do
      
c     Step 7: solve linear system K * phi = rhs
      phi_free = rhs
      call DPOSV('U', nfree, 1, Kmat, nfree, phi_free, nfree, info)
      
      if (info .ne. 0) then
        write(*,*) '  ERROR: DPOSV failed, info =', info
        J_torsion = 0.0d0
        deallocate(bound, map, Kmat, rhs, phi_free, phi)
        return
      end if
      
c     Step 8: reconstruct full phi vector
      do i = 1, nnode
        phi(i) = 0.0d0
        if (bound(i) .eq. 0) then
          phi(i) = phi_free(map(i))
        end if
      end do
      
c     Step 9: compute J = 2 * ∫ phi dA
      J_torsion = 0.0d0
      do k = 1, ntri
        i1 = conn(1, k)
        i2 = conn(2, k)
        i3 = conn(3, k)
        
        y1 = y(i1)
        y2 = y(i2)
        y3 = y(i3)
        z1 = z(i1)
        z2 = z(i2)
        z3 = z(i3)
        
        area = 0.5d0 * abs((y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))
        
        J_torsion = J_torsion + (area / 3.0d0) * (phi(i1) + phi(i2) + phi(i3))
      end do
      
      J_torsion = 2.0d0 * J_torsion
      
      if (verb) then
        call cpu_time(t2)
        write(*,*) '  J = ', J_torsion, ' mm⁴'
        write(*,*) '  Time: ', t2 - t1, ' seconds'
        write(*,*) ''
      end if
      
      deallocate(bound, map, Kmat, rhs, phi_free, phi)
      
      end subroutine compute_torsion_J
      
c=======================================================================
      subroutine find_boundary_nodes(nnode, ntri, conn, bound, verbose)
c=======================================================================
c     Identify boundary nodes (nodes belonging to an edge that appears
c     only once in the mesh)
c=======================================================================
      implicit none
      
      integer, intent(in) :: nnode, ntri
      integer, intent(in) :: conn(3, ntri)
      integer, intent(out) :: bound(nnode)
      logical, intent(in), optional :: verbose
      
      integer :: i, j, k, i1, i2, i3, cnt
      integer, allocatable :: ecount(:,:)
      logical :: verb
      
      verb = .false.
      if (present(verbose)) verb = verbose
      
c     Initialise
      do i = 1, nnode
        bound(i) = 0
      end do
      
c     Allocate edge counter
      allocate(ecount(nnode, nnode))
      do i = 1, nnode
        do j = 1, nnode
          ecount(i, j) = 0
        end do
      end do
      
c     Count edge occurrences
      do k = 1, ntri
        i1 = conn(1, k)
        i2 = conn(2, k)
        i3 = conn(3, k)
        call inc_edge_count(i1, i2, ecount, nnode)
        call inc_edge_count(i2, i3, ecount, nnode)
        call inc_edge_count(i3, i1, ecount, nnode)
      end do
      
c     Mark boundary nodes (edges with count = 1)
      do k = 1, ntri
        i1 = conn(1, k)
        i2 = conn(2, k)
        i3 = conn(3, k)
        if (get_edge_count(i1, i2, ecount, nnode) .eq. 1) then
          bound(i1) = 1
          bound(i2) = 1
        end if
        if (get_edge_count(i2, i3, ecount, nnode) .eq. 1) then
          bound(i2) = 1
          bound(i3) = 1
        end if
        if (get_edge_count(i3, i1, ecount, nnode) .eq. 1) then
          bound(i3) = 1
          bound(i1) = 1
        end if
      end do
      
      if (verb) then
        cnt = 0
        do i = 1, nnode
          if (bound(i) .eq. 1) cnt = cnt + 1
        end do
        write(*,*) '  Boundary nodes identified:', cnt
      end if
      
      deallocate(ecount)
      
      end subroutine find_boundary_nodes
      
c=======================================================================
      subroutine inc_edge_count(i, j, ecount, nnode)
c=======================================================================
      implicit none
      
      integer, intent(in) :: i, j, nnode
      integer, intent(inout) :: ecount(nnode, nnode)
      integer :: imin, imax
      
      imin = min(i, j)
      imax = max(i, j)
      ecount(imin, imax) = ecount(imin, imax) + 1
      
      end subroutine inc_edge_count
      
c=======================================================================
      function get_edge_count(i, j, ecount, nnode) result(cnt)
c=======================================================================
      implicit none
      
      integer, intent(in) :: i, j, nnode
      integer, intent(in) :: ecount(nnode, nnode)
      integer :: cnt
      integer :: imin, imax
      
      imin = min(i, j)
      imax = max(i, j)
      cnt = ecount(imin, imax)
      
      end function get_edge_count
      
c=======================================================================
      subroutine analytical_J_approx(nnode, ntri, conn, y, z, J_torsion,
     &                                verbose)
c=======================================================================
c     Analytical approximation for coarse meshes with no interior nodes
c=======================================================================
      implicit none
      
      integer, intent(in) :: nnode, ntri
      integer, intent(in) :: conn(3, ntri)
      double precision, intent(in) :: y(nnode), z(nnode)
      double precision, intent(out) :: J_torsion
      logical, intent(in), optional :: verbose
      
      double precision :: ymin, ymax, zmin, zmax, width, height
      logical :: verb
      integer :: i
      
      verb = .false.
      if (present(verbose)) verb = verbose
      
      ymin = y(1)
      ymax = y(1)
      zmin = z(1)
      zmax = z(1)
      
      do i = 2, nnode
        if (y(i) .lt. ymin) ymin = y(i)
        if (y(i) .gt. ymax) ymax = y(i)
        if (z(i) .lt. zmin) zmin = z(i)
        if (z(i) .gt. zmax) zmax = z(i)
      end do
      
      width = ymax - ymin
      height = zmax - zmin
      
      if (verb) then
        write(*,*) '  Bounding box: ', width, ' x ', height, ' mm'
      end if
      
c     Approximate J based on shape
      if (abs(width - height) .lt. 1.0d-6) then
c       Square/circle-like
        J_torsion = 3.141592653589793d0 * (width/2.0d0)**4 / 2.0d0
        if (verb) write(*,*) '  Using circular approximation'
      else
c       Rectangle-like
        if (height .lt. width) then
          J_torsion = width * height**3 / 3.0d0
        else
          J_torsion = height * width**3 / 3.0d0
        end if
        if (verb) write(*,*) '  Using rectangular approximation'
      end if
      
      if (verb) write(*,*) '  Approximate J = ', J_torsion, ' mm⁴'
      
      end subroutine analytical_J_approx
      
c==============================================================================
      subroutine read_mesh_wrapper(filename, nnode, ntri, conn, y, z)
c==============================================================================
c     Wrapper to read UNV mesh with allocatable arrays
c==============================================================================
      implicit none
      character(len=*), intent(in) :: filename
      integer, intent(out) :: nnode, ntri
      integer, allocatable, intent(out) :: conn(:,:)
      double precision, allocatable, intent(out) :: y(:), z(:)
      
      integer, parameter :: MAX_NODES = 10000
      integer, parameter :: MAX_TRI = 20000
      double precision :: ytmp(MAX_NODES), ztmp(MAX_NODES)
      integer :: conntmp(3, MAX_TRI)
      
      call read_section_mesh_unv(trim(filename), MAX_NODES, MAX_TRI,
     &                           nnode, ntri, conntmp, ytmp, ztmp)
      
      allocate(y(nnode), z(nnode), conn(3, ntri))
      y(1:nnode) = ytmp(1:nnode)
      z(1:nnode) = ztmp(1:nnode)
      conn(:, 1:ntri) = conntmp(:, 1:ntri)
      
      end subroutine read_mesh_wrapper
      
      end module torsion_j
