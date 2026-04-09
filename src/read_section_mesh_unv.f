c=======================================================================
c read_section_mesh_unv.f
c
c Read UNV mesh file (Universal File Format)
c Supports both simple format and Salome/ASTER format with 2-line nodes
c
c Authors: Bruno Zilli & DeepSeek
c Licence: MIT
c=======================================================================

      subroutine read_section_mesh_unv(filename, max_nodes, 
     &                                 max_triangles, n_nodes, 
     &                                 n_triangles, conn, y, z)

      implicit none
      
c     INPUT
      character(len=*), intent(in) :: filename
      integer, intent(in) :: max_nodes, max_triangles
      
c     OUTPUT
      integer, intent(out) :: n_nodes, n_triangles
      integer, intent(out) :: conn(3, max_triangles)
      double precision, intent(out) :: y(max_nodes), z(max_nodes)
      
c     LOCALS
      integer :: io, node_id, elem_id, elem_type
      integer :: n1, n2, n3, dummy, phys_prop
      character(len=256) :: line
      logical :: in_nodes, in_elements, salome_format
      double precision :: x, yc, zc
      
c     Initialise
      n_nodes = 0
      n_triangles = 0
      in_nodes = .false.
      in_elements = .false.
      salome_format = .false.
      
c     Open file
      open(10, file=filename, status='old', iostat=io)
      if (io .ne. 0) then
        print *, 'ERROR: Cannot open file ', filename
        return
      end if
      
c     Scan file for sections
      do
        read(10, '(A)', iostat=io) line
        if (io .ne. 0) exit
        
c       Detect node section (2411)
        if (index(line, '2411') .ne. 0) then
          in_nodes = .true.
          in_elements = .false.
c         Check format: read next line to see if it has 4 integers
          read(10, '(A)', iostat=io) line
          backspace(10)
          read(line, *, iostat=io) dummy, dummy, dummy, dummy
          if (io .eq. 0) then
            salome_format = .true.      ! Salome/ASTER format
          else
            salome_format = .false.     ! Simple format
          endif
          cycle
        endif
        
c       Detect element section (2412)
        if (index(line, '2412') .ne. 0) then
          in_nodes = .false.
          in_elements = .true.
          cycle
        endif
        
c       End of section marker
        if (line .eq. '    -1' .or. line .eq. '-1') then
          in_nodes = .false.
          in_elements = .false.
          cycle
        endif
        
c       Read nodes
        if (in_nodes) then
          if (salome_format) then
c           Salome format: "id 1 1 11" on first line, coordinates on next line
            read(line, *, iostat=io) node_id, dummy, dummy, dummy
            if (io .eq. 0 .and. node_id .gt. 0) then
              read(10, *, iostat=io) x, yc, zc
              if (io .eq. 0) then
                n_nodes = n_nodes + 1
                if (n_nodes .le. max_nodes) then
                  y(n_nodes) = x
                  z(n_nodes) = yc
                endif
              endif
            endif
          else
c           Simple format: "id x y z" all on one line
            read(line, *, iostat=io) node_id, x, yc, zc
            if (io .eq. 0 .and. node_id .gt. 0) then
              n_nodes = n_nodes + 1
              if (n_nodes .le. max_nodes) then
                y(n_nodes) = x
                z(n_nodes) = yc
              endif
            endif
          endif
          cycle
        endif
        
c       Read triangular elements (type 41 = 3-node triangle)
        if (in_elements) then
c         Element header: "id type phys_prop"
          read(line, *, iostat=io) elem_id, elem_type, phys_prop
          if (io .eq. 0 .and. elem_type .eq. 41) then
            n_triangles = n_triangles + 1
            if (n_triangles .le. max_triangles) then
c             Read connectivity: n1 n2 n3
              read(10, *, iostat=io) n1, n2, n3
              if (io .eq. 0) then
                conn(1, n_triangles) = n1
                conn(2, n_triangles) = n2
                conn(3, n_triangles) = n3
              endif
            endif
          endif
        endif
      end do
      
      close(10)
      
      print *, 'Read mesh from: ', filename
      print *, '   Nodes:       ', n_nodes
      print *, '   Triangles:   ', n_triangles
      
      end subroutine read_section_mesh_unv
