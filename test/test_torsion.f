c test/test_torsion.f
c Test program for torsional constant J computation
c
c Authors: Bruno Zilli & DeepSeek
c Licence: MIT

      program test_torsion
      use torsion_j
      implicit none
      
      integer :: nnode, ntri
      integer, allocatable :: conn(:,:)
      double precision, allocatable :: y(:), z(:)
      double precision :: J_result
      
c     Test 1: Rectangle 10x20 mm (manually built)
      print *, ''
      print *, '========================================'
      print *, 'TEST 1: Rectangle 10x20 mm'
      print *, '========================================'
      
      nnode = 4
      ntri = 2
      allocate(y(nnode), z(nnode), conn(3, ntri))
      
      y = (/0.0d0, 10.0d0, 10.0d0, 0.0d0/)
      z = (/0.0d0, 0.0d0, 20.0d0, 20.0d0/)
      
      conn(1,1) = 1
      conn(2,1) = 2
      conn(3,1) = 3
      conn(1,2) = 1
      conn(2,2) = 3
      conn(3,2) = 4
      
      call compute_torsion_J(nnode, ntri, conn, y, z, J_result,
     &                       .true.)
      
      deallocate(y, z, conn)
      
c     Test 2: Circle diameter 10 mm
      print *, ''
      print *, '========================================'
      print *, 'TEST 2: Circle diameter 10 mm'
      print *, '========================================'
      
      call read_mesh_wrapper('meshes/circle_dia10.unv', nnode,
     &                       ntri, conn, y, z)
      call compute_torsion_J(nnode, ntri, conn, y, z, J_result,
     &                       .true.)
      
      deallocate(y, z, conn)
      
c     Test 3: HEB200
      print *, ''
      print *, '========================================'
      print *, 'TEST 3: HEB200'
      print *, '========================================'
      
      call read_mesh_wrapper('meshes/HEB200_mm.unv', nnode,
     &                       ntri, conn, y, z)
      call compute_torsion_J(nnode, ntri, conn, y, z, J_result,
     &                       .true.)
      
      deallocate(y, z, conn)
      
      end program test_torsion
