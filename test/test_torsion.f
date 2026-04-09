c test/test_torsion.f
c Test program for torsional constant J computation
c
c Authors: Bruno Zilli & DeepSeek
c License: MIT
c Copyright (c) 2025 Bruno Zilli & DeepSeek
c     
c Permission is hereby granted, free of charge, to any person obtaining
c a copy of this software and associated documentation files (the
c "Software"), to deal in the Software without restriction, including
c without limitation the rights to use, copy, modify, merge, publish,
c distribute, sublicense, and/or sell copies of the Software, and to
c permit persons to whom the Software is furnished to do so, subject to
c the following conditions:
c     
c The above copyright notice and this permission notice shall be
c included in all copies or substantial portions of the Software.
c     
c THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
c EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
c MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
c IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
c CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
c TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
c SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

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
