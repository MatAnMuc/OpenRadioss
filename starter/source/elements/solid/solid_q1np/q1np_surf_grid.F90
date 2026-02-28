!Copyright>        OpenRadioss
!Copyright>        Copyright (C) 1986-2024 Altair Engineering Inc.
!Copyright>
!Copyright>        This program is free software: you can redistribute it and/or modify
!Copyright>        it under the terms of the GNU Affero General Public License as published by
!Copyright>        the Free Software Foundation, either version 3 of the License, or
!Copyright>        (at your option) any later version.
!Copyright>
!Copyright>        This program is distributed in the hope that it will be useful,
!Copyright>        but WITHOUT ANY WARRANTY; without even the implied warranty of
!Copyright>        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!Copyright>        GNU Affero General Public License for more details.
!Copyright>
!Copyright>        You should have received a copy of the GNU Affero General Public License
!Copyright>        along with this program.  If not, see <https://www.gnu.org/licenses/>.
!Copyright>
!Copyright>
!Copyright>        Commercial Alternative: Altair Radioss Software
!Copyright>
!Copyright>        As an alternative to this open-source version, Altair also offers Altair Radioss
!Copyright>        software under a commercial license.  Contact Altair to discuss further if the
!Copyright>        commercial version may interest you: https://www.altair.com/radioss/.
!Chd|====================================================================
!Chd|  Q1NP_BUILD_SURF_GRID            source/elements/solid/solid_q1np/q1np_surf_grid.F90
!Chd|====================================================================
!=======================================================================
!   Build surface grid topology from quad segment connectivity.
!   Each segment is a quad with 4 nodes; segments share edges.
!   Uses BFS from segment 1 to assign (I,J) grid indices and to fill
!   GRID_NODE (node IDs at grid points) and GRID_TO_SEG (segment at (I,J)).
!   IERR: 0=ok, 3=rotation mismatch, 4=disconnected, 5=non-rectangular, 7=bounds.
!=======================================================================
      module q1np_surf_grid_mod
        use groupdef_mod
        implicit none
      contains
        subroutine q1np_build_surf_grid(surf, nseg, &
     &                                  nx, ny, seg_i, seg_j, &
     &                                  grid_node, grid_to_seg, ierr)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Implicit none
! ----------------------------------------------------------------------------------------------------------------------
          implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
          type(surf_), intent(in) :: surf
          integer, intent(in)  :: nseg
          integer, intent(out) :: nx, ny, ierr
          integer, intent(inout), dimension(nseg) :: seg_i, seg_j
          integer, intent(inout), dimension(nseg*4, nseg*4) :: grid_node
          integer, intent(inout), dimension(nseg, nseg) :: grid_to_seg
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local variables
! ----------------------------------------------------------------------------------------------------------------------
          integer, parameter :: off = 2048
          integer, allocatable :: neighbor(:, :)
          integer :: nodes(4), n1, n2, na, nb, na2, nb2
          integer :: iseg, jseg, k, kk, dir, idir, d
          integer :: ii, jj, i_min, i_max, j_min, j_max
          integer :: head, tail, nassign
          integer :: idx, i1, i2, next(4)
          integer :: gi(4), gj(4), base(4), cand(4), exist(4)
          integer :: grid_node_tmp(1-off:off*2+1, 1-off:off*2+1)
          integer :: qseg(4*off), qi(4*off), qj(4*off)
!     Direction order for BFS: right(2), top(4), left(1), bottom(3)
          integer :: dir_list(4), delta_i(4), delta_j(4)
          integer :: gioff(4, 4), gjoff(4, 4)
          integer :: jseg_found
          data next / 2, 3, 4, 1 /
          data dir_list / 2, 4, 1, 3 /
          data delta_i / -1, 1, 0, 0 /
          data delta_j / 0, 0, -1, 1 /
!     Grid offsets for neighbor quad corners (dir, corner): 1=left, 2=right, 3=bottom, 4=top
          data gioff / -1, 1, 0, 0,  0, 2, 1, 1,  0, 2, 1, 1,  -1, 1, 0, 0 /
          data gjoff / 0, 0, -1, 1,  0, 0, -1, 1,  1, 1, 0, 2,  1, 1, 0, 2 /
!=======================================================================
          ierr = 0
          nx = 0
          ny = 0
          allocate(neighbor(nseg, 4))

          seg_i(1:nseg) = 0
          seg_j(1:nseg) = 0
          neighbor(1:nseg, 1:4) = 0

!     --- Build neighbor list ---
!     For each segment and each edge (DIR): find segment sharing that edge.
!     DIR 1=left(edge 4-1), 2=right(2-3), 3=bottom(1-2), 4=top(3-4).
          do iseg = 1, nseg
            do k = 1, 4
              nodes(k) = surf%nodes(iseg, k)
            end do
            do dir = 1, 4
              if (dir .eq. 1) then
                na = nodes(4)
                nb = nodes(1)
              else if (dir .eq. 2) then
                na = nodes(2)
                nb = nodes(3)
              else if (dir .eq. 3) then
                na = nodes(1)
                nb = nodes(2)
              else
                na = nodes(3)
                nb = nodes(4)
              end if
              n1 = min(na, nb)
              n2 = max(na, nb)
              jseg_found = 0
              do jseg = 1, nseg
                if (jseg .eq. iseg) cycle
                do kk = 1, 4
                  na2 = surf%nodes(jseg, kk)
                  nb2 = surf%nodes(jseg, next(kk))
                  if (min(na2, nb2) .eq. n1 .and. max(na2, nb2) .eq. n2) then
                    neighbor(iseg, dir) = jseg
                    jseg_found = 1
                    exit
                  end if
                end do
                if (jseg_found .ne. 0) exit
              end do
            end do
          end do

!     Temporary grid for node IDs (offset indexing so (1,1) is valid)
          do ii = 1-off, off*2+1
            do jj = 1-off, off*2+1
              grid_node_tmp(ii, jj) = 0
            end do
          end do

!     --- BFS: assign (I,J) to each segment and fill GRID_NODE_TMP ---
!     Start with segment 1 at (1,1); quad nodes 1,2,3,4 at (1,1),(2,1),(2,2),(1,2).
          head = 1
          tail = 1
          qseg(1) = 1
          qi(1) = 1
          qj(1) = 1
          seg_i(1) = 1
          seg_j(1) = 1
          do k = 1, 4
            nodes(k) = surf%nodes(1, k)
          end do
          grid_node_tmp(1, 1) = nodes(1)
          grid_node_tmp(2, 1) = nodes(2)
          grid_node_tmp(2, 2) = nodes(3)
          grid_node_tmp(1, 2) = nodes(4)
          nassign = 1

          do while (head .le. tail)
            iseg = qseg(head)
            ii = qi(head)
            jj = qj(head)
            head = head + 1

            do k = 1, 4
              nodes(k) = surf%nodes(iseg, k)
            end do

!       Process each direction: if unvisited neighbor exists, enqueue and assign node IDs
            do idir = 1, 4
              d = dir_list(idir)
              jseg = neighbor(iseg, d)
              if (jseg .le. 0 .or. seg_i(jseg) .ne. 0) cycle

              seg_i(jseg) = ii + delta_i(d)
              seg_j(jseg) = jj + delta_j(d)
              nassign = nassign + 1
              tail = tail + 1
              qseg(tail) = jseg
              qi(tail) = seg_i(jseg)
              qj(tail) = seg_j(jseg)

              do k = 1, 4
                base(k) = surf%nodes(jseg, k)
              end do
              do k = 1, 4
                gi(k) = ii + gioff(d, k)
                gj(k) = jj + gjoff(d, k)
                exist(k) = grid_node_tmp(gi(k), gj(k))
              end do

!         Find rotation (0..3) of neighbor's nodes so CAND matches already-assigned EXIST
              idx = -1
              do i1 = 0, 3
                do k = 1, 4
                  i2 = 1 + mod(k - 1 + i1, 4)
                  cand(k) = base(i2)
                  if (exist(k) .gt. 0 .and. exist(k) .ne. cand(k)) exit
                end do
                if (k .gt. 4) then
                  idx = i1
                  exit
                end if
              end do

              if (idx .lt. 0) then
                ierr = 3
                deallocate(neighbor)
                return
              end if

              do k = 1, 4
                grid_node_tmp(gi(k), gj(k)) = cand(k)
              end do
            end do
          end do

!     --- Sanity: all segments must have been assigned (connected surface) ---
          if (nassign .ne. nseg) then
            ierr = 4
            deallocate(neighbor)
            return
          end if

!     --- Normalize (I,J) to 1..NX, 1..NY and compute NX, NY ---
          i_min = seg_i(1)
          i_max = seg_i(1)
          j_min = seg_j(1)
          j_max = seg_j(1)
          do iseg = 2, nseg
            i_min = min(i_min, seg_i(iseg))
            i_max = max(i_max, seg_i(iseg))
            j_min = min(j_min, seg_j(iseg))
            j_max = max(j_max, seg_j(iseg))
          end do
          nx = i_max - i_min + 1
          ny = j_max - j_min + 1

          if (nx*ny .ne. nseg) then
            ierr = 5
            deallocate(neighbor)
            return
          end if

          do iseg = 1, nseg
            seg_i(iseg) = seg_i(iseg) - i_min + 1
            seg_j(iseg) = seg_j(iseg) - j_min + 1
          end do

!     --- Fill GRID_TO_SEG: segment index at each (I,J) ---
          grid_to_seg(1:nx, 1:ny) = 0
          do iseg = 1, nseg
            ii = seg_i(iseg)
            jj = seg_j(iseg)
            grid_to_seg(ii, jj) = iseg
          end do

!     --- Copy GRID_NODE_TMP into output GRID_NODE (1-based, size (NX+1)*(NY+1)) ---
          if (nx+1 .gt. nseg*4 .or. ny+1 .gt. nseg*4) then
            ierr = 7
            deallocate(neighbor)
            return
          end if

          do jj = j_min, j_max + 1
            do ii = i_min, i_max + 1
              grid_node(ii - i_min + 1, jj - j_min + 1) = grid_node_tmp(ii, jj)
            end do
          end do

          deallocate(neighbor)

          return
        end subroutine q1np_build_surf_grid
      end module q1np_surf_grid_mod
