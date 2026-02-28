!Chd|====================================================================
!Chd|  Q1NP_MASS3                    source/elements/solid/solid_q1np/q1np_mass3.F90
!Chd|====================================================================
!=======================================================================
!   Calculate lumped mass for Q1Np enriched elements
!
!   This routine distributes mass for Q1Np elements:
!   - Total mass equals the underlying HEX8 element mass
!   - 50% of mass distributed to 4 bulk nodes (12.5% each)
!   - 50% of mass distributed equally to all control point nodes
!   - HEX8 top nodes (nodes 5-8) explicitly set to zero mass
!=======================================================================
      module q1np_mass3_mod
        use message_mod
        use q1np_restart_mod
        use q1np_volume_mod
        use precision_mod, only : WP
        use constant_mod,  only : ZERO, HALF, FOURTH, SIX
        use element_mod   , only : NIXS
        implicit none
      contains
!
        subroutine q1np_mass3( &
     &      rho, ms, partsav, x, v, &
     &      ipart, mss, volu, &
     &      in, &
     &      vr, rhof, frac, fill, &
     &      kq1np_tab, iq1np_tab, iq1np_bulk_tab, &
     &      ixs, numelq1np_in, npropm, nummat, pm, &
     &      numels, numnod, q1np_ktab_g)
!-----------------------------------------------------------------------
          implicit none
!-----------------------------------------------------------------------
!     Dummy arguments
!-----------------------------------------------------------------------
          integer, intent(in) :: numelq1np_in
          integer, intent(in) :: npropm, nummat
          integer, intent(in) :: kq1np_tab(15, *)
          integer, intent(in) :: iq1np_tab(*)
          integer, intent(in) :: iq1np_bulk_tab(*)
          integer, intent(in) :: ixs(NIXS,*)
          integer, intent(in) :: ipart(*)
          real(kind=WP), intent(inout) :: ms(*), partsav(20,*)
          real(kind=WP), intent(inout) :: x(3,*), v(3,*)
          real(kind=WP), intent(inout) :: mss(8,*)
          real(kind=WP), intent(inout) :: in(*), vr(3,*)
          real(kind=WP), intent(in)    :: rho(*), volu(*), rhof(*), frac(*), fill(*)
          real(kind=WP), intent(in)    :: pm(npropm, nummat)
          integer, intent(in)          :: numels, numnod
          real(kind=WP), intent(in)    :: q1np_ktab_g(*)
!-----------------------------------------------------------------------
!     Local variables
!-----------------------------------------------------------------------
          integer :: iel, iel_hex8, ip
          integer :: i, j, k, n_ctrl, offset_ctrl, offset_bulk
          integer :: node_bulk(4), node_cp, node_hex8_top(4)
          integer :: mid
          real(kind=WP) :: mass_total, mass_bulk, mass_cp, mass_per_cp
          real(kind=WP) :: mass_per_bulk
          real(kind=WP) :: xx, yy, zz, xy, yz, zx
          real(kind=WP) :: volu_q1np
          integer, parameter :: IDEBUG_Q1NP = 0
!=======================================================================
!   Early return if no Q1Np elements
!=======================================================================
          if (numelq1np_in == 0) return
!=======================================================================
!   Loop over all Q1Np elements
!=======================================================================
          do iel = 1, numelq1np_in
!           Get original HEX8 element index
            iel_hex8 = kq1np_tab(10, iel)
            if (iel_hex8 <= 0 .or. iel_hex8 > numels) cycle
!           Get material ID and part ID
            mid = kq1np_tab(1, iel)
            ip  = kq1np_tab(11, iel)
            if (ip <= 0) cycle
!           Compute Q1NP volume
            call q1np_compute_volume_element(iel, &
     &        kq1np_tab, iq1np_tab, iq1np_bulk_tab, &
     &        q1np_ktab_g, x, q1np_nx_g, q1np_ny_g, &
     &        volu_q1np)
!           Optional debug: compare Q1NP volume with original HEX8 volume
            if (IDEBUG_Q1NP >= 2) then
              if (iel_hex8 > 0 .and. iel_hex8 <= numels) then
                write(*,'(A,I6,A,I6,2(A,1P,E12.5),A,1P,E12.5)') &
     &            'Q1NP MASS DBG: IEL_Q1NP=', iel, ' HEX8=', iel_hex8, &
     &            ' VOL_HEX8=', volu(iel_hex8), &
     &            ' VOL_Q1NP=', volu_q1np, &
     &            ' RATIO_Q1NP_HEX8=', volu_q1np / max(volu(iel_hex8), ZERO)
              else
                write(*,'(A,I6,A,I6,A,1P,E12.5)') &
     &            'Q1NP MASS DBG: IEL_Q1NP=', iel, ' HEX8=', iel_hex8, &
     &            ' VOL_Q1NP=', volu_q1np
              end if
            end if
!           Get density
            if (mid > 0) then
              mass_total = pm(1,mid) * volu_q1np
              if (fill(iel_hex8) > ZERO) then
                mass_total = fill(iel_hex8) * mass_total
              end if
            else
              if (fill(iel_hex8) > ZERO) then
                mass_total = fill(iel_hex8) * rho(iel_hex8) * volu_q1np
              else
                mass_total = rho(iel_hex8) * volu_q1np
              end if
            end if
!           Distribute 50% to bulk nodes, 50% to control nodes
            mass_bulk = HALF * mass_total
            mass_cp   = HALF * mass_total
!           Bulk nodes
            offset_bulk = kq1np_tab(14, iel)
            do i = 1, 4
              node_bulk(i) = iq1np_bulk_tab(offset_bulk + i - 1)
            end do
            mass_per_bulk = mass_bulk * FOURTH
            do i = 1, 4
              if (node_bulk(i) > 0 .and. node_bulk(i) <= numnod) then
                ms(node_bulk(i)) = ms(node_bulk(i)) + mass_per_bulk
              end if
            end do
!           Control points
            n_ctrl      = kq1np_tab(3, iel)
            offset_ctrl = kq1np_tab(4, iel)
            if (n_ctrl > 0) then
              mass_per_cp = mass_cp / real(n_ctrl, kind=WP)
              do i = 1, n_ctrl
                node_cp = iq1np_tab(offset_ctrl + i - 1)
                if (node_cp > 0 .and. node_cp <= numnod) then
                  ms(node_cp) = ms(node_cp) + mass_per_cp
                end if
              end do
            end if
!           PARTSAV mass and inertia
            partsav(1,ip) = partsav(1,ip) + mass_total
            do i = 1, 4
              if (node_bulk(i) > 0 .and. node_bulk(i) <= numnod) then
                partsav(2,ip) = partsav(2,ip) + mass_per_bulk * x(1, node_bulk(i))
                partsav(3,ip) = partsav(3,ip) + mass_per_bulk * x(2, node_bulk(i))
                partsav(4,ip) = partsav(4,ip) + mass_per_bulk * x(3, node_bulk(i))
              end if
            end do
            do i = 1, n_ctrl
              node_cp = iq1np_tab(offset_ctrl + i - 1)
              if (node_cp > 0 .and. node_cp <= numnod) then
                partsav(2,ip) = partsav(2,ip) + mass_per_cp * x(1, node_cp)
                partsav(3,ip) = partsav(3,ip) + mass_per_cp * x(2, node_cp)
                partsav(4,ip) = partsav(4,ip) + mass_per_cp * x(3, node_cp)
              end if
            end do
!           Inertia (simplified)
            xx = ZERO; yy = ZERO; zz = ZERO
            xy = ZERO; yz = ZERO; zx = ZERO
            do i = 1, 4
              if (node_bulk(i) > 0 .and. node_bulk(i) <= numnod) then
                xx = xx + mass_per_bulk * x(1, node_bulk(i)) * x(1, node_bulk(i))
                yy = yy + mass_per_bulk * x(2, node_bulk(i)) * x(2, node_bulk(i))
                zz = zz + mass_per_bulk * x(3, node_bulk(i)) * x(3, node_bulk(i))
                xy = xy + mass_per_bulk * x(1, node_bulk(i)) * x(2, node_bulk(i))
                yz = yz + mass_per_bulk * x(2, node_bulk(i)) * x(3, node_bulk(i))
                zx = zx + mass_per_bulk * x(3, node_bulk(i)) * x(1, node_bulk(i))
              end if
            end do
            do i = 1, n_ctrl
              node_cp = iq1np_tab(offset_ctrl + i - 1)
              if (node_cp > 0 .and. node_cp <= numnod) then
                xx = xx + mass_per_cp * x(1, node_cp) * x(1, node_cp)
                yy = yy + mass_per_cp * x(2, node_cp) * x(2, node_cp)
                zz = zz + mass_per_cp * x(3, node_cp) * x(3, node_cp)
                xy = xy + mass_per_cp * x(1, node_cp) * x(2, node_cp)
                yz = yz + mass_per_cp * x(2, node_cp) * x(3, node_cp)
                zx = zx + mass_per_cp * x(3, node_cp) * x(1, node_cp)
              end if
            end do
            partsav(5,ip)  = partsav(5,ip)  + (yy + zz)
            partsav(6,ip)  = partsav(6,ip)  + (zz + xx)
            partsav(7,ip)  = partsav(7,ip)  + (xx + yy)
            partsav(8,ip)  = partsav(8,ip)  - xy
            partsav(9,ip)  = partsav(9,ip)  - yz
            partsav(10,ip) = partsav(10,ip) - zx
!           Momentum
            do i = 1, 4
              if (node_bulk(i) > 0 .and. node_bulk(i) <= numnod) then
                partsav(11,ip) = partsav(11,ip) + mass_per_bulk * v(1, node_bulk(i))
                partsav(12,ip) = partsav(12,ip) + mass_per_bulk * v(2, node_bulk(i))
                partsav(13,ip) = partsav(13,ip) + mass_per_bulk * v(3, node_bulk(i))
              end if
            end do
            do i = 1, n_ctrl
              node_cp = iq1np_tab(offset_ctrl + i - 1)
              if (node_cp > 0 .and. node_cp <= numnod) then
                partsav(11,ip) = partsav(11,ip) + mass_per_cp * v(1, node_cp)
                partsav(12,ip) = partsav(12,ip) + mass_per_cp * v(2, node_cp)
                partsav(13,ip) = partsav(13,ip) + mass_per_cp * v(3, node_cp)
              end if
            end do
!           Kinetic energy
            do i = 1, 4
              if (node_bulk(i) > 0 .and. node_bulk(i) <= numnod) then
                partsav(14,ip) = partsav(14,ip) + HALF * mass_per_bulk * &
     &            ( v(1,node_bulk(i))*v(1,node_bulk(i)) + &
     &              v(2,node_bulk(i))*v(2,node_bulk(i)) + &
     &              v(3,node_bulk(i))*v(3,node_bulk(i)) )
              end if
            end do
            do i = 1, n_ctrl
              node_cp = iq1np_tab(offset_ctrl + i - 1)
              if (node_cp > 0 .and. node_cp <= numnod) then
                partsav(14,ip) = partsav(14,ip) + HALF * mass_per_cp * &
     &            ( v(1,node_cp)*v(1,node_cp) + &
     &              v(2,node_cp)*v(2,node_cp) + &
     &              v(3,node_cp)*v(3,node_cp) )
              end if
            end do
          end do
          return
        end subroutine q1np_mass3
!
      end module q1np_mass3_mod

