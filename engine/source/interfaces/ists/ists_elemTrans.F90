!||====================================================================
!||    sts_handle_element_transition  ../engine/source/interfaces/ists/ists_elemTrans.F90
!||--- called by ------------------------------------------------------
!||    sts_tangentvel_global   ../engine/source/interfaces/ists/ists_tangentvel.F90
!||====================================================================
!-----------------------------------------------
!   Handle element transition for global xi tracking
!-----------------------------------------------
      subroutine sts_handle_element_transition(xi1_local_prev, xi1_local_new, &
     &     xi2_local_prev, xi2_local_new, &
     &     gp_xi1_period, gp_xi2_period, tol)

      implicit none

!-- Input arguments
      real*8  xi1_local_prev, xi1_local_new
      real*8  xi2_local_prev, xi2_local_new
      real*8  tol

!-- In/out arguments
      integer gp_xi1_period
      integer gp_xi2_period

!-----------------------------------------------
!  Detect crossings of +/-1 in xi1 and xi2 and
!  update period counters accordingly.
!-----------------------------------------------

!-- Check transition in xi1 direction
      if (dabs(xi1_local_prev - 1.0d0) < tol) then
!       Was at +1 boundary, likely crossed to next element
        if (xi1_local_new < 0.0d0) then
!         Crossed from +1 to -1 (forward)
          gp_xi1_period = gp_xi1_period + 1
        endif
      else if (dabs(xi1_local_prev + 1.0d0) < tol) then
!       Was at -1 boundary, likely crossed to previous element
        if (xi1_local_new > 0.0d0) then
!         Crossed from -1 to +1 (backward)
          gp_xi1_period = gp_xi1_period - 1
        endif
      endif

!-- Check transition in xi2 direction (same logic)
      if (dabs(xi2_local_prev - 1.0d0) < tol) then
        if (xi2_local_new < 0.0d0) then
          gp_xi2_period = gp_xi2_period + 1
        endif
      else if (dabs(xi2_local_prev + 1.0d0) < tol) then
        if (xi2_local_new > 0.0d0) then
          gp_xi2_period = gp_xi2_period - 1
        endif
      endif

      return
      end