!=======================================================================
SUBROUTINE BldVib ( nseg, ind, mdsqss, psi0, sprlag, eflap, elag, evr, evl, gm )

   !     * Build blade model and compute its rotating natural vibration
   !       characteristics.
   !     * Obtain blade flap hinge location (eflap) and lag hinge location
   !       (elag.)

   ! *** OTHER SYMBOLS:

   !     iartic = 0; hingeless blade
   !            = 1; flap hinge articulated blade
   !            = 2; flap then lag hinge articulated blade
   !            = 3; flap then lag hinge artic. blade with del3 effect
   !            = 4; flap and lag hinge coincident articulated blade
   !            = 3,4; pitch link stiffness for trim analysis
   !     nselt  = no. of spatial elements
   !     nconf  = 0; ......
   !            = 1; ......
   !            = 2; bmr
   !     nndof  = no. of nodal dof
   !     nintn  = no. of internal nodes
   !     nedof  = no. of element dof
   !     ngd    = no. of bounded(?) dof
   !     nefs(i) = nr. of elements, segment i
   !     nefs1  = nefs(1), nr. of elements, segment 1
   !     nefs2  = nefs(2), nr. of elements, segment 2
   !     nefs12 = sum of nr. of elements, segments 1 and 2.
   !     nesh   = total dofs per element, counting only one node and
   !              the internal dofs.
   !     iend   = sum of nr. of elements, segments 1 and 2. (maybe should
   !              equivalence with nefs12?)
   !     xbi    = location of inboard end of an element.
   !     xbi(i) = location of inboard end of i'th element.
   !     xbims  = square of outboard station of an element.
   !     cfei   = centrifugal force at inboard end of an element.
   !     cfe(i) = centrifugal force at inboard end of i'th element.
   !     el(i)  = length of element i.  see MODULE Gauss in modules.f90
   !     rmas(i) = see MODULE Struc in modules.f90
   !     indrns = analysis switch for subroutine asbgbm.
   !     idisbd = 1 for dissimilar blades. See um1.f.
   !     jb     = a counter used to count through one blade, or all blades,
   !              if idisbd = 1. (maybe take this out of common block.)
   !     gk     = global stiffness matrix.
   !     wa     = global displacement vector - in bldvibs is always zero.
   !              Need it to pass to struct.

USE Basic
USE CNPar
USE Conf
USE DatRef
USE DisBld
USE Gauss
USE Gravity
USE NWTC_Library
USE Omg
USE Param
USE pbrng
USE plgeom
USE Pitch
USE plgeom
USE PLink
USE SftLnk
USE Struc
USE Swept
USE TipDat
USE TowWires


   ! Argument declarations.

REAL(ReKi), INTENT(OUT)       :: eflap                            ! Flap hinge offset.
REAL(ReKi), INTENT(OUT)       :: elag                             ! Lag hinge offset.
REAL(ReKi), INTENT(OUT)       :: evl      (ngd)                   ! Vector of eigenvalues.
REAL(ReKi), INTENT(OUT)       :: evr      (ngd,ngd)               ! Matrix of eigenvectors.
REAL(ReKi), INTENT(OUT)       :: gm       (ngd,ngd)               ! Global mass matrix.
REAL(ReKi), INTENT(IN)        :: psi0                             ! Azimuth at which to evaluate blade normal modes.
REAL(ReKi), INTENT(IN)        :: sprlag                           ! Lag damper spring constant.

INTEGER, INTENT(IN)           :: ind      (nedof)                 ! Elemental connectivity vector.
INTEGER, INTENT(IN)           :: mdsqss   (nmodet)                ! Sequence of normal modes for response.
INTEGER, INTENT(IN)           :: nseg                             ! Number of blade segments (bmr).


   ! Local declarations.

REAL(ReKi)                    :: AZBAR
REAL(ReKi)                    :: CFEI
REAL(ReKi)                    :: CFTT
REAL(ReKi)                    :: DEN
REAL(ReKi), ALLOCATABLE       :: dfxc     (:,:,:)
REAL(ReKi), ALLOCATABLE       :: dfxcd    (:,:,:)
REAL(ReKi)                    :: dfxs     (nedof,nedof)
REAL(ReKi)                    :: ec       (nedof,nedof)
REAL(ReKi)                    :: EFFLP
REAL(ReKi)                    :: ek       (nedof,nedof)
REAL(ReKi), ALLOCATABLE       :: elv      (:)
REAL(ReKi)                    :: em       (nedof,nedof)
REAL(ReKi)                    :: eq       (nedof)
REAL(ReKi)                    :: eu       (nedof,3)
REAL(ReKi)                    :: FREQ
REAL(ReKi), ALLOCATABLE       :: gk       (:,:)
REAL(ReKi), ALLOCATABLE       :: gu       (:)
REAL(ReKi), ALLOCATABLE       :: gud      (:)
REAL(ReKi)                    :: HUB_RAD
REAL(ReKi)                    :: JZ
REAL(ReKi)                    :: SPAN_LOC
REAL(ReKi), ALLOCATABLE       :: tcsd     (:,:,:)
REAL(ReKi)                    :: TH0
REAL(ReKi), ALLOCATABLE       :: tht      (:)
REAL(ReKi), ALLOCATABLE       :: tksd     (:,:,:)
REAL(ReKi)                    :: tlam     (3,3)
REAL(ReKi)                    :: tlam2    (nedof,nedof)
REAL(ReKi)                    :: TWST
REAL(ReKi)                    :: VDISP
REAL(ReKi)                    :: VP
REAL(ReKi), ALLOCATABLE       :: vpt      (:)
REAL(ReKi), ALLOCATABLE       :: wa       (:)
REAL(ReKi)                    :: WDISP
REAL(ReKi)                    :: WP
REAL(ReKi), ALLOCATABLE       :: wpt      (:)
REAL(ReKi)                    :: XBI
REAL(ReKi)                    :: XBIMS
REAL(ReKi)                    :: XBJ
REAL(ReKi)                    :: XBLHS
REAL(ReKi)                    :: XBRH
REAL(ReKi)                    :: XBRHS
REAL(ReKi)                    :: XSPAN
REAL(ReKi)                    :: ZBAR

INTEGER                       :: glob_node
INTEGER                       :: I
INTEGER                       :: IE
INTEGER                       :: IEND
INTEGER                       :: IJK
INTEGER                       :: ILG1
INTEGER                       :: ILG2
INTEGER                       :: INODE
INTEGER                       :: J
INTEGER                       :: J1
INTEGER                       :: JF
INTEGER                       :: JJ
INTEGER                       :: JP
INTEGER                       :: JT
INTEGER                       :: JV
INTEGER                       :: JW
INTEGER                       :: KF
INTEGER                       :: L
INTEGER                       :: LBB
INTEGER                       :: LL
INTEGER                       :: LMN
INTEGER                       :: nd_gv2
INTEGER                       :: nd_gw2
INTEGER                       :: NEFS1
INTEGER                       :: NEFS12
INTEGER                       :: NEFS2
INTEGER                       :: NEFS3
INTEGER                       :: NESH
INTEGER                       :: NNBLAD
INTEGER                       :: NSH
INTEGER                       :: Sttus

CHARACTER(200)                :: Fmt
CHARACTER(200)                :: fmt_data1
CHARACTER(200)                :: fmt_data2
CHARACTER(200)                :: fmt_tabs
CHARACTER(200)                :: FmtDashes = "('=================================================================================')"



   ! Allocate some local arrays.

ALLOCATE ( dfxc(nselt,nedof,nedof), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, dfxc, in subroutine bldvib.' )

ALLOCATE ( dfxcd(nselt,nedof,nedof), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, dfxcd, in subroutine bldvib.' )

ALLOCATE ( elv(nselt), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, elv, in subroutine bldvib.' )

ALLOCATE ( gk(ngd,ngd), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, gk, in subroutine bldvib.' )

ALLOCATE ( gu(ngd), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, gu, in subroutine bldvib.' )

ALLOCATE ( gud(ngd), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, gud, in subroutine bldvib.' )

ALLOCATE ( tcsd(nselt,maxhub,nedof), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, tcsd, in subroutine bldvib.' )

ALLOCATE ( tht(nselt), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, tht, in subroutine bldvib.' )

ALLOCATE ( tksd(nselt,maxhub,nedof), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, tksd, in subroutine bldvib.' )

ALLOCATE ( vpt(nselt), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, vpt, in subroutine bldvib.' )

ALLOCATE ( wa(ngd), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, wa, in subroutine bldvib.' )

ALLOCATE ( wpt(nselt), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, wpt, in subroutine bldvib.' )


   !-----------------------------------------------------------------------
   !sv
   ! *** Copy el into elv. elv is what gets sent down to 'struct'.

do lmn = 1,nselt
   elv(lmn) = el(lmn)
end do ! lmn

   !     'j1' seems to remain at 1 for entire subprogram.
j1 = 1

   ! *** Calculate flap and lag offsets:

efflp = 0.0
SELECT CASE ( nconf )
CASE ( 0 )

   DO ijk=1,nselt
      efflp = efflp + el(ijk)
   END DO ! ijk

CASE ( 2 )

   nefs12 = nefs(1) + nefs(2)

   DO ijk=1,nefs12
      efflp = efflp + el(ijk)
   END DO ! ijk

END SELECT

   !check later
eflap = 1.0 - efflp
if (iartic==4) then
   elag = eflap
  else
   if (nconf==0) then
      elag = eflap + el(nselt)
   else if (nconf==2) then
      elag = eflap + el(nefs12)
   end if
end if
   !--


   ! *** Determine global connectivity vector by relating element dof
   !     to global dof. Also: determine centrifugal force term.

   !     Initialize:
nesh = nndof + nintn
nefs1 = nefs(1)
nefs2 = nefs(2)
iend = nefs1 + nefs2
   !     Values at tip:
xbi   = 1.
cfei  = 0.
   !tower
   !       cfei  = 25580.0 * grav/(rm*radius)
   !---

   !     loop on elements, from tip to inboard:

DO i=1,iend

   xbims  = xbi * xbi
   xbi    = xbi - el(i)
   xb(i)  = xbi
   cfe(i) = cfei
   cfei   = cfei + 0.5 * rmas(i) * (xbims - xbi * xbi) ! improve later
      !tower
      !         cfei   = cfei - grav * rmas(i) * el(i)
      !---

      !        For the i'th element, set the connectivity vector indeg:
   nsh = ( i - 1 )*nesh

   DO j=1,nedof
      indeg(j,i) = ind(j) + nsh
   END DO ! j

END DO ! i

   !--------------------------------------------------------------------

   !     compute element boundary locations and centrifugal force
   !     distribution for the bmr blade
   !-------------------------------------------------------------------

   !free   For a free-free blade, skip the following and go to the
   !       end of IF block

if (ifree /= 1) then

   if (nseg == 3) then

      nefs3 = nefs(3)

      xbj = xb(nefs1)

      DO i=1,nefs3
         ll = iend + i
         xbj = xbj  - el(ll)
         xb(ll) = xbj
      END DO ! i

      cftt = 0.

      DO i=1,nefs3
         ll = nselt - (i-1)
         xblhs = xb(ll)**2
         xbrh  = xb(ll) + el(ll)
         xbrhs = xbrh*xbrh
         cfe(ll) = cftt + .5*rmas(ll)*(xblhs - xbrhs)
         cftt = cfe(ll)
      END DO ! i

      DO i=1,nefs2
         ll      = nefs1 + i
         cfe(ll) = cfe(ll) - cftt
      END DO ! i

      nsh = nsh + nesh - nndof
      DO i=1,nefs3

         l = iend + i

         DO j=1,nedof
            indeg(j,l) = ind(j) + nsh
         END DO ! j

         nsh = nsh + nesh

      END DO ! i

   !     impose compatibility conditions at the cuff

      l = iend + 1
      lbb = nefs1

      indeg( 4,l) = indeg(1,lbb)
      indeg( 7,l) = indeg(5,lbb)
      indeg( 8,l) = indeg(6,lbb)
      indeg(11,l) = indeg(9,lbb)
      indeg(12,l) = indeg(10,lbb)
      indeg(15,l) = indeg(13,lbb)

   end if

   !comp fix above for 19dof element bmr
   !--------------------------------------------------------------

   !     cantilevered  boundary   constraints

   iend = nefs(1) + nefs(2)
   !end-----------------------------------------------------------

   ! *** Set all dof numbers at root of blade to zero. They may be reset
   !     to non-zero values later, depending on which if any dofs, if any,
   !     remain unconstrained.

   indeg( 1, iend) = 0
   indeg( 5, iend) = 0
   indeg( 6, iend) = 0
   indeg( 9, iend) = 0
   indeg(10, iend) = 0
   indeg(13, iend) = 0

   !comp
   !     restraint shear displacements at cantilever constraint
   !     set indeg(16, iend) = 0 and indeg(18, iend) = 0 without
   !     explicit assignment statement.  this is required for complilation

   if ( nedof == 19 )  then
      DO i=1,3
         IF ( i == 1 .OR. i == 3 )  indeg(nedof-i,iend) = 0
      END DO ! i
   end if
   !comp

   !gg--------------------------------------------------------

   !comp
   !     fix this from here on down (all bmr calculations) for shear

   !     boundary conditions for the bmr/lag-pin combination

   if( nseg >= 3 ) then
      ll = nefs(1) + nefs(2) + nefs(3)

   !
      if( npin == 1 ) then
         indeg( 1,ll) = ngd - 2 - lsft
         indeg( 5,ll) = 0
         indeg( 6,ll) = ngd - 1 - lsft
         indeg( 9,ll) = 0
         indeg(10,ll) = ngd - lsft
         indeg(13,ll) = ngd * lsft

      else if( npin == 0 ) then
         indeg( 1,ll) = ngd - 4 - lsft
         indeg( 5,ll) = ngd - 3 - lsft
         indeg( 6,ll) = ngd - 2 - lsft
         indeg( 9,ll) = ngd - 1 - lsft
         indeg(10,ll) = ngd - lsft
         indeg(13,ll) = ngd * lsft
      end if

   end if

      !end------------------------------------------------------

      !comp no changes required here. release bending slope constraints

      ! *** Set articulated blade boundary conditions

   if (iartic == 1) then
      indeg(10,nselt) = indeg(2,nselt) + 1

   else if (iartic == 2 .or. iartic == 3) then
      indeg(10,nselt) = indeg(2,nselt) + 1
      indeg(8,nselt)  = indeg(2,nselt) + 2

   else if (iartic == 4) then
      indeg(10,nselt) = indeg(2,nselt) + 1
      indeg(6,nselt)  = indeg(2,nselt) + 2
   end if
      !free
end if

   !nrel
write (UnOu,'(A)')  'Results generated by '//TRIM( ProgName )//TRIM( ProgVer )//' on '//CurDate()//' at '//CurTime()//'.'

write (UnOu,'(A)')  title
write (UnOu,FmtDashes)


   !     End output.


   ! *** Determine blade natural modes and natural frequencies for use in
   !     normal mode equations.


   !jmw
   ! *** If idisbd==1, start loop on all blades; else analyze just one
   !     blade.

nnblad=1

if ( idisbd == 1 )  nnblad = nblade

DO jb=1,nnblad

   if ( idisbd == 1 )  THEN
      write(*   ,*)  ' free vibration in undeformed frame for blade', jb
      write(UnOu,*)  ' free vibration in undeformed frame for blade', jb
   end if

      !Initialize global mass and stiffness matrices. Displacement
      !vector 'wa' is also initialized - need it only to pass to
      !'struct'.
      !jmw
   wa(:  ) = 0.0
   gk(:,:) = 0.0
   gm(:,:) = 0.0

      !        Initialize certain sweep parameters:  tht,wpt,vpt need comment.
      !sweep
   tht(:) = 0.0
   wpt(:) = 0.0
   vpt(:) = 0.0

   DO i=1,nselt

      jt = indeg(13,i)
      jw = indeg(10,i)
      jv = indeg( 6,i)

      IF ( jt /= 0 )  tht(i) = wa(jt)
      IF ( jw /= 0 )  wpt(i) = wa(jw)
      IF ( jv /= 0 )  vpt(i) = wa(jv)

   END DO ! i

      !gsbend

      ! ***    Determine structural mass and stiffness matrices, then assemble
      !        into global eigen matrices for this blade.

      !nrel
      !         open(unit=12,file = 'Element.matrices')
      !--
      ! ***    Start loop on blade elements:
   DO i=1,nselt
      !sweep2
      !           Obtain displacement vector for this element. Except it's by
      !           definition zero. ielrot requires comment.
         if ( ielrot(i) == 1 )  THEN
      !              Element is swept.
      !agb
         th0 = 0.0

         call evfrsw ( wa, eu(1,j1), indeg(1,i), tlam2 )
      else
         call evfrgv ( wa, eu(1,j1), indeg(1,i) )
      end if
         !sweep
         !comp
      if ( nedof == 15 )  then     ! 15 DOF Bernoulli-Euler element.

         call struct ( eu(1,j1), ek, ec, em, eq, dfxs, psi0, xb(i), el(i), i, cfe(i), indrns, indnl, tlam, gu, gud, dfxc, &
                       dfxcd, elv )
      end if

         !----------------------------------------------------------------

         !     modifications for bmr blade:
         !     impose boundary conditions at the torque tube end
         !     modify stifness matrix if pitch link is soft

         !gbmr

      IF (nconf == 2 )  THEN

         if( i == nselt) then

               ! xxa = (ttx0 - ttp) / tta
               ! zxa = xxa * ttz

            xmp = ttx0 - ttp
            ttb = sqrt( (xmp-xpl)**2 + (tta-ypl)**2 + zpl**2 )

            zbar = zpl


               ! Case 1: ------------------------------------------

            if(npin == 1 .and. lsft /= 1) then

               den = ttz*(tta - ypl) + tta*zbar

               b1wpw = ttx0
               b1vpph = (xmp*ypl - tta*xpl) / den
               b1wpph = - xmp*zbar / den

               b1vpv = ttx0 + ttz*b1vpph
               b1wpv = ttz*b1wpph

            end if


               ! Case 2: ------------------------------------------

            if(npin == 0 .and. lsft /= 1) then

               azbar = tta*zbar

               b2vph = (ypl - tta) / azbar
               b2wph = -1. / tta
               b2vpph = (tta*(ttx0-xpl) - ttp*ypl) / azbar
               b2wpph = ttp / tta

            end if


               ! Case 3: ------------------------------------------

            if(npin == 1 .and. lsft == 1) then
                  ! no change in kinematic bcs

               svp  = (tta*xpl - xmp*ypl) / ttb
               swp  = zbar*xmp / ttb
               sphi = (tta*zbar + ttz*(tta-ypl)) / ttb

            end if


               ! Case 4: ------------------------------------------

            if(npin == 0 .and. lsft == 1) then

               tv   = (tta - ypl) / ttb
               tvp  = (ttp*ypl - tta*(ttx0-xpl)) / ttb
               tw   = zbar / ttb
               twp  = -tw * ttp
               tphi = tw * tta

            end if


               !--------------   ----------------------------------

            call bmrmod (em, ek, gk, i)

         end if
      end if

         !end-----------------------------------------------------------------

         ! ***       Assemble this element's matrices into global matrices:

      call asbgmk(gm,gk,em,ek,indeg(1,i))

   END DO ! i     End loop on blade elements

      !print

      !nrel        open(unit=13,file = 'Global_Mass_Matrix')
      !        open(unit=14,file = 'Global_Stiffness_Matrix')
      !        write(13,611) ngd, ngd
      ! 611     format(//1x, 'Global mass matrix (size=',I3,'x',I3,')'/)
      !        write(13,*) gm

      !        write(14,612) ngd, ngd
      ! 612     format(//1x,'Global stiffness matrix (size=',I3,'x',I3,')'/)
      !        write(14,*) gk


   if(beam_type == 2) then  ! tower tension wires contribution
      do i = 1, n_attachments
         glob_node = nselt+2-node_attach(i)
         nd_gv2 = 2+9*(glob_node-1)
         nd_gw2 = 4+9*(glob_node-1)
         gk(nd_gv2,nd_gv2) = gk(nd_gv2,nd_gv2)+k_tower(i)
         gk(nd_gw2,nd_gw2) = gk(nd_gw2,nd_gw2)+k_tower(i)
      end do
   end if


      ! Tip inertia contribution (using inertias per wt sign convention)

   gm(2,2) = gm(2,2) + tip_mass
   gm(4,4) = gm(4,4) + tip_mass
   gm(4,6) = gm(4,6) + tip_mass*cm_loc

   if(beam_type == 1) then
      gm(3,3) = gm(3,3) + ixx_tip         !add lag of inertia
      gm(5,5) = gm(5,5) + iyy_tip         !add flap mom of inertia
      gm(6,6) = gm(6,6) + izz_tip         !add torsion mom of inertia
      gm(3,5) = gm(3,5) - ixy_tip         !add flap-lag cross mom of inertia
      gm(5,6) = gm(5,6) - iyz_tip         !add flap-torsion cross mom of inertia
      gm(6,3) = gm(6,3) + izx_tip         !add torsion-lag cross mom of inertia
      gm(5,3) = gm(5,3) - ixy_tip         !add flap-lag cross mom of inertia
      gm(6,5) = gm(6,5) - iyz_tip         !add flap-torsion cross mom of inertia
      gm(3,6) = gm(3,6) + izx_tip         !add torsion-lag cross mom of inertia
   else
      gm(3,3) = gm(3,3) + iyy_tip         !add f-a mom of inertia
      gm(5,5) = gm(5,5) + ixx_tip         !add s-s mom of inertia
      gm(6,6) = gm(6,6) + izz_tip         !add torsion mom of inertia
      gm(3,5) = gm(3,5) + ixy_tip         !add cross mom of inertia
      gm(5,6) = gm(5,6) + izx_tip         !add cross mom of inertia
      gm(6,3) = gm(6,3) + iyz_tip         !add cross mom of inertia
      gm(5,3) = gm(5,3) + ixy_tip         !add cross mom of inertia
      gm(6,5) = gm(6,5) + izx_tip         !add cross mom of inertia
      gm(3,6) = gm(3,6) + iyz_tip         !add cross mom of inertia
   end if

      !----------------------------------------------------------------------
      !comp fix this for shear element

      ! ***  Include articulated lag spring and pitch link stiffness

   if (iartic /= 0) then
      ilg1 = 9*(nselt - 1) + 3
      ilg2 = 9*nselt + 2
      gk(ilg1,ilg1) = gk(ilg1,ilg1) + sprlag
      gk(ilg2,ilg2) = gk(ilg2,ilg2) + sprlag
      gk(ilg1,ilg2) = gk(ilg1,ilg2) - sprlag
      gk(ilg2,ilg1) = gk(ilg2,ilg1) - sprlag
   end if

      !--------------------------------------------------------------------

      !free
      !     Add tau*[gm] to [gk] for eivl-shift

      !2000      if (ifree == 1) then

   do kf=1,ngd
      do jf=1,ngd
         gk(jf,kf) = gk(jf,kf) + tau*gm(jf,kf)
      end do ! jf
   end do ! kf

      !test      print *, ' tau=', tau

      !2000end      end if

      ! *** Perform the eigenanalysis on the assembled matrices:


   call jacobi(ngd,gk,gm,evl,evr)


      !test
      !      do kk=1, ngd
      !        write(13,*) (gm(kk,jj),jj=1,ngd)
      !        write(23,*) (gk(kk,jj),jj=1,ngd)
      !     enddo

      !free  Undo the shift

      !2000      if (ifree == 1) then

   DO jf=1,ngd
      evl(jf) = evl(jf) - tau
   END DO ! jf


      ! Note: this evl (=freq squared) is normalized wrt actual omega

      !  Note gm may be singular, as long as it is symmetric positive
      !  definite. Loop on i makes the eigenvalues look like they did
      !  with the old JACOBI subroutine. Set matz to positive integer
      !  to get eigenvectors.

      !debug
      !  matz = 1
      !  call rsg (mgd,ngd,gm,gk,evl,matz,evr,fv1,fv2,ierr)
      !  if (ierr/=0) call bomb('bldvib:       rsg returned non-zero error code.')
      !  do i=1,ngd
      !     if (1.d0+evl(i)==1.d0)  call bomb("bldvib: can't handle infinite frequency yet.")
      !     evl(i) = 1.d0/evl(i)
      !  end do ! i

   write(*,*) ' '
   write(*,*)' ******** modal analysis results ', '**********'
   write(*,*) ' '

      !----------------------------------------------------------------

      ! ***    Normalize  eigenvectors: ( norm2 )

   call  normev ( evr, ngd )

      !  --------------------------------------------------------------

      ! *** Output eigenvalues and eigenvectors:

   if (ngd < modepr) then

      modepr = ngd
      Fmt    = "(/,3x,'WARNING: number of output modes requested is', A,'this exceeds the max computed', A," &
             //"/, 'therefore, only', A, 'modes are printed',/)"

      write(*   ,Fmt)  TRIM( Int2LStr( modepr ) ), TRIM( Int2LStr( ngd    ) ), TRIM( Int2LStr( ngd    ) )
      write(UnOu,Fmt)  TRIM( Int2LStr( modepr ) ), TRIM( Int2LStr( ngd    ) ), TRIM( Int2LStr( ngd    ) )

   end if

   if(beam_type == 1) then
      Fmt = "(1x,/,12x,'rotating blade frequencies & mode shapes',/,12x,'--- only first',i4,' modes printed')"
      WRITE (UnOu,Fmt) modepr
   else
      Fmt = "(1x,/,12x,'tower frequencies & mode shapes',/,12x,'--- only first',i4,' modes printed')"
      WRITE (UnOu,Fmt) modepr
   end if

   DO j=1,modepr

      jp = ngd + 1 -j
      freq = sqrt( abs( evl(jp) ) )


         ! Output eivalue and freq normalized wrt ref omega

      Fmt = "(6x,'eigenvalue(',i3,') = ', d13.6, 8x, 'mode', i3, 1x,'frequency = ', f13.6)"
      WRITE (*,Fmt)  jp, evl(jp)*omegar**2, j, freq*omegar*romg/TwoPi

         !nrel           write (UnOu,7200) jp, evl(jp)*omegar**2, freq*omegar
         !nrel           write (UnOu,7300) (i,evr(i,jp),i=1,ngd)
         !new
         !        reformat modal displacement output: show spanwise disp
         !        distribution for mode 'jp'

      if (iartic == 0) then  !modify later for other hub articulations

         Fmt =  "(//,1x,'-------- Mode No.',i4,'  (freq =',e12.5,' Hz)',/)"
         WRITE (UnOu,Fmt)  j, freq*omegar*romg/TwoPi

         if(TabDelim ) then
            if(beam_type == 1) then
              fmt_tabs="('span_loc"//tab//"flap disp"//tab//"flap slope"//tab//"lag disp"//tab//"lag slope"//tab//"twist',/)"
            else  ! for tower
              fmt_tabs="('span_loc"//tab//"s-s disp"//tab//"s-s slope"//tab//"f-a disp"//tab//"f-a slope"//tab//"twist',/)"
            end if
            write(UnOu, fmt_tabs)
         else

            IF ( beam_type == 1 )  THEN
               Fmt = "(1x,'span_loc',t13,'flap disp',t27,'flap slope',t43,'lag disp',t58,'lag slope',t73,'twist',/)"
               WRITE (UnOu,Fmt)
            ELSE
               Fmt = "(1x,'span_loc',t13,'s-s disp',t27,'s-s slope',t43,'f-a disp',t58,'f-a slope',t73,'twist',/)"
               WRITE (UnOu,Fmt) ! for tower
            END IF
         END IF

            !nrel
         xspan = xb(iend)

         if(id_form == 1) then  ! for wt
            hub_rad  = radius - bl_len
            span_loc = (xspan*radius - hub_rad)/bl_len
         end if

         if ( ifree == 1 ) then
            wdisp = evr(ngd-2,jp)
            wp    = evr(ngd-1,jp)
            vdisp = evr(ngd-4,jp)
            vp    = evr(ngd-3,jp)
            twst  = evr(ngd  ,jp)
         else
            wdisp = 0.0
            wp    = 0.0
            vdisp = 0.0
            vp    = 0.0
            twst  = 0.0
         end if

            !             write (UnOu,1213) xspan*radius, wdisp, wp,
         if(TabDelim ) then
            fmt_data1="(f7.4,'"//tab//"',f9.6,'"//tab//"',f9.6,'"//tab//"',f9.6,'"//tab//"',f9.6,'"//tab//"',f9.6)"
            write (UnOu,fmt_data1) span_loc, wdisp, wp, vdisp, vp, twst
         else
            Fmt = "(f8.4,t12,f9.6,t27,f9.6,t42,f9.6,t57,f9.6,t72,f9.6)"
            WRITE (UnOu,Fmt) span_loc, wdisp, wp, vdisp, vp, twst
         end if
         DO ie=iend,1,-1
            DO inode= 1,2

               ! inode=11: mid-node
               ! inode=22: rhs end node

               ! ispan = 2*(iend-ie) + inode
               xspan = xb(ie) + 0.5*el(ie)*inode
               if(id_form == 1) then  ! for wt
                span_loc = (xspan*radius - hub_rad)/bl_len
               end if

               if ( inode == 1 .and. mid_node_tw ) then
                  twst = evr(9*(ie-1)+8,jp)
                     !later                      if(id_form == 1) then  ! wt
                     !later                    twst = -twst
                     !later                   end if

                  if(TabDelim ) then
                     fmt_data2="(f7.4,'"//tab//" "//tab//" "//tab//" "//tab//" "//tab//"',f9.6)"
                     write (UnOu,fmt_data2) span_loc, twst
                  else
                     Fmt = "(1x,f7.4,t72,f9.6)"
                     WRITE (UnOu,Fmt) span_loc, twst
                  end if

               elseif (inode == 2) then
                  wdisp = evr(9*(ie-1)+4,jp)
                  wp    = evr(9*(ie-1)+5,jp)
                  vdisp = evr(9*(ie-1)+2,jp)
                  vp    = evr(9*(ie-1)+3,jp)
                  twst  = evr(9*(ie-1)+6,jp)

                  !later                    if(id_form == 1) then  ! wt
                  !later                  if(beam_type == 1) then
                  !later                    vdisp = -vdisp
                  !later                    vp = -vp
                  !later                    twst = -twst
                  !later                  end if
                  !later                 end if

                  if ( TabDelim )  then
                     fmt_data1="(f7.4,'"//tab//"',f9.6,'"//tab//"',f9.6,'"//tab//"',f9.6,'"//tab//"',f9.6,'"//tab//"',f9.6)"
                     write (UnOu,fmt_data1) span_loc, wdisp, wp, vdisp, vp, twst
                  else
                     Fmt = "(f8.4,t12,f9.6,t27,f9.6,t42,f9.6,t57,f9.6,t72,f9.6)"
                     write (UnOu,Fmt) span_loc, wdisp, wp, vdisp, vp, twst
                  end if

               end if

            END DO ! inode

         END DO ! ie

      end if

   END DO ! j

   write(UnOu,FmtDashes)

      !---------------------------------------------------------------------

      !jmw
   if ( idisbd == 1 )  then
      do jj=1,nmodet
         jz  =  ngd + 1 - mdsqss(jj)
      end do ! jj
   end if

      ! *** End loop on blades:   (just one blade, unless idisbd == 1 )

END DO ! jb


RETURN
END SUBROUTINE BldVib ! ( nseg, ind, mdsqss, psi0, sprlag, eflap, elag, evr, evl, gm )
!=======================================================================
