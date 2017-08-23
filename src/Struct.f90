!=======================================================================
subroutine Struct ( euj1, ek, ec, em, eq, dfx, shi, xbi, eli, l, axfi, indrns, indnl, tlam, gu, gud, dfxc, dfxcd, elv )


   !     this subroutine calculates the structural part of the response
   !     mass, stiffness, and damping matrices and load vector

   !     Some global variables used in this routine:

   !     gqp    = gauss quadrature sampling points (input); size (ngauss)
   !     gqw    = gauss quadrature weighting factors (input); size
   !              (ngauss)
   !     ngauss = number of gauss quadrature sampling points (input)
   !     nedof  = total no. of degrees of freedom for the element (input)

   !--------------------------------------------------------------------

USE CFunc
USE Conf
USE DisBd
USE DisBld
USE Gauss
USE Gravity
USE NWTC_Library
USE Omg ! March-08
USE Param
USE pbrng
USE Pitch
USE PLink
USE PreCone
USE Struc
USE Swept
USE TrimV

IMPLICIT NONE


   ! Argument declarations.

REAL(ReKi), INTENT(IN)        :: axfi                             ! Total centrifugal force at right end of element.
REAL(ReKi), INTENT(OUT)       :: dfx      (nedof,nedof)           !
REAL(ReKi), INTENT(OUT)       :: dfxc     (nselt,nedof,nedof)     !
REAL(ReKi), INTENT(OUT)       :: dfxcd    (nselt,nedof,nedof)     !
REAL(ReKi), INTENT(OUT)       :: ec       (nedof,nedof)           ! Element damping matrix.
REAL(ReKi), INTENT(OUT)       :: ek       (nedof,nedof)           ! Element stiffness matrix.
REAL(ReKi), INTENT(IN)        :: eli                              ! length of the element
REAL(ReKi), INTENT(IN)        :: elv      (nselt)                 !
REAL(ReKi), INTENT(OUT)       :: em       (nedof,nedof)           ! Element mass matrix.
REAL(ReKi), INTENT(OUT)       :: eq       (nedof)                 ! Element load vector.
REAL(ReKi), INTENT(IN)        :: euj1     (nedof)                 !
REAL(ReKi), INTENT(IN)        :: gu       (ngd)                   !
REAL(ReKi), INTENT(IN)        :: gud      (ngd)                   !
REAL(ReKi), INTENT(IN)        :: shi                              !
REAL(ReKi), INTENT(IN)        :: tlam     (3,3)                   !
REAL(ReKi), INTENT(IN)        :: xbi                              ! X coordinate of the left end of the element.

INTEGER, INTENT(IN)           :: l                                ! Element number.
INTEGER, INTENT(IN)           :: indnl                            ! Include nonlinear terms (hardwired to 1).
INTEGER, INTENT(IN)           :: indrns                           ! rotating vacuum freqs about undeformed frame (hardwired to -1).


   ! Local declarations.

REAL(ReKi)                    :: A1
REAL(ReKi)                    :: A11
REAL(ReKi)                    :: A12
REAL(ReKi)                    :: A1DOT
REAL(ReKi)                    :: A1T
REAL(ReKi)                    :: A2
REAL(ReKi)                    :: A22
REAL(ReKi)                    :: A23
REAL(ReKi)                    :: A2DOT
REAL(ReKi)                    :: A2T
REAL(ReKi)                    :: A3
REAL(ReKi)                    :: A31
REAL(ReKi)                    :: A33
REAL(ReKi)                    :: A3DOT
REAL(ReKi)                    :: A3T
REAL(ReKi)                    :: AEI
REAL(ReKi)                    :: AEITHP
REAL(ReKi)                    :: alpha    (3,3)
REAL(ReKi)                    :: AVWPS
REAL(ReKi)                    :: BETA12
REAL(ReKi)                    :: BETA13
REAL(ReKi)                    :: BETA22
REAL(ReKi)                    :: BETA23
REAL(ReKi)                    :: BETA32
REAL(ReKi)                    :: BETA33
REAL(ReKi)                    :: C2T
REAL(ReKi)                    :: CBTP
REAL(ReKi)                    :: CEICSS
REAL(ReKi)                    :: CEISCS
REAL(ReKi)                    :: cf1      (ngauss)
REAL(ReKi)                    :: CF2
REAL(ReKi)                    :: CI
REAL(ReKi)                    :: CIC
REAL(ReKi)                    :: CIN
REAL(ReKi)                    :: CIS
REAL(ReKi)                    :: CS
REAL(ReKi)                    :: CST
REAL(ReKi)                    :: CT
REAL(ReKi)                    :: CTS
REAL(ReKi)                    :: DEI
REAL(ReKi)                    :: DEIC2T
REAL(ReKi)                    :: DEICST
REAL(ReKi)                    :: DEIS2T
REAL(ReKi), ALLOCATABLE       :: delta1   (:)
REAL(ReKi), ALLOCATABLE       :: delta2   (:)
REAL(ReKi), ALLOCATABLE       :: delta3   (:)
REAL(ReKi)                    :: DSK
REAL(ReKi)                    :: DSKC2T
REAL(ReKi)                    :: DSKCST
REAL(ReKi)                    :: DSKS2T
REAL(ReKi)                    :: DUFPFP
REAL(ReKi)                    :: DUUPF
REAL(ReKi)                    :: DUUPFP
REAL(ReKi)                    :: DUUPVS
REAL(ReKi)                    :: DUUPWP
REAL(ReKi)                    :: DUUPWS
REAL(ReKi)                    :: DUVSF
REAL(ReKi)                    :: DUVSFP
REAL(ReKi)                    :: DUVSVS
REAL(ReKi)                    :: DUVSWP
REAL(ReKi)                    :: DUVSWS
REAL(ReKi)                    :: DUWPFP
REAL(ReKi)                    :: DUWSF
REAL(ReKi)                    :: DUWSWS
REAL(ReKi)                    :: EACEA
REAL(ReKi)                    :: EAEACT
REAL(ReKi)                    :: EAEAST
REAL(ReKi)                    :: EB2TCT
REAL(ReKi)                    :: EB2TST
REAL(ReKi)                    :: EGCT
REAL(ReKi)                    :: EGST
REAL(ReKi)                    :: EGXCT
REAL(ReKi)                    :: EGXST
REAL(ReKi)                    :: EL0
REAL(ReKi)                    :: ELL
REAL(ReKi), ALLOCATABLE       :: eqv      (:)
REAL(ReKi), ALLOCATABLE       :: eqvd     (:)
REAL(ReKi), ALLOCATABLE       :: eqvd1    (:)
REAL(ReKi)                    :: eqvd1l   (ngauss)
REAL(ReKi)                    :: eqvdl    (ngauss)
REAL(ReKi)                    :: eqvl     (ngauss)
REAL(ReKi)                    :: eu1      (nedof)
REAL(ReKi)                    :: eu2      (nedof)
REAL(ReKi)                    :: FI
REAL(ReKi)                    :: GJEB1T
REAL(ReKi)                    :: GQWL
REAL(ReKi)                    :: GQWML
REAL(ReKi)                    :: h        (4)
REAL(ReKi), ALLOCATABLE       :: h2m      (:,:)
REAL(ReKi), ALLOCATABLE       :: h2m1     (:,:)
REAL(ReKi)                    :: h2m1l    (ngauss,4)
REAL(ReKi)                    :: h2ml     (ngauss,4)
REAL(ReKi)                    :: hf       (3)
REAL(ReKi)                    :: HFHF
REAL(ReKi)                    :: hfp      (3)
REAL(ReKi)                    :: HFPHFP
REAL(ReKi)                    :: hfs      (3)
REAL(ReKi)                    :: HH
REAL(ReKi)                    :: HHF
REAL(ReKi)                    :: HHP
REAL(ReKi)                    :: hp       (4)
REAL(ReKi)                    :: HPH
REAL(ReKi)                    :: HPHF
REAL(ReKi)                    :: HPHFP
REAL(ReKi)                    :: HPHP
REAL(ReKi)                    :: HPHS
REAL(ReKi)                    :: HPHUP
REAL(ReKi), ALLOCATABLE       :: hpvp     (:,:)
REAL(ReKi), ALLOCATABLE       :: hpvpd    (:,:)
REAL(ReKi)                    :: hpvpdl   (ngauss,4)
REAL(ReKi)                    :: hpvpl    (ngauss,4)
REAL(ReKi), ALLOCATABLE       :: hpwp     (:,:)
REAL(ReKi), ALLOCATABLE       :: hpwpd    (:,:)
REAL(ReKi)                    :: hpwpdl   (ngauss,4)
REAL(ReKi)                    :: hpwpl    (ngauss,4)
REAL(ReKi)                    :: hs       (4)
REAL(ReKi)                    :: HSHF
REAL(ReKi)                    :: HSHFP
REAL(ReKi)                    :: HSHFS
REAL(ReKi)                    :: HSHP
REAL(ReKi)                    :: HSHS
REAL(ReKi)                    :: HSHUP
REAL(ReKi)                    :: hu       (4)
REAL(ReKi)                    :: HUH
REAL(ReKi)                    :: HUHP
REAL(ReKi)                    :: HUHU
REAL(ReKi)                    :: hup      (4)
REAL(ReKi)                    :: HUPHF
REAL(ReKi)                    :: HUPHFP
REAL(ReKi)                    :: HUPHP
REAL(ReKi)                    :: HUPHS
REAL(ReKi)                    :: PH
REAL(ReKi)                    :: PHP
REAL(ReKi)                    :: S2T
REAL(ReKi)                    :: SBTP
REAL(ReKi)                    :: SKM
REAL(ReKi)                    :: SKMCS
REAL(ReKi)                    :: SKMSC
REAL(ReKi)                    :: SKMTDD
REAL(ReKi)                    :: SS
REAL(ReKi)                    :: ST
REAL(ReKi)                    :: STS
REAL(ReKi)                    :: TFF
REAL(ReKi)                    :: TFV
REAL(ReKi)                    :: TFW
REAL(ReKi)                    :: TH0D
REAL(ReKi)                    :: TH0DD
REAL(ReKi)                    :: TH0P
REAL(ReKi), ALLOCATABLE       :: tht      (:)
REAL(ReKi)                    :: tkp      (3,3)
REAL(ReKi)                    :: TKP11
REAL(ReKi)                    :: tl       (3)
REAL(ReKi), ALLOCATABLE       :: tl1      (:)
REAL(ReKi), ALLOCATABLE       :: tl2      (:)
REAL(ReKi), ALLOCATABLE       :: tl3      (:)
REAL(ReKi)                    :: tlam2t   (nedof,nedof)
REAL(ReKi)                    :: tlami    (3,3)
REAL(ReKi)                    :: tlamk    (3,3)
REAL(ReKi)                    :: tlamkt   (3,3)
REAL(ReKi)                    :: tlamp    (3,3)
REAL(ReKi)                    :: tlamt    (3,3)
REAL(ReKi)                    :: tlddot   (3)
REAL(ReKi)                    :: tldot    (3)
REAL(ReKi)                    :: TWST
REAL(ReKi)                    :: TWX
REAL(ReKi)                    :: UNLF
REAL(ReKi)                    :: UNLFP
REAL(ReKi)                    :: UNLUP
REAL(ReKi)                    :: UNLVS
REAL(ReKi)                    :: UNLWP
REAL(ReKi)                    :: UNLWS
REAL(ReKi)                    :: UP
REAL(ReKi)                    :: VD
REAL(ReKi)                    :: VP
REAL(ReKi)                    :: VPD
REAL(ReKi), ALLOCATABLE       :: vpt      (:)
REAL(ReKi)                    :: VS
REAL(ReKi)                    :: VSPH
REAL(ReKi)                    :: VZ
REAL(ReKi)                    :: WD
REAL(ReKi)                    :: WP
REAL(ReKi)                    :: WPD
REAL(ReKi), ALLOCATABLE       :: wpt      (:)
REAL(ReKi)                    :: WS
REAL(ReKi)                    :: WSPH
REAL(ReKi)                    :: WZ
REAL(ReKi)                    :: X
REAL(ReKi)                    :: XBILS
REAL(ReKi)                    :: XI
REAL(ReKi)                    :: XT1
REAL(ReKi)                    :: XT2
REAL(ReKi)                    :: XT3
REAL(ReKi)                    :: XX
REAL(ReKi)                    :: Z1
REAL(ReKi)                    :: Z2
REAL(ReKi)                    :: Z3

INTEGER                       :: I
INTEGER                       :: I1
INTEGER                       :: I12
INTEGER                       :: I4
INTEGER                       :: I8
INTEGER                       :: ID_TW
INTEGER                       :: II
INTEGER                       :: IJK
INTEGER                       :: IK
INTEGER                       :: IKSEG
INTEGER                       :: ILSEG
INTEGER                       :: ISEC
INTEGER                       :: J
INTEGER                       :: J12
INTEGER                       :: J4
INTEGER                       :: J8
INTEGER                       :: JT
INTEGER                       :: JU
INTEGER                       :: JV
INTEGER                       :: JVV
INTEGER                       :: JW
INTEGER                       :: JWW
INTEGER                       :: K
INTEGER                       :: KK
INTEGER                       :: N
INTEGER                       :: N1
INTEGER                       :: NEFS12
INTEGER                       :: Sttus



   ! Allocate some variable-length arrays.

ALLOCATE ( delta1(nselt), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, delta1, in subroutine Struct.' )

ALLOCATE ( delta2(nselt), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, delta2, in subroutine Struct.' )

ALLOCATE ( delta3(nselt), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, delta3, in subroutine Struct.' )

ALLOCATE ( eqv(nselt), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, eqv, in subroutine Struct.' )

ALLOCATE ( eqvd(nselt), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, eqvd, in subroutine Struct.' )

ALLOCATE ( eqvd1(nselt), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, eqvd1, in subroutine Struct.' )

ALLOCATE ( h2m(nselt,4), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, h2m, in subroutine Struct.' )

ALLOCATE ( h2m1(nselt,4), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, h2m1, in subroutine Struct.' )

ALLOCATE ( hpwp(nselt,4), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, hpwp, in subroutine Struct.' )

ALLOCATE ( hpwpd(nselt,4), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, hpwpd, in subroutine Struct.' )

ALLOCATE ( hpvp(nselt,4), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, hpvp, in subroutine Struct.' )

ALLOCATE ( hpvpd(nselt,4), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, hpvpd, in subroutine Struct.' )

ALLOCATE ( tht(nselt), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, tht, in subroutine Struct.' )

ALLOCATE ( tl1(nselt), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, tl1, in subroutine Struct.' )

ALLOCATE ( tl2(nselt), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, tl2, in subroutine Struct.' )

ALLOCATE ( tl3(nselt), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, tl3, in subroutine Struct.' )

ALLOCATE ( vpt(nselt), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, vpt, in subroutine Struct.' )

ALLOCATE ( wpt(nselt), STAT=Sttus )
IF ( Sttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, wpt, in subroutine Struct.' )



   !jmw  begin dissimilar blade structure modification by jw nov. 29, 1991

if ( idisbd == 1 )  THEN

   if ( InitSt )  then

      do ii=1,nblade
         do i=1,nselt
            eiytmp(ii,i) = eiy (i)*eiyfac(ii,i)
            eiztmp(ii,i) = eiz (i)*eizfac(ii,i)
            gjtmp (ii,i) = gj  (i)*gjfac (ii,i)
            rmastm(ii,i) = rmas(i)*rmasfc(ii,i)
         end do ! i
      end do ! ii

      InitSt = .FALSE.

   end if

      !
   if ( indrns==1 .or. indrns==2 .or. indrns==-1 ) then
      eiy (l) = eiytmp(jb,l)
      eiz (l) = eiztmp(jb,l)
      gj  (l) = gjtmp (jb,l)
      rmas(l) = rmastm(jb,l)
   end if

   ! indrns is hardcoded to -1.
   !    if ( indrns==2 )  then
   !
   !       if ( eiyfac(jb,l) /= 1.0 )  then
   !          Fmt = "( ' blade ',i2,' at shi=',f6.2,'  eiy(',i2,')=',"
   !&             //"f10.6  )"
   !          write(9,Fmt) ib, shi*57.3, l,  eiy(l)
   !       end if
   !
   !       if ( eizfac(jb,l) /= 1.0 )  then
   !          Fmt = "( ' blade ',i2,' at shi=',f6.2,'  eiz(',i2,')=',"
   !&             //"f10.6  )"
   !         write(9,Fmt) jb, shi*57.3, l,  eiz(l)
   !       end if
   !
   !       if ( gjfac(jb,l) /= 1.0 )  then
   !          Fmt = "( ' blade ',i2,' at shi=',f6.2,'   gj(',i2,')=',"
   !&             //"f10.6  )"
   !         write(9,Fmt) jb, shi*57.3, l,   gj(l)
   !       end if
   !
   !       if ( rmasfc(jb,l) /= 1.0 )  then
   !          Fmt = "( ' blade ',i2,' at shi=',f6.2,' rmas(',i2,')=',"
   !&             //"f10.6  )"
   !         write(9,Fmt) jb, shi*57.3, l,   gj(l)
   !       end if
   !
   !    END IF

END IF ! if ( idisbd == 1 )

   !
   !jmw---------end dissimilar blade structure modifications--------------

   !
   !     initialize matrices and vectors

ek(:,:) = 0.0
em(:,:) = 0.0

if (indrns==1) then
   eq(:) = 0.0
end if

if (indrns==-2 .or. indrns>=1) then
   dfx(:,:) = 0.0
end if

if (indrns>=1) then
   ec(:,:) = 0.0

   dfxc (:,:,:) = 0.0
   dfxcd(:,:,:) = 0.0

   cs = cos( shi )
   ss = sin( shi )
else
   cs = 0.0
   ss = 0.0
end if

   !
   !--------------------------------------------------------------------

   !     calculate pitch twist control velocity and acceleration

th0d  = -th1c*ss + th1s*cs
th0dd = -th1c*cs - th1s*ss

   !--------------------------------------------------------------------

   !     calculate products of elastic constants

eacea  =  eac(l)*ea(l)
dei    =  eiz(l) - eiy(l)
aei    =  eiz(l) + eiy(l)
skm    =  skm1(l) + skm2(l)
dsk    =  skm2(l) - skm1(l)

   !--------------------------------------------------------------------
   !sweep
xbils  =  (xbi + eli) ** 2

   !
   !     set sweep connectivity parametres

tht(:) = 0.0
wpt(:) = 0.0
vpt(:) = 0.0

do i=1,nselt
   jt = indeg(13,i)
   jw = indeg(10,i)
   jv = indeg( 6,i)
   ju = indeg( 1,i)
   jvv= indeg( 5,i)
   jww= indeg( 9,i)
   !gsb
   !gsb        if(jt /= 0) tht(i) = ggu(jt)
   !gsb        if(jw /= 0) wpt(i) = ggu(jw)
   !gsb        if(jv /= 0) vpt(i) = ggu(jv)

end do ! i

   !agb
   !       calculate coefficients needed for swept blade
   !       centrifugal force

if ( netip /= 0 ) then

   tl1(:) = 0.0
   tl2(:) = 0.0
   tl3(:) = 0.0

   !        offset

   el0 = 1.0 - SUM( elv(:) )

   !        current segment

   k = nselt + 1 - l

   !      transformation for current swept element
   !agb
    Ptch = th75 + th1c*cs + th1s*ss

   !       transformation for current element

   do ikseg=k,nselt

      ik = nselt + 1 -ikseg

      call dist ( ik, tl, tldot, tlddot, elv, cs, ss )

      tl1(ikseg) = tl(1)
      tl2(ikseg) = tl(2)
      tl3(ikseg) = tl(3)

   end do ! ikseg

   !        calculate delta

   Ptch = th75 + th1c*cs + th1s*ss

   do ilseg=k,nselt
      i = nselt + 1 - ilseg
      Ptch = th75 + th1c*cs + th1s*ss
      tlamkt = TRANSPOSE( tlamk )
      alpha  = MATMUL   ( tlami, tlamkt )
      cbtp = cos(btp)
      sbtp = sin(btp)
      z1 = tlamk(1,1)*sbtp + tlamk(1,3)*cbtp
      z2 = tlamk(2,1)*sbtp + tlamk(2,3)*cbtp

      z3 = tlamk(3,1)*sbtp + tlamk(3,3)*cbtp

      beta13 = z1*alpha(1,2) - z2*alpha(1,1)
      beta12 = z3*alpha(1,1) - z1*alpha(1,3)
      beta23 = z1*alpha(2,2) - z2*alpha(2,1)
      beta22 = z3*alpha(2,1) - z1*alpha(2,3)
      beta33 = z1*alpha(3,2) - z2*alpha(3,1)
      beta32 = z3*alpha(3,1) - z1*alpha(3,3)

      delta1(ilseg) = z2*beta13 - z3*beta12
      delta2(ilseg) = z2*beta23 - z3*beta22
      delta3(ilseg) = z2*beta33 - z3*beta32

   end do ! ilseg

   !       calculate componenets of centrifugal force
   !       due to outboard elements

   cf2 = 0.0
   k = nselt + 1 - l
   do ilseg=k+1,nselt
      i = nselt + 1 - ilseg
      do n=1,ngauss
         gqwml = gqw(n)*rmas(i)*elv(i)
         ci = gqp(n)
         cf2 = cf2 - gqwml*( (ci*elv(i) + tl1(ilseg) )*delta1(ilseg) + tl2(ilseg)*delta2(ilseg) + tl3(ilseg)*delta3(ilseg) )
      end do ! n
   end do ! ilseg

   !       calculate components of centrifugal force
   !       due to current elements

   do kk=1,ngauss
      cf1(kk) = 0.0
   END DO ! kk

   !       calculate centrifugal force at gauss point n

   ilseg = nselt + 1 - l

   do n=1,ngauss

   !        distance of gauss point n from left end

      cin = gqp(n)

      do n1=1,ngauss

   !         distance of gauss point location in reduced element

         ci    =  cin + gqp(n1)*(1.0-cin)

         gqwml = gqw(n1)*rmas(l)*eli*(1.0 - cin)

         cf1(n) = cf1(n) - gqwml*( (ci*eli + tl1(ilseg))*delta1(ilseg) + tl2(ilseg)*delta2(ilseg) + tl3(ilseg)*delta3(ilseg) )

      END DO ! n1

   END DO ! n

end if
   !---

   !sv--------------------------------------------------------------------

   !     for response or stability (and if nonlinear terms included):

   !     calculate and store required terms for nonlinear double integrals
   !     these terms arise due to axial foreshortening

   !
IF ( ( indrns >= 1 ) .AND. ( indnl == 1 ) )  THEN
   !gsbmr
   nefs12 = nefs(1) + nefs(2)

   if ( ( l == 1 ) .AND. ( nselt > 1 ) )  THEN

      hpwp (:,:) = 0.0
      hpwpd(:,:) = 0.0
      hpvp (:,:) = 0.0
      hpvpd(:,:) = 0.0
      h2m  (:,:) = 0.0
      h2m1 (:,:) = 0.0

      do i=1,nselt

         if ( ielrot(i) == 1 )  THEN
               !agb
            Ptch = th75 + th1c*cs + th1s*ss
               !agb
            tlamp = TRANSPOSE( tlamt )
            tkp   = MATMUL   ( tlam, tlamp )
            tkp11 = tkp(1,1)

            call evfrsw(gu ,eu1,indeg(1,i),tlam2t)
            call evfrsw(gud,eu2,indeg(1,i),tlam2t)
         else
            call evfrgv(gu ,eu1,indeg(1,i))
            call evfrgv(gud,eu2,indeg(1,i))
         end if ! ( ielrot(i) == 1 )

         if ( ielrot(i) == 1 )  then
            cbtp = cos(btp)
            sbtp = sin(btp)
               !agb
            a1 = tlamt(1,1)*( sbtp + th0d ) + tlamt(1,3)*cbtp
            a2 = tlamt(2,1)*( sbtp + th0d ) + tlamt(2,3)*cbtp
            a3 = tlamt(3,1)*( sbtp + th0d ) + tlamt(3,3)*cbtp

         end if ! ( ielrot(i) == 1 )

         eqv(i)  = 0.0
         eqvd(i) = 0.0
         eqvd1(i)= 0.0

         ell = elv(i)

         !gsbmr
         if ( i <= nefs12 ) then

            DO n=1,ngauss

               gqwl  =  gqw(n)*ell
               gqwml =  gqwl*rmas(i)
               ci    =  gqp(n)
               cis   =  ci*ci
               cic   =  ci*cis

               h(1) = 2.*cic - 3.*cis + 1.
               h(2) = ell* (cic - 2.*cis + ci)
               h(3) =-2.*cic + 3.*cis
               h(4) = ell* (cic - cis)

               hp(1) = 6.*(cis - ci) / ell
               hp(2) = 3.*cis - 4.*ci + 1.
               hp(3) =-hp(1)
               hp(4) = 3.*cis - 2.*ci

               ijk = 1

               vp  = hp(1)*eu1(ijk+4) + hp(2)*eu1(ijk+5) + hp(3)*eu1(ijk+6 ) + hp(4)*eu1(ijk+7 )
               wp  = hp(1)*eu1(ijk+8) + hp(2)*eu1(ijk+9) + hp(3)*eu1(ijk+10) + hp(4)*eu1(ijk+11)
               vpd = hp(1)*eu2(ijk+4) + hp(2)*eu2(ijk+5) + hp(3)*eu2(ijk+6 ) + hp(4)*eu2(ijk+7 )
               wpd = hp(1)*eu2(ijk+8) + hp(2)*eu2(ijk+9) + hp(3)*eu2(ijk+10) + hp(4)*eu2(ijk+11)
               vd  = h (1)*eu2(ijk+4) + h (2)*eu2(ijk+5) + h (3)*eu2(ijk+6 ) + h (4)*eu2(ijk+7 )
               wd  = h (1)*eu2(ijk+8) + h (2)*eu2(ijk+9) + h (3)*eu2(ijk+10) + h (4)*eu2(ijk+11)

               if ( ielrot(i) == 1 ) then
                     !agb
                  eqv(i)  = eqv(i)  + gqwl*(vp*vpd+wp*wpd)*tkp11
                  eqvd(i) = eqvd(i) + gqwml *2.0* a3* vd
                  eqvd1(i) = eqvd1(i) + gqwml* 2.0*a2*wd
               else
                  eqv(i)  = eqv(i)  + gqwl*(vp*vpd+wp*wpd)
                  eqvd(i) = eqvd(i) + gqwml*2.0*vd
               end if ! ( ielrot(i) == 1 )

               DO i1=1,4
                  hpwpd(i,i1) = hpwpd(i,i1) + wpd*hp(i1)*gqwl
                  hpwp (i,i1) = hpwp (i,i1) + wp *hp(i1)*gqwl
                  hpvpd(i,i1) = hpvpd(i,i1) + vpd*hp(i1)*gqwl
                  hpvp (i,i1) = hpvp (i,i1) + vp *hp(i1)*gqwl
               END DO ! i1

               if ( ielrot(i) == 1 )  then
                  DO i1=1,4
                     h2m (i,i1) = h2m (i,i1) + 2.0*gqwml*a3*h(i1)
                     h2m1(i,i1) = h2m1(i,i1) + 2.0*gqwml*a2*h(i1)
                  END DO ! i1
               else
                  DO i1=1,4
                     h2m  (i,i1) = h2m  (i,i1) + 2.0*gqwml*h(i1)
                  END DO ! i1
               end if ! ( ielrot(i) == 1 )

            END DO ! n

         end if ! ( i <= nefs12 )

      end do ! i

   end if ! ( ( l == 1 ) .AND. ( nselt > 1 ) )

   !--------------------------------------------------------------------

   !       nonlinear foreshortening terms for "current" (lth) element

   hpwpl (:,:) = 0.0
   hpwpdl(:,:) = 0.0
   hpvpl (:,:) = 0.0
   hpvpdl(:,:) = 0.0
   h2ml  (:,:) = 0.0
   h2m1l (:,:) = 0.0

   if ( ielrot(l) == 1 )  then
         !agb
      Ptch = th75 + th1c*cs + th1s*ss
         !agb
      tlamp = TRANSPOSE( tlamt )
      tkp   = MATMUL   ( tlam, tlamp )
      tkp11 = tkp(1,1)

      call evfrsw(gu ,eu1,indeg(1,l),tlam2t)
      call evfrsw(gud,eu2,indeg(1,l),tlam2t)
   else
      call evfrgv(gu ,eu1,indeg(1,l))
      call evfrgv(gud,eu2,indeg(1,l))
   end if ! ( ielrot(l) == 1 )

   if ( ielrot(l) == 1 ) then
      cbtp = cos(btp)
      sbtp = sin(btp)
         !agb
      a1 = tlamt(1,1)*( sbtp + th0d ) + tlamt(1,3)*cbtp
      a2 = tlamt(2,1)*( sbtp + th0d ) + tlamt(2,3)*cbtp
      a3 = tlamt(3,1)*( sbtp + th0d ) + tlamt(3,3)*cbtp
   end if ! ( ielrot(l) == 1 )


      ! Gauss loop over current element.

   DO n=1,ngauss

      eqvl(n)  = 0.0
      eqvdl(n) = 0.0
      eqvd1l(n) = 0.0
      cin   =  gqp(n)

         ! gauss loop over "partial current" element
         ! "partial" = length only up to gauss point n

         !gsbmr
      if ( l <= nefs12 )  then

         DO n1=1,ngauss

            gqwl  =  gqw(n1)*eli*cin
            gqwml =  gqwl*rmas(l)
            ci    =  gqp(n1)*cin
            cis   =  ci*ci
            cic   =  ci*cis

            h(1) = 2.*cic - 3.*cis + 1.
            h(2) = eli* (cic - 2.*cis + ci)
            h(3) =-2.*cic + 3.*cis
            h(4) = eli* (cic - cis)

            hp(1) = 6.*(cis - ci) / eli
            hp(2) = 3.*cis - 4.*ci + 1.
            hp(3) =-hp(1)
            hp(4) = 3.*cis - 2.*ci

            ijk = 1

            vp  = hp(1)*eu1(ijk+4) + hp(2)*eu1(ijk+5) + hp(3)*eu1(ijk+6 ) + hp(4)*eu1(ijk+7 )
            wp  = hp(1)*eu1(ijk+8) + hp(2)*eu1(ijk+9) + hp(3)*eu1(ijk+10) + hp(4)*eu1(ijk+11)
            vpd = hp(1)*eu2(ijk+4) + hp(2)*eu2(ijk+5) + hp(3)*eu2(ijk+6 ) + hp(4)*eu2(ijk+7 )
            wpd = hp(1)*eu2(ijk+8) + hp(2)*eu2(ijk+9) + hp(3)*eu2(ijk+10) + hp(4)*eu2(ijk+11)
            vd  = h (1)*eu2(ijk+4) + h (2)*eu2(ijk+5) + h (3)*eu2(ijk+6 ) + h (4)*eu2(ijk+7 )
            wd  = h (1)*eu2(ijk+8) + h (2)*eu2(ijk+9) + h (3)*eu2(ijk+10) + h (4)*eu2(ijk+11)

            eqvl(n)  = eqvl(n) + gqwl*(vp*vpd+wp*wpd)

            DO i1=1,4
               hpwpdl(n,i1) = hpwpdl(n,i1) + wpd*hp(i1)*gqwl
               hpwpl (n,i1) = hpwpl (n,i1) + wp *hp(i1)*gqwl
               hpvpdl(n,i1) = hpvpdl(n,i1) + vpd*hp(i1)*gqwl
               hpvpl (n,i1) = hpvpl (n,i1) + vp *hp(i1)*gqwl
            END DO ! i1

            gqwl  =  gqw(n1)*eli*(1.0-cin)
            gqwml =  gqwl*rmas(l)
            ci    =  cin + gqp(n1)*(1.0-cin)
            cis   =  ci*ci
            cic   =  ci*cis

            h(1) = 2.*cic - 3.*cis + 1.
            h(2) = eli* (cic - 2.*cis + ci)
            h(3) =-2.*cic + 3.*cis
            h(4) = eli* (cic - cis)

            vd = h(1)*eu2(ijk+4) + h(2)*eu2(ijk+5) + h(3)*eu2(ijk+6 ) + h(4)*eu2(ijk+7 )
            wd = h(1)*eu2(ijk+8) + h(2)*eu2(ijk+9) + h(3)*eu2(ijk+10) + h(4)*eu2(ijk+11)

            if(ielrot(l)==1)then
               eqvdl(n) = eqvdl(n) + gqwml*2.0*vd*a3
               eqvd1l(n)= eqvd1l(n)+ gqwml*2.0*wd*a2
            else
               eqvdl(n) = eqvdl(n) + gqwml*2.0*vd
            end if

            if(ielrot(l)==1)then
                DO i1=1,4
                   h2ml (n,i1) = h2ml (n,i1) + 2.0*gqwml*h(i1)*a3
                   h2m1l(n,i1) = h2m1l(n,i1) + 2.0*gqwml*h(i1)*a2
                END DO ! i1
            else
                DO i1=1,4
                   h2ml(n,i1) = h2ml(n,i1) + 2.0*gqwml*h(i1)
                END DO ! i1
            end if

         END DO ! n1

      END IF ! ( l <= nefs12 )

   END DO ! n

END IF ! ( ( indrns >= 1 ) .AND. ( indnl == 1 ) )

   !     end of preliminary foreshortening calculations

   !sv--------------------------------------------------------------------

   !     begin calculations for matrices and load vector for "lth" element

   !
   !sweep
if ( ielrot(l) == 1 )  then
      !agb
    Ptch = th75 + th1c*cs + th1s*ss

   cbtp = cos(btp)
   sbtp = sin(btp)
      !agbx
   a1 = tlam(1,1)*( sbtp + th0d ) + tlam(1,3)*cbtp
   a2 = tlam(2,1)*( sbtp + th0d ) + tlam(2,3)*cbtp
   a3 = tlam(3,1)*( sbtp + th0d ) + tlam(3,3)*cbtp

   a1dot = th0dd*tlam(1,1) + th0d*cbtp*tlam(1,2)
   a2dot = th0dd*tlam(2,1) + th0d*cbtp*tlam(2,2)
   a3dot = th0dd*tlam(3,1) + th0d*cbtp*tlam(3,2)

   a11 = a1*a1
   a22 = a2*a2
   a33 = a3*a3
   a12 = a1*a2
   a23 = a2*a3
   a31 = a3*a1
   a1t = a1*tlam(1,1)
   a2t = a2*tlam(2,1)
   a3t = a3*tlam(3,1)
      !        a12t = a1t*a2t
      !        a23t = a2t*a3t
      !        a31t = a3t*a1t
      !agb
   call dist ( l, tl, tldot, tlddot, elv, cs, ss )

   xt1 = tl(1)
   xt2 = tl(2)
   xt3 = tl(3)
      !        xt1 = xbi*tlam(1,1)
      !        xt2 = xbi*tlam(2,1)
      !        xt3 = xbi*tlam(3,1)
end if

   !--------------------------------------------------------------------

   ! +++++++++++++++++++++++    gauss loop over "current (lth)" element

DO n=1,ngauss

   gqwl  =  gqw(n)*eli
   gqwml =  gqwl*rmas(l)
   ci    =  gqp(n)
   cis   =  ci*ci
   cic   =  ci*cis
   xi    =  ci*eli

   if ( ielrot(l) == 1 )  then
      x  = xt1 + xi
      xx = xbi + xi*tlam(1,1)
   else
      x = xbi + xi
   end if ! ( ielrot(l) == 1 )
      !     xs    =  x*x

      !--------------------------------------------------------------------
      ! Calculate shape functions

   hu(1) =- 4.5*cic +  9.0*cis - 5.5*ci + 1.
   hu(2) = 13.5*cic - 22.5*cis + 9.0*ci
   hu(3) =-13.5*cic + 18.0*cis - 4.5*ci
   hu(4) =  4.5*cic -  4.5*cis + ci

   hup(1) = (-13.5*cis + 18.*ci - 5.5) / eli
   hup(2) = ( 40.5*cis - 45.*ci + 9. ) / eli
   hup(3) = (-40.5*cis + 36.*ci - 4.5) / eli
   hup(4) = ( 13.5*cis -  9.*ci + 1. ) / eli

   h(1) = 2.*cic - 3.*cis + 1.
   h(2) = eli*(cic - 2.*cis + ci)
   h(3) =-2.*cic + 3.*cis
   h(4) = eli*(cic - cis)

   hp(1) = 6.*(cis - ci) / eli
   hp(2) = 3.*cis - 4.*ci + 1.
   hp(3) =-hp(1)
   hp(4) = 3.*cis - 2.*ci

   hs(1) = (12.*ci - 6.) / eli / eli
   hs(2) = ( 6.*ci - 4.) / eli
   hs(3) =  -hs(1)
   hs(4) = ( 6.*ci - 2.) / eli

   hf(1) = 2.*cis - 3.*ci + 1.
   hf(2) =-4.*cis + 4.*ci
   hf(3) = 2.*cis - ci

   hfp(1) = ( 4.*ci - 3.) / eli
   hfp(2) = (-8.*ci + 4.) / eli
   hfp(3) = ( 4.*ci - 1.) / eli

   hfs(1) = ( 4./ eli) / eli
   hfs(2) = (-8./ eli) / eli
   hfs(3) =  (4./ eli) / eli

   !--------------------------------------------------------------------

   !     calculate deflections (for use in non-linear terms)

   ijk = 1

   up  = hup(1)*euj1(ijk   ) + hup(2)*euj1(ijk+1 ) + hup(3)*euj1(ijk+2 ) + hup(4)*euj1(ijk+3 )
   vz  = h  (1)*euj1(ijk+4 ) + h  (2)*euj1(ijk+5 ) + h  (3)*euj1(ijk+6 ) + h  (4)*euj1(ijk+7 )
   wz  = h  (1)*euj1(ijk+8 ) + h  (2)*euj1(ijk+9 ) + h  (3)*euj1(ijk+10) + h  (4)*euj1(ijk+11)
   ph  = hf (1)*euj1(ijk+12) + hf (2)*euj1(ijk+13) + hf (3)*euj1(ijk+14)
   vp  = hp (1)*euj1(ijk+4 ) + hp (2)*euj1(ijk+5 ) + hp (3)*euj1(ijk+6 ) + hp (4)*euj1(ijk+7 )
   wp  = hp (1)*euj1(ijk+8 ) + hp (2)*euj1(ijk+9 ) + hp (3)*euj1(ijk+10) + hp (4)*euj1(ijk+11)
   php = hfp(1)*euj1(ijk+12) + hfp(2)*euj1(ijk+13) + hfp(3)*euj1(ijk+14)
   vs  = hs (1)*euj1(ijk+4 ) + hs (2)*euj1(ijk+5 ) + hs (3)*euj1(ijk+6 ) + hs (4)*euj1(ijk+7 )
   ws  = hs (1)*euj1(ijk+8 ) + hs (2)*euj1(ijk+9 ) + hs (3)*euj1(ijk+10) + hs (4)*euj1(ijk+11)

      !--------------------------------------------------------------------

      !     calculate rigid pitch angle (pretwist + control pitch)

      !nrel if (itwist == 0) then
      !        th0  = th75 + th0p*(x-.75) + th1c*cs + th1s*ss

      !  new computation of blade pre-twist

      !  interpolate for twx(x)
   id_tw = 0

   do isec=1,n_secs-1
      IF ( ( x >= sec_loc(isec)-tiny(x) ) .AND. ( x < sec_loc(isec+1) ) ) THEN
         twx   = str_tw(isec) + ( str_tw(isec+1)-  str_tw(isec) )*( x - sec_loc(isec) )/( sec_loc(isec+1) - sec_loc(isec) )
         id_tw = id_tw+1
         exit
      end if ! ( ( x >= sec_loc(isec)-tiny(x) ) .AND. ( x < sec_loc(isec+1) ) )
   end do ! isec

   if ( id_tw == 0 )  then
      CALL ProgAbort ( 'Interpolation failed for TWX in struct().' )
   end if ! ( id_tw == 0 )

   if ( id_tw >= 2 )  CALL ProgAbort ( 'multiply interpolated' )

   th0p = 0.  ! (modify later)

      !  add pre-twist to pitch setting
   th0 = th75 + twx

      !      else
      !        call twistx(x,th75,thx,th0p,l)
      !       th0 = thx + th1c*cs + th1s*ss
      !--      end if


      !gsb : check the following statment
   if ( iartic == 3 )  then
      th0 = th0 + wp*(vp - del3pb)
   end if

   !     compute transformation due to pitch control

   ct = cos( th0 )
   st = sin( th0 )
   !agb

   !     for swept blade, ct and st are calculated using
   !     blade pitch

   if ( ielrot(l) == 1 ) then

         !     collective  pitch

      Ptch = th75 + th1c*cs + th1s*ss

      !     blade twist

      !           twst  = th0p*(x-.75)
      !         else
      !           call twistx(x,th75,thx,th0p,l)
      !           twst  = thx - th75
      !        end if
      !nrel
      twst = th0   !(modify later for curved blade)
      !--
      ct     =  cos(twst)
      st     =  sin(twst)
      th0d = 0.0  ! visit later
      th0dd = 0.0 ! visit later

   end if ! ( ielrot(l) == 1 )

   !---
   skmtdd =  skm*th0dd
   cst   =  ct*st
   cts   =  ct*ct
   sts   =  st*st
   c2t   =  cts - sts
   s2t   =  2.* cst

   aeithp =  aei*th0p
   gjeb1t =  gj(l) + eb1(l)*th0p*th0p
   eb2tct =  eb2(l)*th0p*ct
   eb2tst =  eb2(l)*th0p*st
   egct   =  eg(l)*ct
   egst   =  eg(l)*st
   egxct  =  egct *x
   egxst  =  egst *x

   tfv    =  egct  + th0dd*egst
   tfw    = -btp*x - th0dd*egct
   tff    = -skmtdd - dsk*cst - btp*egxct

   eaeact =  eacea*ct
   eaeast =  eacea*st

   ceicss =  eiz(l)*cts + eiy(l)*sts
   ceiscs =  eiz(l)*sts + eiy(l)*cts
   deicst =  dei*cst
   deic2t =  dei*c2t
   deis2t =  dei*s2t
   dskc2t =  dsk*c2t
      !agb
   dsks2t =  dsk*s2t
      !---
   dskcst = dsk*cst
   skmsc  = skm1(l)*sts + skm2(l)*cts
   skmcs  = skm1(l)*cts + skm2(l)*sts

   vsph  =  vs*ph
   wsph  =  ws*ph
   avwps =  .5*(vp*vp + wp*wp)

      !--------------------------------------------------------------------
      !sweep
      !agb
   if ( ielrot(l) /= 1 )  then
!      fi = axfi + 0.5*rmas(l)*( xbils - x**2 )           !old
      fi = ( axfi + 0.5*rmas(l)*( xbils - x**2 ) )*omega2 !march-08
      !tower -- for gravity-induced com load (e.g. for a tower)
      !        fi = (axfi - grav *rmas(l)*(eli - xi))*omega2
      !--
   end if
      !agb
      !     calculate centrifugal force for advanced geometry blade

   if ( netip /= 0 )  then
      fi = (cf1(n) +  cf2)*omega2
   end if

      !--------------------------------------------------------------------
      ! Calculate linear contributions to element matrices and vector

   do i=1,4

      i4  = i + 4
      i8  = i + 8

         !sweep
      if ( indrns == 1 )  then

         if ( ielrot(l) == 1 )  then

            eq(i) = eq(i) + gqwml*hu(i)*(x*(a22+a33) - a1*(a2*xt2+a3*xt3)-a1*(a2*egct + a3*egst))  &
               !agb
                          + gqwml*hu(i)*(a3dot*xt2 - a2dot*xt3)

            eq(i4) = eq(i4) + gqwml*h(i)*((a11 + a33)*(xt2+egct) - a2*(a1*x +a3*xt3+a3*egst) + th0dd*egst  &
                   + 2.0*a1*th0d*egct) + gqwml*hp(i)* (-(a22+a33)*egxct + a1*(a2*xt2+a3*xt3)*egct + a12*skmsc +a31*dskcst)  &
               !agb
                   + gqwml*h(i) *(-a3dot*x + a1dot*xt3) + gqwml*hp(i)*(-a2dot*ct*ct*skm1(l) - a3dot*ct*ct*skm1(l)  &
                   - a3dot*skm2(l)*ct*ct + a2dot*skm2(l)*ct*st)

            eq(i8) = eq(i8) +gqwml*h(i)* ( (a11+a22)*(xt3+egst) - a3*(a2 *xt2 +a2*egct +a1*x) - th0dd*egct  &
                   + 2.0*a1*th0d*egst ) + gqwml*hp(i)* (-(a22+a33)*egxst + a1*(a2*xt2+a3*xt3)*egst + a31*skmcs +a12*dskcst)  &
               !agb
                   + gqwml*h(i) *(a2dot*x - a1dot*xt2) + gqwml*hp(i)*(a2dot*skm1(l)*ct*ct + a3dot*skm1(l)*ct*ct  &
                   - a3dot*skm2(l)*ct*st + a2dot*skm2(l)*st*st)

         else

            eq(i)  = eq(i)  + gqwml*hu(i)*x
            eq(i4) = eq(i4) + gqwml*(h(i)*tfv - hp(i)*egxct)
            eq(i8) = eq(i8) + gqwml*(h(i)*tfw - hp(i)*egxst)

         end if ! ( ielrot(l) == 1 )

      end if ! ( indrns == 1 )

      DO j=1,4

         j4  = j + 4
         j8  = j + 8

            !sweep
         huhu = hu(i)*hu(j)
         huh  = hu(i)*h(j)
         huhp = hu(i)*hp(j)
         hph  = hp(i)*h(j)

         huphs = hup(i)*hs(j)
         hh    = h(i)  *h(j)
         hhp   = h(i)  *hp(j)
         hphp  = hp(i) *hp(j)
         hshs  = hs(i) *hs(j)


         if ( ielrot(l) == 1 ) then

            ek(i ,j ) = ek(i ,j ) + gqwl*hup(i)*hup(j)*eac(l) - gqwml*huhu*(a22 + a33)
            ek(i ,j4) = ek(i ,j4) - gqwl*huphs*eaeact + gqwml*huh*a12
            ek(i ,j8) = ek(i ,j8) - gqwl*huphs*eaeast + gqwml*huh*a31
            ek(i4,j4) = ek(i4,j4) - gqwml* hh*(a11 + a33) + gqwl*(hphp*fi + hshs*ceicss) - gqwml*(hhp+hph)*a12*egct  &
                                  + gqwml*(a3*th0d*hhp*egst*2.0)
            ek(i4,j8) = ek(i4,j8) + gqwl*hshs*deicst - gqwml* (-hh*a23 -hhp*a12*egst-a31*hph*egct)

            ek(i8,j8) = ek(i8,j8) + gqwl*(hphp*fi + hshs*ceiscs) - gqwml*hh*(a11+a22) - gqwml*(hhp+hph)*a31*egst

         else

            ek(i ,j ) = ek(i,j)   + gqwl*hup(i)*hup(j)*eac(l)

            ek(i ,j4) = ek(i,j4)  - gqwl*huphs *eaeact
            ek(i ,j8) = ek(i,j8)  - gqwl*huphs *eaeast
!           ek(i4,j4) = ek(i4,j4) - gqwml* hh + gqwl*(hphp*fi + hshs*ceicss)                     ! old
            ek(i4,j4) = ek(i4,j4) - gqwml* omega2*hh + gqwl*(hphp*fi + hshs*ceicss + hh*elm_dist_k(l))  !new- omega2 and additional dist stiffness feb-08
            ek(i4,j8) = ek(i4,j8) + gqwl*hshs  *deicst
!           ek(i8,j8) = ek(i8,j8) + gqwl*(hphp*fi + hshs*ceiscs)                      ! old
            ek(i8,j8) = ek(i8,j8) + gqwl*(hphp*fi + hshs*ceiscs + hh*elm_dist_k(l))   !new- additional dist stiffness feb-08

         end if ! ( ielrot(l) == 1 )

         ek(j4,i)  = ek(i,j4)
         ek(j8,i)  = ek(i,j8)
         ek(j8,i4) = ek(i4,j8)

            !agb
         if ( ielrot(l) == 1 )  then
            ek(i ,j4) = ek(i ,j4) - gqwml*huh*a3dot
            ek(j4,i ) = ek(j4,i ) + gqwml*huh*a3dot
            ek(i ,j8) = ek(i ,j8) + gqwml*huh*a2dot
            ek(j8,i ) = ek(j8,i ) - gqwml*huh*a2dot
            ek(i4,j8) = ek(i4,j8) - gqwml*hh *a1dot
            ek(j8,i4) = ek(j8,i4) + gqwml*hh *a1dot
         end if ! ( ielrot(l) == 1 )

         if ( indrns >= 1 )  then

            if ( ielrot(l) == 1 ) then
               ec(i,j4)  = ec(i,j4) - gqwml*hu(i)*h(j)*2.0*a3
                  !agb
               ec(i,j8)  = ec(i,j8) + gqwml*hu(i)*h(j)*2.0*a2

               ec(i4,j4) = ec(i4,j4)+gqwml *(hp(i)*h(j) - hhp)*2.*a3*egct
               ec(i4,j8) = ec(i4,j8) - gqwml *(hh*a1+a3*hhp*egst)*2.0 - gqwml *hph*egct*a2*2.
            else
               ec(i,j4)  = ec(i,j4)  - gqwml*omega*hu(i)*h(j)*2.0
               ec(i4,j4) = ec(i4,j4) + gqwml*omega *(hp(i)*h(j) - hhp)*2.* egct
               ec(i4,j8) = ec(i4,j8) - gqwml*omega *(hh*btp + hhp*egst)*2.
               
!               ec(i,j4)  = ec(i,j4)  - gqwml*hu(i)*h(j)*2.0 !old
!               ec(i4,j4) = ec(i4,j4) + gqwml *(hp(i)*h(j) - hhp)*2.* egct !old
!               ec(i4,j8) = ec(i4,j8) - gqwml *(hh*btp + hhp*egst)*2. !old
            end if ! ( ielrot(l) == 1 )

            ec(j4,i)  =-ec(i,j4)
             !agb
            ec(j8,i)  = -ec(i,j8)
             !---
            ec(j8,i4) = -ec(i4,j8)

         end if ! ( indrns >= 1 )

         em(i,j)   = em(i,j)   + gqwml*hu(i)*hu(j)
!        em(i4,j4) = em(i4,j4) + gqwml*hh     !old
!        em(i8,j8) = em(i8,j8) + gqwml*hh     !old
         em(i4,j4) = em(i4,j4) + (gqwml + gqwl*elm_dist_m(l))*hh  !new (added mass distr)
         em(i8,j8) = em(i8,j8) + (gqwml + gqwl*elm_dist_m(l))*hh  !new (added mass distr)

      END DO ! J

      DO j=1,3

         j12 = j + 12

         hhf   = h(i)  *hf(j)
         hphf  = hp(i) *hf(j)
         hshf  = hs(i) *hf(j)
         hshfp = hs(i) *hfp(j)
         hshfs = hs(i) *hfs(j)

         if ( ielrot(l) == 1 ) then

            ek(i,j12)   = ek(i,j12)  + gqwl*hup(i)*hfp(j)*aeithp + gqwml*hu(i)*hf(j)*(a31*egct-a12*egst)
            ek(i4,j12)  = ek(i4,j12) + gqwml*(hhf*(a11+a33)*egst - hphf*(a22+a33)*egxst + hhf*a23*egct) &
                        + gqwml*hphf*xbi*a1*(a2t+a3t)*egst - gqwl*(hshfp*eb2tct + hshfs*ec2(l)*st)
            ek(i8,j12)  = ek(i8,j12) + gqwml *hphf*(a22+a33)*egxct - gqwml*hhf*(a11+a22)*egct - gqwml*hhf*a23*egst &
                        - gqwml*hphf*xbi*a1*(a2t+a3t)*egct - gqwl *(hshfp*eb2tst - hshfs*ec2(l)*ct)
         else
            ek(i,j12)   = ek(i,j12)  + gqwl *hup(i)*hfp(j)*aeithp

            ek(i4,j12)  = ek(i4,j12) + gqwml*omega2*(hhf*egst - hphf*egxst) - gqwl *(hshfp*eb2tct + hshfs*ec2(l)*st) !march-08 omega2 
 !           ek(i4,j12)  = ek(i4,j12) + gqwml*(hhf*egst - hphf*egxst) - gqwl *(hshfp*eb2tct + hshfs*ec2(l)*st) !old
            ek(i8,j12)  = ek(i8,j12) + gqwml*omega2*hphf*egxct - gqwl*(hshfp*eb2tst - hshfs*ec2(l)*ct) !march-08 omega2 
 !           ek(i8,j12)  = ek(i8,j12) + gqwml*hphf*egxct - gqwl*(hshfp*eb2tst - hshfs*ec2(l)*ct) !old

         end if ! ( ielrot(l) == 1 )

         ek(j12,i)   =  ek(i,j12)
         ek(j12,i4)  =  ek(i4,j12)
         ek(j12,i8)  =  ek(i8,j12)

         em(i4,j12)  =  em(i4,j12) - gqwml*hhf*egst
         em(i8,j12)  =  em(i8,j12) + gqwml*hhf*egct
         em(j12,i4)  =  em(i4,j12)
         em(j12,i8)  =  em(i8,j12)

      end do ! j

   end do ! i

   do i=1,3

      i12 = i + 12

      if ( indrns == 1 )  then

         if ( ielrot(l) == 1 )  then
           eq(i12) = eq(i12) +gqwml*hf(i)*( dsk*cst*(a11+a22) - (a11 +a33)*xt2*egst + (a11 +a22)*xt3*egct  &
                   - a31*x*egct + a12*x*egst - a23*(xt2*egst-xt3*egst) - skmtdd + a23*(skmcs-skmsc)  )  &
            !agb
                   + gqwml*hf(i)*(-2.0*a1dot*skm1(l)*ct*ct - a1dot*skm2(l))

         else

            eq(i12) = eq(i12) + gqwml*hf(i)*tff

         end if ! ( ielrot(l) == 1 )

      end if ! ( indrns == 1 )

      do j=1,3

         j12  =  j  +  12
         hfhf   =  hf(i)*hf(j)

         if ( ielrot(l) == 1 )  then
            ek(i12,j12) = ek(i12,j12) + gqwml*hfhf*(a33-a22)*dskc2t*omega2 + gqwl*hfp(i)*hfp(j)*gjeb1t + gqwml*hfhf*xt2*((a11+a33)*egct  &
                        - a23*egst) + gqwml*hfhf*xt3*( (a11+a33)*egst - a23*egct ) + gqwl*hfs(i)*hfs(j)*ec1(l)  &
               !agb
                        - gqwml*hfhf*2.0*a2*a3*dsks2t
         else
            ek(i12,j12)  = ek(i12,j12) + gqwml*hfhf*dskc2t*omega2 + gqwl *hfp(i)*hfp(j)*gjeb1t + gqwl *hfs(i)*hfs(j)*ec1(l) !march-08 omega2
!            ek(i12,j12)  = ek(i12,j12) + gqwml*hfhf*dskc2t + gqwl *hfp(i)*hfp(j)*gjeb1t + gqwl *hfs(i)*hfs(j)*ec1(l) !old
         end if ! ( ielrot(l) == 1 )

         em(i12,j12) = em(i12,j12) + gqwml*hfhf*skm

      end do ! j

   end do ! i

      !----------------------------------------------------------------------
      ! Calculate non-linear element load vector and stiffness matrix

   if ( ( indnl /= 0 ) .AND. ( indrns /= -1 ) )  THEN

         !
         ! Non-linear load vector terms

      if ( indrns == 1 )  then

         unlup =  eac(l)*(ea(l)*(vsph*st - wsph*ct)) + 0.5*(aei*php*php) + aeithp*wp*vs

         unlvs = -deis2t*vsph  + deic2t*wsph - eaeact*avwps +  eaeast*up*ph + aeithp*wp*up +  gjeb1t*wp*php

         unlws = deic2t*vsph  + deis2t*wsph - eaeast*avwps - eaeact*up*ph

         unlwp = aeithp*up*vs + gjeb1t*vs*php

         unlf  = deicst*ws*ws + deic2t*vs*ws - deicst*vs*vs

         unlfp = gjeb1t*wp*vs + aei*up*php

         DO i=1,4

            i4  =  i + 4
            i8  =  i + 8
            i12 =  i + 12

            eq(i)  = eq(i)  - gqwl*hup(i)*unlup

            eq(i4) = eq(i4) - gqwl*hs(i) *unlvs
            eq(i8) = eq(i8) - gqwl*hs(i) *unlws - gqwl*hp(i)*unlwp


               !sv--------------------------------------------------------------------

               ! 0 to x double integral term

            if ( l < nselt )  then
               if ( ielrot(l) == 1 )  then
                  do i1=nselt,l+1,-1
                     eq(i4) = eq(i4) + 2.0*gqwml*a3*eqv(i1)*h(i)
                     eq(i8) = eq(i8) - 2.0*gqwml*a2*eqv(i1)*h(i)
                  end do ! i1
               else
                  do i1=nselt,l+1,-1
                     eq(i4) = eq(i4) + 2.0*gqwml*eqv(i1)*h(i)
                  end do ! i1
               end if ! ( ielrot(l) == 1 )
            end if ! ( l < nselt )

               ! x to 1 double integral term

            if ( l > 1 )  then
               if ( ielrot(i) == 1 )  then
                  DO i1=1,l-1
                     eq(i4) = eq(i4) - gqwl*vp*eqvd(i1)*hp(i) + gqwl*vp*eqvd1(i1)*hp(i)
                     eq(i8) = eq(i8) - gqwl*wp*eqvd(i1)*hp(i) + gqwl*wp*eqvd1(i1)*hp(i)
                  END DO ! i1
               else
                  DO i1=1,l-1
                     eq(i4) = eq(i4) - gqwl*vp*eqvd(i1)*hp(i)
                     eq(i8) = eq(i8) - gqwl*wp*eqvd(i1)*hp(i)
                  END DO ! i1
               end if ! ( ielrot(i) == 1 )
            end if ! ( l > 1 )

               ! Foreshortening related load vector

               !gsb      if(ielrot(l)==1)then
               !gsb        eq(i4) = eq(i4) + 2.0*gqwml*a3*eqvl(n)*h(i)
               !gsb        eq(i8) = eq(i8) - 2.0*gqwml*a2*eqvl(n)*h(i)
               !gsb        eq(i4) = eq(i4) - gqwl*vp*eqvdl(n)*hp(i)
               !gsb     &                  + gqwl*vp*eqvd1l(n)*hp(i)
               !gsb        eq(i8) = eq(i8) - gqwl*wp*eqvdl(n)*hp(i)
               !gsb     &                  + gqwl*wp*eqvd1l(n)*hp(i)
               !gsb      else
            eq(i4) = eq(i4)+2.0*gqwml*eqvl(n)*h(i) - gqwl*vp*eqvdl(n)*hp(i)
            eq(i8) = eq(i8)-gqwl*wp*eqvdl(n)*hp(i)
               !gsb      end if

               !sv-------------------------------------------------------------------

            if ( i < 4 )  then
               eq(i12) =  eq(i12)- gqwl *(hf(i)*unlf + hfp(i)*unlfp)
            end if ! ( i < 4 )

         END DO ! i

      end if ! ( indrns == 1 )


         ! Non-linear stiffness matrix terms

      duupwp  =  aeithp*vs
      duupvs  =  eaeast*ph + aeithp*wp
      duupws  = -eaeact*ph
      duupf   =  eac(l)*ea(l)*(vs*st - ws*ct)
      duupfp  =  aei*php

      duvsvs  = -deis2t*ph
      duvswp  =  aeithp*up + gjeb1t*php
      duvsws  =  deic2t*ph
      duvsf   = -dei*(vs*s2t - ws*c2t) + eaeast*up
      duvsfp  =  gjeb1t*wp

      duwpfp  =  gjeb1t*vs
      duwsws  =  deis2t*ph

      duwsf   =  dei*(vs*c2t + ws*s2t) - eaeact*up

      dufpfp  =  aei*up

      do i=1,4

         i4  = i + 4
         i8  = i + 8
         i12 = i + 12

         DO j=1,4

            j4  =  j + 4
            j8  =  j + 8

            huphp = hup(i)*hp (j)
            huphs = hup(i)*hs (j)
            hshup = hs (i)*hup(j)
            hphup = hp (i)*hup(j)
            hshp  = hs (i)*hp (j)
            hshs  = hs (i)*hs (j)
            hphp  = hp (i)*hp (j)
            hphs  = hp (i)*hs (j)

            dfx(i,j4)  = dfx(i,j4) - gqwl *(huphs*duupvs)
            dfx(i,j8)  = dfx(i,j8) - gqwl *(huphs*duupws + huphp*duupwp)
            dfx(j4,i)  = dfx(i,j4)
            dfx(j8,i)  = dfx(i,j8)

            dfx(i4,j4) = dfx(i4,j4) - gqwl *(hshs*duvsvs)
            dfx(i4,j8) = dfx(i4,j8) - gqwl *(hshp*duvswp + hshs*duvsws)
               ! jm, ecs: V1.1: this term deleted from dfx(i4,j8):
               !     &               - gqwl *(hphs*duvpws)
            dfx(j8,i4) = dfx(i4,j8)
            dfx(i8,j8) = dfx(i8,j8) - gqwl *(hshs*duwsws)

               !sv--------------------------------------------------------------------

            if ( indrns /= -2 )  then

                  ! 0 to x double integral term

               if ( l < nselt )  then
                     !gsb       if(ielrot(l)==1) then
                     !gsb       do 111 i1 = mselt, l+1, -1
                     !gsb        dfxc (i1,i4,j4)=dfxc(i1,i4,j4)+gqwml*2.0*h(i)*hpvpd(i1,j)*a3
                     !gsb        dfxc (i1,i8,j4)=dfxc(i1,i8,j4)-gqwml*2.0*h(i)*hpvpd(i1,j)*a2
                     !gsb        dfxcd(i1,i4,j4)=dfxcd(i1,i4,j4)+gqwml*2.0*h(i)*hpvp(i1,j)*a3
                     !gsb        dfxcd(i1,i8,j4)=dfxcd(i1,i8,j4)-gqwml*2.0*h(i)*hpvp(i1,j)*a2
                     !gsbc
                     !gsb        dfxc (i1,i4,j8) = dfxc (i1,i4,j8)+gqwml*2.0*h(i)*hpwpd(i1,j)*a3
                     !gsb        dfxc (i1,i8,j8) = dfxc (i1,i8,j8)-gqwml*2.0*h(i)*hpwpd(i1,j)*a2
                     !gsb        dfxcd(i1,i4,j8) = dfxcd(i1,i4,j8)+gqwml*2.0*h(i)*hpwp(i1,j)*a3
                     !gsb        dfxcd(i1,i8,j8) = dfxcd(i1,i8,j8)-gqwml*2.0*h(i)*hpwp(i1,j)*a2
                     !gsb 111   continue
                     !gsb      else
                  dfxc (l+1:nselt,i4,j4) = dfxc (l+1:nselt,i4,j4) + gqwml*2.0*h(i)*hpvpd(l+1:nselt,j)
                  dfxcd(l+1:nselt,i4,j4) = dfxcd(l+1:nselt,i4,j4) + gqwml*2.0*h(i)*hpvp (l+1:nselt,j)

                  dfxc (l+1:nselt,i4,j8) = dfxc (l+1:nselt,i4,j8) + gqwml*2.0*h(i)*hpwpd(l+1:nselt,j)
                  dfxcd(l+1:nselt,i4,j8) = dfxcd(l+1:nselt,i4,j8) + gqwml*2.0*h(i)*hpwp (l+1:nselt,j)
                     !gsb       end if
               end if ! ( l < nselt )

                  !gsb      if(ielrot(l)==1)then
                  !gsb        dfxc (l,i4,j4)=dfxc(l,i4,j4)+gqwml*2.0*h(i)*hpvpdl(n,j)*a3
                  !gsb        dfxc (l,i8,j4)=dfxc(l,i8,j4)-gqwml*2.0*h(i)*hpvpdl(n,j)*a2
                  !gsb        dfxcd(l,i4,j4)=dfxcd(l,i4,j4)+gqwml*2.0*h(i)*hpvpl(n,j)*a3
                  !gsb        dfxcd(l,i8,j4)=dfxcd(l,i8,j4)-gqwml*2.0*h(i)*hpvpl(n,j)*a2
                  !gsbc
                  !gsb        dfxc (l,i4,j8) = dfxc (l,i4,j8)+gqwml*2.0*h(i)*hpwpdl(n,j)*a3
                  !gsb        dfxc (l,i8,j8) = dfxc (l,i8,j8)-gqwml*2.0*h(i)*hpwpdl(n,j)*a2
                  !gsb        dfxcd(l,i4,j8) = dfxcd(l,i4,j8)+gqwml*2.0*h(i)*hpwpl(n,j)*a3
                  !gsb        dfxcd(l,i8,j8) = dfxcd(l,i8,j8)-gqwml*2.0*h(i)*hpwpl(n,j)*a2
                  !gsb      else
               dfxc (l,i4,j4) = dfxc (l,i4,j4) + gqwml*2.0*h(i)*hpvpdl(n,j)
               dfxcd(l,i4,j4) = dfxcd(l,i4,j4) + gqwml*2.0*h(i)*hpvpl (n,j)
               dfxc (l,i4,j8) = dfxc (l,i4,j8) + gqwml*2.0*h(i)*hpwpdl(n,j)
               dfxcd(l,i4,j8) = dfxcd(l,i4,j8) + gqwml*2.0*h(i)*hpwpl (n,j)
                  !gsb      end if

                  !     x to 1 double integral term

               if ( l > 1 )  then

                     !gsb      if(ielrot(l)==1) then
                     !gsb       do 121 i1 = 1,l-1
                     !gsb        dfxc (l,i4,j4) = dfxc (l ,i4,j4) - gqwl*hp(i)*hp(j)*eqvd(i1)
                     !gsb     &                                     + gqwl*hp(i)*hp(j)*eqvd1(i1)
                     !gsb        dfxc (l,i8,j8) = dfxc (l ,i8,j8) - gqwl*hp(i)*hp(j)*eqvd(i1)
                     !gsb     &                                     + gqwl*hp(i)*hp(j)*eqvd1(i1)
                     !gsb 121   continue
                     !gsbc
                     !gsb       do 122 i1 = 1,l-1
                     !gsb        dfxcd(i1,i4,j4) = dfxcd(i1 ,i4,j4) - gqwl*vp* hp(i)*h2m(i1,j)
                     !gsb     &                                     + gqwl*vp* hp(i)*h2m1(i1,j)
                     !gsb        dfxcd(i1,i8,j4) = dfxcd(i1 ,i8,j4) - gqwl*hp(i)*wp*h2m(i1,j)
                     !gsb     &                                     + gqwl*hp(i)*wp*h2m1(i1,j)
                     !gsb 122   continue
                     !gsb      else


                  DO i1=1,l-1
                     dfxc (l ,i4,j4) = dfxc (l  ,i4,j4) - gqwl*hp(i)*hp(j)*eqvd(i1)
                     dfxc (l ,i8,j8) = dfxc (l  ,i8,j8) - gqwl*hp(i)*hp(j)*eqvd(i1)
                     dfxcd(i1,i4,j4) = dfxcd(i1 ,i4,j4) - gqwl*vp*hp(i)*h2m(i1,j)
                     dfxcd(i1,i8,j4) = dfxcd(i1 ,i8,j4) - gqwl*wp*hp(i)*h2m(i1,j)
                  END DO ! i1

                     !gsb      end if

               end if ! ( l > 1 )

                  !gsb      if(ielrot(l)==1)then
                  !gsb        dfxc (l ,i4,j4) = dfxc (l ,i4,j4) - gqwl*hp(i)*hp(j)*eqvdl(n)
                  !gsb     &                                    + gqwl*hp(i)*hp(j)*eqvd1l(n)
                  !gsb        dfxc (l ,i8,j8) = dfxc (l ,i8,j8) - gqwl*hp(i)*hp(j)*eqvdl(n)
                  !gsb     &                                    + gqwl*hp(i)*hp(j)*eqvd1l(n)
                  !gsb        dfxcd(l,i4,j4) = dfxcd(l ,i4,j4) - gqwl*vp* hp(i)*h2ml(n,j)
                  !gsb     &                                   + gqwl*vp* hp(i)*h2m1l(n,j)
                  !gsb        dfxcd(l,i8,j4) = dfxcd(l ,i8,j4) - gqwl*hp(i)*wp*h2ml(n,j)
                  !gsb     &                                   + gqwl*hp(i)*wp*h2m1l(n,j)
                  !gsb      else

               dfxc (l,i4,j4) = dfxc (l,i4,j4)-gqwl*hp(i)*hp(j)*eqvdl(n)
               dfxc (l,i8,j8) = dfxc (l,i8,j8)-gqwl*hp(i)*hp(j)*eqvdl(n)
               dfxcd(l,i4,j4) = dfxcd(l,i4,j4)-gqwl*vp*hp(i)*h2ml(n,j)
               dfxcd(l,i8,j4) = dfxcd(l,i8,j4)-gqwl*wp*hp(i)*h2ml(n,j)
         !gsb      end if

            end if ! ( indrns /= -2 )

         !sv------------------------------------------------------------------

         END DO ! j

         DO j=1,3

            j12  =  j + 12
            hshf   =  hs(i) *hf(j)
            huphf  =  hup(i)*hf(j)
            huphfp =  hup(i)*hfp(j)
            hphfp  =  hp(i) *hfp(j)
            hshfp  =  hs(i) *hfp(j)

            dfx(i,j12) = dfx(i,j12) - gqwl*(huphf*duupf + huphfp*duupfp)
            dfx(i4,j12) = dfx(i4,j12) - gqwl*(hshf*duvsf + hshfp*duvsfp)
            dfx(i8,j12) = dfx(i8,j12) - gqwl*(hshf*duwsf + hphfp*duwpfp)
            dfx(j12,i)  = dfx(i,j12)
            dfx(j12,i4) = dfx(i4,j12)
            dfx(j12,i8) = dfx(i8,j12)

            if ( i < 4 )  then
               hfphfp = hfp(i)*hfp(j)
               dfx(i12,j12) = dfx(i12,j12) - gqwl*(hfphfp*dufpfp)
            end if ! ( i < 4 )

         end do ! j

      end do ! i

   END IF ! ( ( indnl /= 0 ) .AND. ( indrns /= -1 ) )

END DO ! n


RETURN
END subroutine Struct ! ( euj1, ek, ec, em, eq, dfx, shi, xbi, eli, l, axfi, indrns, indnl, tlam, gu, gud, dfxc, dfxcd, elv )
!=======================================================================
