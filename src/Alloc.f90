MODULE Allocation

   ! This module contains routines used to allocate many arrays.

CONTAINS
!=======================================================================
   SUBROUTINE Alloc1

      ! This routine allocates many of the global arrays used in the program.

   USE Gauss
   USE NWTC_Library
   USE Param
   USE Struc
   USE Swept

   IMPLICIT NONE


      ! Local declarations.

   INTEGER                             :: Sttus



      ! Arrays in module Gauss:

   ALLOCATE ( ea(nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, ea, in Alloc1().' )
   END IF

   ALLOCATE ( eac(nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, eac, in Alloc1().' )
   END IF

   ALLOCATE ( eb1(nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, eb1, in Alloc1().' )
   END IF

   ALLOCATE ( eb2(nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, eb2, in Alloc1().' )
   END IF

   ALLOCATE ( ec1(nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, ec1, in Alloc1().' )
   END IF

   ALLOCATE ( ec2(nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, ec2, in Alloc1().' )
   END IF

   ALLOCATE ( eg(nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, eg, in Alloc1().' )
   END IF

   ALLOCATE ( eiy(nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, eiy, in Alloc1().' )
   END IF

   ALLOCATE ( eiz(nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, eiz, in Alloc1().' )
   END IF

   ALLOCATE ( elm_dist_k(nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, elm_dist_k, in Alloc1().' )
   END IF

   ALLOCATE ( elm_dist_m(nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, elm_dist_m, in Alloc1().' )
   END IF

   ALLOCATE ( gay(nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, gay, in Alloc1().' )
   END IF

   ALLOCATE ( gaz(nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, gaz, in Alloc1().' )
   END IF

   ALLOCATE ( gj(nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, gj, in Alloc1().' )
   END IF

   ALLOCATE ( rmas(nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, rmas, in Alloc1().' )
   END IF

   ALLOCATE ( skm1(nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, skm1, in Alloc1().' )
   END IF

   ALLOCATE ( skm2(nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, skm2, in Alloc1().' )
   END IF


      ! Arrays in module Gauss:

   ALLOCATE ( el(nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, el, in Alloc1().' )
   END IF

   ALLOCATE ( gqp(ngauss), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, gqp, in Alloc1().' )
   END IF

   ALLOCATE ( gqw(ngauss), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, gqw, in Alloc1().' )
   END IF

   ALLOCATE ( xb(nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, xb, in Alloc1().' )
   END IF

   ALLOCATE ( indeg(nedof,nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, xb, in Alloc1().' )
   END IF


      ! Arrays in module Swept:

   ALLOCATE ( ielrot(nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, ielrot, in Alloc1().' )
   END IF


   RETURN
   END SUBROUTINE Alloc1
!=======================================================================
   SUBROUTINE Alloc2

      ! This routine allocates the global arrays in module Struc.

   USE NWTC_Library
   USE Struc

   IMPLICIT NONE


      ! Local declarations.

   INTEGER                             :: Sttus



      ! Arrays in module Struc:

   ALLOCATE ( amass_den(n_secs), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, amass_den, in Alloc2().' )
   END IF

   ALLOCATE ( axial_stff(n_secs), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, axial_stff, in Alloc2().' )
   END IF

   ALLOCATE ( cg_offst(n_secs), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, cg_offst, in Alloc2().' )
   END IF

   ALLOCATE ( edge_stff(n_secs), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, edge_stff, in Alloc2().' )
   END IF

   ALLOCATE ( flp_stff(n_secs), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, flp_stff, in Alloc2().' )
   END IF

   ALLOCATE ( sc_offst(n_secs), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, sc_offst, in Alloc2().' )
   END IF

   ALLOCATE ( sec_loc(n_secs), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, sec_loc, in Alloc2().' )
   END IF

   ALLOCATE ( sq_km1(n_secs), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, sq_km1, in Alloc2().' )
   END IF

   ALLOCATE ( sq_km2(n_secs), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, sq_km2, in Alloc2().' )
   END IF

   ALLOCATE ( str_tw(n_secs), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, str_tw, in Alloc2().' )
   END IF

   ALLOCATE ( tc_offst(n_secs), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, tc_offst, in Alloc2().' )
   END IF

   ALLOCATE ( tor_stff(n_secs), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, tor_stff, in Alloc2().' )
   END IF


   RETURN
   END SUBROUTINE Alloc2
!=======================================================================
   SUBROUTINE Alloc3

      ! This routine allocates many of the global arrays used in the program.

   USE Basic
   USE DisBd
   USE Eig
   USE NWTC_Library
   USE Param

   IMPLICIT NONE


      ! Local declarations.

  INTEGER                             :: Sttus



      ! Arrays in module Basic:

   ALLOCATE ( cfe(nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, cfe, in Alloc3().' )
   END IF



      ! Arrays in module DisBd:

   ALLOCATE ( eiyfac(nblade,nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, eiyfac, in Alloc3().' )
   END IF

   ALLOCATE ( eiytmp(nblade,nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, eiytmp, in Alloc3().' )
   END IF

   ALLOCATE ( eizfac(nblade,nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, eizfac, in Alloc3().' )
   END IF

   ALLOCATE ( eiztmp(nblade,nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, eiztmp, in Alloc3().' )
   END IF

   ALLOCATE ( gjfac(nblade,nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, gjfac, in Alloc3().' )
   END IF

   ALLOCATE ( gjtmp(nblade,nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, gjtmp, in Alloc3().' )
   END IF

   ALLOCATE ( rmasfc(nblade,nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, rmasfc, in Alloc3().' )
   END IF

   ALLOCATE ( rmastm(nblade,nselt), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, rmastm, in Alloc3().' )
   END IF


      ! Arrays in module Eig:

   ALLOCATE ( d(ngd), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, d, in Alloc3().' )
   END IF

   ALLOCATE ( eigv(ngd), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, eigv, in Alloc3().' )
   END IF

   ALLOCATE ( jdd(ngd), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, jdd, in Alloc3().' )
   END IF

   ALLOCATE ( x(ngd,ngd), STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort( 'Cannot allocate array, x, in Alloc3().' )
   END IF


   RETURN
   END SUBROUTINE Alloc3
!=======================================================================

END MODULE Allocation