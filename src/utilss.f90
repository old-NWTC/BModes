!=======================================================================
SUBROUTINE asbgmk ( gm, gk, em, ek, indeg )

USE Param
USE Precision
USE Swept

IMPLICIT NONE


   ! Argument declarations.

REAL(ReKi), INTENT(IN)        :: ek       (nedof,nedof)
REAL(ReKi), INTENT(IN)        :: em       (nedof,nedof)
REAL(ReKi), INTENT(INOUT)     :: gk       (ngd,ngd)
REAL(ReKi), INTENT(INOUT)     :: gm       (ngd,ngd)

INTEGER, INTENT(IN)           :: indeg    (nedof)


   ! Local declarations.

INTEGER                       :: i
INTEGER                       :: ii
INTEGER                       :: j
INTEGER                       :: jj



   !sweep

   !gsb      if (ielrot(l) == 1) then
   !gsbc
   !gsb        do 14 i = 1, nedof
   !gsb          do 15 j = 1,nedof
   !gsb            ekt(i,j)= 0.
   !gsb            emt(i,j)= 0.
   !gsb            do 16 k = 1,nedof
   !gsb              ekt(i,j) = ekt(i,j) + ek(i,k) * tlam2(k,j)
   !gsb              emt(i,j) = emt(i,j) + em(i,k) * tlam2(k,j)
   !gsb 16         continue
   !gsb 15       continue
   !gsb 14     continue
   !gsbc
   !gsbcsweep
   !gsb        do 17 i = 1, nedof
   !gsb          do 18 j = 1, nedof
   !gsb            tekt(i,j) =0.
   !gsb            temt(i,j) = 0.
   !gsb            do 19 k = 1, nedof
   !gsb              tekt(i,j) = tekt(i,j) + tlam2(k,i)*ekt(k,j)
   !gsb              temt(i,j) = temt(i,j) + tlam2(k,i)*emt(k,j)
   !gsb 19         continue
   !gsb 18       continue
   !gsb 17     continue
   !gsbc
   !gsb      end if

   !sweep


do j=1,nedof

   jj = indeg(j)

   if ( jj >= 1 )  then

      do i=1,nedof

         ii = indeg(i)

         if ( ii >= 1 )  then
            gm(ii,jj) = gm(ii,jj) + em(i,j)
            gk(ii,jj) = gk(ii,jj) + ek(i,j)
         end if

      end do ! i

   end if

end do ! j


RETURN
END SUBROUTINE asbgmk ! ( gm, gk, em, ek, indeg )
!=======================================================================
SUBROUTINE evfrgv (gq, eq, indeg)

   !     this subroutine obtains the element vector for element ne from the global vector

   !     nedof = total no. of degrees of freedom for the element

USE Param
USE Precision

IMPLICIT NONE


   ! Argument declarations.

REAL(ReKi), INTENT(OUT)    :: eq       (nedof)                          ! element vector
REAL(ReKi), INTENT(IN)     :: gq       (ngd)                            ! global vector

INTEGER, INTENT(IN)        :: indeg    (nedof)                          ! element dof to global dof connector


   ! Local declarations.

INTEGER                    :: i                                         ! Generic DO counter



eq(:) = 0.0

do i=1,nedof

   if ( indeg(i) >= 1 )  eq(i) = gq(indeg(i))

end do ! i


RETURN
END
!=======================================================================
SUBROUTINE evfrsw (gq, eq, indeg, tlam2)

   !     this subroutine obtains the element vector for element ne from
   !     the global vector for spatial fem (sweep and droop)

   !     nedof = total no. of degrees of freedom for the element

USE Param
USE Precision

IMPLICIT NONE


   ! Argument declarations.

REAL(ReKi), INTENT(OUT)       :: eq       (nedof)                       ! element vector
REAL(ReKi), INTENT(IN)        :: gq       (ngd)                         ! global vector
REAL(ReKi), INTENT(IN)        :: tlam2    (nedof,nedof)                 !

INTEGER, INTENT(IN)           :: indeg    (nedof)                       ! element dof to global dof connector


   ! Local declarations.

REAL(ReKi)                    :: eqs      (nedof)                       !

INTEGER                       :: i                                      ! Generic DO counter



eqs(:) = 0.0

do i=1,nedof
  if ( indeg(i) >= 1 )  eqs(i) = gq(indeg(i))
end do ! i


   !sweep

do i=1,nedof
   eq(i) = SUM( eqs(:)*tlam2(i,:) )
end do ! i


RETURN
END SUBROUTINE evfrsw ! (gq, eq, indeg, tlam2)
!=======================================================================
SUBROUTINE gqi5pt ( gqp, gqw, ngauss )

   !
   !     gauss quadrature sampling points and weighting factors for the interval  (0,1)

USE Param
USE Precision

IMPLICIT NONE


   ! Argument declarations.

INTEGER, INTENT(IN)           :: ngauss                                 ! Number of Gauss points.

REAL(ReKi), INTENT(OUT)       :: gqp      (ngauss)
REAL(ReKi), INTENT(OUT)       :: gqw      (ngauss)


   ! Local declarations.

REAL(ReKi)                    :: gau      (6) = (/ 0.906179845938664, 0.538469310105683, 0.0 , 0.236926885056189, &
                                                   0.478628670499366, 0.568888888888889 /)

INTEGER                       :: i                                      ! Generic DO counter
INTEGER                       :: i3
INTEGER                       :: j



do i=1,2
   j      = 6 - i
   i3     = i + 3
   gqp(i) = 0.5*( 1.0 - gau(i) )
   gqp(j) = 0.5*( 1.0 + gau(i) )
   gqw(i) = 0.5*gau(i3)
   gqw(j) = gqw(i)
end do ! i

gqp(3) = 0.5*( 1.0 + gau(3) )
gqw(3) = 0.5*gau(6)


RETURN
END Subroutine gqi5pt ! ( gqp, gqw, ngauss )
!=======================================================================
SUBROUTINE gqi6pt ( gqp, gqw, ngauss )

   ! Gauss quadrature sampling points and weighting factors for the interval  (0,1)

USE Param
USE Precision

IMPLICIT NONE


   ! Argument declarations.

INTEGER, INTENT(IN)           :: ngauss                                 ! Number of Gauss points.

REAL(ReKi), INTENT(OUT)       :: gqp      (ngauss)
REAL(ReKi), INTENT(OUT)       :: gqw      (ngauss)


   ! Local declarations.

REAL(ReKi)                    :: gau      (6) = (/ 0.93246951420315, 0.66120938646627, 0.23861918608320, 0.17132449237917, &
                                                   0.36076157304814, 0.46791393457269 /)


INTEGER                       :: i                                      ! Generic DO counter
INTEGER                       :: i3
INTEGER                       :: j



do i=1,3
   j      = 7 - i
   i3     = i + 3
   gqp(i) = 0.5*( 1.0 - gau(i) )
   gqp(j) = 0.5*( 1.0 + gau(i) )
   gqw(i) = 0.5*gau(i3)
   gqw(j) = gqw(i)
end do ! i


RETURN
END SUBROUTINE gqi6pt ! ( gqp, gqw, ngauss )
!=======================================================================
SUBROUTINE NormEv ( evr, ngd )


USE CNPar
USE Precision

IMPLICIT NONE



   ! Argument declarations.

INTEGER, INTENT(IN)           :: ngd

REAL(ReKi), INTENT(INOUT)     :: evr(ngd,1)


   ! Local declarations.

REAL(ReKi)                    :: SumEvrJ

INTEGER                       :: j                                      !
INTEGER                       :: jz                                     ! Generic DO counter.



do jz=1,ngd

   j = ngd + 1 - jz

   SumEvrJ  = SUM( evr(:,j)**2 )
   evr(:,j) = evr(:,j)/sqrt( SumEvrJ )

end do ! jz


return
end SUBROUTINE NormEv ! ( evr, ngd )
!=======================================================================
