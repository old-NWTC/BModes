      !=======================================================================
      SUBROUTINE jacobi ( n, a, b, evl, evr )

         ! Initialize eigenvalue and eigenvector matrices.

      USE Eig
      USE NWTC_Library
      USE Precision

      IMPLICIT NONE


         ! Argument declarations.

      REAL(ReKi), INTENT(INOUT)    :: a   (n,n)
      REAL(ReKi), INTENT(INOUT)    :: b   (n,n)
      REAL(ReKi), INTENT(  OUT)    :: evl (n)
      REAL(ReKi), INTENT(  OUT)    :: evr (n,n)


         ! Local declarations.

      REAL(ReKi)                   :: ab
      REAL(ReKi)                   :: aj
      REAL(ReKi)                   :: ajj
      REAL(ReKi)                   :: ak
      REAL(ReKi)                   :: akk
      REAL(ReKi)                   :: bb
      REAL(ReKi)                   :: bj
      REAL(ReKi)                   :: bk
      REAL(ReKi)                   :: ca
      REAL(ReKi)                   :: cg
      REAL(ReKi)                   :: check
      REAL(ReKi)                   :: d1
      REAL(ReKi)                   :: d2
      REAL(ReKi)                   :: den
      REAL(ReKi)                   :: dif
      REAL(ReKi)                   :: eps
      REAL(ReKi)                   :: epsa
      REAL(ReKi)                   :: epsb
      REAL(ReKi)                   :: eptola
      REAL(ReKi)                   :: eptolb
      REAL(ReKi), PARAMETER        :: rtol    = 1.0e-8                  ! Convergence tolerance.
      REAL(ReKi)                   :: sqch
      REAL(ReKi)                   :: tol
      REAL(ReKi)                   :: xj
      REAL(ReKi)                   :: xk

      INTEGER                      :: i
      INTEGER                      :: j
      INTEGER                      :: jj
      INTEGER                      :: jm1
      INTEGER                      :: jp1
      INTEGER                      :: k
      INTEGER                      :: km1
      INTEGER                      :: kp1
      INTEGER                      :: n
      INTEGER                      :: nr
      INTEGER, PARAMETER           :: nsmax   = 10                      ! Maximum number of sweeps allowed.
      INTEGER                      :: nsweep

      do i=1,n

         if ( ( a(i,i) <= 0.0 ) .OR. ( b(i,i) <= 0.0 ) )  then
         !nrel
            write(*,*) i
            CALL ProgAbort ( 'Matrices not positive definite in'
     &                 //' jacobi() at check #1.' )
         end if

         d(i)    = a(i,i)/b(i,i)
         eigv(i) = d(i)

      end do ! i

      x(:,:) = 0.0

      do i=1,n
         x(i,i) = 1.0
      end do ! i

      if(n==1) return

         !     initialize sweep counter and begin iteration

      nsweep = 0
      nr = n-1
   40 nsweep = nsweep+1

         !     check if present off-diagonal element is large enough to require z

      eps = (.01**nsweep)**2

      do j=1,nr
         jj = j+1
         do k=jj,n
            eptola = (a(j,k)*a(j,k))/(a(j,j)*a(k,k))
            eptolb = (b(j,k)*b(j,k))/(b(j,j)*b(k,k))
            if((eptola<eps).and.(eptolb<eps)) CYCLE

         !           if zeroing is required, calculate the rotation matrix element ca a

            akk = a(k,k)*b(j,k)-b(k,k)*a(j,k)
            ajj = a(j,j)*b(j,k)-b(j,j)*a(j,k)
            ab = a(j,j)*b(k,k)-a(k,k)*b(j,j)

            check = (ab*ab + 4.*akk*ajj)/4.
            if(check) 50,60,60
         !nrel
            write(*,*) check
   50       CALL ProgAbort ( 'Matrices not positive definite in'
     &                 //' jacobi() at check #2.' )
   60       sqch = sqrt( check )
            d1 = ab/2.+sqch
            d2 = ab/2.-sqch
            den = d1
            if ( abs(d2) > abs(d1) )  den = d2
            if(den)80,70,80
   70       ca = 0.
            cg = -a(j,k)/a(k,k)
            go to 90
   80       ca = akk/den
            cg = -ajj/den

         !           perform the generalized rotation to zero the present off-diagonal

   90       if(n-2)100,190,100
  100       jp1 = j+1
            jm1 = j-1
            kp1 = k+1
            km1 = k-1
            if(jm1-1)130,110,110

  110       do i=1,jm1
               aj     = a(i,j)
               bj     = b(i,j)
               ak     = a(i,k)
               bk     = b(i,k)
               a(i,j) = aj + cg*ak
               b(i,j) = bj + cg*bk
               a(i,k) = ak + ca*aj
               b(i,k) = bk + ca*bj
            end do ! i

  130       if(kp1-n)140,140,160

  140       do i=kp1,n
               aj     = a(j,i)
               bj     = b(j,i)
               ak     = a(k,i)
               bk     = b(k,i)
               a(j,i) = aj + cg*ak
               b(j,i) = bj + cg*bk
               a(k,i) = ak + ca*aj
               b(k,i) = bk + ca*bj
            end do ! i

  160       if(jp1-km1) 170,170,190

  170       do i=jp1,km1
               aj     = a(j,i)
               bj     = b(j,i)
               ak     = a(i,k)
               bk     = b(i,k)
               a(j,i) = aj + cg*ak
               b(j,i) = bj + cg*bk
               a(i,k) = ak + ca*aj
               b(i,k) = bk + ca*bj
            end do ! i

  190       ak = a(k,k)
            bk = b(k,k)
            a(k,k)  =  ak+2.*ca*a(j,k)+ca*ca*a(j,j)
            b(k,k) = bk+2.*ca*b(j,k)+ca*ca*b(j,j)
            a(j,j) = a(j,j)+2.*cg*a(j,k)+cg*cg*ak
            b(j,j) = b(j,j)+2.*cg*b(j,k)+cg*cg*bk
            a(j,k) = 0.
            b(j,k) = 0.

         !     update the eigenvector matrix after each rotation

            do i=1,n
               xj     = x(i,j)
               xk     = x(i,k)
               x(i,j) = xj + cg*xk
               x(i,k) = xk + ca*xj
            end do ! i


         end do ! k

      end do ! j

         ! update the eigenvalues after each sweep

      do i=1,n
         !nrel  modify later
         if ( b(i,i) <= 0.0 )  b(i,i) = 0.1e-7
         !--
         if ( ( a(i,i) <= 0.0 ) .OR. ( b(i,i) <= 0.0 ) )  then
            write(*,*) i, a(i,i), b(i,i)
            CALL ProgAbort ( 'Matrices not positive definite in'
     &                 //' jacobi() at check #3.' )
         end if

         eigv(i) = a(i,i)/b(i,i)

      end do ! i

         !     check for convergence

  230 do i=1,n
         tol = rtol*d(i)
         dif = abs( eigv(i) - d(i) )
         if ( dif > tol )  go to 280
      END DO ! i

         !     check all off-diagonal elements to see if another sweep is require

      eps = rtol**2

      do j=1,nr

         jj = j + 1

         do k=jj,n

            if ( a(j,k) < TINY( epsa ) ) then
               epsa  =  0.0
            else
               epsa = (a(j,k)/a(j,j))*(a(j,k)/a(k,k))
            end if

            if ( b(j,k) < TINY( epsb ) ) then
               epsb  =  0.0
            else
               epsb = (b(j,k)/b(j,j))*(b(j,k)/b(k,k))
            end if

            if ( ( epsa >= eps ) .OR. ( epsb >= eps ) )  go to 280
         end do ! k
      end do ! j

         !     fill out bottom triangle of resultant matrices and scale eigenvect

  255 do i=1,n
         do j=1,n
            a(j,i) = a(i,j)
            b(j,i) = b(i,j)
         end do ! j
      end do ! i

      do j=1,n

         bb = sqrt( b(j,j) )

         do k=1,n
            x(k,j) = x(k,j)/bb
         end do ! k

      end do ! j

      call ordrch (n,evl,evr)
      return

         !     update  d  matrix and start new sweep, if allowed
  280 do i=1,n
         d(i) = eigv(i)
      end do ! i
      if(nsweep<nsmax) go to 40
      go to 255


      RETURN
      END SUBROUTINE jacobi ! ( n, a, b, evl, evr )
      !=======================================================================
      SUBROUTINE ordrch (n, evl, evr )

      USE Eig
      USE Precision
      IMPLICIT NONE


         ! Argument declarations.


      REAL(ReKi), INTENT(OUT)      :: evl (n)
      REAL(ReKi), INTENT(OUT)      :: evr (n,n)

      INTEGER, INTENT(IN)          :: n


         ! Local declarations.

      INTEGER                      :: i
      INTEGER                      :: ii
      INTEGER                      :: j
      INTEGER                      :: jj
      INTEGER                      :: k
      INTEGER                      :: l
      INTEGER                      :: ll
      INTEGER                      :: n1



      do k=1,n

         l  = 1
         ll = 0

         do i=1,n

            if ( ( k /= i ) .AND. ( eigv(k) >= eigv(i) ) )  then
               l = l + 1
               go to 650
            else

               ll = ll + 1

               if ( ll /= n )  then
                  if (i-n)  600,640,600
               end if

               l = 1
            end if

  650       if ( ( l > n ) .OR. ( i /= n ) )  CYCLE

  640       jdd(l) = k

  600       continue
         end do ! i
      end do ! k

      n1 = n+1

      DO i=1,n
         jj     = n1-i
         j      = jdd(jj)
         evl(i) = eigv(j)
      END DO ! i

      do i =1,n
         do ii=1,n
            jj        = n1-ii
            j         = jdd(jj)
            evr(i,ii) = x(i,j)
         end do ! ii
      end do ! i

      return
      end SUBROUTINE ordrch ! (n, evl, evr )
      !=======================================================================
