!*****************************************************
!*  Sorts an array RA of length N in ascending order *
!*                by the Heapsort method             *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*	    N	  size of table RA                   *
!*          RA	  table to be sorted                 *
!* OUTPUT:                                           *
!*	    RA    table sorted in ascending order    *
!*                                                   *
!* NOTE: The Heapsort method is a N Log2 N routine,  *
!*       and can be used for very large arrays.      *
!*****************************************************
SUBROUTINE HPSORT(N,RA,IND)

  INTEGER        :: N
  double precision  RA(N)


  INTEGER        :: L, IR, J, TIND, IND(N)
  double precision  RRA

  L=N/2+1
  IR=N

  !The index L will is decremented from its initial value during the
  !"hiring" (heap creation) phase. Once it reaches 1, the index IR
  !is be decremented from its initial value down to 1 during the
  !"retirement-and-promotion" (heap selection) phase.

10 continue
  if(L > 1)then
    L=L-1
    RRA=RA(L)
    TIND = IND(L)
  else
    RRA=RA(IR)
    TIND = IND(IR)
    RA(IR)=RA(1)
    IND(IR) = IND(1)
    IR=IR-1
    if(IR.eq.1)then
      RA(1)=RRA
      IND(1) = TIND
      return
    end if
  end if
  I=L
  J=L+L
20 if(J.le.IR)then
  if(J < IR)then
    if(RA(J) < RA(J+1))  J=J+1
  end if
  if(RRA < RA(J))then
    RA(I)=RA(J)
    IND(I) = IND(J)
    I=J; J=J+J
  else
    J=IR+1
  end if
  goto 20
  end if
  RA(I)=RRA
  IND(I) = TIND
  goto 10
END

