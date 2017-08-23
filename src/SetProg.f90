!=======================================================================
SUBROUTINE SetProg


   ! This routine sets the version number.  By doing it this way instead
   !   of the old way of initializing it in a module, we will no longer
   !   have to recompile everything every time we change versions.


USE                                NWTC_Library

IMPLICIT                           NONE


   ! Local Variables:

CHARACTER(27)                   :: Version = 'v1.03.01, 25-Sept-2007'             ! String containing the current version.



ProgName = 'BModes'

IF ( ReKi == 4 )  THEN     ! Single precision
   ProgVer = ' ('//TRIM( Version )//')'
ELSEIF ( ReKi == 8 )  THEN     ! Double precision
   ProgVer = ' ('//TRIM( Version )//', compiled using double precision)'
ELSE                       ! Unknown precision - it should be impossible to compile using a KIND that is not 4 or 8, but I'll put this check here just in case.
   ProgVer = ' ('//TRIM( Version )//', compiled using '//TRIM( Int2LStr( ReKi ) )//'-byte precision)'
end if


RETURN
END SUBROUTINE SetProg
