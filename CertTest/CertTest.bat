@ECHO OFF
@ECHO.

REM **********************************
REM **  Test batch file for BModes  **
REM **********************************


REM  Set up environment variables.  You will probably have to change these.

@SET Compare=FC
@SET DateTime=DateTime.exe
@SET Editor="C:\Program Files (x86)\Notepad++\notepad++.exe"

@SET ProgName=BModes
@SET Program=..\bin\BModes_Win32.exe
@SET InExt=bmi
@SET OutExt=out


REM  Program test sequence definition:

@SET  TEST01=Test #01: Non-uniform blade (Space-delimited).
@SET  TEST02=Test #02: Non-uniform blade with tip inertia (Tab-delimited).
@SET  TEST03=Test #03: Non-uniform tower with tip mass (Space-delimited).
@SET  TEST04=Test #04: Non-uniform tower with tip mass and wires-support (Tab-delimited).

@SET  DASHES=---------------------------------------------------------------------------------------------
@SET  POUNDS=#############################################################################################

@IF EXIST CertTest.out  DEL CertTest.out

ECHO.                                                 >> CertTest.out
ECHO           ************************************** >> CertTest.out
ECHO           **  BModes Acceptance Test Results  ** >> CertTest.out
ECHO           ************************************** >> CertTest.out

ECHO.                                                                             >> CertTest.out
ECHO ############################################################################ >> CertTest.out
ECHO # Inspect this file for any differences between your results and the saved # >> CertTest.out
ECHO # results.  Any differing lines and the two lines surrounding them will be # >> CertTest.out
ECHO # listed.  The only differences should be the time stamps at the start of  # >> CertTest.out
ECHO # each file.                                                               # >> CertTest.out
ECHO #                                                                          # >> CertTest.out
ECHO # If you are running on something other than a PC, you may see differences # >> CertTest.out
ECHO # in the last significant digit of many of the numbers.                    # >> CertTest.out
ECHO ############################################################################ >> CertTest.out

ECHO.                                            >> CertTest.out
ECHO Date and time this acceptance test was run: >> CertTest.out
%DateTime%                                       >> CertTest.out
ECHO.                                            >> CertTest.out


rem *******************************************************

@echo %ProgName% %TEST01%

@SET TEST=01_nonunif_blade

rem Run %ProgName%.

%Program% Test%TEST%.%InExt% > NUL

IF ERRORLEVEL 1  GOTO ERROR

@IF NOT EXIST Test%TEST%.%OutExt%  GOTO ERROR

echo.          >> CertTest.out
echo %POUNDS%  >> CertTest.out
echo.          >> CertTest.out
echo %TEST01%  >> CertTest.out
echo %DASHES%  >> CertTest.out

%Compare% Test%TEST%.%OutExt% TestFiles\Test%TEST%.%OutExt% >> CertTest.out
TortoiseGitMerge Test%TEST%.%OutExt% TestFiles\Test%TEST%.%OutExt%


rem *******************************************************

@echo %ProgName% %TEST02%

@SET TEST=02_blade_with_tip_mass

rem Run %ProgName%.

%Program% Test%TEST%.%InExt% > NUL

IF ERRORLEVEL 1  GOTO ERROR

@IF NOT EXIST Test%TEST%.%OutExt%  GOTO ERROR

echo.          >> CertTest.out
echo %POUNDS%  >> CertTest.out
echo.          >> CertTest.out
echo %TEST02%  >> CertTest.out
echo %DASHES%  >> CertTest.out
%Compare% Test%TEST%.%OutExt% TestFiles\Test%TEST%.%OutExt% >> CertTest.out
TortoiseGitMerge Test%TEST%.%OutExt% TestFiles\Test%TEST%.%OutExt%

rem *******************************************************

@echo %ProgName% %TEST03%

@SET TEST=03_tower

rem Run %ProgName%.

%Program% Test%TEST%.%InExt% > NUL

IF ERRORLEVEL 1  GOTO ERROR

@IF NOT EXIST Test%TEST%.%OutExt%  GOTO ERROR

echo.          >> CertTest.out
echo %POUNDS%  >> CertTest.out
echo.          >> CertTest.out
echo %TEST03%  >> CertTest.out
echo %DASHES%  >> CertTest.out
%Compare% Test%TEST%.%OutExt% TestFiles\Test%TEST%.%OutExt% >> CertTest.out
TortoiseGitMerge Test%TEST%.%OutExt% TestFiles\Test%TEST%.%OutExt%

rem *******************************************************

@echo %ProgName% %TEST04%

@SET TEST=04_wires_supported_tower

rem Run %ProgName%.

%Program% Test%TEST%.%InExt% > NUL

IF ERRORLEVEL 1  GOTO ERROR

@IF NOT EXIST Test%TEST%.%OutExt%  GOTO ERROR

echo.          >> CertTest.out
echo %POUNDS%  >> CertTest.out
echo.          >> CertTest.out
echo %TEST04%  >> CertTest.out
echo %DASHES%  >> CertTest.out
%Compare% Test%TEST%.%OutExt% TestFiles\Test%TEST%.%OutExt% >> CertTest.out
TortoiseGitMerge Test%TEST%.%OutExt% TestFiles\Test%TEST%.%OutExt%

rem ******************************************************
rem  Let's look at the comparisons.

%Editor% CertTest.out
goto END

:ERROR
@echo ** An error has occured in Test #%TEST% **

:END

@SET Compare=
@SET DASHES=
@SET DateTime=
@SET Editor=
@SET InExt=
@SET WT_Perf=
@SET OutExt=
@SET POUNDS=
@SET ProgName=
@SET Program=
@SET TEST=
@SET TEST01=
@SET TEST02=
@SET TEST03=
@SET TEST04=

@echo 
@echo Processing complete.
