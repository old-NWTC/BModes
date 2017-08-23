   !=======================================================================
   !                                                                      c
   !                   BModes v3.00.00 (OFFSHORE)                         c
   !                                                                      c
   !          ( Property of National Renewable Energy Laboratory )        c
   !                                                                      c
   !                   Developed by: Gunjit S. Bir                        c
   !                        (303) 384-6953                                c
   !                                                                      c
   !        Partially Coverted to FORTRAN 2003 by: Marshall Buhl          c
   !                                                                      c
   !     Note: This code is still in the development/validation stage.    c
   !     BModes is derived from UMARC, a code deveoped earlier            c
   !     for rotrcraft.  Remanants of the old code are still there; these c
   !     eventually need to cleaned out. Ignore current comments and      c
   !     variable definitions.  Comments pertaining to upgrades and       c
   !     comprehensive input checks will be introduced later.             c
   !                                                                      c
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !                                                                      c
   !    ***************************************************************   c
   !    *                                                             *   c
   !    *    Finite element program to compute modes for a rotating   *   c
   !    *    blade or a tower (with or without tension wires)         *   c
   !    *                                                             *   c
   !    ***************************************************************   c
   !                                                                      c
   !                                                                      c
   !     main program ===>                                         exec   c
   !                                                                      c
   !     index of subroutines:                                            c

   !     --- stand-alone files -------------------------------------------c

   !             build blade model and compute blade vib chs     bldvib   c
   !             structural matrix calculation                   struct   c
   !                                                                      c
   !     --- BModes file contains ----------------------------------------c
   !                                                                      c
   !             exec                                            exec     c
   !                   (calls all other subroutines)                      c
   !                                                                      c
   !             blade sweep transformation matrices             sweep    c
   !            element matrix/vector modifications                       c
   !                   impose bmr kinematics                     bmrmod   c
   !                   linearly combine the rows & columns       modmat   c

   !     --- utils file contains -----------------------------------------c
   !                                                                      c
   !            assembly subroutines                                      c
   !                   space:    (global mass & stiffness)       asbgmk   c
   !                                                                      c
   !            utility subroutines                                       c
   !                   change "normalized" vector to geometric   cnvtgv   c
   !                   element vector from global vector         evfrgv   c
   !                   gauss pts/wts -- 5 point                  gqi5pt   c
   !                   gauss pts/wts -- 6 point                  gqi6pt   c
   !                   initialize routine                        init     c
   !                   normalize eigenvectors                    normev   c
   !                   time connectivity vector                  tcntvt   c
   !                   matrix multiplication                     mult     c
   !                                                                      c
   !                                                                      c
   !     --- math file contains ------------------------------------------c
   !                                                                      c
   !            eigen analysis subroutines                                c
   !                    real non-symmetric matrices &            jacobi   c
   !                    complex eigenvalues --->                 ordrch   c
   !                                                                      c
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !                                                                      c
   !           d e s c r i p t i o n  of  v a r i a b l e s               c
   !                                                                      c

   !     UnIn     = I/O unit for the primary input file.
   !     iartic = 0; hingeless or f-f blade
   !            = 1; flap hinge articulated blade
   !            = 2; flap then lag hinge articulated blade
   !            = 3; flap then lag hinge artic. blade with del3 effect
   !            = 4; flap and lag hinge coincident articulated blade

   !     nselt  = no. of spatial elements
   !     ntelt  = no. of time elements
   !     nblade = no. of blades
   !     nseg   = no. of blade segments
   !     nconf  = 0; ......
   !            = 1; ......
   !            = 2; bmr
   !     npin   = 0; no shear lag pin
   !            = 1; shear lag pin exists
   !     ttk    = pitch link stiffness, nondimensionaized
   !     tta    = pith link arm length / rotor radius
   !     ttz    = vertical offset of the lag shear pin / radius
   !     ttp    = pitch link bend / radius
   !     ttx0   = horizontal offset of the lag shear pin / radius
   !     infus  = 0 : hub motion not included
   !              1 : hub motion included

   !     nflap  = no. of flap modes used in analysis
   !     nlag   = no. of lag modes used in analysis
   !     ntorsn = no. of torsion modes used in analysis
   !     naxial = no. of axial modes used in response analysis
   !       *** note: axial modes not used in stability analysis ***
   !       ***    or determination of vehicle jacobian matrix   ***
   !     nndof  = no. of nodal dof
   !     nintn  = no. of internal nodes
   !     ncon   = no. of zero boundary conditions at root
   !     negf   = no. of modes used in response analysis
   !     negpt  = no. of modes used in stability analysis
   !     ndt    = no. of total dof
   !     nedof  = no. of element dof
   !     ngd    = no. of bounded dof
   !     nsect  = order of time integration within time elements

   !     negpt2,negf2,nfxd,nfxd2 are derived dimension parameters

   !     UnOu     = write control
   !     icpltr = 0; uncoupled trim using rigid trim controls (vtrim)
   !            = 1; coupled trim/initial controls from rigid trim (vtrim)
   !            = 2; uncoupled trim using inputted controls
   !                 but inflow is calculated from vtrim
   !            = 3; coupled trim/initial controls inputted
   !            = 4; uncoupled trim/input + 5% increase in alph
   !            = 5; uncoupled trim/input + 5% increase in phis
   !            = 6; uncoupled trim/input + 5% increase in th0
   !            = 7; uncoupled trim/input + 5% increase in th1c
   !            = 8; uncoupled trim/input + 5% increase in th1s
   !            = 9; uncoupled trim/input + 5% increase in t0tail
   !            =10; uncoupled hover trim (vtrim controls)
   !            =11; coupled hover trim (thrust - collective trim)
   !            =12; coupled wind tunnel trim (zero tip flapping)
   !     itrmcv = 1; coupled trim will stop automatically when the f values
   !                 converge to smaller than the chosen criterion.
   !            = 0; coupled trim will stop when iflagn is reached.
   !     tswch  = -n; n iterations before jacobian calculation
   !             n iterations are at fixed controls
   !     indnl  = 1 ; include non-linear terms  (else = 0)
   !     istab  = 0 ; stability analysis not performed
   !            = 1 ; quasi-steady floquet stability (rotating frame for
   !                  isolated blade)
   !            = 3 ; quasi-steady const coefficient
   !            = 4 ; quasi-steady floquet stab
   !            = 5 ; linear unsteady floquet stab (rotating frame for
   !                  isolated blade)
   !            = 6 ; linear unsteady stability analysis (transient)
   !            = 7 ; non-linear unsteady stability analysis (transient)
   !              8 ; quasi-steady transient analysis for all the
   !                  blades in the rotating frame, and the option of
   !                  with or without body motion.
   !     iflagn = total number of iterations to be performed
   !     mdsqss = mode sequence -- response
   !     nefs(i)= no of elements for segment i
   !     indrns = *** used in subroutine asbgbm ***
   !     indrns =-2; calculate only global mass and stiffness
   !                 matrices for stability natural modes
   !                 and frequencies
   !            =-1; calculate only global mass and stiffness
   !                 matrices for response natural modes and
   !                 frequencies
   !            = 1; calculate mass, damping, stiffness, and load
   !                 matrices for response analysis & accelerations
   !            = 2; calculate mass, damping, stiffness, and load
   !                 matrices for stability analysis
   !     irunge = 0; uses dgear routine for floquet integration
   !              1; uses 4th order runge kutta for floquet integration
   !     tol    =    tolerance for dgear integration method
   !     idisbd = 1; calculate natural frequency and blade response for
   !                 each blade and sum the blade loads to get the fixed
   !                 frame hub loads. the blade properties can also be
   !                 diffferent if desired.
   !     th75i  = pitch setting at 75% blade radius

   !        nbm     = number of flap and lag bending modes (max = 10)
   !        rroot   = blade root cutout

   !        notes: six harmonics of each bending mode are used

   !     TO RUN PRESCRIBED/FREE-WAKE...

   !     *  Set mpsi,rgmax,mra & ra(i) in MAIN

   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !                                                                      c
   !     parameter dictionary                                             c

   !     user defined parameters:
   !           ngaust= no. of time gaussian integration points per elem.
   !                   can use 5 or 6, but 5 is sufficient for 15 dof elem
   !           lsft  = 1: pitch link is soft
   !                 = 0: pitch link is rigid
   !           nmodet= no. of normal modes used for trim
   !           nmodes= no. of normal modes used for stability
   !           maxmod= maximum no. of modes for dimensioning
   !                   (currently set to 30)
   !           ninteg= integer multiplicant which determines npsi
   !                   ( npsi = nblade*ntelt*ninteg)
   !                   suggestion:  choose 1 if ickout = 0
   !                                choose 3 if ickout > 0

   !     fixed parameters:
   !           maxhub= max no. of hub dof used for aeromech stab (5)
   !           mfust = no. of rigid fuselage dof used for trim (6)
   !           ngauss= no. of spatial gaussian integration points per elem
   !                   need to use 6

   !
   !     derived parameters:
   !           ndt   = total number of spatial dof for the blade
   !           ngd   =   "     "     "    "     "     "     "  minus nbcs
   !           ndt   = total number of spatial dof for the blade
   !           mfxd  = no. of states for stability analysis in fixed frame
   !           npsi  = no. of azimuth locations =  blade* telt* integ

   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

program BModes

USE Allocation
USE Basic
USE CNPar
USE Conf
USE DatRef
USE DisBd
USE DisBld
USE Gauss
USE NWTC_Library
USE Omg
USE Param
USE pbrng
USE PLink
USE PreCone
USE SftLnk
USE Struc
USE Swept
USE TipDat
USE TowDat
USE TowWires
USE TrimV

IMPLICIT NONE

REAL(ReKi)                    :: avg_iner
REAL(ReKi)                    :: axial_stff_mult
REAL(ReKi)                    :: bl_thp
REAL(ReKi)                    :: c2tg
REAL(ReKi)                    :: cg_offst_mult
REAL(ReKi)                    :: cm_pform
REAL(ReKi)                    :: cthw
REAL(ReKi), ALLOCATABLE       :: distr_k (:)
REAL(ReKi), ALLOCATABLE       :: distr_m (:)
REAL(ReKi)                    :: dmplag
REAL(ReKi)                    :: dpsi
REAL(ReKi)                    :: draft = 0.
REAL(ReKi)                    :: dumm (6,6)
REAL(ReKi)                    :: edge_iner
REAL(ReKi)                    :: edge_iner_r
REAL(ReKi)                    :: edge_stff_mult
REAL(ReKi)                    :: eflap
REAL(ReKi), ALLOCATABLE       :: el_loc  (:)
REAL(ReKi)                    :: elag
REAL(ReKi)                    :: elt
REAL(ReKi), PARAMETER         :: eps1    = 1.0e-10
REAL(ReKi), ALLOCATABLE       :: evl     (:)
REAL(ReKi), ALLOCATABLE       :: evr     (:,:)
REAL(ReKi)                    :: flp_iner
REAL(ReKi)                    :: flp_iner_mult
REAL(ReKi)                    :: flp_iner_r
REAL(ReKi)                    :: flp_stff_mult
REAL(ReKi), ALLOCATABLE       :: gm      (:,:)
REAL(ReKi)                    :: gqpt    (ngaust)
REAL(ReKi)                    :: gqwt    (ngaust)

REAL(ReKi)                    :: ixx_tp !jan-08
REAL(ReKi)                    :: iyy_tp
REAL(ReKi)                    :: izz_tp
REAL(ReKi)                    :: ixy_tp 
REAL(ReKi)                    :: iyz_tp
REAL(ReKi)                    :: izx_tp

REAL(ReKi)                    :: lag_iner_mult
REAL(ReKi)                    :: mass_pform
REAL(ReKi)                    :: om2
REAL(ReKi)                    :: om2r
REAL(ReKi)                    :: p_base
REAL(ReKi)                    :: psi0
REAL(ReKi)                    :: rad_iner
REAL(ReKi)                    :: ref1
REAL(ReKi)                    :: ref2
REAL(ReKi)                    :: ref3
REAL(ReKi)                    :: ref4
REAL(ReKi)                    :: ref_mr
REAL(ReKi)                    :: ref_mr3
REAL(ReKi)                    :: ref_msl
REAL(ReKi)                    :: rmom2 ! feb-08
REAL(ReKi)                    :: rot_rpm  !march-08
REAL(ReKi)                    :: rpm_mult !march-08
REAL(ReKi)                    :: rroot
REAL(ReKi)                    :: sc_offst_mult
REAL(ReKi)                    :: sec_mass_mult
REAL(ReKi)                    :: sprlag
REAL(ReKi)                    :: tc_offst_mult
REAL(ReKi)                    :: tfinal
REAL(ReKi)                    :: th75
REAL(ReKi)                    :: th75i
REAL(ReKi)                    :: thetg
REAL(ReKi)                    :: tol
REAL(ReKi)                    :: tor_stff_mult
REAL(ReKi), Dimension(6,6)    :: T_pform = 0.
REAL(ReKi)                    :: tw_iner
REAL(ReKi)                    :: xfact
REAL(ReKi)                    :: xjfct1
REAL(ReKi)                    :: xjfct2
REAL(ReKi)                    :: xlhs
REAL(ReKi)                    :: xmid
REAL(ReKi), ALLOCATABLE             :: z_distr_k (:)
REAL(ReKi), ALLOCATABLE             :: z_distr_m (:)

INTEGER                       :: hub_conn
INTEGER                       :: i
INTEGER                       :: ibd
INTEGER                       :: icpltr
INTEGER                       :: id_intp
INTEGER                       :: id_mat
INTEGER                       :: iflagn
INTEGER                       :: ifrwke
INTEGER                       :: ind     (nedof)
INTEGER                       :: inewke
INTEGER                       :: inpjac
INTEGER                       :: irunge
INTEGER                       :: ise
INTEGER                       :: isttus
INTEGER                       :: iterfw
INTEGER                       :: itrmcv
INTEGER                       :: itrtab
INTEGER                       :: itswch
INTEGER                       :: ivecbe(nedfbe) = (/ 10,  9, 7, 1, 11, 12, 2, 3, 13, 14, 4, 5, 15,  8, 6  /)               ! See Bernoilli-EulerBeamVectors.txt for description.
INTEGER                       :: ivects(nedfts) = (/ 12, 11, 9, 1, 13, 14, 2, 3, 15, 16, 4, 5, 17, 10, 6, 18, 7, 19, 8 /)  ! See Bernoilli-EulerBeamVectors.txt for description.
INTEGER                       :: j
INTEGER                       :: lnksft
INTEGER                       :: mdsqss  (nmodet)
INTEGER                       :: n_secs_k_distr
INTEGER                       :: n_secs_m_distr
INTEGER                       :: naxial
INTEGER                       :: nazpsi
INTEGER                       :: nbm
INTEGER                       :: ncon
INTEGER                       :: negf
INTEGER                       :: negf2
INTEGER                       :: negpt
INTEGER                       :: negpt2
INTEGER                       :: nflap
INTEGER                       :: nlag
INTEGER                       :: npsicc
INTEGER                       :: nseg
INTEGER                       :: ntorsn

CHARACTER(200)                :: EchoFile                                     ! Name of the echo file.
CHARACTER(33)                 :: EchoFrmt  = "(2X,L11,2X,A, T30,' - ',A)"     ! Output format for logical parameters.
CHARACTER(33)                 :: EchoSPFmt = "(F7.5,5F10.3,4ES10.2,3F10.3)"   ! Output format for echoing the section properties.
CHARACTER(200)                :: InFile                                       ! Name of the primary input file.
CHARACTER(200)                :: OutFile                                      ! Name of the output file.
CHARACTER(200)                :: sec_props_file                               ! Name of the section-properties input file.


CALL NWTC_Init()

!bjj replaced with NWTC_Init: call OpenCon   ! opens the console !for_unix?
call SetProg   ! retrieves Prog_name and Prog_version
call DispNVD   ! displays the name of the program, version number, and date of release on console screen

!--- !for_unix?
call GET_COMMAND_ARGUMENT(1,InFile)

IF ( LEN_TRIM( InFile ) == 0 )  THEN
   CALL WrScr     ( ' ' )
   CALL ProgAbort ( '  You must specify the name of the input file on the command line.' )
END IF

call OpenFInpFile (UnIn, InFile)     !opens a formatted input file
call NameOFile (1, 'out', OutFile)   !name OutFile to be *.out
!---

   ! Compute some common constants using Pi.
!bjj replaced with NWTC_Init: CALL PiConsts  !for_unix?

   !hardwired vlaues
id_form = 1  ! 2 for h/c !GBBB
rm   = 10.0 !GBBB
romg = 10.0 !march-08


   ! Process the header from the primary input file.

CALL ReadCom ( UnIn, InFile, 'BModes main header line' )
CALL ReadStr ( UnIn, InFile, title, 'title', 'the title line for the input file' )

   ! Process the general parameters from the primary input file.

CALL ReadCom ( UnIn, InFile, 'blank line' )
CALL ReadCom ( UnIn, InFile, 'General Parameters comment line' )

CALL ReadVar ( UnIn, InFile, Echo       , 'Echo'  , 'Echo input file contents to *.echo file if true'   )

IF ( Echo )  THEN
   CALL NameOFile ( 1, 'echo', EchoFile )
   CALL OpenFOutFile ( UnEc, EchoFile )
   WRITE (UnEc,EchoFrmt)  Echo, 'Echo', 'Echo input file contents to *.echo file if true'
END IF ! ( Echo )

CALL ReadVar ( UnIn, InFile, beam_type  , 'beam_type'  , 'beam type, 1: blade; 2: tower'                     )
CALL ReadVar ( UnIn, InFile, rot_rpm    , 'rot_rpm'    , 'rotor speed (rpm)'                                 )
CALL ReadVar ( UnIn, InFile, rpm_mult   , 'rpm_mult'   , 'rotor speed multiplicative factor'                 )
CALL ReadVar ( UnIn, InFile, radius     , 'radius'     , 'rotor tip radius or tower height above ground (m)' )
CALL ReadVar ( UnIn, InFile, rroot      , 'rroot'      , 'hub radius or tower-base height (m)'               )
CALL ReadVar ( UnIn, InFile, btp        , 'btp'        , 'precone (deg)'                                     )
CALL ReadVar ( UnIn, InFile, bl_thp     , 'bl_thp'     , 'blade pitch setting (deg)'                         )
CALL ReadVar ( UnIn, InFile, hub_conn   , 'hub_conn'   , 'hub-to-blade connectivity identifier'              ) !GBBB
CALL ReadVar ( UnIn, InFile, modepr     , 'modepr'     , 'number of modes to be printed'                     )
CALL ReadVar ( UnIn, InFile, TabDelim   , 'TabDelim'   , 'output format (t: std; f: tab-delimited)'          )
CALL ReadVar ( UnIn, InFile, mid_node_tw, 'mid_node_tw', 't: output twist at mid nodes; f: do otherwise'     )


   !input checks (expand later)
if ( rot_rpm  <    0.0 )  CALL ProgAbort ( 'Rotor speed, rot_rpm, must be positive.' )
if ( rpm_mult <    0.0 )  CALL ProgAbort ( 'Muliplicative factor, rpm_mult, must be positive.' )
if ( radius   <    0.0 )  CALL ProgAbort ( 'Rotor radius or tower height must be positive.' )
if ( rroot    <    0.0 )  CALL ProgAbort ( 'Hub radius or tower rigid base height, rroot, must be positive.' )
if (beam_type == 1) then
  if ( radius   <= rroot )  CALL ProgAbort ( 'Rotor radius must be greater than hub radius.' )
endif
if ( hub_conn < 1 .or. hub_conn > 3 )  CALL ProgAbort ( 'hub_conn must be 1, 2 or 3' ) !GBBB
if ( modepr   <      1 )  CALL ProgAbort ( 'modepr must be a positive integer.' )

rot_rpm = rot_rpm*rpm_mult

if ( beam_type == 2 )  THEN
   rot_rpm  = 0.0
   btp      = 0.0
!jan-08   hub_conn = 1 !GBBB
end if

   ! Process the blade-tip or tower-top parameters from the primary input file.

CALL ReadCom ( UnIn, InFile, 'blank line' )
CALL ReadCom ( UnIn, InFile, 'blade-tip or tower-top mass properties comment line' )

CALL ReadVar ( UnIn, InFile, tip_mass, 'tip_mass', 'tip mass'                                          )
CALL ReadVar ( UnIn, InFile, cm_loc  , 'cm_loc'  , 'tip-mass c.m. location wrt the reference axis'     )
CALL ReadVar ( UnIn, InFile, cm_loc  , 'cm_axial'  , 'tip-mass c.m. axial offset wrt tip'               )
CALL ReadVar ( UnIn, InFile, ixx_tip , 'ixx_tip' , 'mass moment of inertia about x axis (wt-specific)' )
CALL ReadVar ( UnIn, InFile, iyy_tip , 'iyy_tip' , 'mass moment of inertia about y axis'               )
CALL ReadVar ( UnIn, InFile, izz_tip , 'izz_tip' , 'mass moment of inertia about z axis'               )
CALL ReadVar ( UnIn, InFile, ixy_tip , 'ixy_tip' , 'cross product of inertia'                          )
CALL ReadVar ( UnIn, InFile, izx_tip , 'izx_tip' , 'cross product of inertia'                          )
CALL ReadVar ( UnIn, InFile, iyz_tip , 'iyz_tip' , 'cross product of inertia'                          )


   ! Process the distributed-property parameters from the primary input file.

CALL ReadCom ( UnIn, InFile, 'blank line' )
CALL ReadCom ( UnIn, InFile, 'Distributed-property identifiers comment line' )

CALL ReadVar ( UnIn, InFile, id_mat        , 'id_mat'        , 'material isotropy identifier (not used; use later)' ) ! LATER
CALL ReadVar ( UnIn, InFile, sec_props_file, 'sec_props_file', 'name of beam section properties file (-)'           )

if ( id_mat /= 1 )  CALL ProgAbort ( 'Material identifier, id_mat, must be 1.' ) !LATER


   ! Process the property scaling factors from the primary input file.

CALL ReadCom ( UnIn, InFile, 'blank line' )
CALL ReadCom ( UnIn, InFile, 'Property scaling factors comment line' )

CALL ReadVar ( UnIn, InFile, sec_mass_mult  , 'sec_mass_mult'  , 'mass density multiplier (-)'                              )
CALL ReadVar ( UnIn, InFile, flp_iner_mult  , 'flp_iner_mult'  , 'blade flap or tower f-a inertia multiplier (-)'           )
CALL ReadVar ( UnIn, InFile, lag_iner_mult  , 'lag_iner_mult'  , 'blade lag or tower s-s inertia multiplier (-)'            )
CALL ReadVar ( UnIn, InFile, flp_stff_mult  , 'flp_stff_mult'  , 'blade flap or tower f-a bending stiffness multiplier (-)' )
CALL ReadVar ( UnIn, InFile, edge_stff_mult , 'edge_stff_mult' , 'blade lag or tower s-s bending stiffness multiplier (-)'  )
CALL ReadVar ( UnIn, InFile, tor_stff_mult  , 'tor_stff_mult'  , 'torsion stiffness multiplier (-)'                         )
CALL ReadVar ( UnIn, InFile, axial_stff_mult, 'axial_stff_mult', 'axial stiffness multiplier (-)'                           )
CALL ReadVar ( UnIn, InFile, cg_offst_mult  , 'cg_offst_mult'  , 'cg offset multiplier (-)'                                 )
CALL ReadVar ( UnIn, InFile, sc_offst_mult  , 'sc_offst_mult'  , 'shear center multiplier (-)'                              )
CALL ReadVar ( UnIn, InFile, tc_offst_mult  , 'tc_offst_mult'  , 'tension center multiplier (-)'                            )


   ! Process the finite-element-discretization parameters from the primary input file.

CALL ReadCom ( UnIn, InFile, 'blank line' )
CALL ReadCom ( UnIn, InFile, 'Finite-element discretization comment line' )

CALL ReadVar ( UnIn, InFile, nselt, 'nselt', 'number of blade or tower elements (-)' )

ALLOCATE ( el_loc(nselt+1), STAT=isttus )
if ( isttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, el_loc, in main program.' )

CALL Alloc1

CALL ReadCom ( UnIn, InFile, 'normalized element locations comment line' )
CALL ReadAry ( UnIn, InFile, el_loc, nselt+1, 'el_loc', 'array of normalized element locations (-)' )

   ! If  tower has support subsystem (platform, monopile, or tension wires)

if ( beam_type == 2 ) then

   CALL ReadCom ( UnIn, InFile, 'blank line' )
   CALL ReadCom ( UnIn, InFile, 'Properties of tower support subsystem -- comment line' )
   CALL ReadVar ( UnIn, InFile, tow_support, 'tow_support', 'aditional tower support (-)' )
   if (tow_support < 0 .or. tow_support > 2) CALL ProgAbort ( ' tow_support cannot be less than 0 or greater than 2' )
   
   if ( tow_support == 2 ) then
    
 ! platform data
 
     CALL ReadVar ( UnIn, InFile, draft, 'draft', 'depth of tower base from the MSL (mean sea level) (m)'  )
     CALL ReadVar ( UnIn, InFile, cm_pform, 'cm_pform', 'distance of platform c.m. below the MSL (m)'      )
     CALL ReadVar ( UnIn, InFile, mass_pform, 'mass_pform', 'platform mass (kg)' )
     
     do i =1,3
       i_matrix_pform(i,i)=mass_pform
     enddo
     
     CALL ReadCom ( UnIn, InFile, 'Platform inertia matrix comment line' )
     CALL ReadAry ( UnIn, InFile, i_matrix_pform(4,4:6), 3, 'i_matrix_pform', 'platform inertia matrix-row1' )
     CALL ReadAry ( UnIn, InFile, i_matrix_pform(5,4:6), 3, 'i_matrix_pform', 'platform inertia matrix-row2' )
     CALL ReadAry ( UnIn, InFile, i_matrix_pform(6,4:6), 3, 'i_matrix_pform', 'platform inertia matrix-row3' )
          
     CALL ReadVar ( UnIn, InFile, ref_msl, 'ref_msl', 'distance of platform reference point below the MSL (m)' )
     
     CALL ReadCom ( UnIn, InFile, 'Platform-reference-point-referred hydrodynamic (added-mass) matrix - comment line' )   
     do i=1,6
       CALL ReadAry ( UnIn, InFile, hydro_M(i,:), 6, 'hydro_M', 'platform added-mass inertia matrix' )
     enddo !i
     
     CALL ReadCom ( UnIn, InFile, 'Platform-reference-point-referred hydrodynamic stiffness matrix - comment line' )    
     do i=1,6
       CALL ReadAry ( UnIn, InFile, hydro_K(i,:), 6, 'hydro_K', 'platform hydrodynamic stiffness matrix' )
     enddo !i
     
!     CALL ReadCom ( UnIn, InFile, 'Platform-reference-point-referred mooring-system mass matrix - comment line' )    
!     do i=1,6
!       CALL ReadAry ( UnIn, InFile, mooring_M(i,:), 6, 'mooring_M', 'platform mooring-system mass matrix' )
!     enddo !i
     
     CALL ReadCom ( UnIn, InFile, 'Platform-reference-point-referred mooring-system stiffness matrix - comment line' )    
     do i=1,6
       CALL ReadAry ( UnIn, InFile, mooring_K(i,:), 6, 'mooring_K', 'platform mooring-system stiffness matrix' )
     enddo !i
     
! distributed added-mass data
     
     CALL ReadCom ( UnIn, InFile, 'blank line' )
     CALL ReadCom ( UnIn, InFile, 'Distributed added-mass per unit length along the tower base' )
     CALL ReadVar ( UnIn, InFile, n_secs_m_distr, 'n_secs_m_distr', 'number of points at which added mass per unit length is specified (-)'      )
     ALLOCATE (z_distr_m(n_secs_m_distr), STAT=isttus )
     if ( isttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, z_distr_m, in main program.' )
     ALLOCATE (distr_m(n_secs_m_distr), STAT=isttus )
     if ( isttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, distr_m, in main program.' )
     CALL ReadAry ( UnIn, InFile, z_distr_m, n_secs_m_distr, 'z_distr_m', 'section heights on flexible tower at which added mass is specified (m)' )
     CALL ReadAry ( UnIn, InFile, distr_m, n_secs_m_distr, 'distr_m', 'added distributed masses per unit length (kg/m)' )
   
! distributed stiffness data
     
     CALL ReadCom ( UnIn, InFile, 'blank line' )
     CALL ReadCom ( UnIn, InFile, 'Distributed elastic stiffness per unit length along the tower base' )
     CALL ReadVar ( UnIn, InFile, n_secs_k_distr, 'n_secs_k_distr', 'number of points at which distributed stiffness per unit length is specified (-)'      )
     ALLOCATE (z_distr_k(n_secs_k_distr), STAT=isttus )
     if ( isttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, z_distr_k, in main program.' )
     ALLOCATE (distr_k(n_secs_k_distr), STAT=isttus )
     if ( isttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, distr_k, in main program.' )
     CALL ReadAry ( UnIn, InFile, z_distr_k, n_secs_k_distr, 'z_distr_k', 'section heights on flexible tower at which distributed stiffness are specified (m)' )
     CALL ReadAry ( UnIn, InFile, distr_k, n_secs_k_distr, 'distr_k', 'distributed stiffnesses per unit length (N/m^2)' )

   elseif (tow_support == 1) then  ! tension wires
  
     CALL ReadCom ( UnIn, InFile, 'Tension wires data -- comment line' )
     CALL ReadVar ( UnIn, InFile, n_attachments, 'n_attachments', 'number of wire-attachment locations on tower (-)' )

     if (n_attachments < 0 ) CALL ProgAbort ( '   n_attachments cannot be -ve' )
     
     ALLOCATE (n_wires(n_attachments), STAT=isttus )
     if ( isttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, n_wires, in main program.' )
     ALLOCATE (node_attach(n_attachments), STAT=isttus )
     if ( isttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, node_attach, in main program.' )
     ALLOCATE (wire_stfness(n_attachments), STAT=isttus )
     if ( isttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, wire_stfness, in main program.' )
     ALLOCATE (th_wire(n_attachments), STAT=isttus )
     if ( isttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, th_wire, in main program.' )
     ALLOCATE (k_tower(n_attachments), STAT=isttus )
     if ( isttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, k_tower, in main program.' )
       
     CALL ReadAry ( UnIn, InFile, n_wires     , n_attachments, 'n_wires'     , 'no of wires attached at each location (-)'       )
     CALL ReadAry ( UnIn, InFile, node_attach , n_attachments, 'node_attach' , 'node numbers of attacments location (-)'         )
     CALL ReadAry ( UnIn, InFile, wire_stfness, n_attachments, 'wire_stfness', 'wire spring constant in each set (N/m)'          )
     CALL ReadAry ( UnIn, InFile, th_wire     , n_attachments, 'th_wire'     , 'angle of tension wires wrt the tower axis (deg)' )

     do i=1,n_attachments

       if ( n_wires    (i) <       3 )  CALL ProgAbort ( '   n_wires must be 3 or higher' )
       if ( node_attach(i) <=      1 )  CALL ProgAbort ( '   node_attach must be more than 1' ) !GBBB
       if ( node_attach(i) > nselt+1 )  CALL ProgAbort ( '   node_attach must be less than nselt+1' ) !GBBB
       if ( ( th_wire(i) > 90.0 ) .or. ( th_wire(i) < 0.0 ) )  then
          CALL ProgAbort ( ' *WARNING: th_wire angle is outside 0-90 deg range.' )
       end if

     enddo ! i

     do i=1,n_attachments
       th_wire(i) = th_wire(i)*d2r
       cthw       = cos(th_wire(i))
       k_tower(i) = n_wires(i)*wire_stfness(i)*cthw*cthw !GBBB
     end do ! i
   endif !tow_support
end if   ! beam_type

!-------------    
radius = radius + draft
bl_len = radius - rroot ! flexible-blade or flexible-tower length
if (bl_len < 0.) CALL ProgAbort ( 'rroot must be smaller than radius + dfaft' )

do i=1,nselt+1
   el_loc(i) = ( el_loc(i)*bl_len + rroot )/radius  ! non-d location wrt the axis of rotation
enddo ! i

do i=1,nselt
   el(i) = el_loc(nselt+2-i) - el_loc(nselt+1-i)   !elm numbering from tip
enddo ! i

if (tow_support == 2) then
  do i=1,n_secs_m_distr
    z_distr_m(i) = ( z_distr_m(i) + rroot ) / radius ! non-d
  enddo ! i

  do i=1,n_secs_k_distr
    z_distr_k(i) = ( z_distr_k(i) + rroot ) / radius ! non-d
!ctest
! write(55,*) i, z_distr_k(i)
!end test
  enddo ! i
endif

!-------------    

   ! Read the section-properties file

CALL OpenFInpFile ( UnSP, sec_props_file )

CALL ReadCom ( UnSP, sec_props_file, 'section properties header line' )
CALL ReadVar ( UnSP, sec_props_file, n_secs, 'n_secs', 'number of blade or tower sections (-)' )
CALL ReadCom ( UnSP, sec_props_file, 'blank line' )
CALL ReadCom ( UnSP, sec_props_file, 'section-properties column-titles line' )
CALL ReadCom ( UnSP, sec_props_file, 'section-properties column-units line' )

if ( n_secs < 2 )  CALL ProgAbort ( 'n_secs must be at least 2' )

   ! Allocate arrays for section properties.
CALL Alloc2

do i= 1, n_secs
   read (UnSP,*) sec_loc(i), str_tw(i), tw_iner, amass_den(i), flp_iner, edge_iner, flp_stff(i), edge_stff(i), tor_stff(i), &
                 axial_stff(i), cg_offst(i), sc_offst(i), tc_offst(i)

   IF ( Echo )  THEN
      WRITE (UnEc,EchoSPFmt) sec_loc(i), str_tw(i), tw_iner, amass_den(i), flp_iner, edge_iner, flp_stff(i), edge_stff(i), &
                 tor_stff(i), axial_stff(i), cg_offst(i), sc_offst(i), tc_offst(i)
   END IF ! ( Echo )

  if(beam_type == 2) then  ! if beam is a tower
    str_tw(i) = 0.
    tw_iner = 0.
    cg_offst(i) = 0.
    sc_offst(i) = 0.
    tc_offst(i) = 0.
    edge_iner = flp_iner
    edge_stff(i) = flp_stff(i)
!    edge_stff(i) = flp_stff(i)*(1.+tau*0.1e-030)+0.1e-7
  end if

  if(beam_type == 1) then   ! bld
     if (id_form == 1) then          ! wt bld
      cm_loc = -cm_loc   
      thetg = -tw_iner + str_tw(i) 
      cg_offst(i) = -cg_offst(i)
      tc_offst(i) = -tc_offst(i)
      sc_offst(i) = -sc_offst(i)
     else
      thetg = tw_iner - str_tw(i)   ! h/c bld
     endif
  else
     thetg = 0.
  end if

  if(abs(thetg) > 45.0) then   !GBBB
    write(*,*) ' **ERROR: thetg exceeds 45 deg. Contact NREL'
  end if

  rad_iner = 0.5*(edge_iner - flp_iner)
  avg_iner = 0.5*(edge_iner + flp_iner)
  c2tg = cos(thetg*pi/90.)
  flp_iner_r = avg_iner - rad_iner*c2tg  ! about chord ref axis  !GBBB
  edge_iner_r = avg_iner + rad_iner*c2tg

!  sq_km1(i) = sqrt(flp_iner_r/amass_den(i))      !jan-08
!  sq_km2(i) = sqrt(edge_iner_r/amass_den(i))     !jan-08
  sq_km1(i) = flp_iner_r/amass_den(i)      !jan-08
  sq_km2(i) = edge_iner_r/amass_den(i)     !jan-08
  tc_offst(i) = tc_offst(i) - sc_offst(i) !modify later  !GBBB

  sec_loc(i) = ( sec_loc(i)*bl_len + rroot ) / radius
enddo

close ( UnSP )

   !     input reading completed   ----
   
 !  tip mass related xtions
 
  if (id_form == 1) then         !wt
    if(beam_type == 1) then      !wt bld
     ixx_tp = izz_tip   !jan-08
     iyy_tp = iyy_tip   !jan-08
     izz_tp = ixx_tip   !jan-08
     ixy_tp = -iyz_tip  !jan-08
     iyz_tp = -ixy_tip  !jan-08
     izx_tp = izx_tip   !jan-08
    else                         !wt tow
     ixx_tp = izz_tip   !jan-08
     iyy_tp = ixx_tip   !jan-08
     izz_tp = iyy_tip   !jan-08
     ixy_tp = izx_tip   !jan-08
     iyz_tp = ixy_tip   !jan-08
     izx_tp = iyz_tip   !jan-08
    endif
  else                           !h/c bld
     ixx_tp = ixx_tip   !jan-08
     iyy_tp = iyy_tip   !jan-08
     izz_tp = izz_tip   !jan-08
     ixy_tp = ixy_tip   !jan-08
     iyz_tp = iyz_tip   !jan-08
     izx_tp = izx_tip   !jan-08
  end if


 !  platform related xtions
  
  if (tow_support == 2) then     !aditional tower support
  
  ! ---  i_matrix_pform of platform
  
     p_base = cm_pform - draft + rroot
     
     T_pform(1,2) =  1.
     T_pform(1,3) = -p_base
     T_pform(2,4) =  1.
     T_pform(2,5) = -p_base
     T_pform(3,1) =  1.
     T_pform(4,5) = -1.
     T_pform(5,3) =  1.
     T_pform(6,6) =  1.

     dumm = matmul( i_matrix_pform, T_pform )
     i_matrix_pform = matmul( transpose(T_pform), dumm )
     
     
  ! ---  hydro_M of platform
    
    p_base = ref_msl - draft + rroot 
     
     T_pform(1,3) = -p_base
     T_pform(2,5) = -p_base
     
     dumm = matmul( hydro_M, T_pform )
     hydro_M = matmul( transpose(T_pform), dumm )
     
   ! ---  hydro_K of platform  
        
     dumm = matmul( hydro_K, T_pform )
     hydro_K = matmul( transpose(T_pform), dumm )
   
   ! --- mooring_K
         
     dumm = matmul( mooring_K, T_pform )
     mooring_K = matmul( transpose(T_pform), dumm )

  endif


   ! normalizations and conversions ----------------------

omega = rot_rpm*pi/30. !rad/sec
btp   = btp*d2r        !rad
th75  = -bl_thp*d2r    !rad  !GBBB

omega = omega/romg   !March-08
omega2 = omega*omega !March-08

om2 = romg*romg
om2r = om2*radius
rmom2 = rm*om2
ref1 = rm*om2r
ref2 = ref1*radius
ref3 = ref2*radius
ref4 = ref3*radius

ref_mr = rm*radius
ref_mr3 = ref_mr*radius*radius

tip_mass = tip_mass/ref_mr
cm_loc  = cm_loc/radius
cm_axial = cm_axial/radius ! jan-08
ixx_tip = ixx_tp/ref_mr3 ! jan-08
iyy_tip = iyy_tp/ref_mr3 ! jan-08
izz_tip = izz_tp/ref_mr3 ! jan-08
ixy_tip = ixy_tp/ref_mr3 ! jan-08
iyz_tip = iyz_tp/ref_mr3 ! jan-08
izx_tip = izx_tp/ref_mr3 ! jan-08

if (tow_support == 2) then

! i_matrix_pform
  i_matrix_pform(1,:) = i_matrix_pform(1,:)/ref2
  i_matrix_pform(2,:) = i_matrix_pform(2,:)/ref2
  i_matrix_pform(4,:) = i_matrix_pform(4,:)/ref2
  i_matrix_pform(3,:) = i_matrix_pform(3,:)/ref3
  i_matrix_pform(5,:) = i_matrix_pform(5,:)/ref3
  i_matrix_pform(6,:) = i_matrix_pform(6,:)/ref3
  
  i_matrix_pform(:,1) = i_matrix_pform(:,1)*om2r
  i_matrix_pform(:,2) = i_matrix_pform(:,2)*om2r
  i_matrix_pform(:,4) = i_matrix_pform(:,4)*om2r
  i_matrix_pform(:,3) = i_matrix_pform(:,3)*om2
  i_matrix_pform(:,5) = i_matrix_pform(:,5)*om2
  i_matrix_pform(:,6) = i_matrix_pform(:,6)*om2
  

! hydro_M of platform
  hydro_M(1,:) = hydro_M(1,:)/ref2
  hydro_M(2,:) = hydro_M(2,:)/ref2
  hydro_M(4,:) = hydro_M(4,:)/ref2
  hydro_M(3,:) = hydro_M(3,:)/ref3
  hydro_M(5,:) = hydro_M(5,:)/ref3
  hydro_M(6,:) = hydro_M(6,:)/ref3
  
  hydro_M(:,1) = hydro_M(:,1)*om2r
  hydro_M(:,2) = hydro_M(:,2)*om2r
  hydro_M(:,4) = hydro_M(:,4)*om2r
  hydro_M(:,3) = hydro_M(:,3)*om2
  hydro_M(:,5) = hydro_M(:,5)*om2
  hydro_M(:,6) = hydro_M(:,6)*om2
  
! hydro_K of platform
  hydro_K(1,:) = hydro_K(1,:)/ref2
  hydro_K(2,:) = hydro_K(2,:)/ref2
  hydro_K(4,:) = hydro_K(4,:)/ref2
  hydro_K(3,:) = hydro_K(3,:)/ref3
  hydro_K(5,:) = hydro_K(5,:)/ref3
  hydro_K(6,:) = hydro_K(6,:)/ref3
  
  hydro_K(:,1) = hydro_K(:,1)*radius
  hydro_K(:,2) = hydro_K(:,2)*radius
  hydro_K(:,4) = hydro_K(:,4)*radius
  
! mooring_K of platform
  mooring_K(1,:) = mooring_K(1,:)/ref2
  mooring_K(2,:) = mooring_K(2,:)/ref2
  mooring_K(4,:) = mooring_K(4,:)/ref2
  mooring_K(3,:) = mooring_K(3,:)/ref3
  mooring_K(5,:) = mooring_K(5,:)/ref3
  mooring_K(6,:) = mooring_K(6,:)/ref3
  
  mooring_K(:,1) = mooring_K(:,1)*radius
  mooring_K(:,2) = mooring_K(:,2)*radius
  mooring_K(:,4) = mooring_K(:,4)*radius
  
endif

if (tow_support == 1) then
! tower tension wires
  do i = 1,n_attachments
    k_tower(i) = k_tower(i)/ref1
  enddo
endif

do i= 1, n_secs

  if (id_form == 1) then          ! wt
    str_tw(i) = -str_tw(i)*d2r     !rad
  else
    str_tw(i) = str_tw(i)*d2r     !rad
  endif
  amass_den(i) = sec_mass_mult*amass_den(i)/rm
  flp_stff(i) = flp_stff_mult*flp_stff(i)/ref4
  edge_stff(i) = edge_stff_mult*edge_stff(i)/ref4
  tor_stff(i) = tor_stff_mult*tor_stff(i)/ref4
  axial_stff(i) = axial_stff_mult*axial_stff(i)/ref2
  cg_offst(i) = cg_offst_mult*cg_offst(i)/radius
  tc_offst(i) = tc_offst_mult*sc_offst_mult*tc_offst(i)/radius !later

  sq_km1(i) = (flp_iner_mult*sq_km1(i)/sec_mass_mult)/radius**2 !jan-08
  sq_km2(i) = (lag_iner_mult*sq_km2(i)/sec_mass_mult)/radius**2 !jan-08

enddo

!ctest
!    write(*,*) amass_den(2), sq_km1(2), sq_km2(2)
!end test

   !     end of normalizations and conversions ---------------------------

   iartic = 0 !GBBB
   nconf  = 0 !GBBB

   SELECT CASE ( hub_conn ) !GBBB
   CASE ( 1 )
     ifree = 0   ! cantilevered root
   CASE ( 2 )
     ifree = 1   ! free root
   CASE ( 3 )    !feb-08
     ifree = 0   
     iartic = 10 ! partially free root (axial and torsion constraints) !feb-08
   END SELECT !hub_conn

nseg   = 1 !GBBB
npin   = 0 !GBBB
lnksft = 0 !GBBB

nflap  = 2 !GBBB
nlag   = 2 !GBBB
ntorsn = 2 !GBBB
naxial = 0 !GBBB

   !     check for input/parameter compatability

if ( lsft /= lnksft)  CALL ProgAbort ( 'lsft must equal lnksft' )

if ( nmodet /= ( nflap + nlag + ntorsn + naxial) )  CALL ProgAbort ( 'nmodet must equal nflap + nlag + ntorsn + naxial' )

if ( nmodes /= ( nflap + nlag + ntorsn ) )  CALL ProgAbort ( 'nmodes must equal nflap + nlag + ntorsn' )

ndt = 9*nselt + 6
!xx ngd = ndt - nbcs !GBBB

   !----------------------------------------------------------------------
   !comp
   !     Copy either 15 or 19 dof elemental connectivity vector into 'ind':

if ( nedof == 15 )  then
   !       bernoulli-euler elements
   do i=1,nedof
      ind(i) = ivecbe(i)
   end do ! i
else if ( nedof == 19 )  then
   !       transverse shear flexible elements
   do i=1,nedof
      ind(i) = ivects(i)
   end do ! i
end if
   !comp

   !----------------------------------------------------------------------
   !     determine global and local dof parameters
   !comp

if ( nedof == 15 )  then
      ! bernoulli-euler elements
   nndof = 6
   nintn = 3
   ncon  = 6 !GBBB

   if(nseg>2) then !later
      ncon = 6*(nseg-1)
   end if

   if(nconf==2) then !later
      ncon = ncon - 3 - lnksft
      if (npin==0) then
         ncon = ncon - 2
      end if
   end if
else if ( nedof == 19 )  then

   !       transverse shear flexible elements
  nndof = 8
  nintn = 3
  ncon  = 8

      !comp: fix  here down from here

   if ( nseg > 2 )  then
      ncon = 6*( nseg - 1 )
   end if

   if ( nconf == 2 )  then
      ncon = ncon - 3 - lnksft
      if ( npin == 0 )  then
         ncon = ncon - 2
      end if
   end if

end if

   !----------------------------------------------------------------------
   !     relieve boundary constraints for articulated blade configurations

IF ( ifree == 1 )  THEN !GBBB
   ncon = 0
ELSE
   SELECT CASE ( iartic ) !GBBB
   CASE ( 1 )
      ncon = ncon - 1
   CASE ( 2, 3, 4 )
      ncon = ncon - 2
   CASE ( 10 )         ! feb-08  axial and torsion constraints
      ncon = ncon - 4
   END SELECT
END IF

negf  = nmodet !GBBB
negpt = nmodes !GBBB
   !um   ndt   = (nselt+1) * nndof + nselt * nintn
   !um   nedof = nndof * 2 + nintn
ngd   = ndt - ncon !GBBB

!xx if ( nbcs /= ncon  )  CALL ProgAbort ( 'nbcs must equal ncon' )
   !xx Allocate global arrays whose sizes are based upon input parameters

CALL Alloc3
   ! Allocate local arrays whose sizes are based upon input parameters.

ALLOCATE ( evl(ngd), STAT=isttus )
if ( isttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, evl, in main program.' )

ALLOCATE ( evr(ngd,ngd), STAT=isttus )
if ( isttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, evr, in main program.' )

ALLOCATE ( gm(ngd,ngd), STAT=isttus )
if ( isttus /= 0 )  CALL ProgAbort ( 'Unable to allocate array, gm, in main program.' )

   !----------------------------------------------------------------------
  write(6,*)'   '
  write(6,*)' all parameters checked : o.k.'
   !----------------------------------------------------------------------
   !comp fix down from here

negpt2= negpt*2 !GBBB
negf2 = 2*negf !GBBB

   !comp
   !-------------------------- input -------------------------------

   !nrel hardwired inputs continued (handle/remove later)

iflagn = 50
npsicc = 2

itrmcv = 1
ifrwke = 0
inewke = 0

iterfw = 0
icpltr = 12
itswch = 0

xjfct1 = 1.0
xjfct2 = 0.0

inpjac = 0

sprlag = 0.0
dmplag = 0.0
del3pb = 0.0

nefs(1) = nselt      ! no of elements for seg 1 !LATER

   !     do i = 1, negf
mdsqss(1) = 1
mdsqss(2) = 2
mdsqss(3) = 3
mdsqss(4) = 4
mdsqss(5) = 5
mdsqss(6) = 8
   !     enddo

   !gsb      do i = 1, nselt
   !gsb        ielrot(i) = 0
   !gsb        swrad(i) = 0.
   !gsb        drrad(i) = 0.
   !gsb        prrad(i) = 0.
   !gsb      enddo

ielrot(:) = 0 !LATER

if ( SUM( nefs(1:nseg) ) /= nselt )  then
   write (UnOu,*) ' the number of elements did not add up to nselt.'
   CALL ProgAbort ( 'The number of elements did not add up to nselt.' )
end if

do i=1,nselt
   if ( ielrot(i) /= 0 )  netip = i
end do ! i

open(unit=UnOu, file=TRIM( OutFile ))

   ! Compute Gaussian weights and abscissas.
call gqi6pt ( gqp , gqw , ngauss )
call gqi5pt ( gqpt, gqwt, ngaust )

   !  interpolate for tower or blade element props

 elm_dist_k = 0. ! corrp to distr stiffness
 elm_dist_m = 0. ! corrp to distr mass
 
 xlhs = 1.0

do i=1,nselt

   xlhs = xlhs - el(i)
   xmid = xlhs + el(i)*0.5
   id_intp = 0

   do j=1,n_secs-1

!      if(xmid >= sec_loc(j)-eps1 .and. xmid < sec_loc(j+1)) then
      if(xmid >= sec_loc(j) .and. xmid <= sec_loc(j+1)) then
         xfact = (xmid-sec_loc(j)) / (sec_loc(j+1)-sec_loc(j))
         rmas(i) = amass_den(j)+(amass_den(j+1)-amass_den(j))*xfact
         skm1(i) = sq_km1(j) + (sq_km1(j+1)-sq_km1(j))*xfact
         skm2(i) = sq_km2(j) + (sq_km2(j+1)-sq_km2(j))*xfact
         eiy(i) = flp_stff(j) + (flp_stff(j+1)-flp_stff(j))*xfact
         eiz(i) = edge_stff(j) + (edge_stff(j+1)-edge_stff(j))*xfact
         gj(i) = tor_stff(j) + (tor_stff(j+1)-tor_stff(j))*xfact
         eac(i) = axial_stff(j)+(axial_stff(j+1)-axial_stff(j))*xfact
         eg(i) = cg_offst(j) + (cg_offst(j+1)-cg_offst(j))*xfact
         ea(i) = tc_offst(j) + (tc_offst(j+1)-tc_offst(j))*xfact
         id_intp = id_intp+1

         exit ! do j loop

      end if
      
   enddo ! j
   
   if (tow_support == 2) then
   
     do j = 1, n_secs_m_distr-1

       if(xmid >= z_distr_m(j)-eps1 .and. xmid < z_distr_m(j+1)) then

         xfact = (xmid-z_distr_m(j)) / (z_distr_m(j+1)-z_distr_m(j))
         elm_dist_m(i) = (distr_m(j)+(distr_m(j+1)-distr_m(j))*xfact)/rm ! non-d

         exit ! do j loop

       end if
      
     enddo ! j
   
   
     do j = 1, n_secs_k_distr-1

       if(xmid >= z_distr_k(j)-eps1 .and. xmid < z_distr_k(j+1)) then

         xfact = (xmid-z_distr_k(j)) / (z_distr_k(j+1)-z_distr_k(j))
         elm_dist_k(i) = (distr_k(j)+(distr_k(j+1)-distr_k(j))*xfact)/rmom2 ! non-d
!ctest
      if(xmid == sec_loc(j) .and. xmid == sec_loc(j+1)) then
       write(*,*) j, xmid, z_distr_k(j), z_distr_k(j+1), xfact
      endif 
!end test

         exit ! do j loop

       end if
      
     enddo ! j
   
   endif ! tow_support=2

   
!ctest
!      write(55,5551) i, elm_dist_m(i),elm_dist_k(i)
!5551  format (1x, i3, 3x, 2e17.6)
!end test

enddo ! i

if(id_intp == 0) call ProgAbort('stop 3333: interpolation failed for sections')
if(id_intp >= 2) call ProgAbort('stop 4444: multiple interpolation detected for sections')


   !comp
   !nrel      if (medof == 19) then
   !       read (UnIn,*   ) (gay(i),i=1,nselt)
   !       read (UnIn,*   ) (gaz(i),i=1,nselt)
   !--      end if
   !comp

   !      read (UnIn,*) (eg(i),i=1,nselt)
   !      read (UnIn,*) (ea(i),i=1,nselt)

   !      read (UnIn,*) (rmas(i),i=1,nselt)
   !      read (UnIn,*) (skm1(i),i=1,nselt)
   !      read (UnIn,*) (skm2(i),i=1,nselt)

      eb1(:) = 0.0 !GBBB
      eb2(:) = 0.0
      ec1(:) = 0.0
      ec2(:) = 0.0

         !     normalize convert chordwise offsets by radius (input wrt chord)
         !     calculate radius of gyration from bending and axial properties

itrtab = 0 !LATER
th75i = 0. !LATER

   !nrel introduced recently (how were these handled earlier???)
th1c = 0. !LATER
th1s = 0. !LATER
   !--
   !     input propulsive trim controls  (if icpltr is in the range 2 to 9)


   !jmw  inputs for runge kutta integration flag used for floquet

  irunge = 0 !LATER

  tol = 0 !LATER

   !     inputs for dissimilar blades

  idisbd = 0 !LATER

   !--
if ( idisbd == 1 )  then
   do ibd=1,nblade
      read(UnIn,*)(eiyfac(ibd,ise),ise=1,nselt)
      read(UnIn,*)(eizfac(ibd,ise),ise=1,nselt)
      read(UnIn,*)(gjfac (ibd,ise),ise=1,nselt)
      read(UnIn,*)(rmasfc(ibd,ise),ise=1,nselt)
   end do ! ibd
end if


   !  ------------- add rotor speed variation for ground -----------------
   !                      and air resonance analyses


!! ttk = ttk/omegar**2

!!if ( iartic /= 0 )  then
!!  sprlag = sprlag/omegar**2
!!  dmplag = dmplag/omegar
!!end if

   !     define initial logic and stepsize parameters

tfinal  =  2*pi  !LATER
elt  =  tfinal/ntelt  !LATER
psi0  =  0.  !LATER

nazpsi = npsi/ntelt  !LATER
dpsi   =  2.0*pi/npsi  !LATER

nbm = nflap + nlag  !LATER

   !---------------------------------------------------------------------
   !      write(*,*) 'begin free vibration analysis'
   !      write(*,*) ' '

call bldvib ( nseg, ind, mdsqss, psi0, sprlag, eflap, elag, evr, evl, gm )

CALL NormStop

END

!=======================================================================
SUBROUTINE dist ( l, tl, tldot, tlddot, elv, cs, ss )  !LATER

   ! calculates distance of left edge of element
   ! from the root for swept blade

USE Param
USE PLink
USE PreCone
USE Struc
USE Swept
USE TrimV

IMPLICIT NONE


   ! Argument declarations.

REAL(ReKi), INTENT(IN)        :: CS
REAL(ReKi), INTENT(IN)        :: elv      (nselt)
REAL(ReKi), INTENT(OUT)       :: TL       (3)
REAL(ReKi), INTENT(OUT)       :: TLDot    (3)
REAL(ReKi), INTENT(IN)        :: SS

INTEGER, INTENT(IN)           :: L


   ! Local declarations.

REAL(ReKi)                    :: DDTLAI   (3,3)
REAL(ReKi)                    :: DDTLAK   (3,3)
REAL(ReKi)                    :: DDTPK1
REAL(ReKi)                    :: DDTPK2
REAL(ReKi)                    :: DDTPK3
REAL(ReKi)                    :: DTLAMI   (3,3)
REAL(ReKi)                    :: DTLAMK   (3,3)
REAL(ReKi)                    :: DTPK1
REAL(ReKi)                    :: DTPK2
REAL(ReKi)                    :: DTPK3
REAL(ReKi)                    :: EL0
REAL(ReKi)                    :: TH0D
REAL(ReKi)                    :: TH0DD
REAL(ReKi)                    :: TH75
REAL(ReKi)                    :: TLAMI    (3,3)
REAL(ReKi)                    :: TLAMK    (3,3)
REAL(ReKi)                    :: TLDDOT   (3)
REAL(ReKi)                    :: TPK1
REAL(ReKi)                    :: TPK2
REAL(ReKi)                    :: TPK3

INTEGER                       :: I
INTEGER                       :: IKSEG
INTEGER                       :: IPSEG



   !     collective pitch and time derivatives

th0 = th75 + th1c * cs + th1s * ss
th0d =  -th1c * ss + th1s * cs
th0dd =  -th1c * cs - th1s * ss

tl    (:) = 0.0
tldot (:) = 0.0
tlddot(:) = 0.0


   ! offset

el0 = 1.0 - SUM( elv )

   !        calculate transformation for current swept element

   !bsb        call swder(l,th0,th0d,th0dd,swrad(l),drrad(l),
   !bsb     &             tlamk,dtlamk,ddtlak)

   ! initialize

tl    (:) = 0.0
tldot (:) = 0.0
tlddot(:) = 0.0

   !          calculate location of left end of elements

   !
   !          segment for current element

ikseg = nselt + 1 - l

do ipseg=0,ikseg-1

   if ( ipseg /= 0 ) then

      i = nselt + 1 - ipseg

         !gsb           call swder(i,th0,th0d,th0dd,swrad(i),drrad(i),
         !gsb     &             tlami,dtlami,ddtlai)

      tpk1 = tlami(1,1)*tlamk(1,1) + tlami(1,2)*tlamk(1,2) + tlami(1,3)*tlamk(1,3)
      tpk2 = tlami(1,1)*tlamk(2,1) + tlami(1,2)*tlamk(2,2) + tlami(1,3)*tlamk(2,3)
      tpk3 = tlami(1,1)*tlamk(3,1) + tlami(1,2)*tlamk(3,2) + tlami(1,3)*tlamk(3,3)

      dtpk1 = tlami(1,1)*dtlamk(1,1) + dtlami(1,1)*tlamk(1,1) + tlami(1,2)*dtlamk(1,2) + dtlami(1,2)*tlamk(1,2) &
            + tlami(1,3)*dtlamk(1,3) + dtlami(1,3)*tlamk(1,3)

      dtpk2 = tlami(1,1)*dtlamk(2,1) + dtlami(1,1)*tlamk(2,1) + tlami(1,2)*dtlamk(2,2) + dtlami(1,2)*tlamk(2,2) &
            + tlami(1,3)*dtlamk(2,3) + dtlami(1,3)*tlamk(2,3)

      dtpk3 = tlami(1,1)*dtlamk(3,1) + dtlami(1,1)*tlamk(3,1) + tlami(1,2)*dtlamk(3,2) + dtlami(1,2)*tlamk(3,2) &
            + tlami(1,3)*dtlamk(3,3) + dtlami(1,3)*tlamk(3,3)

      ddtpk1 = tlami(1,1)*ddtlak(1,1) + ddtlai(1,1)*tlamk(1,1) + 2.0*dtlami(1,1)*dtlamk(1,1) &
             + tlami(1,2)*ddtlak(1,2) + ddtlai(1,2)*tlamk(1,2) + 2.0*dtlami(1,2)*dtlamk(1,2) &
             + tlami(1,3)*ddtlak(1,3) + ddtlai(1,3)*tlamk(1,3) + 2.0*dtlami(1,3)*dtlamk(1,3)

      ddtpk2 = tlami(1,1)*ddtlak(2,1) + ddtlai(1,1)*tlamk(2,1) + 2.0*dtlami(1,1)*dtlamk(2,1) &
             + tlami(1,2)*ddtlak(2,2) + ddtlai(1,2)*tlamk(2,2) + 2.0*dtlami(1,2)*dtlamk(2,2) &
             + tlami(1,3)*ddtlak(2,3) + ddtlai(1,3)*tlamk(2,3) + 2.0*dtlami(1,3)*dtlamk(2,3)

      ddtpk3 = tlami(1,1)*ddtlak(3,1) + ddtlai(1,1)*tlamk(3,1) + 2.0*dtlami(1,1)*dtlamk(3,1) &
             + tlami(1,2)*ddtlak(3,2) + ddtlai(1,2)*tlamk(3,2) + 2.0*dtlami(1,2)*dtlamk(3,2) &
             + tlami(1,3)*ddtlak(3,3) + ddtlai(1,3)*tlamk(3,3) + 2.0*dtlami(1,3)*dtlamk(3,3)

      tl(1) = tl(1) + elv(i)*tpk1
      tl(2) = tl(2) + elv(i)*tpk2
      tl(3) = tl(3) + elv(i)*tpk3

      tldot(1) = tldot(1) + elv(i)*dtpk1
      tldot(2) = tldot(2) + elv(i)*dtpk2
      tldot(3) = tldot(3) + elv(i)*dtpk3

      tlddot(1) = tlddot(1) + elv(i)*ddtpk1
      tlddot(2) = tlddot(2) + elv(i)*ddtpk2
      tlddot(3) = tlddot(3) + elv(i)*ddtpk3

   else

         !           add effect of offset from hub

         !gsb           call swder(1,th0,th0d,th0dd,0.d0,0.d0,
         !gsb     &             tlami,dtlami,ddtlai)

      tpk1 = tlami(1,1)*tlamk(1,1) + tlami(1,2)*tlamk(1,2) + tlami(1,3)*tlamk(1,3)
      tpk2 = tlami(1,1)*tlamk(2,1) + tlami(1,2)*tlamk(2,2) + tlami(1,3)*tlamk(2,3)
      tpk3 = tlami(1,1)*tlamk(3,1) + tlami(1,2)*tlamk(3,2) + tlami(1,3)*tlamk(3,3)

      dtpk1 = tlami(1,1)*dtlamk(1,1) + dtlami(1,1)*tlamk(1,1) &
            + tlami(1,2)*dtlamk(1,2) + dtlami(1,2)*tlamk(1,2) &
            + tlami(1,3)*dtlamk(1,3) + dtlami(1,3)*tlamk(1,3)

      dtpk2 = tlami(1,1)*dtlamk(2,1) + dtlami(1,1)*tlamk(2,1) &
            + tlami(1,2)*dtlamk(2,2) + dtlami(1,2)*tlamk(2,2) &
            + tlami(1,3)*dtlamk(2,3) + dtlami(1,3)*tlamk(2,3)

      dtpk3 = tlami(1,1)*dtlamk(3,1) + dtlami(1,1)*tlamk(3,1) &
            + tlami(1,2)*dtlamk(3,2) + dtlami(1,2)*tlamk(3,2) &
            + tlami(1,3)*dtlamk(3,3) + dtlami(1,3)*tlamk(3,3)

      ddtpk1 = tlami(1,1)*ddtlak(1,1) + ddtlai(1,1)*tlamk(1,1) + 2.0*dtlami(1,1)*dtlamk(1,1) &
             + tlami(1,2)*ddtlak(1,2) + ddtlai(1,2)*tlamk(1,2) + 2.0*dtlami(1,2)*dtlamk(1,2) &
             + tlami(1,3)*ddtlak(1,3) + ddtlai(1,3)*tlamk(1,3) + 2.0*dtlami(1,3)*dtlamk(1,3)

      ddtpk2 = tlami(1,1)*ddtlak(2,1) + ddtlai(1,1)*tlamk(2,1) + 2.0*dtlami(1,1)*dtlamk(2,1) &
             + tlami(1,2)*ddtlak(2,2) + ddtlai(1,2)*tlamk(2,2) + 2.0*dtlami(1,2)*dtlamk(2,2) &
             + tlami(1,3)*ddtlak(2,3) + ddtlai(1,3)*tlamk(2,3) + 2.0*dtlami(1,3)*dtlamk(2,3)

      ddtpk3 = tlami(1,1)*ddtlak(3,1) + ddtlai(1,1)*tlamk(3,1) + 2.0*dtlami(1,1)*dtlamk(3,1) &
             + tlami(1,2)*ddtlak(3,2) + ddtlai(1,2)*tlamk(3,2) + 2.0*dtlami(1,2)*dtlamk(3,2) &
             + tlami(1,3)*ddtlak(3,3) + ddtlai(1,3)*tlamk(3,3) + 2.0*dtlami(1,3)*dtlamk(3,3)

      tl(1) = tl(1) + el0*tpk1
      tl(2) = tl(2) + el0*tpk2
      tl(3) = tl(3) + el0*tpk3

      tldot(1) = tldot(1) + elv(i)*dtpk1
      tldot(2) = tldot(2) + elv(i)*dtpk2
      tldot(3) = tldot(3) + elv(i)*dtpk3

      tlddot(1) = tlddot(1) + elv(i)*ddtpk1
      tlddot(2) = tlddot(2) + elv(i)*ddtpk2
      tlddot(3) = tlddot(3) + elv(i)*ddtpk3

   end if

END DO ! ipseg

RETURN
END SUBROUTINE dist ! ( l, tl, tldot, tlddot, elv, cs, ss )
!=======================================================================
SUBROUTINE gztwst ( x, th75, th, thp )  !LATER

   ! Calculate gazelle's non-linear blade twist.

USE     Precision

IMPLICIT NONE


   ! Argument declarations.

REAL(ReKi), INTENT(OUT)       :: TH
REAL(ReKi), INTENT(IN)        :: TH75
REAL(ReKi), INTENT(OUT)       :: THP
REAL(ReKi), INTENT(IN)        :: X


   ! Local declarations.

REAL(ReKi)                    :: R
REAL(ReKi)                    :: THBU
REAL(ReKi)                    :: THC



r = x*5.25


   !     thbu = built in twist at r = 3.9375 meter (@ 75% radius)

thbu = .135238
thc  = th75 - thbu

if (r>4.855) then
    th = thc + .112944
    thp = 0.0
    return
end if

if (r>1.235) then
    th = thc + .2009572 - .024299*(r-1.235)
    thp = -.1276437
    return
end if

if (r>0.8244) then
    th = thc + 0.4894233*(r-0.8244)
    thp = 2.5694753
    return
end if

th = thc
thp = 0.0


RETURN
END SUBROUTINE gztwst ! ( x, th75, th, thp )
!=======================================================================
SUBROUTINE bmrmod (em, ek, gk, i)  !LATER

   !        this routine modifies the last torque element mass and
   !        stiffness matrices in accordance with the bmr bcs.
   !        for soft pitch link, the global stiffness matrix is
   !        also modified.
   !        ( caled by bldvib & stab routines )

USE Gauss
USE Param
USE plgeom
USE Precision
USE SftLnk

IMPLICIT NONE


   ! Argument declarations.

REAL(ReKi), INTENT(INOUT)  :: em       (nedof,nedof)
REAL(ReKi), INTENT(INOUT)  :: ek       (nedof,nedof)
REAL(ReKi), INTENT(INOUT)  :: gk       (ngd,ngd)

INTEGER, INTENT(IN)        :: i


   ! Local declarations.

INTEGER                    :: iphi
INTEGER                    :: iv
INTEGER                    :: ivp
INTEGER                    :: iw
INTEGER                    :: iwp
INTEGER                    :: jphi
INTEGER                    :: jv
INTEGER                    :: jvp
INTEGER                    :: jw
INTEGER                    :: jwp

iv = 5
ivp = 6
iw = 9
iwp = 10
iphi = 13

   !       note:  i = nselt
   !           u, vp, wp are independent dof
   !           phi is independent dof if pitch link is soft

   !
   !         apply kinematic boundary conditions

call modmat ( em )
call modmat ( ek )


   !           augment gm & gk matrices to include contribution of pitch
   !           link stiffness

if ( lsft == 1 )  then
   !             ttxp = ttx0 - ttp
   jv  = indeg(iv,i)
   jvp = indeg(ivp,i)
   jw  = indeg(iw,i)
   jwp = indeg(iwp,i)
   jphi= indeg(iphi,i)

   if ( npin == 1 )  then
      gk(jvp,jvp) = gk(jvp,jvp) + ttk*svp*svp
      gk(jvp,jwp) = gk(jvp,jwp) + ttk*svp*swp
      gk(jvp,jphi) = gk(jvp,jphi) + ttk*svp*sphi

      gk(jwp,jvp) = gk(jwp,jvp) + ttk*swp*svp
      gk(jwp,jwp) = gk(jwp,jwp) + ttk*swp*swp
      gk(jwp,jphi) = gk(jwp,jphi) + ttk*swp*sphi

      gk(jphi,jphi) = gk(jphi,jphi) + ttk*sphi*sphi
      gk(jphi,jwp) = gk(jphi,jwp) + ttk*sphi*swp
      gk(jphi,jvp) = gk(jphi,jvp) + ttk*sphi*svp

   else if ( npin == 0 )  then

      gk(jv ,jv ) = gk(jv ,jv ) + ttk*tv*tv
      gk(jv ,jvp ) = gk(jv ,jvp ) + ttk*tv*tvp
      gk(jv ,jw ) = gk(jv ,jw ) + ttk*tv*tw
      gk(jv ,jwp ) = gk(jv ,jwp ) + ttk*tv*twp
      gk(jv  ,jphi) = gk(jv  ,jphi) + ttk*tv*tphi

      gk(jvp,jv ) = gk(jv ,jv ) + ttk*tvp*tv
      gk(jvp,jvp ) = gk(jv ,jvp ) + ttk*tvp*tvp
      gk(jvp,jw ) = gk(jv ,jw ) + ttk*tvp*tw
      gk(jvp,jwp ) = gk(jv ,jwp ) + ttk*tvp*twp
      gk(jvp ,jphi) = gk(jv  ,jphi) + ttk*tvp*tphi

      gk(jw ,jv ) = gk(jv ,jv ) + ttk*tw*tv
      gk(jw ,jvp ) = gk(jv ,jvp ) + ttk*tw*tvp
      gk(jw ,jw ) = gk(jv ,jw ) + ttk*tw*tw
      gk(jw ,jwp ) = gk(jv ,jwp ) + ttk*tw*twp
      gk(jw  ,jphi) = gk(jv  ,jphi) + ttk*tw*tphi

      gk(jwp,jv ) = gk(jv ,jv ) + ttk*twp*tv
      gk(jwp,jvp ) = gk(jv ,jvp ) + ttk*twp*tvp
      gk(jwp,jw ) = gk(jv ,jw ) + ttk*twp*tw
      gk(jwp,jwp ) = gk(jv ,jwp ) + ttk*twp*twp
      gk(jwp ,jphi) = gk(jv  ,jphi) + ttk*twp*tphi

      gk(jphi,jv ) = gk(jv ,jv ) + ttk*tphi*tv
      gk(jphi,jvp ) = gk(jv ,jvp ) + ttk*tphi*tvp
      gk(jphi,jw ) = gk(jv ,jw ) + ttk*tphi*tw
      gk(jphi,jwp ) = gk(jv ,jwp ) + ttk*tphi*twp
      gk(jphi ,jphi) = gk(jv  ,jphi) + ttk*tphi*tphi

   end if

end if


RETURN
END SUBROUTINE bmrmod ! (em, ek, gk, i)
!=======================================================================
SUBROUTINE modmat ( emp )  !LATER

   !----------------------------------------------------------------------

   !
   !     modifies the mass and stiffness matrices in accordance with the
   !     boundary conditions at the torque tube end

   !----------------------------------------------------------------------

USE Param
USE plgeom
USE Precision
USE SftLnk

IMPLICIT NONE


   ! Argument declarations.

REAL(ReKi), INTENT(INOUT)  :: emp   (nedof,nedof)                 ! Mass or stiffness matrix.


   ! Local declarations.

INTEGER                    :: j
INTEGER                    :: iphi
INTEGER                    :: iv
INTEGER                    :: ivp
INTEGER                    :: iw
INTEGER                    :: iwp



iphi = 13
iv   =  5
ivp  =  6
iw   =  9
iwp  = 10

if ( npin == 1 )  then
   !           u, vp, wp are independent dof
   !           phi is independent dof if pitch link is soft


   do j=1,nedof

      if(lsft /= 1) then
         emp(j,iwp) = emp(j,iwp) + b1wpw*emp(j,iw) + b1wpph*emp(j,iphi) + b1wpv*emp(j,iv)
         emp(j,ivp) = emp(j,ivp) + b1vpv*emp(j,iv) + b1vpph*emp(j,iphi)
      end if

      if(lsft == 1) then
         emp(j,iwp ) = emp(j,iwp ) + ttx0*emp(j,iw)
         emp(j,ivp ) = emp(j,ivp ) + ttx0*emp(j,iv)
         emp(j,iphi) = emp(j,iphi) + ttz *emp(j,iv)
      end if

   end do ! j

   do j=1,nedof

      if ( lsft /= 1 )  then
         emp(iwp,j) = emp(iwp,j) + b1wpw*emp(iw,j) + b1wpph*emp(iphi,j) + b1wpv*emp(iv,j)
         emp(ivp,j) = emp(ivp,j) + b1vpv*emp(iv,j) + b1vpph*emp(iphi,j)
      else if ( lsft == 1 )  then
         emp(iwp ,j) = emp(iwp ,j) + ttx0*emp(iw,j)
         emp(ivp ,j) = emp(ivp ,j) + ttx0*emp(iv,j)
         emp(iphi,j) = emp(iphi,j) + ttz *emp(iv,j)
      end if

   end do ! j

else if ( npin == 0 )  then

      ! phi is   dependent dof if pitch link is rigid

   do j=1,nedof

      if ( lsft /= 1 )  then
         emp(j,iwp) = emp(j,iwp) + b2wpph*emp(j,iphi)
         emp(j,iw ) = emp(j,iw ) + b2wph *emp(j,iphi)
         emp(j,ivp) = emp(j,ivp) + b2vpph*emp(j,iphi)
         emp(j,iv ) = emp(j,iv ) + b2vph *emp(j,iphi)
      end if

   end do ! j

   do j=1,nedof

      if ( lsft /= 1 )  then
         emp(iwp,j) = emp(iwp,j) + b2wpph*emp(iphi,j)
         emp(iw ,j) = emp(iw ,j) + b2wph *emp(iphi,j)
         emp(ivp,j) = emp(ivp,j) + b2vpph*emp(iphi,j)
         emp(iv ,j) = emp(iv ,j) + b2vph *emp(iphi,j)
      end if

   end do ! j

end if


RETURN
END SUBROUTINE modmat ! ( emp )
!=======================================================================
