!=======================================================================
MODULE Basic

USE Precision

REAL(ReKi), ALLOCATABLE       :: cfe      (:)

END MODULE Basic
!=======================================================================
MODULE Conf

INTEGER                       :: nconf

END MODULE Conf
!=======================================================================
MODULE CFunc

USE Precision

REAL(ReKi)                    :: th75


END MODULE CFunc
!=======================================================================
MODULE CNPar

INTEGER                       :: indnl  = 1                             ! 1: include nonlinear terms, 0: d not include nonlinear terms
INTEGER                       :: indrns = -1                            ! Use -1 so asbgbm calculates only global mass and stiffness matrices for response natural modes and frequencies.

END MODULE CNPar
!=======================================================================
MODULE DatRef

USE Precision

REAL(ReKi)                    :: rm
REAL(ReKi)                    :: romg
REAL(ReKi)                    :: radius
REAL(ReKi)                    :: bl_len

INTEGER                       :: id_form
INTEGER                       :: modepr

LOGICAL                       :: mid_node_tw
LOGICAL                       :: TabDelim

END MODULE DatRef
!=======================================================================
MODULE DisBd

USE Precision

REAL(ReKi), ALLOCATABLE       :: eiyfac   (:,:)
REAL(ReKi), ALLOCATABLE       :: eiytmp   (:,:)
REAL(ReKi), ALLOCATABLE       :: eizfac   (:,:)
REAL(ReKi), ALLOCATABLE       :: eiztmp   (:,:)
REAL(ReKi), ALLOCATABLE       :: gjfac    (:,:)
REAL(ReKi), ALLOCATABLE       :: gjtmp    (:,:)
REAL(ReKi), ALLOCATABLE       :: rmasfc   (:,:)
REAL(ReKi), ALLOCATABLE       :: rmastm   (:,:)

LOGICAL                       :: InitSt   = .TRUE.

END MODULE DisBd
!=======================================================================
MODULE DisBld

INTEGER                       :: idisbd
INTEGER                       :: jb

END MODULE DisBld
!=======================================================================
MODULE Eig

USE Precision

REAL(ReKi), ALLOCATABLE       :: d        (:)
REAL(ReKi), ALLOCATABLE       :: eigv     (:)
REAL(ReKi), ALLOCATABLE       :: x        (:,:)

INTEGER, ALLOCATABLE          :: jdd      (:)

END MODULE Eig
!=======================================================================
MODULE Gauss

USE Precision

REAL(ReKi), ALLOCATABLE       :: el       (:)                           ! Length of element i.
REAL(ReKi), ALLOCATABLE       :: gqp      (:)                           ! i'th Gauss quadrature point.
REAL(ReKi), ALLOCATABLE       :: gqw      (:)                           ! i'th Gauss quadrature weight.
REAL(ReKi), ALLOCATABLE       :: xb       (:)                           ! Location of inboard end of i'th element.

INTEGER, ALLOCATABLE          :: indeg    (:,:)                         ! Connectivity vector. The global dof number of the i'th elemental dof of the j'th element.
INTEGER, PARAMETER            :: ngauss   = 6                           ! Number of gauss points per spatial element.

END MODULE Gauss
!=======================================================================
MODULE Gravity

USE Precision

REAL(ReKi)                    :: Grav  = 0.0                            ! Acceleration due to gravity.

END MODULE Gravity
!=======================================================================
MODULE Omg

USE Precision

REAL(ReKi)                    :: omegar

END MODULE Omg
!=======================================================================
MODULE Param


   ! Pure parameters:

INTEGER, PARAMETER            :: ifree    =  0                          !
INTEGER, PARAMETER            :: lsft     =  0                          ! 1: pitch link is soft, 0: pitch link is rigid
INTEGER, PARAMETER            :: maxhub   =  5                          !
INTEGER, PARAMETER            :: maxmd    = 10                          !
INTEGER, PARAMETER            :: nbcs     = 6                           ! Number of boundary conditions.  NOTE: This will become an input parameter later.
INTEGER, PARAMETER            :: nblade   =  1                          ! no. of blades
INTEGER                       :: ndt                                    ! Total number of DOFs.
INTEGER, PARAMETER            :: nedfbe   = 15                                        !
INTEGER, PARAMETER            :: nedfts   = 19                          !
INTEGER, PARAMETER            :: nedof    = 15                          ! Number of element DOFs.
INTEGER, PARAMETER            :: ngaust   =  5                          ! Number of time gaussian integration points per elem; can use 5 or 6, but 5 is sufficient for 15 dof elem.
INTEGER                       :: ngd                                    ! Number of bounded(?) DOFs.     GUNJIT NOTE: Check later.
INTEGER, PARAMETER            :: ninteg   =  3                          ! Integer multiplicant which determines npsi ( npsi = nblade*ntelt*ninteg)
INTEGER                       :: nintn                                  ! Number of internal nodes.
INTEGER, PARAMETER            :: nmodes   =  6                          ! Number of normal modes used for stability.
INTEGER, PARAMETER            :: nmodet   =  6                          ! Number of normal modes used for trim.
INTEGER                       :: nndof                                  ! Number of nodal DOFs.
INTEGER, PARAMETER            :: nsect    =  3                          ! Order of time integration within each time element.
INTEGER                       :: nselt                                  ! Number of spatial elements.
INTEGER, PARAMETER            :: ntelt    = 1                           ! Number of time elements.
INTEGER, PARAMETER            :: tau      =  5.5                        !
INTEGER, PARAMETER            :: UnIn     = 1                           ! I/O unit for primary input file.
INTEGER, PARAMETER            :: UnOu     = 3                           ! I/O unit for output file.
INTEGER, PARAMETER            :: UnSP     = 2                           ! I/O unit for section properties input file.


   ! Derived parameters:

INTEGER, PARAMETER            :: npsi     = nblade*ntelt*ninteg         ! Number of azimuth locations = nblade*ntelt*ninteg

END MODULE Param
!=======================================================================
MODULE pbrng

USE Precision

REAL(ReKi)                    :: del3pb

INTEGER                       :: iartic

END MODULE pbrng
!=======================================================================
MODULE Pitch

USE Precision

REAL(ReKi)                    :: Ptch

END MODULE Pitch
!=======================================================================
MODULE plgeom

USE Precision

REAL(ReKi)                    :: b1vpph
REAL(ReKi)                    :: b1vpv
REAL(ReKi)                    :: b1wpph
REAL(ReKi)                    :: b1wpv
REAL(ReKi)                    :: b1wpw
REAL(ReKi)                    :: b2vph
REAL(ReKi)                    :: b2vpph
REAL(ReKi)                    :: b2wph
REAL(ReKi)                    :: b2wpph
REAL(ReKi)                    :: sphi
REAL(ReKi)                    :: svp
REAL(ReKi)                    :: swp
REAL(ReKi)                    :: tphi
REAL(ReKi)                    :: tv
REAL(ReKi)                    :: tvp
REAL(ReKi)                    :: tw
REAL(ReKi)                    :: twp

END MODULE plgeom
!=======================================================================
MODULE PLink

USE Precision

REAL(ReKi)                    :: nefs(4)

END MODULE PLink
!=======================================================================
MODULE PreCone

USE Precision

REAL(ReKi)                    :: btp

END MODULE PreCone
!=======================================================================
MODULE SftLnk

USE Precision

REAL(ReKi)                    :: tta  = 0.0                             ! pith link arm length / rotor radius
REAL(ReKi)                    :: ttb                                    !
REAL(ReKi)                    :: ttk  = 0.0                             ! pitch link stiffness, nondimensionaized
REAL(ReKi)                    :: ttp  = 0.0                             ! pitch link bend / radius
REAL(ReKi)                    :: ttx0 = 0.0                             ! horizontal offset of the lag shear pin / radius
REAL(ReKi)                    :: ttz  = 0.0                             ! vertical offset of the lag shear pin / radius
REAL(ReKi)                    :: xmp                                    !
REAL(ReKi)                    :: xpl  = 0.0                             !
REAL(ReKi)                    :: ypl  = 0.0                             !
REAL(ReKi)                    :: zpl  = 0.0                             !

INTEGER                       :: npin = 0                               ! 0: no shear lag pin, 1: shear lag pin exists

END MODULE SftLnk
!=======================================================================
MODULE Struc

USE Precision

REAL(ReKi), ALLOCATABLE       :: amass_den   (:)
REAL(ReKi), ALLOCATABLE       :: axial_stff  (:)
REAL(ReKi), ALLOCATABLE       :: cg_offst    (:)
REAL(ReKi), ALLOCATABLE       :: ea          (:)
REAL(ReKi), ALLOCATABLE       :: eac         (:)
REAL(ReKi), ALLOCATABLE       :: eb1         (:)
REAL(ReKi), ALLOCATABLE       :: eb2         (:)
REAL(ReKi), ALLOCATABLE       :: ec1         (:)
REAL(ReKi), ALLOCATABLE       :: ec2         (:)
REAL(ReKi), ALLOCATABLE       :: edge_stff   (:)
REAL(ReKi), ALLOCATABLE       :: eg          (:)
REAL(ReKi), ALLOCATABLE       :: eiy         (:)
REAL(ReKi), ALLOCATABLE       :: eiz         (:)
REAL(ReKi), ALLOCATABLE       :: flp_stff    (:)
REAL(ReKi), ALLOCATABLE       :: gay         (:)
REAL(ReKi), ALLOCATABLE       :: gaz         (:)
REAL(ReKi), ALLOCATABLE       :: gj          (:)
REAL(ReKi), ALLOCATABLE       :: rmas        (:)
REAL(ReKi), ALLOCATABLE       :: sc_offst    (:)
REAL(ReKi), ALLOCATABLE       :: sec_loc     (:)
REAL(ReKi), ALLOCATABLE       :: skm1        (:)
REAL(ReKi), ALLOCATABLE       :: skm2        (:)
REAL(ReKi), ALLOCATABLE       :: sq_km1      (:)
REAL(ReKi), ALLOCATABLE       :: sq_km2      (:)
REAL(ReKi), ALLOCATABLE       :: str_tw      (:)
REAL(ReKi), ALLOCATABLE       :: tc_offst    (:)
REAL(ReKi), ALLOCATABLE       :: tor_stff    (:)

INTEGER, PARAMETER            :: max_secs    = 50
INTEGER                       :: n_secs

CHARACTER(99)                 :: title

END MODULE Struc
!=======================================================================
MODULE Swept

INTEGER, ALLOCATABLE          :: ielrot   (:)                           ! Connectivity vector. The global dof number of the i'th elemental dof of the j'th element.
INTEGER                       :: netip    = 0

END MODULE Swept
!=======================================================================
MODULE TipDat

USE Precision

REAL(ReKi)                    :: cm_loc
REAL(ReKi)                    :: ixx_tip
REAL(ReKi)                    :: ixy_tip
REAL(ReKi)                    :: iyy_tip
REAL(ReKi)                    :: iyz_tip
REAL(ReKi)                    :: izx_tip
REAL(ReKi)                    :: izz_tip
REAL(ReKi)                    :: tip_mass

END MODULE TipDat
!=======================================================================
MODULE TowWires

USE Precision

REAL(ReKi)                    :: k_tower(2)
REAL(ReKi)                    :: th_wire(2)
REAL(ReKi)                    :: wire_stfness(2)

INTEGER                       :: beam_type
INTEGER                       :: n_attachments
INTEGER                       :: n_wires(2)
INTEGER                       :: node_attach(2)

END MODULE TowWires
!=======================================================================
MODULE TrimV

USE Precision

REAL(ReKi)                    :: th0
REAL(ReKi)                    :: th1c
REAL(ReKi)                    :: th1s

END MODULE TrimV
!=======================================================================
