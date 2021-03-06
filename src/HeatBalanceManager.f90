MODULE HeatBalanceManager

          ! Module containing the heat balance simulation routines
          ! calculation (initialization) routines

          ! MODULE INFORMATION:
          !       AUTHOR         Richard J. Liesen
          !       DATE WRITTEN   February 1998
          !       MODIFIED       November 1998, FW
          !       MODIFIED       April 1999, LKL
          !       MODIFIED       Dec 2006 DJS of PSU for ecoroof
          !       Added          Dec 2008 TH for thermochromic windows:
          !                       new subroutine CreateTCConstructions called by GetHeatBalanceInput
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS MODULE:
          ! To encapsulate the data and algorithms required to
          ! manage the heat balance simulation on the building.

          ! METHODOLOGY EMPLOYED:
          !

          ! REFERENCES:
          ! The heat balance method is outlined in the "Tarp Algorithms Manual"
          ! The methods are also summarized in many BSO Theses and papers.



          ! OTHER NOTES:
          ! This module was created from IBLAST subroutines
          !

          ! USE STATEMENTS:
          ! Use statements for data only modules
USE DataPrecisionGlobals
USE DataGlobals
USE DataEnvironment
USE DataHeatBalFanSys
USE DataHeatBalance
USE DataHeatBalSurface
USE DataRoomAirModel
USE DataIPShortCuts
USE DataInterfaces
          ! Use statements for access to subroutines in other modules
USE InputProcessor, ONLY: GetNumObjectsFound,GetObjectItem,SameString,FindItemInList,VerifyName,GetObjectItemNum
USE ScheduleManager, ONLY: GetScheduleIndex
USE DataSurfaces, ONLY: TotSurfaces,FrameDivider,FrameDividerProperties,CalcSolRefl,SurfaceWindow,StormWindow,TotStormWin,  &
                        DividedLite,Suspended,ShadingTransmittanceVaries
USE WindowManager, ONLY:W5LsqFit
USE DataContaminantBalance, ONLY: Contaminant, ZoneAirCO2, ZoneAirCO2Temp, ZoneAirCO2Avg, OutdoorCO2, &
                                  ZoneAirGC, ZoneAirGCTemp, ZoneAirGCAvg, OutdoorGC
USE ScheduleManager, ONLY: GetCurrentScheduleValue

IMPLICIT NONE         ! Enforce explicit typing of all variables

PRIVATE


  ! MODULE PARAMETER DEFINITIONS
CHARACTER(len=*), PARAMETER :: Blank=' '
CHARACTER(len=*), PARAMETER :: fmtA="(A)"

CHARACTER(len=*), DIMENSION(2), PARAMETER :: PassFail=(/'Fail','Pass'/)

  ! DERIVED TYPE DEFINITIONS
TYPE WarmupConvergence
  INTEGER,DIMENSION(4) :: PassFlag   = 2     ! one flag (1=Fail), (2=Pass) for each of the 4 conditions of convergence from
                                             ! warmup (PassFlag(1)=Max Temp, PassFlag(2)=Min Temp, PassFlag(3)=Max Heat Load
                                             ! PassFlag(4)=Max Cool Load)
  ! Following are stored test values for temperature and loads convergence
  REAL(r64) :: TestMaxTempValue      =0.0d0  ! Max Temperature convergence value=ABS(MaxTempPrevDay(ZoneNum)-MaxTempZone(ZoneNum))
  REAL(r64) :: TestMinTempValue      =0.0d0  ! Min Temperature convergence value=ABS(MinTempPrevDay(ZoneNum)-MinTempZone(ZoneNum))
  REAL(r64) :: TestMaxHeatLoadValue  =0.0d0  ! Max Heat Load convergence value=
                                             !  ABS((MaxHeatLoadZone(ZoneNum)-MaxHeatLoadPrevDay(ZoneNum))/MaxHeatLoadZone(ZoneNum))
  REAL(r64) :: TestMaxCoolLoadValue  =0.0d0  ! Max Cool Load convergence value=
                                             !  ABS((MaxCoolLoadZone(ZoneNum)-MaxCoolLoadPrevDay(ZoneNum))/MaxCoolLoadZone(ZoneNum))
END TYPE


  ! MODULE VARIABLE DECLARATIONS:

  !Real Variables for the Heat Balance Simulation
  !Variables used to determine warmup convergence
REAL(r64), ALLOCATABLE, DIMENSION(:) :: MaxCoolLoadPrevDay   !Max cooling load from the previous day
REAL(r64), ALLOCATABLE, DIMENSION(:) :: MaxCoolLoadZone      !Maximum zone cooling load from the current day
REAL(r64), ALLOCATABLE, DIMENSION(:) :: MaxHeatLoadPrevDay   !Max heating load from the previous day
REAL(r64), ALLOCATABLE, DIMENSION(:) :: MaxHeatLoadZone      !Maximum zone heating load from the current day
REAL(r64), ALLOCATABLE, DIMENSION(:) :: MaxTempPrevDay       !Max temperature from the previous day
REAL(r64), ALLOCATABLE, DIMENSION(:) :: MaxTempZone          !Maximum zone temperature from the current day
REAL(r64), ALLOCATABLE, DIMENSION(:) :: MinTempPrevDay       !Min temperature from the previous day
REAL(r64), ALLOCATABLE, DIMENSION(:) :: MinTempZone          !Minimum zone temperature from the current day

  !Variables used to report difference in temperature and load from the last two warmup days
REAL(r64), ALLOCATABLE, DIMENSION(:) :: WarmupTempDiff       !Temperature difference between the last two warmup days
REAL(r64), ALLOCATABLE, DIMENSION(:) :: WarmupLoadDiff       !Zone load differences between the last two warmup days
REAL(r64), ALLOCATABLE, DIMENSION(:) :: TempZoneSecPrevDay   !Zone air temperature from the second last warmup day
REAL(r64), ALLOCATABLE, DIMENSION(:) :: LoadZoneSecPrevDay   !Zone load from the second last warmup day
REAL(r64), ALLOCATABLE, DIMENSION(:) :: TempZonePrevDay      !Zone air temperature from the previous day
REAL(r64), ALLOCATABLE, DIMENSION(:) :: LoadZonePrevDay      !Zone load from the previuos day
REAL(r64), ALLOCATABLE, DIMENSION(:) :: TempZone             !Zone air temperature from the current warmup day
REAL(r64), ALLOCATABLE, DIMENSION(:) :: LoadZone             !Zone load from the current warmup day

REAL(r64), ALLOCATABLE, DIMENSION(:,:) :: TempZoneRpt        !Zone air temperature to report (average over all warmup days)
REAL(r64), ALLOCATABLE, DIMENSION(:)   :: TempZoneRptStdDev  !Zone air temperature to report (std dev over all warmup days)
REAL(r64), ALLOCATABLE, DIMENSION(:,:) :: LoadZoneRpt        !Zone load to report (average over all warmup days)
REAL(r64), ALLOCATABLE, DIMENSION(:)   :: LoadZoneRptStdDev  !Zone load to report (std dev over all warmup days)
REAL(r64), ALLOCATABLE, DIMENSION(:,:) :: MaxLoadZoneRpt     !Maximum zone load for reporting calcs
INTEGER :: CountWarmupDayPoints                              !Count of warmup timesteps (to achieve warmup)
TYPE(WarmUpConvergence), ALLOCATABLE, DIMENSION(:)  :: WarmupConvergenceValues

CHARACTER(len=MaxNameLength) :: CurrentModuleObject ! to assist in getting input

          ! Subroutine Specifications for the Heat Balance Module
          ! Driver Routines
PUBLIC  ManageHeatBalance

          ! Input reader routines for the module
PRIVATE GetHeatBalanceInput
PRIVATE GetProjectControlData
PRIVATE GetSiteAtmosphereData
PRIVATE GetMaterialData
PRIVATE GetWindowGlassSpectralData
PRIVATE ValidateMaterialRoughness
PRIVATE GetConstructData
PRIVATE GetBuildingData
PRIVATE CheckValidSimulationObjects
PRIVATE GetZoneData
PRIVATE GetFrameAndDividerData
PRIVATE SearchWindow5DataFile
PRIVATE SetupSimpleWindowGlazingSystem

          ! Initialization routines for module
PRIVATE InitHeatBalance
PRIVATE AllocateHeatBalArrays
PRIVATE SetStormWindowControl

          ! Record Keeping/Utility Routines for Module
PRIVATE RecKeepHeatBalance
PRIVATE CheckWarmupConvergence
PRIVATE ReportWarmupConvergence

          ! Reporting routines for module
PRIVATE ReportHeatBalance

CONTAINS

! MODULE SUBROUTINES:
!*************************************************************************

SUBROUTINE ManageHeatBalance

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Rick Strand
          !       DATE WRITTEN   January 1997
          !       MODIFIED       February 1998 Richard Liesen
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! This subroutine manages the heat balance method of calculating
          ! building thermal loads.  It is called from the SimulationManager
          ! at the time step level.  This driver manages the calls to all of
          ! the other modules, drivers, and simulation algorithms.

          ! METHODOLOGY EMPLOYED:
          ! The order of this routine was taken from HeatBalanceModule with routine
          !  and Data Structuring

          ! REFERENCES:
          ! Legacy code from (I)BLAST, subroutine SIMZGD.

          ! USE STATEMENTS:
  USE HeatBalanceSurfaceManager
  USE EMSManager , ONLY: ManageEMS, UpdateEMSTrendVariables
  USE DataGlobals, ONLY: emsCallFromEndZoneTimestepBeforeZoneReporting, emsCallFromEndZoneTimestepAfterZoneReporting, &
                          emsCallFromBeginNewEvironmentAfterWarmUp
  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
          ! na

          ! SUBROUTINE PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  LOGICAL,SAVE :: GetInputFlag=.true.

          ! FLOW:

    ! Get the heat balance input at the beginning of the simulation only
  IF (GetInputFlag) THEN
    CALL GetHeatBalanceInput  ! Obtains heat balance related parameters from input file
    GetInputFlag=.false.
  ENDIF

    ! These Inits will still have to be looked at as the routines are re-engineered further
  CALL InitHeatBalance  ! Initialize all heat balance related parameters

          ! Solve the zone heat balance by first calling the Surface Heat Balance Manager
          ! and then the Air Heat Balance Manager is called by the Surface Heat Balance
          ! Manager.  The order of execution is still important and the zone cannot
          ! go through any record keeping before the HVAC system has run because there
          ! may be a radiant system in the building which will require iteration between
          ! the HVAC system (called from the Air Heat Balance) and the zone (simulated
          ! in the Surface Heat Balance Manager).  In the future, this may be improved.
  CALL ManageSurfaceHeatBalance
  CALL ManageEMS(emsCallFromEndZoneTimestepBeforeZoneReporting)  ! EMS calling point
  CALL RecKeepHeatBalance   ! Do any heat balance related record keeping

 ! This call has been moved to the FanSystemModule and does effect the output file
 !   You do get a shift in the Air Handling System Summary for the building electric loads
 ! IF ((.NOT.WarmupFlag).AND.(DayOfSim.GT.0)) CALL RCKEEP  ! Do fan system accounting (to be moved later)

  CALL ReportHeatBalance    ! Manage heat balance reporting until the new reporting is in place

  CALL ManageEMS(emsCallFromEndZoneTimestepAfterZoneReporting) ! EMS calling point

  CALL UpdateEMSTrendVariables

  IF (WarmupFlag.AND.EndDayFlag) THEN

    CALL CheckWarmupConvergence
    IF (.NOT.WarmupFlag) THEN
      DayOfSim = 0   ! Reset DayOfSim if Warmup converged
      DayOfSimChr='0'

      CALL ManageEMS(emsCallFromBeginNewEvironmentAfterWarmUp) ! calling point
    END IF

  END IF

  IF (.not. WarmupFlag .and. EndDayFlag .and. DayOfSim == 1 .and. .not. DoingSizing) THEN
    CALL ReportWarmupConvergence
  ENDIF

  RETURN

END SUBROUTINE ManageHeatBalance


! Get Input Section of the Module
!******************************************************************************
SUBROUTINE GetHeatBalanceInput  ! Heat Balance Input Reader Manager

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Rick Strand
          !       DATE WRITTEN   September 1997
          !       MODIFIED       February 1998 Richard Liesen
          !                      November 1998 FW
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! This subroutine is the main driver for initializations within the
          ! heat balance.

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InternalHeatGains,      ONLY: ManageInternalHeatGains
  USE DataSystemVariables,    ONLY: DetailedSkyDiffuseAlgorithm

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
          ! na

          ! SUBROUTINE PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
    LOGICAL :: ErrorsFound = .false.   ! If errors detected in input
    LOGICAL :: ValidSimulationWithNoZones

          ! FLOW:


    CALL GetProjectControlData(ErrorsFound)

    CALL GetSiteAtmosphereData(ErrorsFound)

    CALL GetWindowGlassSpectralData(ErrorsFound)

    CALL GetMaterialData(ErrorsFound)    ! Read materials from input file/transfer from legacy data structure

    CALL GetFrameAndDividerData(ErrorsFound)

    CALL GetConstructData(ErrorsFound)   ! Read constructs from input file/transfer from legacy data structure

    CALL GetBuildingData(ErrorsFound)    ! Read building data from input file

    ! Added TH 1/9/2009 to create thermochromic window constructions
    CALL CreateTCConstructions(ErrorsFound)

    IF (TotSurfaces > 0 .and. NumOfZones == 0) THEN
      ValidSimulationWithNoZones=CheckValidSimulationObjects()
      IF (.not. ValidSimulationWithNoZones) THEN
        CALL ShowSevereError('GetHeatBalanceInput: There are surfaces in input but no zones found.  Invalid simulation.')
        ErrorsFound=.true.
      ENDIF
    ENDIF

    CALL CheckUsedConstructions(ErrorsFound)

    IF (ErrorsFound) THEN
      CALL ShowFatalError('Errors found in Building Input, Program Stopped')
    ENDIF

    ! following is done to "get internal heat gains" input so that lights are gotten before
    ! daylighting input
    CALL ManageInternalHeatGains(InitOnly=.true.)

  RETURN

END SUBROUTINE GetHeatBalanceInput

SUBROUTINE CheckUsedConstructions(ErrorsFound)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Linda Lawrie
          !       DATE WRITTEN   August 2011
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! Counts or details unused constructions.

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE General, ONLY: RoundSigDigits
  USE InputProcessor, ONLY: GetNumObjectsFound,GetObjectItem
  USE DataIPShortCuts

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  LOGICAL, INTENT(INOUT) :: ErrorsFound

          ! SUBROUTINE PARAMETER DEFINITIONS:
  INTEGER,PARAMETER :: NumConstrObjects=5
  CHARACTER(len=*), PARAMETER, DIMENSION(NumConstrObjects) :: ConstrObjects=  &
    (/'Pipe:Indoor                ',  &
      'Pipe:Outdoor               ',  &
      'Pipe:Underground           ',  &
      'GroundHeatExchanger:Surface',  &
      'DaylightingDevice:Tubular  '/)

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER :: Unused
  INTEGER :: Loop
  INTEGER :: NumObjects
  INTEGER :: NumAlphas
  INTEGER :: NumNumbers
  INTEGER :: Status
  INTEGER :: CNum
  INTEGER :: ONum
  LOGICAL :: InErrFlag   ! Preserve (no current use) the input status of ErrorsFound

  InErrFlag=ErrorsFound

  ! Needs to account for Pipe:HeatTransfer/indoor, etc constructions.
  DO ONum=1,NumConstrObjects
    NumObjects=GetNumObjectsFound(ConstrObjects(ONum))
    DO Loop=1,NumObjects
      CALL GetObjectItem(ConstrObjects(ONum),Loop,cAlphaArgs,NumAlphas,rNumericArgs,NumNumbers,Status)
      IF (ONum /= 5) THEN
        CNum=FindItemInList(cAlphaArgs(2),Construct%Name,TotConstructs)
      ELSE
        CNum=FindItemInList(cAlphaArgs(4),Construct%Name,TotConstructs)
      ENDIF
      IF (CNum == 0) CYCLE
      Construct(CNum)%IsUsed=.true.
    ENDDO
  ENDDO
  Unused=TotConstructs-Count(Construct%IsUsed)
  IF (Unused > 0) THEN
    IF (.not. DisplayExtraWarnings) THEN
      CALL ShowWarningError('CheckUsedConstructions: There are '//trim(RoundSigDigits(Unused))//  &
         ' nominally unused constructions in input.')
      CALL ShowContinueError('For explicit details on each unused construction, use Output:Diagnostics,DisplayExtraWarnings;')
    ELSE
      CALL ShowWarningError('CheckUsedConstructions: There are '//trim(RoundSigDigits(Unused))//  &
         ' nominally unused constructions in input.')
      CALL ShowContinueError('Each Unused construction is shown.')
      DO Loop=1,TotConstructs
        IF (Construct(Loop)%IsUsed) CYCLE
        CALL ShowMessage('Construction='//trim(Construct(Loop)%Name))
      ENDDO
    ENDIF
  ENDIF

  RETURN

END SUBROUTINE CheckUsedConstructions

FUNCTION CheckValidSimulationObjects() RESULT (ValidSimulation)

          ! FUNCTION INFORMATION:
          !       AUTHOR         Linda Lawrie
          !       DATE WRITTEN   July 2008
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS FUNCTION:
          ! If an input file presents with surfaces but no zones, there are certain objects
          ! that must be present for the simulation to be valid.  This check was necessitated by
          ! an input file that was entirely detached shading surfaces but no zones (and nothing else).
          ! Other objects include Solar Collectors, PV arrays.

          ! METHODOLOGY EMPLOYED:
          ! Check for specific objects that must be present for such a simulation to be valid.

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor, ONLY: GetNumObjectsFound

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

          ! FUNCTION ARGUMENT DEFINITIONS:
  LOGICAL :: ValidSimulation   ! True is other objects appear to make this a valid simulation.

          ! FUNCTION PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! FUNCTION LOCAL VARIABLE DECLARATIONS:
          ! na
  ValidSimulation = .false.
  IF (GetNumObjectsFound('SolarCollector:FlatPlate:Water') > 0) THEN
    ValidSimulation=.true.
  ELSEIF (GetNumObjectsFound('Generator:Photovoltaic') > 0) THEN
    ValidSimulation=.true.
  ELSEIF (GetNumObjectsFound('Generator:InternalCombustionEngine') > 0) THEN
    ValidSimulation=.true.
  ELSEIF (GetNumObjectsFound('Generator:CombustionTurbine') > 0) THEN
    ValidSimulation=.true.
  ELSEIF (GetNumObjectsFound('Generator:FuelCell') > 0) THEN
    ValidSimulation=.true.
  ELSEIF (GetNumObjectsFound('Generator:MicroCHP') > 0) THEN
    ValidSimulation=.true.
  ELSEIF (GetNumObjectsFound('Generator:MicroTurbine') > 0) THEN
    ValidSimulation=.true.
  ELSEIF (GetNumObjectsFound('Generator:WindTurbine') > 0) THEN
    ValidSimulation=.true.
  ENDIF

  RETURN

END FUNCTION CheckValidSimulationObjects

SUBROUTINE GetProjectControlData(ErrorsFound)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Linda Lawrie
          !       DATE WRITTEN   October 2004
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! This subroutine gets the project control data before the rest of the building data (such as
          ! materials) is obtained.

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! This routine gets the following objects:
          ! BUILDING
          ! INSIDE CONVECTION ALGORITHM
          ! OUTSIDE CONVECTION ALGORITHM
          ! SOLUTION ALGORITHM
          ! ASHRAE Handbook of Fundamentals, Chap 16, for the setting of Site Atmospheric defaults based
          !   on terrain.
          ! ZoneAirHeatBalanceAlgorithm, Added by L. Gu, 12/09
          ! ZoneAirContaminantBalance, Added by L. Gu, 06/10

          ! USE STATEMENTS:
  USE General, ONLY: RoundSigDigits
  USE DataSystemVariables,    ONLY: lMinimalShadowing

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  LOGICAL, INTENT(INOUT) :: ErrorsFound   ! Set to true if errors detected during getting data

          ! SUBROUTINE PARAMETER DEFINITIONS:
  CHARACTER(len=*), PARAMETER :: RoutineName='GetProjectControlData: '

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
    CHARACTER(len=MaxNameLength), DIMENSION(4) :: AlphaName
    REAL(r64), DIMENSION(5)         :: BuildingNumbers
    INTEGER                         :: NumAlpha, NumNumber
    INTEGER                         :: IOStat
    INTEGER                         :: NumObjects
    INTEGER                         :: TMP
    INTEGER :: NumEMPDMat
    INTEGER :: NumPCMat
    INTEGER :: NumVTCMat
    INTEGER :: NumHAMTMat1
    INTEGER :: NumHAMTMat2
    INTEGER :: NumHAMTMat3
    INTEGER :: NumHAMTMat4
    INTEGER :: NumHAMTMat5
    INTEGER :: NumHAMTMat6
    INTEGER :: SumHAMTMat
    LOGICAL :: msgneeded
    
    INTEGER :: DebugFile       =150 !RS: Debugging file denotion, hopefully this works.
    
    OPEN(unit=DebugFile,file='Debug.txt')    !RS: Debugging

   !Assign the values to the building data

   CurrentModuleObject='Building'
   NumObjects=GetNumObjectsFound(TRIM(CurrentModuleObject))

   IF (NumObjects > 0) THEN
     CALL GetObjectItem(TRIM(CurrentModuleObject),1,AlphaName,NumAlpha,BuildingNumbers,NumNumber,IOStat,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)
        ! Building Name (remove certain characters)
     BuildingName=AlphaName(1)
     TMP=INDEX(BuildingName,CHAR(1))
     DO WHILE (TMP /= 0)
       BuildingName(TMP:TMP)=','
       TMP=INDEX(BuildingName,CHAR(1))
     END DO
     TMP=INDEX(BuildingName,CHAR(2))
     DO WHILE (TMP /= 0)
       BuildingName(TMP:TMP)='!'
       TMP=INDEX(BuildingName,CHAR(2))
     END DO
     TMP=INDEX(BuildingName,CHAR(3))
     DO WHILE (TMP /= 0)
       BuildingName(TMP:TMP)='\'
       TMP=INDEX(BuildingName,CHAR(3))
     END DO
        ! Building Azimuth (no validation)
     BuildingAzimuth=MOD(BuildingNumbers(1),360.d0)
        ! Terrain
     IF (AlphaName(2) == 'COUNTRY' .or. AlphaName(2) == '1') THEN
       SiteWindExp = 0.14d0
       SiteWindBLHeight = 270.d0
       AlphaName(2)='Country'
     ELSEIF (AlphaName(2) == 'SUBURBS' .or. AlphaName(2) == '2' .or. AlphaName(2) == 'SUBURB') THEN
       SiteWindExp = 0.22d0
       SiteWindBLHeight = 370.d0
       AlphaName(2)='Suburbs'
     ELSEIF (AlphaName(2) == 'CITY' .or. AlphaName(2) == '3') THEN
       SiteWindExp = 0.33d0
       SiteWindBLHeight = 460.d0
       AlphaName(2)='City'
     ELSEIF (AlphaName(2) == 'OCEAN') THEN
       SiteWindExp = 0.10d0
       SiteWindBLHeight = 210.d0
       AlphaName(2)='Ocean'
     ELSEIF (AlphaName(2) == 'URBAN') THEN
       SiteWindExp = 0.22d0
       SiteWindBLHeight = 370.d0
       AlphaName(2)='Urban'
     ELSE
       CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//': '//TRIM(cAlphaFieldNames(2))//  &
          ' invalid='//TRIM(AlphaName(2)))
       SiteWindExp = 0.14d0
       SiteWindBLHeight = 270.d0
       AlphaName(2)=TRIM(AlphaName(2))//'-invalid'
       ErrorsFound=.true.
     ENDIF
        ! Loads Convergence Tolerance Value
     LoadsConvergTol=BuildingNumbers(2)
     IF (LoadsConvergTol <= 0.0) THEN
       CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//': '//TRIM(cNumericFieldNames(2))//  &
          ' value invalid, ['//  &
          TRIM(RoundSigDigits(LoadsConvergTol,3))//']')
       ErrorsFound=.true.
     ENDIF
        ! Temperature Convergence Tolerance Value
     TempConvergTol=BuildingNumbers(3)
     IF (TempConvergTol <= 0.0) THEN
       CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//': '//TRIM(cNumericFieldNames(2))//  &
          ' value invalid, ['//  &
          TRIM(RoundSigDigits(TempConvergTol,3))//']')
       ErrorsFound=.true.
     ENDIF
        ! Solar Distribution
     IF (AlphaName(3)(1:3) == 'MIN' .or. AlphaName(3) == '-1' .or. lMinimalShadowing) THEN
       SolarDistribution=MinimalShadowing
       AlphaName(3)='MinimalShadowing'
       CalcSolRefl = .FALSE.
     ELSEIF (AlphaName(3) == 'FULLEXTERIOR' .or. AlphaName(3) == '0') THEN
       SolarDistribution=FullExterior
       AlphaName(3)='FullExterior'
       CalcSolRefl = .FALSE.
     ELSEIF (AlphaName(3) == 'FULLINTERIORANDEXTERIOR' .or. AlphaName(3) == '1') THEN
       SolarDistribution=FullInteriorExterior
       AlphaName(3)='FullInteriorAndExterior'
       CalcSolRefl = .FALSE.
     ELSEIF (AlphaName(3) == 'FULLEXTERIORWITHREFLECTIONS') THEN
       SolarDistribution=FullExterior
       AlphaName(3)='FullExteriorWithReflectionsFromExteriorSurfaces'
       CalcSolRefl = .TRUE.
     ELSEIF (AlphaName(3) == 'FULLINTERIORANDEXTERIORWITHREFLECTIONS') THEN
       SolarDistribution=FullInteriorExterior
       AlphaName(3)='FullInteriorAndExteriorWithReflectionsFromExteriorSurfaces'
       CalcSolRefl = .TRUE.
     ELSE
       CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//': '//TRIM(cAlphaFieldNames(3))//  &
          ' invalid='//TRIM(AlphaName(3)))
       ErrorsFound=.true.
       AlphaName(3)=TRIM(AlphaName(3))//'-invalid'
     ENDIF
        ! Maximum Number of Warmup Days
     IF (.not. lNumericFieldBlanks(4)) THEN
       MaxNumberOfWarmupDays=BuildingNumbers(4)
       IF (MaxNumberOfWarmupDays <= 0) THEN
          CALL ShowSevereError(RoutineName//TRIM(CurrentModuleObject)//': '//TRIM(cNumericFieldNames(2))//  &
            ' invalid, ['//  &
            TRIM(RoundSigDigits(MaxNumberOfWarmupDays))//'], '//  &
            trim(RoundSigDIgits(DefaultMaxNumberOfWarmupDays))//' will be used')
         MaxNumberOfWarmupDays=DefaultMaxNumberOfWarmupDays
       ENDIF
     ELSE
       MaxNumberOfWarmupDays=DefaultMaxNumberOfWarmupDays
     ENDIF
        ! Minimum Number of Warmup Days
     IF (.not. lNumericFieldBlanks(5)) THEN
       MinNumberOfWarmupDays=BuildingNumbers(5)
       IF (MinNumberOfWarmupDays <= 0) THEN
         CALL ShowWarningError(RoutineName//TRIM(CurrentModuleObject)//': '//TRIM(cNumericFieldNames(3))//  &
            ' invalid, ['//  &
            TRIM(RoundSigDigits(MinNumberOfWarmupDays))//'], '//  &
            trim(RoundSigDIgits(DefaultMinNumberOfWarmupDays))//' will be used')
         MinNumberOfWarmupDays=DefaultMinNumberOfWarmupDays
       ENDIF
     ELSE
       MinNumberOfWarmupDays=DefaultMinNumberOfWarmupDays
     ENDIF
     IF (MinNumberOfWarmupDays > MaxNumberOfWarmupDays) THEN
       CALL ShowWarningError(RoutineName//TRIM(CurrentModuleObject)//': '//TRIM(cNumericFieldNames(2))//  &
          ' is greater than ['//  &
          TRIM(RoundSigDigits(MaxNumberOfWarmupDays))//'], '//TRIM(RoundSigDigits(MinNumberOfWarmupDays))// &
          ' will be used.')
       MaxNumberOfWarmupDays=MinNumberOfWarmupDays
     ENDIF
     IF (MinNumberOfWarmupDays < 6) THEN
       CALL ShowWarningError(RoutineName//TRIM(CurrentModuleObject)//': '//TRIM(cNumericFieldNames(2))//  &
          ' potentially invalid. '//  &
          'Experience has shown that most files will converge within '//TRIM(RoundSigDigits(DefaultMaxNumberOfWarmupDays))//  &
          ' warmup days. ')
       CALL ShowContinueError('...Choosing less than '//TRIM(RoundSigDigits(DefaultMinNumberOfWarmupDays))//  &
          ' warmup days may have adverse effects on the simulation results, '//  &
          'particularly design day simulations. ')
       CALL ShowContinueError('...Users should only alter this default if they are certain that '// &
                 'less than '//TRIM(RoundSigDigits(DefaultMinNumberOfWarmupDays))//  &
                ' warmup days is appropriate for a particular file. ')
       CALL ShowContinueError('...Verify that convergence to desired results are achieved. You can report values'//  &
                 ' during warmup days to ascertain convergence.')
     ENDIF
   ELSE
     CALL ShowSevereError(RoutineName//' A '//TRIM(CurrentModuleObject)//' Object must be entered.')
     ErrorsFound=.true.
     BuildingName='NOT ENTERED'
     AlphaName(2)='NOT ENTERED'
     AlphaName(3)='NOT ENTERED'
     MaxNumberOfWarmupDays=DefaultMaxNumberOfWarmupDays
     MinNumberOfWarmupDays=DefaultMinNumberOfWarmupDays
   ENDIF

! Write Building Information to the initialization output file
   Write(OutputFileInits,721)

721 Format('! <Building Information>, Building Name,North Axis {deg},Terrain, ', &
          ' Loads Convergence Tolerance Value,Temperature Convergence Tolerance Value, ', &
          ' Solar Distribution,Maximum Number of Warmup Days,Minimum Number of Warmup Days')

   Write(OutputFileInits,720) TRIM(BuildingName),trim(RoundSigDigits(BuildingAzimuth,3)),TRIM(AlphaName(2)),  &
          trim(RoundSigDigits(LoadsConvergTol,5)),trim(RoundSigDigits(TempConvergTol,5)),  &
          TRIM(AlphaName(3)),trim(RoundSigDigits(MaxNumberOfWarmupDays)),trim(RoundSigDigits(MinNumberOfWarmupDays))
720 Format(' Building Information',8(',',A))
   ! Above should be validated...

   CurrentModuleObject='SurfaceConvectionAlgorithm:Inside'
   NumObjects=GetNumObjectsFound(TRIM(CurrentModuleObject))
   IF (NumObjects > 0) THEN
     CALL GetObjectItem(TRIM(CurrentModuleObject),1,AlphaName,NumAlpha,BuildingNumbers,NumNumber,IOStat,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

     SELECT CASE (AlphaName(1))

       CASE ('SIMPLE')
         DefaultInsideConvectionAlgo=ASHRAESimple
         AlphaName(1)='Simple'

       CASE ('TARP','DETAILED')
         DefaultInsideConvectionAlgo=ASHRAETARP
         IF (AlphaName(1) == 'DETAILED') THEN
           CALL ShowSevereError('GetInsideConvectionAlgorithm: Deprecated value for '//TRIM(CurrentModuleObject)//', '// &
             'defaulting to TARP, entered value='//TRIM(AlphaName(1)))
         ENDIF
         AlphaName(1)='TARP'

       CASE ('CEILINGDIFFUSER')
         DefaultInsideConvectionAlgo=CeilingDiffuser
         AlphaName(1)='CeilingDiffuser'

       CASE ('TROMBEWALL')
         DefaultInsideConvectionAlgo=TrombeWall
         CALL ShowSevereError('GetInsideConvectionAlgorithm: TrombeWall has been used as a global definition.'//  &
                              ' This is a zone oriented value.  Will be illegal in the future.')
         AlphaName(1)='TrombeWall'

       CASE ('ADAPTIVECONVECTIONALGORITHM')
         DefaultInsideConvectionAlgo=AdaptiveConvectionAlgorithm
          AlphaName(1)='AdaptiveConvectionAlgorithm'

       CASE DEFAULT
         CALL ShowWarningError('GetInsideConvectionAlgorithm: Invalid value for '//TRIM(CurrentModuleObject)//', '// &
           'defaulting to TARP, invalid value='//TRIM(AlphaName(1)))
         DefaultInsideConvectionAlgo=ASHRAETARP
         AlphaName(1)='TARP'

     END SELECT
   ELSE
     ! default value, if not specified
      DefaultInsideConvectionAlgo=ASHRAETARP
      AlphaName(1)='TARP'

   ENDIF
   Write(OutputFileInits,722) TRIM(AlphaName(1))
722 Format('! <Inside Convection Algorithm>, Algorithm {Simple | TARP | CeilingDiffuser | AdaptiveConvectionAlgorithm}',/,  &
           'Inside Convection Algorithm,',A)

   !Get only the first (if more were input)
   CurrentModuleObject='SurfaceConvectionAlgorithm:Outside'
   NumObjects=GetNumObjectsFound(TRIM(CurrentModuleObject))
   IF (NumObjects > 0) THEN
     CALL GetObjectItem('SurfaceConvectionAlgorithm:Outside',1,AlphaName,NumAlpha,BuildingNumbers,NumNumber,IOStat,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)
     SELECT CASE (AlphaName(1))

       CASE ('SIMPLECOMBINED', 'SIMPLE')
         DefaultOutsideConvectionAlgo=ASHRAESimple
         IF (AlphaName(1) == 'SIMPLE') THEN
           CALL ShowSevereError('GetOutsideConvectionAlgorithm: Deprecated value for '//TRIM(CurrentModuleObject)//', '// &
             'defaulting to SimpleCombined, entered value='//TRIM(AlphaName(1)))
         ENDIF
         AlphaName(1)='SimpleCombined'

       CASE ('TARP', 'DETAILED', 'BLAST')
         DefaultOutsideConvectionAlgo=ASHRAETARP
         IF (AlphaName(1) == 'DETAILED') THEN
           CALL ShowSevereError('GetOutsideConvectionAlgorithm: Deprecated value for '//TRIM(CurrentModuleObject)//', '// &
             'defaulting to TARP, entered value='//TRIM(AlphaName(1)))
         ENDIF
         IF (AlphaName(1) == 'BLAST') THEN
           CALL ShowSevereError('GetOutsideConvectionAlgorithm: Deprecated value for '//TRIM(CurrentModuleObject)//', '// &
             'defaulting to TARP, entered value='//TRIM(AlphaName(1)))
         ENDIF
         AlphaName(1)='TARP'

       CASE ('MOWITT')
         DefaultOutsideConvectionAlgo=MoWittHcOutside
         AlphaName(1)='MoWitt'

       CASE ('DOE-2','DOE2')
         DefaultOutsideConvectionAlgo=DOE2HcOutside
         IF (AlphaName(1) == 'DOE2') THEN
           CALL ShowSevereError('GetOutsideConvectionAlgorithm: Deprecated value for '//TRIM(CurrentModuleObject)//', '// &
             'defaulting to DOE-2, entered value='//TRIM(AlphaName(1)))
         ENDIF
         AlphaName(1)='DOE-2'

       CASE ('ADAPTIVECONVECTIONALGORITHM')
         DefaultOutsideConvectionAlgo=AdaptiveConvectionAlgorithm
         AlphaName(1)='AdaptiveConvectionAlgorithm'

       CASE DEFAULT
         CALL ShowWarningError('GetOutsideConvectionAlgorithm: Invalid value for '//TRIM(CurrentModuleObject)//', '// &
           'defaulting to DOE-2, invalid value='//TRIM(AlphaName(1)))
         DefaultOutsideConvectionAlgo=DOE2HcOutside
         AlphaName(1)='DOE-2'

     END SELECT
   ELSE
     ! default value, if not specified
     DefaultOutsideConvectionAlgo=DOE2HcOutside
     AlphaName(1)='DOE-2'

   ENDIF

   Write(OutputFileInits,723) TRIM(AlphaName(1))
723 Format('! <Outside Convection Algorithm>, ',  &
       'Algorithm {SimpleCombined | TARP | MoWitt | DOE-2 | AdaptiveConvectionAlgorithm}', /,  &
           'Outside Convection Algorithm,',A)

   CurrentModuleObject='HeatBalanceAlgorithm'
   NumObjects=GetNumObjectsFound(TRIM(CurrentModuleObject))
   IF (NumObjects > 0) THEN
     CALL GetObjectItem(TRIM(CurrentModuleObject),1,AlphaName,NumAlpha,BuildingNumbers,NumNumber,IOStat,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)
     SELECT CASE (AlphaName(1))
       !The default is CTF = 0.  Then the moisture solution is EMPD =2
       CASE ('CONDUCTIONTRANSFERFUNCTION','DEFAULT','CTF')
         SolutionAlgo = UseCTF
         AlphaName(1)='CTF - Conduction Transfer Function'

       CASE ('MOISTUREPENETRATIONDEPTHCONDUCTIONTRANSFERFUNCTION','EMPD')
         SolutionAlgo = UseEMPD
         AlphaName(1)='EMPD - Effective Moisture Penetration Depth'

       CASE ('CONDUCTIONFINITEDIFFERENCE','CONDFD','CONDUCTIONFINITEDIFFERENCEDETAILED')
         SolutionAlgo = UseCondFD
         AlphaName(1)='CONDFD - Conduction Finite Difference'
         IF (NumOfTimeStepInHour < 20) THEN
           CALL ShowSevereError('GetSolutionAlgorithm: '//TRIM(CurrentModuleObject)//' '//TRIM(cAlphaFieldNames(1))//  &
              ' is Conduction Finite Difference but Number of TimeSteps in Hour < 20, Value is '//   &
              TRIM(RoundSigDigits(NumOfTimeStepInHour))//'.')
           CALL ShowContinueError('...Suggested minimum number of time steps in hour for '//  &
                  'Conduction Finite Difference solutions is 20.'//  &
                  ' Errors or inaccurate calculations may occur.')
         ENDIF

       CASE ('COMBINEDHEATANDMOISTUREFINITEELEMENT','HAMT')
         SolutionAlgo = UseHAMT
         AlphaName(1)='HAMT - Combined Heat and Moisture Transfer Finite Element'

       CASE ('CONDUCTIONFINITEDIFFERENCESIMPLIFIED')
         SolutionAlgo = UseCondFDSimple
         AlphaName(1)='CONDFDSimple - Conduction Finite Difference Simplified'
         IF (NumOfTimeStepInHour < 20) THEN
           CALL ShowWarningError('GetSolutionAlgorithm: '//TRIM(CurrentModuleObject)//' '//TRIM(cAlphaFieldNames(1))//  &
              ' is Conduction Finite Difference (Simplified)  but Number of TimeSteps in Hour < 20, Value is '//   &
              TRIM(RoundSigDigits(NumOfTimeStepInHour))//'.')
           CALL ShowContinueError('...Suggested minimum number of time steps in hour for '//  &
                  'Conduction Finite Difference solutions is 20.'//  &
                  ' Errors may occur.')
         ENDIF

       CASE DEFAULT
         SolutionAlgo = UseCTF
         AlphaName(1)='CTF - Conduction Transfer Function'
     END SELECT

     IF (NumNumber > 0) THEN
       MaxSurfaceTempLimit=BuildingNumbers(1)
       MaxSurfaceTempLimitBeforeFatal=MaxSurfaceTempLimit*2.5d0
       IF (MaxSurfaceTempLimit < MinSurfaceTempLimit) THEN
       ELSEIF (MaxSurfaceTempLimit < 0.0d0) THEN
         MaxSurfaceTempLimit=DefaultSurfaceTempLimit
         MaxSurfaceTempLimitBeforeFatal=MaxSurfaceTempLimit*2.5d0
       ENDIF
     ENDIF

     IF ( .NOT. lNumericFieldBlanks(2)) THEN
       LowHConvLimit = BuildingNumbers(2)
     ENDIF
     IF ( .NOT. lNumericFieldBlanks(3)) THEN
       HighHConvLimit = BuildingNumbers(3)
     ENDIF

   ELSE
     SolutionAlgo = UseCTF
     AlphaName(1)='ConductionTransferFunction'
     MaxSurfaceTempLimit=DefaultSurfaceTempLimit
     MaxSurfaceTempLimitBeforeFatal=MaxSurfaceTempLimit*2.5d0
   ENDIF

   NumEMPDMat=GetNumObjectsFound('MaterialProperty:MoisturePenetrationDepth:Settings')
   NumPCMat=GetNumObjectsFound('MaterialProperty:PhaseChange') ! needs detailed algo
   NumVTCMat=GetNumObjectsFound('MaterialProperty:VariableThermalConductivity')
   NumHAMTMat1=GetNumObjectsFound('MaterialProperty:HeatAndMoistureTransfer:Settings')
   NumHAMTMat2=GetNumObjectsFound('MaterialProperty:HeatAndMoistureTransfer:SorptionIsotherm')
   NumHAMTMat3=GetNumObjectsFound('MaterialProperty:HeatAndMoistureTransfer:Suction')
   NumHAMTMat4=GetNumObjectsFound('MaterialProperty:HeatAndMoistureTransfer:Redistribution')
   NumHAMTMat5=GetNumObjectsFound('MaterialProperty:HeatAndMoistureTransfer:Diffusion')
   NumHAMTMat6=GetNumObjectsFound('MaterialProperty:HeatAndMoistureTransfer:ThermalConductivity')
   SumHAMTMat=NumHAMTMat1+NumHAMTMat2+NumHAMTMat3+NumHAMTMat4+NumHAMTMat5+NumHAMTMat6
   msgneeded=.false.

   IF (SolutionAlgo == UseCTF) THEN
     IF (NumEMPDMat > 0) THEN
       !CALL ShowWarningError(trim(CurrentModuleObject)//'="'//trim(AlphaName(1))//'" but input file includes '// &
       !   trim(RoundSigDigits(NumEMPDMat))//' MaterialProperty:MoisturePenetrationDepth:Settings objects.') !RS: Secret Search String
       WRITE(DebugFile,*) TRIM(CurrentModuleObject)//'="'//TRIM(AlphaName(1))//'" but input file includes '// &
        TRIM(RoundSigDigits(NumEMPDMat))//' MaterialProperty:MoisturePenetrationDepth:Settings objects.'
       msgneeded=.true.
     ENDIF
     IF (NumPCMat > 0 .or. NumVTCMat > 0) THEN
       IF (NumPCMat > 0) THEN
         !CALL ShowWarningError(trim(CurrentModuleObject)//'="'//trim(AlphaName(1))//'" but input file includes '// &
         !   trim(RoundSigDigits(NumPCMat))//' MaterialProperty:PhaseChange objects.')   !RS: Secret Search String
         WRITE(DebugFile,*) TRIM(CurrentModuleObject)//'="'//TRIM(AlphaName(1))//'" but input file includes '// &
            TRIM(RoundSigDigits(NumPCMat))//' MaterialProperty:PhaseChange objects.'
         msgneeded=.true.
       ENDIF
       IF (NumVTCMat > 0) THEN
         !CALL ShowWarningError(trim(CurrentModuleObject)//'="'//trim(AlphaName(1))//'" but input file includes '// &
         !   trim(RoundSigDigits(NumVTCMat))//' MaterialProperty:VariableThermalConductivity objects.')  !RS: Secret Search String
         WRITE(DebugFile,*) TRIM(CurrentModuleObject)//'="'//TRIM(AlphaName(1))//'" but input file includes '// &
            TRIM(RoundSigDigits(NumVTCMat))//' MaterialProperty:VariableThermalConductivity objects.'
         msgneeded=.true.
       ENDIF
     ENDIF
     IF (SumHAMTMat > 0) THEN
       !CALL ShowWarningError(trim(CurrentModuleObject)//'="'//trim(AlphaName(1))//'" but input file includes '// &
       !   trim(RoundSigDigits(SumHAMTMat))//' MaterialProperty:HeatAndMoistureTransfer:* objects.') !RS: Secret Search String
       WRITE(DebugFile,*) TRIM(CurrentModuleObject)//'="'//TRIM(AlphaName(1))//'" but input file includes '// &
        TRIM(RoundSigDigits(SumHAMTMat))//' MaterialProperty:HeatAndMoistureTransfer:* objects.'
       msgneeded=.true.
     ENDIF
   ELSEIF (SolutionAlgo == UseEMPD) THEN
     IF (NumPCMat > 0 .or. NumVTCMat > 0) THEN
       IF (NumPCMat > 0) THEN
         CALL ShowWarningError(trim(CurrentModuleObject)//'="'//trim(AlphaName(1))//'" but input file includes '// &
            trim(RoundSigDigits(NumPCMat))//' MaterialProperty:PhaseChange objects.')
         msgneeded=.true.
       ENDIF
       IF (NumVTCMat > 0) THEN
         CALL ShowWarningError(trim(CurrentModuleObject)//'="'//trim(AlphaName(1))//'" but input file includes '// &
            trim(RoundSigDigits(NumVTCMat))//' MaterialProperty:VariableThermalConductivity objects.')
         msgneeded=.true.
       ENDIF
     ENDIF
     IF (SumHAMTMat > 0) THEN
       CALL ShowWarningError(trim(CurrentModuleObject)//'="'//trim(AlphaName(1))//'" but input file includes '// &
          trim(RoundSigDigits(SumHAMTMat))//' MaterialProperty:HeatAndMoistureTransfer:* objects.')
       msgneeded=.true.
     ENDIF
   ELSEIF (SolutionAlgo == UseHAMT) THEN
     IF (NumEMPDMat > 0) THEN
       CALL ShowWarningError(trim(CurrentModuleObject)//'="'//trim(AlphaName(1))//'" but input file includes '// &
          trim(RoundSigDigits(NumEMPDMat))//' MaterialProperty:MoisturePenetrationDepth:Settings objects.')
       msgneeded=.true.
     ENDIF
     IF (NumPCMat > 0 .or. NumVTCMat > 0) THEN
       IF (NumPCMat > 0) THEN
         CALL ShowWarningError(trim(CurrentModuleObject)//'="'//trim(AlphaName(1))//'" but input file includes '// &
            trim(RoundSigDigits(NumPCMat))//' MaterialProperty:PhaseChange objects.')
         msgneeded=.true.
       ENDIF
       IF (NumVTCMat > 0) THEN
         CALL ShowWarningError(trim(CurrentModuleObject)//'="'//trim(AlphaName(1))//'" but input file includes '// &
            trim(RoundSigDigits(NumVTCMat))//' MaterialProperty:VariableThermalConductivity objects.')
         msgneeded=.true.
       ENDIF
     ENDIF
   ELSEIF (SolutionAlgo == UseCondFD .or. SolutionAlgo == UseCondFDSimple) THEN
     IF (NumEMPDMat > 0) THEN
       CALL ShowWarningError(trim(CurrentModuleObject)//'="'//trim(AlphaName(1))//'" but input file includes '// &
          trim(RoundSigDigits(NumEMPDMat))//' MaterialProperty:MoisturePenetrationDepth:Settings objects.')
       msgneeded=.true.
     ENDIF
   ENDIF
   IF (msgneeded) THEN
     !CALL ShowContinueError('Previous materials will be ignored due to HeatBalanceAlgorithm choice.')   !RS: Secret Search String
     WRITE(DebugFile,*) 'Previous materials will be ignored due to HeatBalanceAlgorithm choice.'
   ENDIF

   IF (SolutionAlgo == UseEMPD) THEN
     IF (NumEMPDMat == 0) THEN
       CALL ShowWarningError(trim(CurrentModuleObject)//'="'//trim(AlphaName(1))//'" but input file includes '// &
          trim(RoundSigDigits(NumEMPDMat))//' MaterialProperty:MoisturePenetrationDepth:Settings objects.')
       msgneeded=.true.
     ENDIF
   ELSEIF (SolutionAlgo == UseHAMT) THEN
     IF (SumHAMTMat == 0) THEN
       CALL ShowWarningError(trim(CurrentModuleObject)//'="'//trim(AlphaName(1))//'" but input file includes '// &
          trim(RoundSigDigits(NumEMPDMat))//' MaterialProperty:HeatAndMoistureTransfer:* objects.')
       msgneeded=.true.
     ENDIF
   ENDIF
   IF (msgneeded) THEN
     !CALL ShowContinueError('Certain materials are necessary to achieve proper results with the HeatBalanceAlgorithm choice.')  !RS: Secret Search String
     WRITE(DebugFile,*) 'Certain materials are necessary to achieve proper results with the HeatBalanceAlgorithm choice.'
   ENDIF

     ! Write Solution Algorithm to the initialization output file for User Verification
   WRITE(OutputFileInits,fmtA) '! <Solution Algorithm>, Value {CTF - ConductionTransferFunction | '//  &
           'EMPD - MoisturePenetrationDepthConductionTransferFunction | '// &
           'CONDFD - ConductionFiniteDifference | CONDFDSimple - Conduction Finite Difference Simplified | '// &
           'HAMT - CombinedHeatAndMoistureFiniteElement} - Description,Inside Surface Max Temperature Limit{C}'
   WRITE(OutputFileInits,725) TRIM(AlphaName(1)),TRIM(RoundSigDigits(MaxSurfaceTempLimit,0))
725 FORMAT('Solution Algorithm, ',A,',',A)


   WRITE(OutputFileInits,724)
724 FORMAT('! <Sky Radiance Distribution>, Value {Anisotropic}',/,  &
           'Sky Radiance Distribution,Anisotropic')


   CurrentModuleObject='Compliance:Building'
   NumObjects=GetNumObjectsFound(TRIM(CurrentModuleObject))

   IF (NumObjects > 0) THEN
     CALL GetObjectItem(TRIM(CurrentModuleObject),1,AlphaName,NumAlpha,BuildingNumbers,NumNumber,IOStat,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)
        ! Building Rotation for Appendix G
     BuildingRotationAppendixG = MOD(BuildingNumbers(1),360.d0)
   END IF

   ! A new object is added by L. Gu, 12/09
   CurrentModuleObject='ZoneAirHeatBalanceAlgorithm'
   NumObjects=GetNumObjectsFound(TRIM(CurrentModuleObject))
   IF (NumObjects > 0) THEN
     CALL GetObjectItem(TRIM(CurrentModuleObject),1,AlphaName,NumAlpha,BuildingNumbers,NumNumber,IOStat,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)
     If (NumAlpha > 0) Then
       SELECT CASE (AlphaName(1))
         CASE ('3RDORDERBACKWARDDIFFERENCE','THIRDORDERBACKWARDDIFFERENCE')
           ZoneAirSolutionAlgo = Use3rdOrder
           AlphaName(1)='ThirdOrderBackwardDifference'
         CASE ('ANALYTICALSOLUTION')
           ZoneAirSolutionAlgo = UseAnalyticalSolution
           AlphaName(1)='AnalyticalSolution'
         CASE ('EULERMETHOD')
           ZoneAirSolutionAlgo = UseEulerMethod
           AlphaName(1)='EulerMethod'
         CASE DEFAULT
           ZoneAirSolutionAlgo = Use3rdOrder
           AlphaName(1)='ThirdOrderBackwardDifference'
           CALL ShowWarningError(TRIM(CurrentModuleObject)//': Invalid input of '//TRIM(cAlphaFieldNames(1))//  &
              '. The default choice is assigned = '//TRIM(AlphaName(1)))
           CALL ShowContinueError('Valid choices are: ThirdOrderBackwardDifference, AnalyticalSolution, or EulerMethod.')
       END SELECT
     End If
   ELSE
     ZoneAirSolutionAlgo = Use3rdOrder
     AlphaName(1)='ThirdOrderBackwardDifference'
   ENDIF

     ! Write Solution Algorithm to the initialization output file for User Verification
   WRITE(OutputFileInits,726)
   WRITE(OutputFileInits,727) TRIM(AlphaName(1))
726 FORMAT('! <Zone Air Solution Algorithm>, Value {ThirdOrderBackwardDifference | AnalyticalSolution | EulerMethod}')
727 FORMAT(' Zone Air Solution Algorithm, ',A)

   ! A new object is added by L. Gu, 06/10
   CurrentModuleObject='ZoneAirContaminantBalance'
   NumObjects=GetNumObjectsFound(TRIM(CurrentModuleObject))
   IF (NumObjects > 0) THEN
     CALL GetObjectItem(TRIM(CurrentModuleObject),1,AlphaName,NumAlpha,BuildingNumbers,NumNumber,IOStat,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)
     IF (NumAlpha > 0) THEN
       SELECT CASE (AlphaName(1))
         CASE ('YES')
           Contaminant%CO2Simulation = .TRUE.
           Contaminant%SimulateContaminants = .TRUE.
         CASE ('NO')
           Contaminant%CO2Simulation = .FALSE.
         CASE DEFAULT
           Contaminant%CO2Simulation = .FALSE.
           AlphaName(1)='NO'
           CALL ShowWarningError(TRIM(CurrentModuleObject)//': Invalid input of '//TRIM(cAlphaFieldNames(1))//  &
              '. The default choice is assigned = NO')
       END SELECT
     END IF
     IF (NumAlpha .EQ. 1 .AND. Contaminant%CO2Simulation) THEN
       IF (Contaminant%CO2Simulation) THEN
         CALL ShowSevereError(TRIM(CurrentModuleObject)//', '//TRIM(cAlphaFieldNames(2))//' is required and not given.')
         ErrorsFound=.true.
       END IF
     ELSEIF (NumAlpha > 1 .AND. Contaminant%CO2Simulation) THEN
       Contaminant%CO2OutdoorSchedPtr = GetScheduleIndex(AlphaName(2))
       IF (Contaminant%CO2OutdoorSchedPtr == 0) THEN
         CALL ShowSevereError(TRIM(CurrentModuleObject)//', '//TRIM(cAlphaFieldNames(2))//' not found: '//TRIM(AlphaName(2)))
         ErrorsFound=.true.
       ENDIF
     END IF
     If (NumAlpha > 2) THEN
       SELECT CASE (AlphaName(3))
         CASE ('YES')
           Contaminant%GenericContamSimulation = .TRUE.
           IF (.NOT. Contaminant%CO2Simulation) Contaminant%SimulateContaminants = .TRUE.
         CASE ('NO')
           Contaminant%GenericContamSimulation = .FALSE.
         CASE DEFAULT
           Contaminant%GenericContamSimulation = .FALSE.
           AlphaName(3)='NO'
           CALL ShowWarningError(TRIM(CurrentModuleObject)//': Invalid input of '//TRIM(cAlphaFieldNames(3))//  &
              '. The default choice is assigned = NO')
       END SELECT
       IF (NumAlpha .EQ. 3 .AND. Contaminant%GenericContamSimulation) THEN
         IF (Contaminant%GenericContamSimulation) THEN
           CALL ShowSevereError(TRIM(CurrentModuleObject)//', '//TRIM(cAlphaFieldNames(4))//' is required and not given.')
           ErrorsFound=.true.
         END IF
       ELSEIF (NumAlpha > 3 .AND. Contaminant%GenericContamSimulation) THEN
         Contaminant%GenericContamOutdoorSchedPtr = GetScheduleIndex(AlphaName(4))
         IF (Contaminant%GenericContamOutdoorSchedPtr == 0) THEN
           CALL ShowSevereError(TRIM(CurrentModuleObject)//', '//TRIM(cAlphaFieldNames(4))//' not found: '//TRIM(AlphaName(4)))
           ErrorsFound=.true.
         ENDIF
       END IF
     END IF
   ELSE
     Contaminant%SimulateContaminants = .FALSE.
     Contaminant%CO2Simulation = .FALSE.
     Contaminant%GenericContamSimulation = .FALSE.
     AlphaName(1)='NO'
     AlphaName(3)='NO'
   ENDIF

   WRITE(OutputFileInits,728)
   IF (Contaminant%SimulateContaminants .AND. Contaminant%CO2Simulation) THEN
     WRITE(OutputFileInits,730) 'Yes',TRIM(AlphaName(1))
   ELSE
     WRITE(OutputFileInits,730) 'No','N/A'
   END IF
728 FORMAT('! <Zone Air Contaminant Balance Simulation>, Simulation {Yes/No}, Carbon Dioxide Concentration')
730 FORMAT(' Zone Air Carbon Dioxide Balance Simulation, ',A,',',A)

   WRITE(OutputFileInits,729)
   IF (Contaminant%SimulateContaminants .AND. Contaminant%GenericContamSimulation) THEN
     WRITE(OutputFileInits,731) 'Yes',TRIM(AlphaName(3))
   ELSE
     WRITE(OutputFileInits,731) 'No','N/A'
   END IF
729 FORMAT('! <Zone Air Contaminant Balance Simulation>, Simulation {Yes/No}, Generic Contaminant Concentration')
731 FORMAT(' Zone Air Generic Contaminant Balance Simulation, ',A,',',A)


  RETURN

END SUBROUTINE GetProjectControlData


SUBROUTINE GetSiteAtmosphereData(ErrorsFound)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Peter Graham Ellis
          !       DATE WRITTEN   January 2006
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! Reads the input data for the SITE ATMOSPHERIC VARIATION object.

          ! METHODOLOGY EMPLOYED:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor, ONLY: GetNumObjectsFound, GetObjectItem
  USE General, ONLY: RoundSigDigits

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  LOGICAL, INTENT(INOUT)                    :: ErrorsFound

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER                                   :: NumObjects
  INTEGER                                   :: NumAlphas ! Number of elements in the alpha array
  INTEGER                                   :: NumNums   ! Number of elements in the numeric array
  INTEGER                                   :: IOStat    ! IO Status when calling get input subroutine
  CHARACTER(len=MaxNameLength),DIMENSION(1) :: AlphArray ! Character string data
  REAL(r64), DIMENSION(3)                   :: NumArray  ! Numeric data

     ! FLOW:
  CurrentModuleObject='Site:HeightVariation'
  NumObjects = GetNumObjectsFound(TRIM(CurrentModuleObject))

  IF (NumObjects == 1) THEN
    CALL GetObjectItem(TRIM(CurrentModuleObject),1,AlphArray,NumAlphas,NumArray,NumNums,IOStat,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

    IF (NumNums > 0) SiteWindExp = NumArray(1)
    IF (NumNums > 1) SiteWindBLHeight = NumArray(2)
    IF (NumNums > 2) SiteTempGradient = NumArray(3)

  ELSE IF (NumObjects > 1) THEN
    CALL ShowSevereError('Too many '//TRIM(CurrentModuleObject)//' objects, only 1 allowed.')
    ErrorsFound = .TRUE.
  ELSE  !  None entered
  ! IDD defaults would have this:
  ! Building object defaults use Terrain to set SiteWindExp and SiteWindBLHeight but would
  ! be overridden by a Site Atmospheric Variation Object.
           !SiteWindExp = 0.22
           !SiteWindBLHeight = 370.0
    SiteTempGradient = 0.0065d0
  END IF

  ! Write to the initialization output file
  WRITE(OutputFileInits,'(A)') '! <Environment:Site Atmospheric Variation>,Wind Speed Profile Exponent {},'// &
    'Wind Speed Profile Boundary Layer Thickness {m},Air Temperature Gradient Coefficient {K/m}'

  WRITE(OutputFileInits,720) TRIM(RoundSigDigits(SiteWindExp,3)), TRIM(RoundSigDigits(SiteWindBLHeight,3)),   &
                    TRIM(RoundSigDigits(SiteTempGradient,6))

720 FORMAT('Environment:Site Atmospheric Variation',3(',',A))

  RETURN

END SUBROUTINE GetSiteAtmosphereData


SUBROUTINE GetMaterialData(ErrorsFound)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Richard Liesen
          !       DATE WRITTEN   September 1997
          !       MODIFIED       April 1999; L.Lawrie
          !                      Sept 1999, FCW, Window5 modifications
          !                      Mar 2001, FCW, WindowShade mods
          !                      Sep 2001, FCW, add Material:WindowGasMixture
          !                      Oct 2001, FCW, add Material:WindowBlind
          !                      Dec 2003, FCW, add glass solar/visible transmittance dirt factor
          !                      Feb 2009, TH, added WindowMaterial:GlazingGroup:Thermochromic

          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! The purpose of this subroutine is to serve as a transfer agent
          ! between the input file and the material derived type.  The new input
          ! file is working, and this file reads the material data directly
          ! from the input file and transfer that information to the new data
          ! structure.  Data read in this routine is stored in a
          ! derived type (Material) defined in the DataHeatBalance module.

          ! In April 1999, a new set of material definitions replaced the one "all-purpose"
          ! material definition.  There are now 10 flavors of materials.  Definitions from
          ! the IDD appear below before their counterpart "gets".

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE General, ONLY: RoundSigDigits, ScanForReports

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  LOGICAL, INTENT(INOUT) :: ErrorsFound  ! set to true if errors found in input

          ! SUBROUTINE PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:

  INTEGER :: IOStat           ! IO Status when calling get input subroutine
  CHARACTER(len=MaxNameLength),DIMENSION(5) &
          :: MaterialNames ! Number of Material Alpha names defined
  INTEGER :: MaterNum         ! Counter to keep track of the material number
  INTEGER :: MaterialNumAlpha ! Number of material alpha names being passed
  INTEGER :: MaterialNumProp  ! Number of material properties being passed
  REAL(r64), DIMENSION(27) :: MaterialProps !Temporary array to transfer material properties
  INTEGER :: RegMat   ! Regular Materials -- full property definition
  INTEGER :: RegRMat  ! Regular Materials -- R only property definition
  INTEGER :: AirMat   ! Air space materias in opaque constructions
  INTEGER :: IRTMat   ! Infrared Transmitting Materials -- R only property definition

  INTEGER :: EcoRoofMat !Materials for ecoRoof
  INTEGER :: NumGas   ! Index for loop over gap gases in a mixture
  INTEGER :: NumGases ! Number of gasses in a mixture
  INTEGER :: GasType  ! Gas type index: 1=air, 2=argon, 3=krypton, 4=xenon
  INTEGER :: Loop
  INTEGER :: ICoeff   ! Gas property coefficient index
  LOGICAL :: ErrorInName
  LOGICAL :: IsBlank
  CHARACTER(len=MaxNameLength) :: TypeOfGas  ! Type of window gas fill (Air, Argon, Krypton, &
                                                 ! Xenon, or Custom
  REAL(r64)    :: MinSlatAngGeom, MaxSlatAngGeom ! Minimum and maximum slat angle allowed by slat geometry (deg)
  REAL(r64)    :: ReflectivitySol   ! Glass reflectivity, solar
  REAL(r64)    :: ReflectivityVis   ! Glass reflectivity, visible
  REAL(r64)    :: TransmittivitySol ! Glass transmittivity, solar
  REAL(r64)    :: TransmittivityVis ! Glass transmittivity, visible
  LOGICAL      :: DoReport=.false.
  REAL(r64)    :: DenomRGas         ! Denominator for WindowGas calculations of NominalR

  ! Added TH 1/9/2009 to read the thermochromic glazings
  INTEGER :: iTC = 0
  INTEGER :: iMat = 0

  ! Added TH 7/27/2009 for constructions defined with F or C factro method
  INTEGER :: TotFfactorConstructs  ! Number of slabs-on-grade or underground floor constructions defined with F factors
  INTEGER :: TotCfactorConstructs  ! Number of underground wall constructions defined with C factors

  INTEGER :: DebugFile       =150 !RS: Debugging file denotion, hopefully this works.
    
  OPEN(unit=DebugFile,file='Debug.txt')    !RS: Debugging
  
        ! FLOW:

  RegMat=GetNumObjectsFound('Material')
  RegRMat=GetNumObjectsFound('Material:NoMass')
  IRTMat=GetNumObjectsFound('Material:InfraredTransparent')
  AirMat=GetNumObjectsFound('Material:AirGap')
  W5GlsMat=GetNumObjectsFound('WindowMaterial:Glazing')
  W5GlsMatAlt=GetNumObjectsFound('WindowMaterial:Glazing:RefractionExtinctionMethod')
  W5GasMat=GetNumObjectsFound('WindowMaterial:Gas')
  W5GasMatMixture=GetNumObjectsFound('WindowMaterial:GasMixture')
  TotShades=GetNumObjectsFound('WindowMaterial:Shade')
  TotScreens=GetNumObjectsFound('WindowMaterial:Screen')
  TotBlinds=GetNumObjectsFound('WindowMaterial:Blind')
  EcoRoofMat=GetNumObjectsFound('Material:RoofVegetation')
  TotSimpleWindow = GetNumObjectsFound('WindowMaterial:SimpleGlazingSystem')

  TotMaterials=RegMat+RegRMat+AirMat+W5GlsMat+W5GlsMatAlt+W5GasMat+W5GasMatMixture+ &
       TotShades+TotScreens+TotBlinds+EcoRoofMat+IRTMat+TotSimpleWindow

  TotFfactorConstructs = GetNumObjectsFound('Construction:FfactorGroundFloor')
  TotCfactorConstructs = GetNumObjectsFound('Construction:CfactorUndergroundWall')
  IF (TotFfactorConstructs + TotCfactorConstructs >=1 ) THEN
    ! Add a new fictitious insulation layer and a thermal mass layer for each F or C factor defined construction
    TotMaterials = TotMaterials + 1 + TotFfactorConstructs + TotCfactorConstructs
  ENDIF

  ALLOCATE (Material(TotMaterials))! Allocate the array Size to the number of materials

  ALLOCATE(NominalR(TotMaterials))
  NominalR=0.0

  MaterNum=0

     ! Regular Materials

  CurrentModuleObject='Material'
  DO Loop=1,RegMat

    !Call Input Get routine to retrieve material data
    CALL GetObjectItem(TRIM(CurrentModuleObject),Loop,MaterialNames,MaterialNumAlpha,MaterialProps,MaterialNumProp,IOSTAT,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

    ErrorInName=.false.
    IsBlank=.false.
    CALL VerifyName(MaterialNames(1),Material%Name,MaterNum,ErrorInName,IsBlank,TRIM(CurrentModuleObject)//' Name')
    IF (ErrorInName) THEN
      CALL ShowContinueError('...All Material names must be unique regardless of subtype.')
      ErrorsFound=.true.
      CYCLE
    ENDIF

    !Load the material derived type from the input data.
    MaterNum=MaterNum+1
    Material(MaterNum)%Group=RegularMaterial
    Material(MaterNum)%Name = MaterialNames(1)

    CALL ValidateMaterialRoughness(MaterNum,MaterialNames(2),ErrorsFound)

    Material(MaterNum)%Thickness     = MaterialProps(1)
    Material(MaterNum)%Conductivity  = MaterialProps(2)
    Material(MaterNum)%Density       = MaterialProps(3)
    Material(MaterNum)%SpecHeat      = MaterialProps(4)
    ! min fields is 6 -- previous four will be there
    IF (MaterialNumProp >= 5) THEN
      Material(MaterNum)%AbsorpThermal = MaterialProps(5)
      Material(MaterNum)%AbsorpThermalInput = MaterialProps(5)
    ELSE
      Material(MaterNum)%AbsorpThermal = .9d0
      Material(MaterNum)%AbsorpThermalInput = .9d0
    ENDIF
    IF (MaterialNumProp >= 6) THEN
      Material(MaterNum)%AbsorpSolar   = MaterialProps(6)
      Material(MaterNum)%AbsorpSolarInput = MaterialProps(6)
    ELSE
      Material(MaterNum)%AbsorpSolar   = .7d0
      Material(MaterNum)%AbsorpSolarInput = .7d0
    ENDIF
    IF (MaterialNumProp >= 7) THEN
      Material(MaterNum)%AbsorpVisible = MaterialProps(7)
      Material(MaterNum)%AbsorpVisibleInput = MaterialProps(7)
    ELSE
      Material(MaterNum)%AbsorpVisible = .7d0
      Material(MaterNum)%AbsorpVisibleInput = .7d0
    ENDIF

    IF (Material(MaterNum)%Conductivity > 0.0) THEN
      NominalR(MaterNum)            = Material(MaterNum)%Thickness/Material(MaterNum)%Conductivity
      Material(MaterNum)%Resistance = NominalR(MaterNum)
    ELSE
      CALL ShowSevereError('Positive thermal conductivity required for material '//TRIM(Material(MaterNum)%Name))
      ErrorsFound = .TRUE.
    END IF

  ENDDO

  ! Add the 6" heavy concrete for constructions defined with F or C factor method
  IF (TotFfactorConstructs + TotCfactorConstructs >=1 ) THEN
    MaterNum = MaterNum + 1

    Material(MaterNum)%Group=RegularMaterial
    Material(MaterNum)%Name = '~FC_Concrete'
    Material(MaterNum)%Thickness     = 0.15d0
    Material(MaterNum)%Conductivity  = 1.95d0
    Material(MaterNum)%Density       = 2240.0d0
    Material(MaterNum)%SpecHeat      = 900.0d0
    Material(MaterNum)%Roughness = MediumRough
    Material(MaterNum)%AbsorpSolar = 0.7d0
    Material(MaterNum)%AbsorpThermal = 0.9d0
    Material(MaterNum)%AbsorpVisible = 0.7d0
    NominalR(MaterNum) = Material(MaterNum)%Thickness / Material(MaterNum)%Conductivity
    Material(MaterNum)%Resistance = NominalR(MaterNum)

    RegMat = RegMat + 1
  ENDIF

  CurrentModuleObject='Material:NoMass'
  DO Loop=1,RegRMat

    !Call Input Get routine to retrieve material data
    CALL GetObjectItem(TRIM(CurrentModuleObject),Loop,MaterialNames,MaterialNumAlpha,MaterialProps,MaterialNumProp,IOSTAT,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

    ErrorInName=.false.
    IsBlank=.false.
    CALL VerifyName(MaterialNames(1),Material%Name,MaterNum,ErrorInName,IsBlank,TRIM(CurrentModuleObject)//' Name')
    IF (ErrorInName) THEN
      CALL ShowContinueError('...All Material names must be unique regardless of subtype.')
      ErrorsFound=.true.
      CYCLE
    ENDIF

    !Load the material derived type from the input data.
    MaterNum=MaterNum+1
    Material(MaterNum)%Group=RegularMaterial
    Material(MaterNum)%Name = MaterialNames(1)

    CALL ValidateMaterialRoughness(MaterNum,MaterialNames(2),ErrorsFound)

    Material(MaterNum)%Resistance    = MaterialProps(1)
    Material(MaterNum)%ROnly         = .true.
    IF (MaterialNumProp >= 2) THEN
      Material(MaterNum)%AbsorpThermal = MaterialProps(2)
      Material(MaterNum)%AbsorpThermalInput = MaterialProps(2)
    ELSE
      Material(MaterNum)%AbsorpThermal = .9d0
      Material(MaterNum)%AbsorpThermalInput = .9d0
    ENDIF
    IF (MaterialNumProp >= 3) THEN
      Material(MaterNum)%AbsorpSolar   = MaterialProps(3)
      Material(MaterNum)%AbsorpSolarInput   = MaterialProps(3)
    ELSE
      Material(MaterNum)%AbsorpSolar   = .7d0
      Material(MaterNum)%AbsorpSolarInput = .7d0
    ENDIF
    IF (MaterialNumProp >= 4) THEN
      Material(MaterNum)%AbsorpVisible = MaterialProps(4)
      Material(MaterNum)%AbsorpVisibleInput = MaterialProps(4)
    ELSE
      Material(MaterNum)%AbsorpVisible = .7d0
      Material(MaterNum)%AbsorpVisibleInput = .7d0
    ENDIF

    NominalR(MaterNum)=Material(MaterNum)%Resistance

  ENDDO

  ! Add a fictitious insulation layer for each construction defined with F or C factor method
  IF (TotFfactorConstructs + TotCfactorConstructs >= 1 ) THEN
    DO Loop = 1, TotFfactorConstructs + TotCfactorConstructs
      MaterNum = MaterNum + 1
      Material(MaterNum)%Group = RegularMaterial
      Material(MaterNum)%Name = '~FC_Insulation_' // RoundSigDigits(Loop,0)
      Material(MaterNum)%ROnly = .true.
      Material(MaterNum)%Roughness = MediumRough
      Material(MaterNum)%AbsorpSolar = 0.0
      Material(MaterNum)%AbsorpThermal = 0.0
      Material(MaterNum)%AbsorpVisible = 0.0
    ENDDO
    RegRMat = RegRMat + TotFfactorConstructs + TotCfactorConstructs
  ENDIF


  ! Air Materials (for air spaces in opaque constructions)
  CurrentModuleObject='Material:AirGap'
  DO Loop=1,AirMat

    !Call Input Get routine to retrieve material data
    CALL GetObjectItem(TRIM(CurrentModuleObject),Loop,MaterialNames,MaterialNumAlpha,MaterialProps,MaterialNumProp,IOSTAT,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

    ErrorInName=.false.
    IsBlank=.false.
    CALL VerifyName(MaterialNames(1),Material%Name,MaterNum,ErrorInName,IsBlank,TRIM(CurrentModuleObject)//' Name')
    IF (ErrorInName) THEN
      CALL ShowContinueError('...All Material names must be unique regardless of subtype.')
      ErrorsFound=.true.
      CYCLE
    ENDIF

    !Load the material derived type from the input data.
    MaterNum=MaterNum+1
    Material(MaterNum)%Group=Air
    Material(MaterNum)%Name = MaterialNames(1)

    Material(MaterNum)%Roughness=MediumRough

    Material(MaterNum)%Resistance      = MaterialProps(1)
    Material(MaterNum)%ROnly         = .true.

    NominalR(MaterNum)=Material(MaterNum)%Resistance

  ENDDO

  CurrentModuleObject='Material:InfraredTransparent'
  DO Loop=1,IRTMat

    !Call Input Get routine to retrieve material data
    CALL GetObjectItem(TRIM(CurrentModuleObject),Loop,MaterialNames,MaterialNumAlpha,MaterialProps,MaterialNumProp,IOSTAT,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

    ErrorInName=.false.
    IsBlank=.false.
    CALL VerifyName(MaterialNames(1),Material%Name,MaterNum,ErrorInName,IsBlank,TRIM(CurrentModuleObject)//' Name')
    IF (ErrorInName) THEN
      CALL ShowContinueError('...All Material names must be unique regardless of subtype.')
      ErrorsFound=.true.
      CYCLE
    ENDIF

    MaterNum=MaterNum+1
    Material(MaterNum)%Group=IRTMaterial

    !Load the material derived type from the input data.
    Material(MaterNum)%Name = MaterialNames(1)

    IF (MaterialNumProp >= 1) THEN
      Material(MaterNum)%Resistance    = MaterialProps(1)
      Material(MaterNum)%ROnly         = .true.
    ELSE
      Material(MaterNum)%Resistance = .01d0
    ENDIF
    IF (MaterialNumProp >= 2) THEN
      Material(MaterNum)%AbsorpThermal = MaterialProps(2)
      Material(MaterNum)%AbsorpThermalInput =  MaterialProps(2)
    ELSE
      Material(MaterNum)%AbsorpThermal = 0.9999d0
      Material(MaterNum)%AbsorpThermalInput = 0.9999d0
    ENDIF
    IF (MaterialNumProp >= 3) THEN
      Material(MaterNum)%AbsorpSolar   = MaterialProps(3)
      Material(MaterNum)%AbsorpSolarInput = MaterialProps(3)
    ELSE
      Material(MaterNum)%AbsorpSolar   = 1.d0
      Material(MaterNum)%AbsorpSolarInput = 1.d0
    ENDIF
    IF (MaterialNumProp >= 4) THEN
      Material(MaterNum)%AbsorpVisible = MaterialProps(4)
      Material(MaterNum)%AbsorpVisibleInput = MaterialProps(4)
    ELSE
      Material(MaterNum)%AbsorpVisible = 1.d0
      Material(MaterNum)%AbsorpVisibleInput = 1.d0
    ENDIF

    NominalR(MaterNum)=Material(MaterNum)%Resistance

  ENDDO

  ! Glass materials, regular input: transmittance and front/back reflectance

  CurrentModuleObject='WindowMaterial:Glazing'
  DO Loop=1,W5GlsMat

    !Call Input Get routine to retrieve material data
    CALL GetObjectItem(TRIM(CurrentModuleObject),Loop,MaterialNames,MaterialNumAlpha,MaterialProps,MaterialNumProp,IOSTAT,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

    ErrorInName=.false.
    IsBlank=.false.
    CALL VerifyName(MaterialNames(1),Material%Name,MaterNum,ErrorInName,IsBlank,TRIM(CurrentModuleObject)//' Name')
    IF (ErrorInName) THEN
      CALL ShowContinueError('...All Material names must be unique regardless of subtype.')
      ErrorsFound=.true.
      CYCLE
    ENDIF

    MaterNum=MaterNum+1
    Material(MaterNum)%Group=WindowGlass

    !Load the material derived type from the input data.

    Material(MaterNum)%Name = MaterialNames(1)
    Material(MaterNum)%Roughness=VerySmooth
    Material(MaterNum)%ROnly         = .true.
    Material(MaterNum)%Thickness           = MaterialProps(1)
    Material(MaterNum)%Trans               = MaterialProps(2)
    Material(MaterNum)%ReflectSolBeamFront = MaterialProps(3)
    Material(MaterNum)%ReflectSolBeamBack  = MaterialProps(4)
    Material(MaterNum)%TransVis            = MaterialProps(5)
    Material(MaterNum)%ReflectVisBeamFront = MaterialProps(6)
    Material(MaterNum)%ReflectVisBeamBack  = MaterialProps(7)
    Material(MaterNum)%TransThermal        = MaterialProps(8)
    Material(MaterNum)%AbsorpThermalFront  = MaterialProps(9)
    Material(MaterNum)%AbsorpThermalBack   = MaterialProps(10)
    Material(MaterNum)%Conductivity        = MaterialProps(11)
    Material(MaterNum)%GlassTransDirtFactor= MaterialProps(12)
    IF(MaterialProps(12) == 0.0) Material(MaterNum)%GlassTransDirtFactor = 1.0
    Material(MaterNum)%AbsorpThermal       = Material(MaterNum)%AbsorpThermalBack

    IF (Material(MaterNum)%Conductivity > 0.0) THEN
      NominalR(MaterNum)=Material(MaterNum)%Thickness/Material(MaterNum)%Conductivity
    ELSE
      ErrorsFound = .true.
      CALL ShowSevereError('Window glass material ' //Trim(Material(MaterNum)%Name)// &
            ' has Conductivity = 0.0, must be >0.0, default = .9')
    END IF

    Material(MaterNum)%GlassSpectralDataPtr = 0
    IF (TotSpectralData > 0 .and. .not. lAlphaFieldBlanks(3)) THEN
    Material(MaterNum)%GlassSpectralDataPtr = FindIteminList(MaterialNames(3),SpectralData%Name,TotSpectralData)
    ENDIF
    IF(SameString(MaterialNames(2),'SpectralAverage')) Material(MaterNum)%GlassSpectralDataPtr = 0
    IF(Material(MaterNum)%GlassSpectralDataPtr == 0 .AND. SameString(MaterialNames(2),'Spectral')) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(trim(CurrentModuleObject)//'="'//Trim(Material(MaterNum)%Name)// &
            '" has '//TRIM(cAlphaFieldNames(2))//' = Spectral but has no matching MaterialProperty:GlazingSpectralData set')
      IF (lAlphaFieldBlanks(3)) THEN
        CALL ShowContinueError('...'//trim(cAlphaFieldNames(3))//' is blank.')
      ELSE
        CALL ShowContinueError('...'//trim(cAlphaFieldNames(3))//'="'//trim(MaterialNames(3))//  &
           '" not found as item in MaterialProperty:GlazingSpectralData objects.')
    END IF
    END IF

    IF(.not. SameString(MaterialNames(2),'SpectralAverage') .AND. .not. SameString(MaterialNames(2),'Spectral')) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(trim(CurrentModuleObject)//'="'//Trim(Material(MaterNum)%Name)//'", invalid specification.')
      CALL ShowContinueError(TRIM(cAlphaFieldNames(2))//' must be SpectralAverage or Spectral, value='//TRIM(MaterialNames(2)))
    END IF

    ! TH 8/24/2011, allow glazing properties MaterialProps(2 to 10) to equal 0 or 1: 0.0 =< Prop <= 1.0
    ! Fixed CR 8413 - modeling spandrel panels as glazing systems
    IF(SameString(MaterialNames(2),'SpectralAverage')) THEN

      IF(MaterialProps(2)+MaterialProps(3) > 1.0) THEN
        ErrorsFound = .true.
        CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
        CALL ShowContinueError(TRIM(cNumericFieldNames(2))//' + '//TRIM(cNumericFieldNames(3))//' not <= 1.0')
      END IF

      IF(MaterialProps(2)+MaterialProps(4) > 1.0) THEN
        ErrorsFound = .true.
        CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
        CALL ShowContinueError(TRIM(cNumericFieldNames(2))//' + '//TRIM(cNumericFieldNames(4))//' not <= 1.0')
      END IF

      IF(MaterialProps(5)+MaterialProps(6) > 1.0) THEN
        ErrorsFound = .true.
        CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
        CALL ShowContinueError(TRIM(cNumericFieldNames(5))//' + '//TRIM(cNumericFieldNames(6))//' not <= 1.0')
      END IF

      IF(MaterialProps(5)+MaterialProps(7) > 1.0) THEN
        ErrorsFound = .true.
        CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
        CALL ShowContinueError(TRIM(cNumericFieldNames(5))//' + '//TRIM(cNumericFieldNames(7))//' not <= 1.0')
      END IF

      IF(MaterialProps(8)+MaterialProps(9) > 1.0) THEN
        ErrorsFound = .true.
        CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
        CALL ShowContinueError(TRIM(cNumericFieldNames(8))//' + '//TRIM(cNumericFieldNames(9))//' not <= 1.0')
      END IF

      IF(MaterialProps(8)+MaterialProps(10) > 1.0) THEN
        ErrorsFound = .true.
        CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
        CALL ShowContinueError(TRIM(cNumericFieldNames(8))//' + '//TRIM(cNumericFieldNames(10))//' not <= 1.0')
      END IF

      IF(MaterialProps(2) < 0.0) THEN
        CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
        CALL ShowContinueError(TRIM(cNumericFieldNames(2))//' not >= 0.0')
        ErrorsFound=.true.
      END IF

      IF(MaterialProps(2) > 1.0) THEN
        ErrorsFound = .true.
        CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
        CALL ShowContinueError(TRIM(cNumericFieldNames(2))//' not <= 1.0')
      END IF

      IF(MaterialProps(3) < 0.0 .or. MaterialProps(3) > 1.0) THEN
        ErrorsFound = .true.
        CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
        CALL ShowContinueError(TRIM(cNumericFieldNames(3))//' not >= 0.0 and <= 1.0')
      END IF

      IF(MaterialProps(4) < 0.0 .or. MaterialProps(4) > 1.0) THEN
        ErrorsFound = .true.
        CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
        CALL ShowContinueError(TRIM(cNumericFieldNames(4))//' not >= 0.0 and <= 1.0')
      END IF

      IF(MaterialProps(5) < 0.0) THEN
        CALL ShowWarningError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", minimal value.')
        CALL ShowWarningError(TRIM(cNumericFieldNames(5))//' not >= 0.0')
      END IF

      IF(MaterialProps(5) > 1.0) THEN
        ErrorsFound = .true.
        CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
        CALL ShowContinueError(TRIM(cNumericFieldNames(5))//' not <= 1.0')
      END IF

      IF(MaterialProps(6) < 0.0 .or. MaterialProps(6) > 1.0) THEN
        ErrorsFound = .true.
        CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
        CALL ShowContinueError(TRIM(cNumericFieldNames(6))//' not >= 0.0 and <= 1.0')
      END IF

      IF(MaterialProps(7) < 0.0 .or. MaterialProps(7) > 1.0) THEN
        ErrorsFound = .true.
        CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
        CALL ShowContinueError(TRIM(cNumericFieldNames(7))//' not >= 0.0 and <= 1.0')
      END IF

    END IF

    IF(MaterialProps(8) > 1.0) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(8))//' not <= 1.0')
    END IF

    IF(MaterialProps(9) <= 0.0 .or. MaterialProps(9) >= 1.0) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(9))//' not > 0.0 and < 1.0')
    END IF

    IF(MaterialProps(10) <= 0.0 .or. MaterialProps(10) >= 1.0) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(10))//' not > 0.0 and < 1.0')
    END IF

    IF(MaterialNames(4) == ' ') THEN
      Material(MaterNum)%SolarDiffusing = .false.
    ELSE IF(MaterialNames(4) == 'YES') THEN
      Material(MaterNum)%SolarDiffusing = .true.
    ELSE IF(MaterialNames(4) == 'NO') THEN
      Material(MaterNum)%SolarDiffusing = .false.
    ELSE
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(4))//' must be Yes or No, entered value='//TRIM(MaterialNames(4)))
    END IF

  ENDDO

! Glass materials, alternative input: index of refraction and extinction coefficient

  CurrentModuleObject='WindowMaterial:Glazing:RefractionExtinctionMethod'
  DO Loop=1,W5GlsMatAlt

    !Call Input Get routine to retrieve material data
    CALL GetObjectItem(TRIM(CurrentModuleObject),Loop,MaterialNames,MaterialNumAlpha,MaterialProps,MaterialNumProp,IOSTAT,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

    ErrorInName=.false.
    IsBlank=.false.
    CALL VerifyName(MaterialNames(1),Material%Name,MaterNum,ErrorInName,IsBlank,TRIM(CurrentModuleObject)//' Name')
    IF (ErrorInName) THEN
      CALL ShowContinueError('...All Material names must be unique regardless of subtype.')
      ErrorsFound=.true.
      CYCLE
    ENDIF

    MaterNum=MaterNum+1
    Material(MaterNum)%Group=WindowGlass

    !Load the material derived type from the input data.

    Material(MaterNum)%Name = MaterialNames(1)
    Material(MaterNum)%Roughness = VerySmooth
    Material(MaterNum)%Thickness = MaterialProps(1)
    Material(MaterNum)%ROnly         = .true.

    ! Calculate solar and visible transmittance and reflectance at normal incidence from thickness,
    ! index of refraction and extinction coefficient. With the alternative input the front and back
    ! properties are assumed to be the same.

    ReflectivitySol = ((MaterialProps(2)-1.d0)/(MaterialProps(2)+1.d0))**2
    ReflectivityVis = ((MaterialProps(4)-1.d0)/(MaterialProps(4)+1.d0))**2
    TransmittivitySol = EXP(-MaterialProps(3)*MaterialProps(1))
    TransmittivityVis = EXP(-MaterialProps(5)*MaterialProps(1))
    Material(MaterNum)%Trans               = TransmittivitySol * ((1.d0-ReflectivitySol)**2) / &
                                                (1.d0-(ReflectivitySol*TransmittivitySol)**2)
    Material(MaterNum)%ReflectSolBeamFront = ReflectivitySol * (1.d0 +  &
                                                ((1.d0-ReflectivitySol)**2)*(TransmittivitySol**2) / &
                                                (1.d0-(ReflectivitySol*TransmittivitySol)**2) )
    Material(MaterNum)%ReflectSolBeamBack  = Material(MaterNum)%ReflectSolBeamFront
    Material(MaterNum)%TransVis            = TransmittivityVis * ((1.d0-ReflectivityVis)**2) / &
                                                (1.d0-(ReflectivityVis*TransmittivityVis)**2)

    Material(MaterNum)%ReflectVisBeamFront = ReflectivityVis * (1.d0 +  &
                                                ((1.d0-ReflectivityVis)**2)*(TransmittivityVis**2) / &
                                                (1.d0-(ReflectivityVis*TransmittivityVis)**2) )
    Material(MaterNum)%ReflectVisBeamBack  = Material(MaterNum)%ReflectSolBeamFront
    Material(MaterNum)%TransThermal        = MaterialProps(6)
    Material(MaterNum)%AbsorpThermalFront  = MaterialProps(7)
    Material(MaterNum)%AbsorpThermalBack   = MaterialProps(7)
    Material(MaterNum)%Conductivity        = MaterialProps(8)
    Material(MaterNum)%GlassTransDirtFactor= MaterialProps(9)
    IF(MaterialProps(9) == 0.0) Material(MaterNum)%GlassTransDirtFactor = 1.0
    Material(MaterNum)%AbsorpThermal       = Material(MaterNum)%AbsorpThermalBack

    IF (Material(MaterNum)%Conductivity > 0.0) NominalR(MaterNum)=Material(MaterNum)%Thickness/Material(MaterNum)%Conductivity

    Material(MaterNum)%GlassSpectralDataPtr = 0

    IF(MaterialProps(6)+MaterialProps(7) >= 1.0) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(6))//' + '//TRIM(cNumericFieldNames(7))//' not < 1.0')
    END IF

    IF(MaterialNames(2) == ' ') THEN
      Material(MaterNum)%SolarDiffusing = .false.
    ELSE IF(MaterialNames(2) == 'YES') THEN
      Material(MaterNum)%SolarDiffusing = .true.
    ELSE IF(MaterialNames(2) == 'NO') THEN
      Material(MaterNum)%SolarDiffusing = .false.
    ELSE
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(2))//' must be Yes or No, entered value='//TRIM(MaterialNames(4)))
    END IF

  ENDDO

  ! Window gas materials (for gaps with a single gas)

  CurrentModuleObject='WindowMaterial:Gas'
  DO Loop=1,W5GasMat

    !Call Input Get routine to retrieve material data
    CALL GetObjectItem(TRIM(CurrentModuleObject),Loop,MaterialNames,MaterialNumAlpha,MaterialProps,MaterialNumProp,IOSTAT,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

    ErrorInName=.false.
    IsBlank=.false.
    CALL VerifyName(MaterialNames(1),Material%Name,MaterNum,ErrorInName,IsBlank,TRIM(CurrentModuleObject)//' Name')
    IF (ErrorInName) THEN
      CALL ShowContinueError('...All Material names must be unique regardless of subtype.')
      ErrorsFound=.true.
      CYCLE
    ENDIF

    MaterNum=MaterNum+1
    Material(MaterNum)%Group=WindowGas
    Material(MaterNum)%GasType(1) = -1
    Material(MaterNum)%NumberOfGasesInMixture = 1
    Material(MaterNum)%GasFract(1) = 1.0

    !Load the material derived type from the input data.

    Material(MaterNum)%Name = MaterialNames(1)
    Material(MaterNum)%NumberOfGasesInMixture = 1
    TypeOfGas = TRIM(MaterialNames(2))
    IF(TypeOfGas == 'AIR')     Material(MaterNum)%GasType(1) = 1
    IF(TypeOfGas == 'ARGON')   Material(MaterNum)%GasType(1) = 2
    IF(TypeOfGas == 'KRYPTON') Material(MaterNum)%GasType(1) = 3
    IF(TypeOfGas == 'XENON')   Material(MaterNum)%GasType(1) = 4
    IF(TypeOfGas == 'CUSTOM')  Material(MaterNum)%GasType(1) = 0

    IF(Material(MaterNum)%GasType(1) == -1) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
      CALL ShowContinueError(trim(cAlphaFieldNames(2))//' entered value="'//TRIM(TypeOfGas)//  &
         '" should be Air, Argon, Krypton, Xenon or Custom.')
    END IF

    Material(MaterNum)%Roughness=MediumRough

    Material(MaterNum)%Thickness              = MaterialProps(1)
    Material(MaterNum)%ROnly         = .true.

    GasType = Material(MaterNum)%GasType(1)
    IF(GasType >= 1 .AND. GasType <= 4) THEN
      Material(MaterNum)%GasWght(1) = GasWght(GasType)
      DO ICoeff = 1,3
        Material(MaterNum)%GasCon(1,ICoeff)   = GasCoeffsCon(GasType,ICoeff)
        Material(MaterNum)%GasVis(1,ICoeff)   = GasCoeffsVis(GasType,ICoeff)
        Material(MaterNum)%GasCp (1,ICoeff)   = GasCoeffsCp (GasType,ICoeff)
      END DO
    END IF

    ! Custom gas

    IF(GasType == 0) THEN
      DO ICoeff = 1,2
        Material(MaterNum)%GasCon(1,ICoeff)   = MaterialProps(1+ICoeff)
        Material(MaterNum)%GasVis(1,ICoeff)   = MaterialProps(3+ICoeff)
        Material(MaterNum)%GasCp (1,ICoeff)   = MaterialProps(5+ICoeff)
      END DO
      Material(MaterNum)%GasWght(1)           = MaterialProps(8)

      ! Check for errors in custom gas properties
!      IF(Material(MaterNum)%GasCon(1,1) <= 0.0) THEN
!        ErrorsFound = .true.
!        CALL ShowSevereError('Conductivity Coefficient A for custom window gas='&
!                 //TRIM(MaterialNames(1))//' should be > 0.')
!      END IF

      IF(Material(MaterNum)%GasVis(1,1) <= 0.0) THEN
        ErrorsFound = .true.
        CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
        CALL ShowContinueError(TRIM(cNumericFieldNames(3+ICoeff))//' not > 0.0')
      END IF
      IF(Material(MaterNum)%GasCp(1,1) <= 0.0) THEN
        ErrorsFound = .true.
        CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
        CALL ShowContinueError(TRIM(cNumericFieldNames(5+ICoeff))//' not > 0.0')
      END IF
      IF(Material(MaterNum)%GasWght(1) <= 0.0) THEN
        ErrorsFound = .true.
        CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
        CALL ShowContinueError(TRIM(cNumericFieldNames(8))//' not > 0.0')
      END IF
    END IF

    ! Nominal resistance of gap at room temperature
    IF(.not.ErrorsFound) THEN
      DenomRGas=(Material(MaterNum)%GasCon(1,1) + Material(MaterNum)%GasCon(1,2)*300.0d0 + Material(MaterNum)%GasCon(1,3)*90000.0d0)
      IF (DenomRGas > 0.0) THEN
        NominalR(MaterNum)=Material(MaterNum)%Thickness/DenomRGas
      ELSE
        CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
        CALL ShowContinueError('Nominal resistance of gap at room temperature calculated at a negative Conductivity=['//  &
                         trim(RoundSigDigits(DenomRGas,3))//'].')
        ErrorsFound=.true.
      ENDIF
    ENDIF

  ENDDO

! Window gas mixtures (for gaps with two or more gases)

  CurrentModuleObject='WindowMaterial:GasMixture'
  DO Loop=1,W5GasMatMixture

    !Call Input Get routine to retrieve material data
    CALL GetObjectItem(TRIM(CurrentModuleObject),Loop,MaterialNames,MaterialNumAlpha,MaterialProps,MaterialNumProp,IOSTAT,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

    ErrorInName=.false.
    IsBlank=.false.
    CALL VerifyName(MaterialNames(1),Material%Name,MaterNum,ErrorInName,IsBlank,TRIM(CurrentModuleObject)//' Name')
    IF (ErrorInName) THEN
      CALL ShowContinueError('...All Material names must be unique regardless of subtype.')
      ErrorsFound=.true.
      CYCLE
    ENDIF

    MaterNum=MaterNum+1
    Material(MaterNum)%Group=WindowGasMixture
    Material(MaterNum)%GasType = -1

    !Load the material derived type from the input data.

    Material(MaterNum)%Name = MaterialNames(1)
    NumGases = MaterialProps(2)
    Material(MaterNum)%NumberOfGasesInMixture = NumGases
    DO NumGas = 1,NumGases
      TypeOfGas = TRIM(MaterialNames(1+NumGas))
      IF(TypeOfGas == 'AIR')     Material(MaterNum)%GasType(NumGas) = 1
      IF(TypeOfGas == 'ARGON')   Material(MaterNum)%GasType(NumGas) = 2
      IF(TypeOfGas == 'KRYPTON') Material(MaterNum)%GasType(NumGas) = 3
      IF(TypeOfGas == 'XENON')   Material(MaterNum)%GasType(NumGas) = 4
      IF(Material(MaterNum)%GasType(NumGas) == -1) THEN
        ErrorsFound = .true.
        CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
        CALL ShowContinueError(trim(cAlphaFieldNames(1+NumGas))//' entered value="'//TRIM(TypeOfGas)//  &
           '" should be Air, Argon, Krypton, or Xenon.')
      END IF
    END DO

    Material(MaterNum)%Roughness=MediumRough  ! Unused

    Material(MaterNum)%Thickness              = MaterialProps(1)
    Material(MaterNum)%ROnly         = .true.

    DO NumGas = 1,NumGases
      GasType = Material(MaterNum)%GasType(NumGas)
      IF(GasType >= 1 .AND. GasType <= 4) THEN
        Material(MaterNum)%GasWght(NumGas) = GasWght(GasType)
        Material(MaterNum)%GasFract(NumGas) = MaterialProps(2+NumGas)
        DO ICoeff = 1,3
          Material(MaterNum)%GasCon(NumGas,ICoeff)   = GasCoeffsCon(GasType,ICoeff)
          Material(MaterNum)%GasVis(NumGas,ICoeff)   = GasCoeffsVis(GasType,ICoeff)
          Material(MaterNum)%GasCp (NumGas,ICoeff)   = GasCoeffsCp (GasType,ICoeff)
        END DO
      END IF
    END DO

    ! Nominal resistance of gap at room temperature (based on first gas in mixture)
    NominalR(MaterNum)=Material(MaterNum)%Thickness/(Material(MaterNum)%GasCon(1,1) + &
       Material(MaterNum)%GasCon(1,2)*300.0 + Material(MaterNum)%GasCon(1,3)*90000.0)

  ENDDO

  ! Window Shade Materials

  CurrentModuleObject='WindowMaterial:Shade'
  DO Loop=1,TotShades

    !Call Input Get routine to retrieve material data
    CALL GetObjectItem(TRIM(CurrentModuleObject),Loop,MaterialNames,MaterialNumAlpha,MaterialProps,MaterialNumProp,IOSTAT,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

    ErrorInName=.false.
    IsBlank=.false.
    CALL VerifyName(MaterialNames(1),Material%Name,MaterNum,ErrorInName,IsBlank,TRIM(CurrentModuleObject)//' Name')
    IF (ErrorInName) THEN
      CALL ShowContinueError('...All Material names must be unique regardless of subtype.')
      ErrorsFound=.true.
      CYCLE
    ENDIF

    MaterNum=MaterNum+1
    Material(MaterNum)%Group=Shade

    !Load the material derived type from the input data.

    Material(MaterNum)%Name            = MaterialNames(1)
    Material(MaterNum)%Roughness       = MediumRough
    Material(MaterNum)%Trans           = MaterialProps(1)
    Material(MaterNum)%ReflectShade    = MaterialProps(2)
    Material(MaterNum)%TransVis        = MaterialProps(3)
    Material(MaterNum)%ReflectShadeVis = MaterialProps(4)
    Material(MaterNum)%AbsorpThermal   = MaterialProps(5)
    Material(MaterNum)%AbsorpThermalInput = MaterialProps(5)
    Material(MaterNum)%TransThermal    = MaterialProps(6)
    Material(MaterNum)%Thickness       = MaterialProps(7)
    Material(MaterNum)%Conductivity    = MaterialProps(8)
    Material(MaterNum)%AbsorpSolar     = MAX(0.d0,1.d0- Material(MaterNum)%Trans - Material(MaterNum)%ReflectShade)
    Material(MaterNum)%AbsorpSolarInput = Material(MaterNum)%AbsorpSolar
    Material(MaterNum)%WinShadeToGlassDist         = MaterialProps(9)
    Material(MaterNum)%WinShadeTopOpeningMult      = MaterialProps(10)
    Material(MaterNum)%WinShadeBottomOpeningMult   = MaterialProps(11)
    Material(MaterNum)%WinShadeLeftOpeningMult     = MaterialProps(12)
    Material(MaterNum)%WinShadeRightOpeningMult    = MaterialProps(13)
    Material(MaterNum)%WinShadeAirFlowPermeability = MaterialProps(14)
    Material(MaterNum)%ROnly         = .true.

    IF (Material(MaterNum)%Conductivity > 0.0) THEN
      NominalR(MaterNum)=Material(MaterNum)%Thickness/Material(MaterNum)%Conductivity
    ELSE
      NominalR(MaterNum)=1.0
    ENDIF

    IF(MaterialProps(1)+MaterialProps(2) >= 1.0) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(1))//' + '//TRIM(cNumericFieldNames(2))//' not < 1.0')
    END IF

    IF(MaterialProps(3)+MaterialProps(4) >= 1.0) THEN
      ErrorsFound = .true.
      !CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      !CALL ShowContinueError(TRIM(cNumericFieldNames(3))//' + '//TRIM(cNumericFieldNames(4))//' not < 1.0') !RS: Secret Search String
      WRITE(DebugFile,*) TRIM(CurrentModuleObject)//'="'//TRIM(MaterialNames(1))//'", Illegal value combination.'
      WRITE(DebugFile,*) TRIM(cNumericFieldNames(3))//' + '//TRIM(cNumericFieldNames(4))//' not < 1.0'
    END IF

    IF(MaterialProps(5)+MaterialProps(6) >= 1.0) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(5))//' + '//TRIM(cNumericFieldNames(6))//' not < 1.0')
    END IF

  ENDDO

  ! Window Screen Materials

  CurrentModuleObject='WindowMaterial:Screen'
  DO Loop=1,TotScreens

    !Call GetObjectItem routine to retrieve material data
    CALL GetObjectItem(TRIM(CurrentModuleObject),Loop,MaterialNames,MaterialNumAlpha,MaterialProps,MaterialNumProp,IOSTAT,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

    ErrorInName=.false.
    IsBlank=.false.
    CALL VerifyName(MaterialNames(1),Material%Name,MaterNum,ErrorInName,IsBlank,TRIM(CurrentModuleObject)//' Name')
    IF (ErrorInName) THEN
      CALL ShowContinueError('...All Material names must be unique regardless of subtype.')
      ErrorsFound=.true.
      CYCLE
    ENDIF

    MaterNum=MaterNum+1
    Material(MaterNum)%Group=Screen

    !Load the material derived type from the input data.

    Material(MaterNum)%Name               = MaterialNames(1)
    Material(MaterNum)%ReflectanceModeling= MaterialNames(2)
    IF(.NOT. (SameString(MaterialNames(2),'DoNotModel') .OR. &
              SameString(MaterialNames(2),'ModelAsDirectBeam') .OR. &
              SameString(MaterialNames(2),'ModelAsDiffuse')))THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
      CALL ShowContinueError(TRIM(cAlphaFieldNames(2))//'="'//Trim(MaterialNames(2))//  &
          '", must be one of DoNotModel, ModelAsDirectBeam or ModelAsDiffuse.')
    END IF
    Material(MaterNum)%Roughness          = MediumRough
    Material(MaterNum)%ReflectShade       = MaterialProps(1)
    IF(Material(MaterNum)%ReflectShade .LT. 0.0 .OR. Material(MaterNum)%ReflectShade .GT. 1.0)THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(1))//' must be >= 0 and <= 1')
    END IF
    Material(MaterNum)%ReflectShadeVis    = MaterialProps(2)
    IF(Material(MaterNum)%ReflectShadeVis .LT. 0.0 .OR. Material(MaterNum)%ReflectShadeVis .GT. 1.0)THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(2))//' must be >= 0 and <= 1 for material ' &
                            //TRIM(Material(MaterNum)%Name)//'.')
    END IF
    Material(MaterNum)%AbsorpThermal      = MaterialProps(3)
    Material(MaterNum)%AbsorpThermalInput = MaterialProps(3)
    IF(Material(MaterNum)%AbsorpThermal .LT. 0.0 .OR. Material(MaterNum)%AbsorpThermal .GT. 1.0)THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(3))//' must be >= 0 and <= 1')
    END IF
    Material(MaterNum)%Conductivity       = MaterialProps(4)
    Material(MaterNum)%Thickness          = MaterialProps(6) ! thickness = diameter

    IF(MaterialProps(5) .GT. 0.0)THEN
!      SurfaceScreens(ScNum)%ScreenDiameterToSpacingRatio = MaterialProps(6)/MaterialProps(5) or 1-SQRT(Material(MaterNum)%Trans
      IF(MaterialProps(6)/MaterialProps(5) .GE. 1.0)THEN
        ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
        CALL ShowContinueError(TRIM(cNumericFieldNames(6))//' must be less than '//TRIM(cNumericFieldNames(5)))
      ELSE
!       Calculate direct normal transmittance (open area fraction)
        Material(MaterNum)%Trans = (1.d0 - MaterialProps(6)/MaterialProps(5))**2
      END IF
    ELSE
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(5))//' must be > 0.')
      MaterialProps(5) = 0.000000001d0
    END IF

    IF(MaterialProps(6) .LE. 0.0)THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(6))//' must be > 0.')
    END IF

!   Modify reflectance to account for the open area in the screen assembly
    Material(MaterNum)%ReflectShade       = Material(MaterNum)%ReflectShade * (1.d0 - Material(MaterNum)%Trans)
    Material(MaterNum)%ReflectShadeVis    = Material(MaterNum)%ReflectShadeVis * (1.d0 - Material(MaterNum)%Trans)

    Material(MaterNum)%WinShadeToGlassDist         = MaterialProps(7)
    IF(Material(MaterNum)%WinShadeToGlassDist .LT. 0.001d0 .OR. Material(MaterNum)%WinShadeToGlassDist .GT. 1.0d0)THEN
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(7))//' must be greater than or equal to 0.001 and less than or equal to 1.')
    ENDIF

    Material(MaterNum)%WinShadeTopOpeningMult      = MaterialProps(8)
    IF(Material(MaterNum)%WinShadeTopOpeningMult .LT. 0.0d0 .OR. Material(MaterNum)%WinShadeTopOpeningMult .GT. 1.0d0)THEN
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(8))//' must be greater than or equal to 0 and less than or equal to 1.')
    ENDIF

    Material(MaterNum)%WinShadeBottomOpeningMult   = MaterialProps(9)
    IF(Material(MaterNum)%WinShadeBottomOpeningMult .LT. 0.0d0 .OR. Material(MaterNum)%WinShadeBottomOpeningMult .GT. 1.0d0)THEN
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(9))//' must be greater than or equal to 0 and less than or equal to 1.')
    ENDIF

    Material(MaterNum)%WinShadeLeftOpeningMult     = MaterialProps(10)
    IF(Material(MaterNum)%WinShadeLeftOpeningMult .LT. 0.0d0 .OR. Material(MaterNum)%WinShadeLeftOpeningMult .GT. 1.0d0)THEN
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(10))//' must be greater than or equal to 0 and less than or equal to 1.')
    ENDIF

    Material(MaterNum)%WinShadeRightOpeningMult    = MaterialProps(11)
    IF(Material(MaterNum)%WinShadeRightOpeningMult .LT. 0.0d0 .OR. Material(MaterNum)%WinShadeRightOpeningMult .GT. 1.0d0)THEN
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(11))//' must be greater than or equal to 0 and less than or equal to 1.')
    ENDIF

    Material(MaterNum)%ScreenMapResolution         = MaterialProps(12)
    IF(Material(MaterNum)%ScreenMapResolution .LT. 0 .OR. Material(MaterNum)%ScreenMapResolution .GT. 5 .OR. &
       Material(MaterNum)%ScreenMapResolution .EQ. 4)THEN
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(12))//' must be 0, 1, 2, 3, or 5.')
      ErrorsFound = .true.
    ENDIF

!   Default air flow permeability to open area fraction
    Material(MaterNum)%WinShadeAirFlowPermeability = Material(MaterNum)%Trans
    Material(MaterNum)%TransThermal                = Material(MaterNum)%Trans
    Material(MaterNum)%TransVis                    = Material(MaterNum)%Trans

    Material(MaterNum)%ROnly           = .true.

!   Calculate absorptance accounting for the open area in the screen assembly (used only in CreateShadedWindowConstruction)
    Material(MaterNum)%AbsorpSolar        = MAX(0.d0,1.d0- Material(MaterNum)%Trans - Material(MaterNum)%ReflectShade)
    Material(MaterNum)%AbsorpSolarInput   = Material(MaterNum)%AbsorpSolar
    Material(MaterNum)%AbsorpVisible      = MAX(0.d0,1.d0- Material(MaterNum)%TransVis - Material(MaterNum)%ReflectShadeVis)
    Material(MaterNum)%AbsorpVisibleInput = Material(MaterNum)%AbsorpVisible
    Material(MaterNum)%AbsorpThermal      = Material(MaterNum)%AbsorpThermal * (1.0d0 - Material(MaterNum)%Trans)
    Material(MaterNum)%AbsorpThermalInput = Material(MaterNum)%AbsorpThermal

    IF (Material(MaterNum)%Conductivity > 0.0) THEN
      NominalR(MaterNum)=(1.d0-Material(MaterNum)%Trans)*Material(MaterNum)%Thickness/Material(MaterNum)%Conductivity
    ELSE
      NominalR(MaterNum)=1.0
      CALL ShowWarningError('Conductivity for material="'//TRIM(Material(MaterNum)%Name)//'" must be greater than' &
                          //' 0 for calculating Nominal R-value, Nominal R is defaulted to 1 and the simulation continues.')
    ENDIF

    IF(Material(MaterNum)%Trans+Material(MaterNum)%ReflectShade >= 1.0) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      CALL ShowContinueError('Calculated solar transmittance + solar reflectance not < 1.0')
      CALL ShowContinueError('See Engineering Reference for calculation procedure for solar transmittance.')
    END IF

    IF(Material(MaterNum)%TransVis+Material(MaterNum)%ReflectShadeVis >= 1.0) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      CALL ShowContinueError('Calculated visible transmittance + visible reflectance not < 1.0')
      CALL ShowContinueError('See Engineering Reference for calculation procedure for visible solar transmittance.')
    END IF

    IF(Material(MaterNum)%TransThermal+Material(MaterNum)%AbsorpThermal >= 1.0) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      CALL ShowSevereError('Thermal hemispherical emissivity plus open area fraction (1-diameter/spacing)**2' &
                           //' not < 1.0')
    END IF

  END DO

! Window Blind Materials

  IF(TotBlinds > 0) THEN
    ALLOCATE (Blind(TotBlinds))! Allocate the array Size to the number of blinds
  ENDIF

  CurrentModuleObject='WindowMaterial:Blind'
  DO Loop=1,TotBlinds

    !Call Input Get routine to retrieve material data
    CALL GetObjectItem(TRIM(CurrentModuleObject),Loop,MaterialNames,MaterialNumAlpha,MaterialProps,MaterialNumProp,IOSTAT,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

    ErrorInName=.false.
    IsBlank=.false.
    CALL VerifyName(MaterialNames(1),Material%Name,MaterNum,ErrorInName,IsBlank,TRIM(CurrentModuleObject)//' Name')
    IF (ErrorInName) THEN
      CALL ShowContinueError('...All Material names must be unique regardless of subtype.')
      ErrorsFound=.true.
      CYCLE
    ENDIF

    MaterNum=MaterNum+1
    Material(MaterNum)%Group=WindowBlind

    !Load the material derived type from the input data.

    Material(MaterNum)%Name              = MaterialNames(1)
    Blind(Loop)%Name                     = MaterialNames(1)
    Material(MaterNum)%Roughness         = Rough
    Material(MaterNum)%BlindDataPtr      = Loop
    Material(MaterNum)%ROnly         = .true.

    Blind(Loop)%MaterialNumber           = MaterNum
    IF (SameString(MaterialNames(2),'Horizontal')) THEN
      Blind(Loop)%SlatOrientation        = Horizontal
    ELSEIF (SameString(MaterialNames(2),'Vertical')) THEN
      Blind(Loop)%SlatOrientation        = Vertical
    ENDIF
    Blind(Loop)%SlatWidth                = MaterialProps(1)
    Blind(Loop)%SlatSeparation           = MaterialProps(2)
    Blind(Loop)%SlatThickness            = MaterialProps(3)
    Blind(Loop)%SlatAngle                = MaterialProps(4)
    Blind(Loop)%SlatConductivity         = MaterialProps(5)
    Blind(Loop)%SlatTransSolBeamDiff     = MaterialProps(6)
    Blind(Loop)%SlatFrontReflSolBeamDiff = MaterialProps(7)
    Blind(Loop)%SlatBackReflSolBeamDiff  = MaterialProps(8)
    Blind(Loop)%SlatTransSolDiffDiff     = MaterialProps(9)
    Blind(Loop)%SlatFrontReflSolDiffDiff = MaterialProps(10)
    Blind(Loop)%SlatBackReflSolDiffDiff  = MaterialProps(11)
    Blind(Loop)%SlatTransVisBeamDiff     = MaterialProps(12)
    Blind(Loop)%SlatFrontReflVisBeamDiff = MaterialProps(13)
    Blind(Loop)%SlatBackReflVisBeamDiff  = MaterialProps(14)
    Blind(Loop)%SlatTransVisDiffDiff     = MaterialProps(15)
    Blind(Loop)%SlatFrontReflVisDiffDiff = MaterialProps(16)
    Blind(Loop)%SlatBackReflVisDiffDiff  = MaterialProps(17)
    Blind(Loop)%SlatTransIR              = MaterialProps(18)
    Blind(Loop)%SlatFrontEmissIR         = MaterialProps(19)
    Blind(Loop)%SlatBackEmissIR          = MaterialProps(20)
    Blind(Loop)%BlindToGlassDist         = MaterialProps(21)
    Blind(Loop)%BlindTopOpeningMult      = MaterialProps(22)
    Blind(Loop)%BlindBottomOpeningMult   = MaterialProps(23)
    Blind(Loop)%BlindLeftOpeningMult     = MaterialProps(24)
    Blind(Loop)%BlindRightOpeningMult    = MaterialProps(25)
    Blind(Loop)%MinSlatAngle             = MaterialProps(26)
    Blind(Loop)%MaxSlatAngle             = MaterialProps(27)

    ! TH 2/11/2010. For CR 8010
    ! By default all blinds have fixed slat angle, new blinds with variable slat angle are created if
    !  they are used with window shading controls that adjust slat angles like ScheduledSlatAngle or BlockBeamSolar
    Blind(Loop)%SlatAngleType            = FixedSlats

    IF (Blind(Loop)%SlatWidth < Blind(Loop)%SlatSeparation) THEN
      CALL ShowWarningError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Slat Angles/Widths')
      CALL ShowContinueError(TRIM(cNumericFieldNames(1))//' ['//TRIM(RoundSigDigits(Blind(Loop)%SlatWidth,2))//  &
         '] is less than '//TRIM(cNumericFieldNames(2))//' ['//TRIM(RoundSigDigits(Blind(Loop)%SlatSeparation,2))//'].')
      CALL ShowContinueError('This will allow direct beam to be transmitted when Slat angle = 0.')
    END IF

    IF(.not. SameString(MaterialNames(2),'Horizontal') .AND. .not. SameString(MaterialNames(2),'Vertical')) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value')
      CALL ShowContinueError(TRIM(cAlphaFieldNames(2))//'="'//trim(MaterialNames(2))//'", must be '//  &
         ' Horizontal or Vertical.')
    END IF

    IF((MaterialProps(6)+MaterialProps(7) >= 1.0)) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(6))//' + '//TRIM(cNumericFieldNames(7))//' not < 1.0')
    END IF
    IF((MaterialProps(6)+MaterialProps(8) >= 1.0)) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(6))//' + '//TRIM(cNumericFieldNames(8))//' not < 1.0')
    END IF

    IF((MaterialProps(9)+MaterialProps(10) >= 1.0)) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(9))//' + '//TRIM(cNumericFieldNames(10))//' not < 1.0')
    END IF
    IF((MaterialProps(9)+MaterialProps(11) >= 1.0)) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(9))//' + '//TRIM(cNumericFieldNames(11))//' not < 1.0')
    END IF

    IF((MaterialProps(12)+MaterialProps(13) >= 1.0).OR.(MaterialProps(12)+MaterialProps(14) >= 1.0)) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(12))//' + '//TRIM(cNumericFieldNames(13))//' not < 1.0 OR')
      CALL ShowContinueError(TRIM(cNumericFieldNames(12))//' + '//TRIM(cNumericFieldNames(14))//' not < 1.0')
    END IF

    IF((MaterialProps(12)+MaterialProps(13) >= 1.0)) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(12))//' + '//TRIM(cNumericFieldNames(13))//' not < 1.0')
    END IF
    IF((MaterialProps(12)+MaterialProps(14) >= 1.0)) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(12))//' + '//TRIM(cNumericFieldNames(14))//' not < 1.0')
    END IF

    IF((MaterialProps(15)+MaterialProps(16) >= 1.0)) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(15))//' + '//TRIM(cNumericFieldNames(16))//' not < 1.0')
    END IF
    IF((MaterialProps(15)+MaterialProps(17) >= 1.0)) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(15))//' + '//TRIM(cNumericFieldNames(17))//' not < 1.0')
    END IF

    ! Require that beam and diffuse properties be the same
    IF(ABS(MaterialProps(9)-MaterialProps(6)) > 1.d-5) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(6))//' must equal '//TRIM(cNumericFieldNames(9)))
    END IF

    IF(ABS(MaterialProps(10)-MaterialProps(7)) > 1.d-5) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(7))//' must equal '//TRIM(cNumericFieldNames(10)))
    END IF

    IF(ABS(MaterialProps(11)-MaterialProps(8)) > 1.d-5) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(8))//' must equal '//TRIM(cNumericFieldNames(11)))
    END IF

    IF(ABS(MaterialProps(15)-MaterialProps(12)) > 1.d-5) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(12))//' must equal '//TRIM(cNumericFieldNames(15)))
    END IF

    IF(ABS(MaterialProps(16)-MaterialProps(13)) > 1.d-5) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(13))//' must equal '//TRIM(cNumericFieldNames(16)))
    END IF

    IF(ABS(MaterialProps(17)-MaterialProps(14)) > 1.d-5) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(14))//' must equal '//TRIM(cNumericFieldNames(17)))
    END IF

    IF((MaterialProps(18)+MaterialProps(19) >= 1.0)) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(18))//' + '//TRIM(cNumericFieldNames(19))//' not < 1.0')
    END IF
    IF((MaterialProps(18)+MaterialProps(20) >= 1.0)) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(18))//' + '//TRIM(cNumericFieldNames(20))//' not < 1.0')
    END IF

    IF(Blind(Loop)%BlindToGlassDist < 0.5d0*Blind(Loop)%SlatWidth) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
      CALL ShowContinueError(TRIM(cNumericFieldNames(21))//' is less than half of the '//    &
         trim(cNumericFieldNames(1)))
    END IF

    ! Minimum and maximum slat angles allowed by slat geometry
    IF(Blind(Loop)%SlatWidth > Blind(Loop)%SlatSeparation) THEN
      MinSlatAngGeom = ASIN(Blind(Loop)%SlatThickness/(Blind(Loop)%SlatThickness + Blind(Loop)%SlatSeparation))/DegToRadians
    ELSE
      MinSlatAngGeom = 0.0
    END IF
    MaxSlatAngGeom = 180.d0- MinSlatAngGeom

    ! Error if input slat angle not in range allowed by slat geometry
    IF((Blind(Loop)%SlatSeparation + Blind(Loop)%SlatThickness) < Blind(Loop)%SlatWidth) THEN
      IF(Blind(Loop)%SlatAngle < MinSlatAngGeom) THEN
        ErrorsFound = .true.
        CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
        CALL ShowContinueError(TRIM(cNumericFieldNames(4))//'=['//TRIM(RoundSigDigits(Blind(Loop)%SlatAngle,1))//  &
           '], is less than smallest allowed by slat dimensions and spacing, ['//  &
           TRIM(RoundSigDigits(MinSlatAngGeom,1))//'] deg.')
      ELSE IF(Blind(Loop)%SlatAngle > MaxSlatAngGeom) THEN
        ErrorsFound = .true.
        CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
        CALL ShowContinueError(TRIM(cNumericFieldNames(4))//'=['//TRIM(RoundSigDigits(Blind(Loop)%SlatAngle,1))//  &
          '], is greater than largest allowed by slat dimensions and spacing, ['//  &
          TRIM(RoundSigDigits(MinSlatAngGeom,1))//'] deg.')
      END IF
    END IF

    ! By default all Blinds are "fixed" slats.  Only with Shading Control is one considered variable and this check
    ! is now done when that happens.  9.3.2009 LKL

!    IF(Blind(Loop)%SlatAngleType == VariableSlats) THEN
!
!      ! Error if maximum slat angle less than minimum
!
!      IF(Blind(Loop)%MaxSlatAngle < Blind(Loop)%MinSlatAngle) THEN
!        ErrorsFound = .true.
!        CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
!        CALL ShowContinueError(TRIM(cNumericFieldNames(26))//'=['//TRIM(RoundSigDigits(Blind(Loop)%MinSlatAngle,1))//  &
!           '], is greater than '//trim(cNumericFieldNames(27))//'=['//  &
!           TRIM(RoundSigDigits(Blind(Loop)%MaxSlatAngle,1))//'] deg.')
!      END IF
!
!      ! Error if input slat angle not in input min/max range
!
!      IF(Blind(Loop)%MaxSlatAngle > Blind(Loop)%MinSlatAngle .AND. (Blind(Loop)%SlatAngle < Blind(Loop)%MinSlatAngle &
!          .OR. Blind(Loop)%SlatAngle > Blind(Loop)%MaxSlatAngle)) THEN
!        ErrorsFound = .true.
!        CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
!        CALL ShowContinueError(TRIM(cNumericFieldNames(4))//'=['//TRIM(RoundSigDigits(Blind(Loop)%SlatAngle,1))//  &
!           '] is outside of the input min/max range, min=['//TRIM(RoundSigDigits(Blind(Loop)%MinSlatAngle,1))//  &
!           '], max=['//TRIM(RoundSigDigits(Blind(Loop)%MaxSlatAngle,1))//'] deg.')
!      END IF
!
!      ! Error if input minimum slat angle is less than that allowed by slat geometry
!
!      IF(Blind(Loop)%MinSlatAngle < MinSlatAngGeom) THEN
!        CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
!        CALL ShowContinueError(TRIM(cNumericFieldNames(26))//'=['//TRIM(RoundSigDigits(Blind(Loop)%MinSlatAngle,1))//  &
!           '] is less than the smallest allowed by slat dimensions and spacing, min=['//  &
!           TRIM(RoundSigDigits(MinSlatAngGeom,1))//'] deg.')
!        CALL ShowContinueError('Minimum Slat Angle will be set to '//TRIM(RoundSigDigits(MinSlatAngGeom,1))//' deg.')
!        Blind(Loop)%MinSlatAngle = MinSlatAngGeom
!      END IF
!
!      ! Error if input maximum slat angle is greater than that allowed by slat geometry
!
!      IF(Blind(Loop)%MaxSlatAngle > MaxSlatAngGeom) THEN
!        CALL ShowWarningError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value combination.')
!        CALL ShowContinueError(TRIM(cNumericFieldNames(27))//'=['//TRIM(RoundSigDigits(Blind(Loop)%MaxSlatAngle,1))//  &
!           '] is greater than the largest allowed by slat dimensions and spacing, ['//  &
!           TRIM(RoundSigDigits(MaxSlatAngGeom,1))//'] deg.')
!        CALL ShowContinueError('Maximum Slat Angle will be set to '//TRIM(RoundSigDigits(MaxSlatAngGeom,1))//' deg.')
!        Blind(Loop)%MaxSlatAngle = MaxSlatAngGeom
!      END IF
!
!    END IF  ! End of check if slat angle is variable

  ENDDO

  ! EcoRoof Materials
  !PSU 2006
  CurrentModuleObject='Material:RoofVegetation'
  DO Loop=1,EcoRoofMat
    !Call Input Get Routine to retrieve material data from ecoroof

    CALL GetObjectItem(TRIM(CurrentModuleObject),Loop,MaterialNames,MaterialNumAlpha,MaterialProps,MaterialNumProp,IOSTAT,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

    ErrorInName=.false.
    IsBlank=.false.
    CALL VerifyName(MaterialNames(1),Material%Name,MaterNum,ErrorInName,IsBlank,TRIM(CurrentModuleObject)//' Name')
    IF (ErrorInName) THEN
      CALL ShowContinueError('...All Material names must be unique regardless of subtype.')
      ErrorsFound=.true.
      CYCLE
    ENDIF

    !this part is similar to the regular material
        !Load the material derived type from the input data.
    MaterNum=MaterNum+1
    Material(MaterNum)%Group = EcoRoof

       !this part is new for Ecoroof properties,
        !especially for the Plant Layer of the ecoroof
    Material(MaterNum)%HeightOfPlants    = MaterialProps(1)
    Material(MaterNum)%LAI               = MaterialProps(2)
    Material(MaterNum)%Lreflectivity     = MaterialProps(3) ! Albedo
    Material(MaterNum)%Lemissitivity     = MaterialProps(4)
    Material(MaterNum)%RStomata          = MaterialProps(5)



    Material(MaterNum)%Name = MaterialNames(1)
    !need to treat the A2 with is just the name of the soil(it is
    ! not important)
    CALL ValidateMaterialRoughness(MaterNum,MaterialNames(3),ErrorsFound)
    IF (SameString(MaterialNames(4),'Simple')) THEN
      Material(MaterNum)%EcoRoofCalculationMethod = 1
    ELSEIF (SameString(MaterialNames(4),'Advanced') .or. lAlphaFieldBlanks(4)) THEN
      Material(MaterNum)%EcoRoofCalculationMethod = 2
    ELSE
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(MaterialNames(1))//'", Illegal value')
      CALL ShowContinueError(trim(cAlphaFieldNames(4))//'="'//trim(MaterialNames(4))//'".')
      CALL ShowContinueError('...Valid values are "Simple" or "Advanced".')
      ErrorsFound=.true.
    ENDIF

    Material(MaterNum)%Thickness     = MaterialProps(6)
    Material(MaterNum)%Conductivity  = MaterialProps(7)
    Material(MaterNum)%Density       = MaterialProps(8)
    Material(MaterNum)%SpecHeat      = MaterialProps(9)
    Material(MaterNum)%AbsorpThermal = MaterialProps(10) ! emissivity
    Material(MaterNum)%AbsorpSolar   = MaterialProps(11) ! (1 - Albedo)
    Material(MaterNum)%AbsorpVisible = MaterialProps(12)
    Material(MaterNum)%Porosity      = MaterialProps(13)
    Material(MaterNum)%MinMoisture   = MaterialProps(14)
    Material(MaterNum)%InitMoisture  = MaterialProps(15)


    IF (Material(MaterNum)%Conductivity > 0.0) THEN
      NominalR(MaterNum)            = Material(MaterNum)%Thickness/Material(MaterNum)%Conductivity
      Material(MaterNum)%Resistance = NominalR(MaterNum)
    ELSE
      CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(cAlphaArgs(1))//'" is not defined correctly.')
      CALL ShowContinueError(trim(cNumericFieldNames(7))//' is <=0.')
      ErrorsFound = .TRUE.
    END IF

  END DO

  ! Thermochromic glazing group
  ! get the number of WindowMaterial:GlazingGroup:Thermochromic objects in the idf file
  CurrentModuleObject='WindowMaterial:GlazingGroup:Thermochromic'
  TotTCGlazings = GetNumObjectsFound(Trim(CurrentModuleObject))
  IF (TotTCGlazings >=1) THEN
    ! Read TC glazings
    ALLOCATE (TCGlazings(TotTCGlazings))

    DO Loop = 1, TotTCGlazings
      !Get each TCGlazings from the input processor
      CALL GetObjectItem(TRIM(CurrentModuleObject),Loop,cAlphaArgs,MaterialNumAlpha,rNumericArgs,MaterialNumProp,IOSTAT, &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

      ErrorInName=.false.
      IsBlank=.false.

      ! Verify unique names
      CALL VerifyName(cAlphaArgs(1),TCGlazings%Name,Loop-1,ErrorInName,IsBlank,TRIM(CurrentModuleObject)//' Name')
      IF (ErrorInName) THEN
        CALL ShowContinueError('...All Thermochromic Glazing names must be unique regardless of subtype.')
        ErrorsFound=.true.
        CYCLE
      ENDIF

      IF ( MaterialNumProp + 1 /= MaterialNumAlpha ) THEN
        CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(cAlphaArgs(1))//'" is not defined correctly.')
        CALL ShowContinueError('Check number of '//Trim(cAlphaFieldNames(2))// &
                               ' compared to number of '//Trim(cNumericFieldNames(1)))
        ErrorsFound=.true.
        CYCLE
      ENDIF

      !Allocate arrays
      ALLOCATE (TCGlazings(Loop)%SpecTemp(MaterialNumProp))
      ALLOCATE (TCGlazings(Loop)%LayerName(MaterialNumProp))
      ALLOCATE (TCGlazings(Loop)%LayerPoint(MaterialNumProp))
      TCGlazings(Loop)%SpecTemp = 0.0
      TCGlazings(Loop)%LayerName = ' '
      TCGlazings(Loop)%LayerPoint = 0

      TCGlazings(Loop)%Name = cAlphaArgs(1)
      TCGlazings(Loop)%NumGlzMat = MaterialNumProp

      DO iTC = 1, MaterialNumProp
        TCGlazings(Loop)%SpecTemp(iTC)= rNumericArgs(iTC)
        TCGlazings(Loop)%LayerName(iTC)= cAlphaArgs(1+iTC)

        ! Find this glazing material in the material list
        iMat = FindIteminList(cAlphaArgs(1+iTC),Material%Name,TotMaterials)
        IF (iMat /= 0) THEN
          !TC glazing
          Material(iMat)%SpecTemp = rNumericArgs(iTC)
          Material(iMat)%TCParent = Loop
          TCGlazings(Loop)%LayerPoint(iTC) = iMat

          !test that named material is of the right type
          IF( Material(iMat)%Group /= WindowGlass) THEN
            !CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(cAlphaArgs(1))//'" is not defined correctly.')
            !CALL ShowContinueError('Material named: '//Trim(cAlphaArgs(1+iTC))//' is not a window glazing ')    !RS: Secret Search String
            WRITE(DebugFile,*) TRIM(CurrentModuleObject)//'="'//TRIM(cAlphaArgs(1))//'" is not defined correctly.'
            WRITE(DebugFile,*) 'Material named: '//TRIM(cAlphaArgs(1+iTC))//' is not a window glazing '
            ErrorsFound=.true.
          ENDIF

        ELSE ! thow error because not found
          !CALL ShowSevereError(TRIM(CurrentModuleObject)//'="'//trim(cAlphaArgs(1))//'" is not defined correctly.')
          !CALL ShowContinueError('Material named: '//Trim(cAlphaArgs(1+iTC))//' was not found ')    !RS: Secret Search String
          WRITE(DebugFile,*) TRIM(CurrentModuleObject)//'="'//TRIM(cAlphaArgs(1))//'" is not defined correctly.'
          WRITE(DebugFile,*) 'Material name: '//TRIM(cAlphaArgs(1+iTC))//' was not found '
          ErrorsFound=.true.
        ENDIF
      ENDDO
    ENDDO
  ENDIF


  cCurrentModuleObject='WindowMaterial:SimpleGlazingSystem'
  DO Loop=1, TotSimpleWindow

    CALL GetObjectItem(TRIM(cCurrentModuleObject), Loop, cAlphaArgs, MaterialNumAlpha, &
                      rNumericArgs, MaterialNumProp,IOSTAT,     &
                      AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                      AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)
    ErrorInName=.false.
    IsBlank=.false.
    CALL VerifyName(cAlphaArgs(1), Material%Name,MaterNum,ErrorInName,IsBlank,TRIM(cCurrentModuleObject)//' Name')
    IF (ErrorInName) THEN
      CALL ShowContinueError('...All Material names must be unique regardless of subtype.')
      ErrorsFound=.true. 
      CYCLE
    ENDIF
    MaterNum=MaterNum+1
    Material(MaterNum)%Group   = WindowSimpleGlazing
    Material(MaterNum)%Name    = cAlphaArgs(1)
    Material(MaterNum)%SimpleWindowUfactor   = rNumericArgs(1)
    Material(MaterNum)%SimpleWindowSHGC      = rNumericArgs(2)
    IF (.not. lNumericFieldBlanks(3)) THEN
      Material(MaterNum)%SimpleWindowVisTran = rNumericArgs(3)
      Material(MaterNum)%SimpleWindowVTinputByUser = .TRUE.
    ENDIF

    CALL SetupSimpleWindowGlazingSystem(MaterNum)

  ENDDO

  IF (SolutionAlgo /= UseCTF .and. SolutionAlgo /= UseCondFD .and. SolutionAlgo /= UseCondFDSimple) THEN
    ! If CTF and the Constructions report is requested, then this info has already been written.
    CALL ScanForReports('Constructions',DoReport,'Constructions')
    IF (DoReport) THEN

      WRITE(OutputFileInits,108)

      Do MaterNum=1,TotMaterials

         WRITE(OutputFileInits,111) Trim(Material(MaterNum)%Name),TRIM(RoundSigDigits(NominalR(MaterNum),4))

      end do

    End If
  End If

  CALL ScanForReports('Constructions',DoReport,'Materials')

  IF (DoReport) THEN

    Write(OutputFileInits,'(A)') '! <Material>,Material Name,ThermalResistance {m2-K/w},Roughness,Thickness {m},'//  &
       'Conductivity {w/m-K},Density {kg/m3},Specific Heat {J/kg-K},Absorptance:Thermal,Absorptance:Solar,Absorptance:Visible'

    Write(OutputFileInits,'(A)') '! <Material:Air>,Material Name,ThermalResistance {m2-K/w}'

    Do MaterNum=1,TotMaterials

      SELECT CASE (Material(MaterNum)%Group)
        CASE (Air)
          Write(OutputFileInits,702) TRIM(Material(MaterNum)%Name),  &
                                       TRIM(RoundSigDigits(Material(MaterNum)%Resistance,4))
        CASE DEFAULT
          Write(OutputFileInits,701) TRIM(Material(MaterNum)%Name),TRIM(RoundSigDigits(Material(MaterNum)%Resistance,4)),  &
                                     TRIM(DisplayMaterialRoughness(Material(MaterNum)%Roughness)),  &
                                     TRIM(RoundSigDigits(Material(MaterNum)%Thickness,4)),   &
                                     TRIM(RoundSigDigits(Material(MaterNum)%Conductivity,3)),  &
                                     TRIM(RoundSigDigits(Material(MaterNum)%Density,3)),  &
                                     TRIM(RoundSigDigits(Material(MaterNum)%SpecHeat,3)),  &
                                     TRIM(RoundSigDigits(Material(MaterNum)%AbsorpThermal,4)),  &
                                     TRIM(RoundSigDigits(Material(MaterNum)%AbsorpSolar,4)),  &
                                     TRIM(RoundSigDigits(Material(MaterNum)%AbsorpVisible,4))

      END SELECT

    end do

  End If

!  FORMATS.
  108 FORMAT('! <Material Nominal Resistance>, Material Name,  Nominal R')
  111 FORMAT('Material Nominal Resistance',2(',',A))
 701  FORMAT(' Material',10(',',A))
 702  FORMAT(' Material:AirGap',2(',',A))

  IF (AnyEnergyManagementSystemInModel) THEN ! setup surface property EMS actuators

    DO MaterNum=1,TotMaterials
      IF (Material(MaterNum)%Group /= RegularMaterial) CYCLE
      CALL SetupEMSActuator('Material', Material(MaterNum)%Name, &
                            'Surface Property Solar Absorptance', '[ ]', &
                            Material(MaterNum)%AbsorpSolarEMSOverrideOn , &
                            Material(MaterNum)%AbsorpSolarEMSOverride )
      CALL SetupEMSActuator('Material', Material(MaterNum)%Name, &
                            'Surface Property Thermal Absorptance', '[ ]', &
                            Material(MaterNum)%AbsorpThermalEMSOverrideOn , &
                            Material(MaterNum)%AbsorpThermalEMSOverride )
      CALL SetupEMSActuator('Material', Material(MaterNum)%Name, &
                            'Surface Property Visible Absorptance', '[ ]', &
                            Material(MaterNum)%AbsorpVisibleEMSOverrideOn , &
                            Material(MaterNum)%AbsorpVisibleEMSOverride )
    ENDDO
  ENDIF

  RETURN

END SUBROUTINE GetMaterialData

SUBROUTINE GetWindowGlassSpectralData(ErrorsFound)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Fred Winkelmann
          !       DATE WRITTEN   May 2000
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! Gets spectral data (transmittance, front reflectance, and back
          ! reflectance at normal incidence vs. wavelength) for glass

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE General, ONLY: TrimSigDigits

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  LOGICAL, INTENT(INOUT) :: ErrorsFound  ! set to true if errors found in input

          ! SUBROUTINE PARAMETER DEFINITIONS:
  CHARACTER(len=*), PARAMETER :: RoutineName='GetWindowGlassSpectralData: '

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:

  INTEGER :: IOStat               ! IO Status when calling get input subroutine
  CHARACTER(len=MaxNameLength),DIMENSION(1) &
          :: SpecDataNames        ! Spectral data alpha names
  INTEGER :: SpecDataNumAlpha ! Number of spectral data alpha names being passed
  INTEGER :: SpecDataNumProp  ! Number of spectral data properties being passed
  REAL(r64), ALLOCATABLE, DIMENSION(:) :: SpecDataProps !Temporary array to transfer spectal data properties
  INTEGER :: Loop
  LOGICAL :: ErrorInName
  LOGICAL :: IsBlank
  INTEGER :: LamNum            ! Wavelength number
  INTEGER :: TotLam            ! Total wavelengths
  REAL(r64)    :: Lam               ! Wavelength (microns)
  REAL(r64)    :: Tau,RhoF,RhoB     ! Transmittance, front reflectance, back reflectance

  INTEGER :: DebugFile       =150 !RS: Debugging file denotion, hopefully this works.
    
  OPEN(unit=DebugFile,file='Debug.txt')    !RS: Debugging

  CurrentModuleObject='MaterialProperty:GlazingSpectralData'
  TotSpectralData=GetNumObjectsFound(TRIM(CurrentModuleObject))
  ALLOCATE (SpectralData(TotSpectralData))
  IF (TotSpectralData > 0) ALLOCATE(SpecDataProps(MaxSpectralDataElements*4))

  DO Loop=1,TotSpectralData

    ! Call Input Get routine to retrieve spectral data
    ! Name is followed by up to 450 sets of normal-incidence measured values of
    ! [wavelength (microns), transmittance, front reflectance, back reflectance] for
    ! wavelengths covering the short-wave solar spectrum (from about 0.25 to 2.5 microns)
    CALL GetObjectItem(TRIM(CurrentModuleObject),Loop,SpecDataNames,SpecDataNumAlpha,SpecDataProps,SpecDataNumProp,IOStat,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

    ErrorInName=.false.
    IsBlank=.false.
    CALL VerifyName(SpecDataNames(1),SpectralData%Name,Loop,ErrorInName,IsBlank,TRIM(CurrentModuleObject)//' Name')
    IF (ErrorInName) THEN
      ErrorsFound=.true.
      CYCLE
    ENDIF

    !Load the spectral data derived type from the input data.
    SpectralData(Loop)%Name = SpecDataNames(1)
    TotLam = SpecDataNumProp/4
    IF (MOD(SpecDataNumProp,4) /= 0) THEN
      CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(SpecDataNames(1))//'" invalid set.')
      CALL ShowContinueError('... set not even multiple of 4 items (Wavelength,Trans,ReflFront,ReflBack),'//  &
                             'number of items in dataset = '//trim(TrimSigDigits(SpecDataNumProp)))
      CALL ShowContinueError('... remainder after div by 4 = '//trim(TrimSigDigits(MOD(SpecDataNumProp,4)))//  &
                               ', remainder items will be set to 0.0')
      SpecDataProps(SpecDataNumProp+1:MIN(SpecDataNumProp+4,MaxSpectralDataElements*4))=0.0
    ENDIF
    IF(TotLam > MaxSpectralDataElements) THEN
      ErrorsFound = .true.
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(SpecDataNames(1))//'" invalid set.')
      CALL ShowContinueError('... More than max ['//trim(TrimSigDigits(MaxSpectralDataElements))//   &
                 '] (Wavelength,Trans,ReflFront,ReflBack) entries in set.')
      CYCLE
    END IF
    SpectralData(Loop)%NumOfWavelengths = TotLam

    ALLOCATE(SpectralData(Loop)%WaveLength(TotLam))   ! Wavelength (microns)
    ALLOCATE(SpectralData(Loop)%Trans(TotLam))        ! Transmittance at normal incidence
    ALLOCATE(SpectralData(Loop)%ReflFront(TotLam))    ! Front reflectance at normal incidence
    ALLOCATE(SpectralData(Loop)%ReflBack(TotLam))     ! Back reflectance at normal incidence

    DO LamNum = 1,TotLam
      SpectralData(Loop)%WaveLength(LamNum) = SpecDataProps(4*LamNum-3)
      SpectralData(Loop)%Trans(LamNum)      = SpecDataProps(4*LamNum-2)
      ! Following is needed since angular calculation in subr TransAndReflAtPhi
      ! fails for Trans = 0.0
      IF(SpectralData(Loop)%Trans(LamNum) < 0.001d0) SpectralData(Loop)%Trans(LamNum) = 0.001d0
      SpectralData(Loop)%ReflFront(LamNum)  = SpecDataProps(4*LamNum-1)
      SpectralData(Loop)%ReflBack(LamNum)   = SpecDataProps(4*LamNum)
    END DO

    ! Check integrity of the spectral data
    DO LamNum = 1,TotLam
      Lam = SpectralData(Loop)%WaveLength(LamNum)
      Tau = SpectralData(Loop)%Trans(LamNum)
      RhoF = SpectralData(Loop)%ReflFront(LamNum)
      RhoB = SpectralData(Loop)%ReflBack(LamNum)
      IF(LamNum < TotLam) THEN
        IF (SpectralData(Loop)%WaveLength(LamNum+1) <= Lam) THEN
          ErrorsFound = .true.
          !CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(SpecDataNames(1))//'" invalid set.')
          !CALL ShowContinueError('... Wavelengths not in increasing order. '//  &
          !              'at wavelength#='//trim(TrimSigDigits(LamNum))//', value=['//trim(TrimSigDigits(Lam,4))//  &
          !              '], next is ['//trim(TrimSigDigits(SpectralData(Loop)%WaveLength(LamNum+1),4))//'].')  !RS: Secret Search String
          WRITE(DebugFile,*) RoutineName//TRIM(CurrentModuleObject)//'="'//TRIM(SpecDataNames(1))//'" invalid set.'
          WRITE(DebugFile,*) '... Wavelengths not in increasing order. At wavelength#='//TRIM(TrimSigDigits(LamNum))// &
            ', value=['//TRIM(TrimSigDigits(Lam,4))//'], next is ['//TRIM(TrimSigDigits(SpectralData(Loop)%WaveLength(LamNum+1),4))//'].'
        END IF
      END IF

      IF(Lam < 0.1d0 .OR. Lam > 4.0d0) THEN
        ErrorsFound = .true.
        !CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(SpecDataNames(1))//'" invalid value.')
        !CALL ShowContinueError('... A wavelength is not in the range 0.1 to 4.0 microns; '//  &
        !                  'at wavelength#='//trim(TrimSigDigits(LamNum))//', value=['//trim(TrimSigDigits(Lam,4))//  &
        !         '].')  !RS: Secret Search String
        WRITE(DebugFile,*) RoutineName//TRIM(CurrentModuleObject)//'="'//TRIM(SpecDataNames(1))//'" invalid value.'
        WRITE(DebugFile,*) '... A wavelength is not in the range 0.1 to 4.0 microns; at wavelength#='// &
            TRIM(TrimSigDigits(LamNum))//', value=['//TRIM(TrimSigDigits(Lam,4))//'].'
      END IF

    ! TH 2/15/2011. CR 8343
    ! IGDB (International Glazing Database) does not meet the above strict restrictions.
    !  Relax rules to allow directly use of spectral data from IGDB
      IF(Tau > 1.01d0) THEN
        ErrorsFound = .true.
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(SpecDataNames(1))//'" invalid value.')
        CALL ShowContinueError('... A transmittance is > 1.0; '//  &
           'at wavelength#='//trim(TrimSigDigits(LamNum))//', value=['//trim(TrimSigDigits(Tau,4))//'].')
      END IF

      IF(RhoF < 0.0d0 .OR. RhoF > 1.01d0 .OR. RhoB < 0.0d0 .OR. RhoB > 1.01d0) THEN
        ErrorsFound = .true.
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(SpecDataNames(1))//'" invalid value.')
        CALL ShowContinueError('... A reflectance is < 0.0 or > 1.0; '//  &
           'at wavelength#='//trim(TrimSigDigits(LamNum))//', value=['//trim(TrimSigDigits(RhoF,4))//'].')
      END IF

      IF((Tau + RhoF) > 1.01d0 .OR. (Tau + RhoB) > 1.01d0) THEN
        ErrorsFound = .true.
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(SpecDataNames(1))//'" invalid value.')
        CALL ShowContinueError('... Transmittance + reflectance) > 1.0 for an entry; '//  &
              'at wavelength#='//trim(TrimSigDigits(LamNum))//', value(Tau+RhoF)=['//trim(TrimSigDigits((Tau + RhoF),4))//  &
              '], value(Tau+RhoB)=['//trim(TrimSigDigits((Tau + RhoB),4))//'].')
      END IF

    END DO

  END DO

  IF (TotSpectralData > 0) DEALLOCATE(SpecDataProps)

  RETURN
END SUBROUTINE GetWindowGlassSpectralData

SUBROUTINE ValidateMaterialRoughness(MaterNum,Roughness,ErrorsFound)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Linda K. Lawrie
          !       DATE WRITTEN   April 1999
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! This subroutine compares the input Roughness value against the
          ! valid values and sets the correct value in the Material Data Structure.

          ! METHODOLOGY EMPLOYED:
          ! Error message provided if not valid.

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
          ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  INTEGER, INTENT(IN)            :: MaterNum   ! Which Material number being validated.
  CHARACTER(len=*), INTENT(IN)   :: Roughness   ! Roughness String
  LOGICAL, INTENT(INOUT)         :: ErrorsFound ! If errors found

          ! SUBROUTINE PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:


    !Select the correct Number for the associated ascii name for the roughness type
    IF (SameString(Roughness,'VeryRough'))    Material(MaterNum)%Roughness=VeryRough
    IF (SameString(Roughness,'Rough'))        Material(MaterNum)%Roughness=Rough
    IF (SameString(Roughness,'MediumRough'))  Material(MaterNum)%Roughness=MediumRough
    IF (SameString(Roughness,'MediumSmooth')) Material(MaterNum)%Roughness=MediumSmooth
    IF (SameString(Roughness,'Smooth'))       Material(MaterNum)%Roughness=Smooth
    IF (SameString(Roughness,'VerySmooth'))   Material(MaterNum)%Roughness=VerySmooth

    ! Was it set?
    IF (Material(MaterNum)%Roughness == 0) THEN
      CALL ShowSevereError('Material='//TRIM(Material(MaterNum)%Name)//',Illegal Roughness='//TRIM(Roughness))
      ErrorsFound=.true.
    ENDIF

  RETURN

END SUBROUTINE ValidateMaterialRoughness

SUBROUTINE GetConstructData(ErrorsFound)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Richard Liesen
          !       DATE WRITTEN   September 1997
          !       MODIFIED       January 2003, FCW: accommodate between-glass shading device
          !                      July 2009, TH: added constructions defined with F and C factors
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! This file reads the input through the input processor for Constructions.
          ! Data read in this routine is stored in a derived type (Construct)
          ! defined in the DataHeatBalance module.
          ! This subroutine only sets those parameters which must be obtained
          ! from the input file--all other portions of the Construct derived
          ! type are set during the initializations.

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE DataStringGlobals

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  LOGICAL, INTENT(INOUT) :: ErrorsFound ! If errors found in input

          ! SUBROUTINE PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:

  INTEGER :: ConstrNum        ! Counter to keep track of the construction number
  INTEGER :: Layer            ! loop index for each of the construction layers
  INTEGER :: ConstructNumAlpha ! Number of construction alpha names being passed
  INTEGER :: DummyNumProp      ! dummy variable for properties being passed
  INTEGER :: IOStat            ! IO Status when calling get input subroutine
  CHARACTER(len=MaxNameLength),DIMENSION(0:MaxLayersInConstruct) &
          :: ConstructAlphas ! Construction Alpha names defined
  REAL(r64), DIMENSION(4) :: DummyProps !Temporary array to transfer construction properties
  LOGICAL :: ErrorInName
  LOGICAL :: IsBlank
  INTEGER :: Loop
  INTEGER :: TotRegConstructs   ! Number of "regular" constructions (no embedded sources or sinks and

  INTEGER :: TotFfactorConstructs  ! Number of slabs-on-grade or underground floor constructions defined with F factors
  INTEGER :: TotCfactorConstructs  ! Number of underground wall constructions defined with C factors

  INTEGER :: TotSourceConstructs   ! Number of constructions with embedded sources or sinks
  INTEGER :: TotWindow5Constructs  ! Number of constructions from Window5 data file
  LOGICAL :: ConstructionFound ! True if input window construction name is found in the
                               !  Window5 data file
  LOGICAL :: EOFonW5File       ! True if EOF encountered reading Window5 data file
  LOGICAL :: NoRegularMaterialsUsed=.true.

  INTEGER :: iMatGlass         ! number of glass layers
  CHARACTER(len=MaxNameLength), ALLOCATABLE, DIMENSION(:) :: WConstructNames

  INTEGER :: DebugFile       =150 !RS: Debugging file denotion, hopefully this works.
    
  OPEN(unit=DebugFile,file='Debug.txt')    !RS: Debugging

       ! FLOW:

     !Get the Total number of Constructions from the input
  TotRegConstructs    = GetNumObjectsFound('Construction')
  TotSourceConstructs = GetNumObjectsFound('Construction:InternalSource')

  TotFfactorConstructs = GetNumObjectsFound('Construction:FfactorGroundFloor')
  TotCfactorConstructs = GetNumObjectsFound('Construction:CfactorUndergroundWall')

  TotWindow5Constructs = GetNumObjectsFound('Construction:WindowDataFile')
  ALLOCATE(WConstructNames(TotWindow5Constructs))
  WConstructNames=' '

  TotConstructs       = TotRegConstructs + TotFfactorConstructs + TotCfactorConstructs + TotSourceConstructs
  ALLOCATE(NominalU(TotConstructs))
  NominalU=0.0

     !Allocate the array to the number of constructions/initialize selected variables
  ALLOCATE(Construct(TotConstructs))
     !Note: If TotWindow5Constructs > 0, additional constructions are created in
     !subr. SearchWindow5DataFile corresponding to those found on the data file.
     !Initialize CTF and History terms.
  Construct%NumCTFTerms    = 0
  Construct%NumHistories   = 0

     !Initialize some heat source/sink variables
  Construct%SourceSinkPresent  = .FALSE. ! "default" is no source or sink present
  Construct%SolutionDimensions = 1       ! "default" is 1-D heat transfer
  Construct%SourceAfterLayer   = 0       ! this has no meaning if a source/sink is not present
  Construct%TempAfterLayer     = 0       ! this has no meaning if a source/sink is not present
  Construct%ThicknessPerpend   = 0.0     ! this has no meaning if a source/sink is not present

  Construct%W5FrameDivider = 0
  Construct%FromWindow5DataFile = .FALSE.

  ConstrNum=0

  CurrentModuleObject='Construction'
  DO Loop = 1, TotRegConstructs ! Loop through all constructs in the input...

      !Get the object names for each construction from the input processor
    CALL GetObjectItem(TRIM(CurrentModuleObject),Loop,ConstructAlphas,ConstructNumAlpha,DummyProps,DummyNumProp,IOSTAT,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

    ErrorInName=.false.
    IsBlank=.false.
    CALL VerifyName(ConstructAlphas(0),Construct%Name,ConstrNum,ErrorInName,IsBlank,TRIM(CurrentModuleObject)//' Name')
    IF (ErrorInName) THEN
      ErrorsFound=.true.
      CYCLE
    ENDIF

    !Glass layer counter
    iMatGlass = 0

    ConstrNum=ConstrNum+1
      !Assign Construction name to the Derived Type using the zeroth position of the array
    Construct(ConstrNum)%Name = ConstructAlphas(0)

      !Set the total number of layers for the construction
    Construct(ConstrNum)%TotLayers = ConstructNumAlpha-1

      ! Loop through all of the layers of the construct to match the material names.
      ! The loop index is the number minus 1
    DO Layer = 1, ConstructNumAlpha-1

        !Find the material in the list of materials

      Construct(ConstrNum)%LayerPoint(Layer) = FindIteminList(ConstructAlphas(Layer),Material%Name,TotMaterials)

      ! count number of glass layers
      IF (Construct(ConstrNum)%LayerPoint(Layer)>0) THEN
        IF (Material(Construct(ConstrNum)%LayerPoint(Layer))%Group == WindowGlass) iMatGlass = iMatGlass + 1
      ENDIF

      IF (Construct(ConstrNum)%LayerPoint(Layer) == 0) THEN
        !This may be a TC GlazingGroup
        Construct(ConstrNum)%LayerPoint(Layer) = FindIteminList(ConstructAlphas(Layer),TCGlazings%Name,TotTCGlazings)

        IF (Construct(ConstrNum)%LayerPoint(Layer)>0) THEN
          !reset layer pointer to the first glazing in the TC GlazingGroup
          Construct(ConstrNum)%LayerPoint(Layer)=TCGlazings(Construct(ConstrNum)%LayerPoint(Layer))%LayerPoint(1)
          Construct(ConstrNum)%TCLayer = Construct(ConstrNum)%LayerPoint(Layer)
          IF (Material(Construct(ConstrNum)%LayerPoint(Layer))%Group == WindowGlass) iMatGlass = iMatGlass + 1
          Construct(ConstrNum)%TCFlag = 1
          Construct(ConstrNum)%TCMasterConst = ConstrNum
          Construct(ConstrNum)%TCGlassID = iMatGlass  ! the TC glass layer ID
          Construct(ConstrNum)%TCLayerID = Layer
          Construct(ConstrNum)%TypeIsWindow = .True.
        ENDIF
      ENDIF

      IF (Construct(ConstrNum)%LayerPoint(Layer) == 0) THEN
        !CALL ShowSevereError('Did not find matching material for '//TRIM(CurrentModuleObject)//' '//  &
        !   TRIM(Construct(ConstrNum)%Name)//', missing material = '//TRIM(ConstructAlphas(Layer)))  !RS: Secret Search String
        WRITE(DebugFile,*) 'Did not find matching material for '//TRIM(CurrentModuleObject)//' '// &
            TRIM(Construct(ConstrNum)%Name)//', missing material for '//TRIM(ConstructAlphas(Layer))
        ErrorsFound=.true.
      ELSE
        NominalU(ConstrNum)=NominalU(ConstrNum)+NominalR(Construct(ConstrNum)%LayerPoint(Layer))
        IF (Material(Construct(ConstrNum)%LayerPoint(Layer))%Group ==  RegularMaterial  &
            .and. .not. Material(Construct(ConstrNum)%LayerPoint(Layer))%ROnly) THEN
          NoRegularMaterialsUsed=.false.
        ENDIF
      ENDIF


    END DO  ! ...end of the Layer DO loop

  END DO  ! ...end of Regular Construction DO loop

  TotRegConstructs = ConstrNum

  ! Added TH 7/2009 for underground walls and floors constructions
  IF (TotFfactorConstructs + TotCfactorConstructs >= 1) THEN
    CALL CreateFCfactorConstructions(ConstrNum,ErrorsFound)
    IF (ErrorsFound) THEN
        !CALL ShowSevereError('Errors found in creating the constructions defined with Ffactor or Cfactor method')   !RS: Secret Search String
        WRITE(DebugFile,*) 'Errors found in creating the constructions defined with Ffactor or Cfactor method'
    ENDIF
    TotRegConstructs = TotRegConstructs + TotFfactorConstructs + TotCfactorConstructs
  ENDIF

  ConstrNum=0

  CurrentModuleObject='Construction:InternalSource'
  DO Loop = 1, TotSourceConstructs  ! Loop through all constructs with sources in the input...

      !Get the object names for each construction from the input processor
    CALL GetObjectItem(TRIM(CurrentModuleObject),Loop,ConstructAlphas,ConstructNumAlpha,DummyProps,DummyNumProp,IOSTAT,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

    ErrorInName=.false.
    IsBlank=.false.
    CALL VerifyName(ConstructAlphas(0),Construct%Name,TotRegConstructs+ConstrNum,ErrorInName,  &
               IsBlank,TRIM(CurrentModuleObject)//' Name')
    IF (ErrorInName) THEN
      ErrorsFound=.true.
      CYCLE
    ENDIF

    ConstrNum=ConstrNum+1
      !Assign Construction name to the Derived Type using the zeroth position of the array
    Construct(TotRegConstructs+ConstrNum)%Name = ConstructAlphas(0)

          ! Obtain the source/sink data
    IF (DummyNumProp /= 4) THEN
      CALL ShowSevereError(TRIM(CurrentModuleObject)//': Wrong number of numerical inputs for '//TRIM(Construct(ConstrNum)%Name))
      ErrorsFound = .TRUE.
    END IF
    Construct(TotRegConstructs+ConstrNum)%SourceSinkPresent  = .TRUE.
    Construct(TotRegConstructs+ConstrNum)%SourceAfterLayer   = INT(DummyProps(1))
    Construct(TotRegConstructs+ConstrNum)%TempAfterLayer     = INT(DummyProps(2))
    Construct(TotRegConstructs+ConstrNum)%SolutionDimensions = INT(DummyProps(3))
    IF ( (Construct(TotRegConstructs+ConstrNum)%SolutionDimensions < 1) .OR. &
         (Construct(TotRegConstructs+ConstrNum)%SolutionDimensions > 2) ) THEN
      CALL ShowWarningError('Construction:InternalSource must be either 1- or 2-D.  Reset to 1-D solution.')
      CALL ShowContinueError('Construction='//TRIM(Construct(TotRegConstructs+ConstrNum)%Name)//' is affected.')
      Construct(TotRegConstructs+ConstrNum)%SolutionDimensions=1
    ENDIF
    Construct(TotRegConstructs+ConstrNum)%ThicknessPerpend   = DummyProps(4)

      !Set the total number of layers for the construction
    Construct(TotRegConstructs+ConstrNum)%TotLayers = ConstructNumAlpha-1
    IF (Construct(TotRegConstructs+ConstrNum)%TotLayers <= 1) THEN
      CALL ShowSevereError('Construction '//TRIM(Construct(TotRegConstructs+ConstrNum)%Name)// &
                           ' has an internal source or sink and thus must have more than a single layer')
      ErrorsFound=.true.
    END IF
    IF ( (Construct(TotRegConstructs+ConstrNum)%SourceAfterLayer >= Construct(TotRegConstructs+ConstrNum)%TotLayers) .OR. &
         (Construct(TotRegConstructs+ConstrNum)%SourceAfterLayer <= 0) ) THEN
      CALL ShowWarningError('Construction '//TRIM(Construct(TotRegConstructs+ConstrNum)%Name)// &
                            ' must have a source that is between two layers')
      CALL ShowContinueError('The source after layer parameter has been set to one less than the number of layers.')
      Construct(TotRegConstructs+ConstrNum)%SourceAfterLayer = Construct(TotRegConstructs+ConstrNum)%TotLayers - 1
    END IF
    IF ( (Construct(TotRegConstructs+ConstrNum)%TempAfterLayer >= Construct(TotRegConstructs+ConstrNum)%TotLayers) .OR. &
         (Construct(TotRegConstructs+ConstrNum)%TempAfterLayer <= 0) ) THEN
      CALL ShowWarningError('Construction '//TRIM(Construct(TotRegConstructs+ConstrNum)%Name)// &
                            ' must have a temperature calculation that is between two layers')
      CALL ShowContinueError('The temperature calculation after layer parameter has been set '//  &
         'to one less than the number of layers.')
      Construct(TotRegConstructs+ConstrNum)%TempAfterLayer = Construct(TotRegConstructs+ConstrNum)%TotLayers - 1
    END IF

      ! Loop through all of the layers of the construct to match the material names.
      ! The loop index is the number minus 1
    DO Layer = 1, ConstructNumAlpha-1

        !Find the material in the list of materials

      Construct(TotRegConstructs+ConstrNum)%LayerPoint(Layer) = FindIteminList(ConstructAlphas(Layer),Material%Name,TotMaterials)

      IF (Construct(TotRegConstructs+ConstrNum)%LayerPoint(Layer) == 0) THEN
        !CALL ShowSevereError('Did not find matching material for '//TRIM(CurrentModuleObject)//' '//  &
        !   TRIM(Construct(ConstrNum)%Name)// &
        !   ', missing material = '//TRIM(ConstructAlphas(Layer)))   !RS: Secret Search String
        WRITE(DebugFile,*) 'Did not find matching material for '//TRIM(CurrentModuleObject)//' '// &
            TRIM(Construct(ConstrNum)%Name)//', missing material = '//TRIM(ConstructAlphas(Layer))
        ErrorsFound=.true.
      ELSE
        NominalU(TotRegConstructs+ConstrNum)=NominalU(TotRegConstructs+ConstrNum)+ &
                                             NominalR(Construct(TotRegConstructs+ConstrNum)%LayerPoint(Layer))
        IF (Material(Construct(TotRegConstructs+ConstrNum)%LayerPoint(Layer))%Group ==  RegularMaterial  &
            .and. .not. Material(Construct(TotRegConstructs+ConstrNum)%LayerPoint(Layer))%ROnly) THEN
          NoRegularMaterialsUsed=.false.
        ENDIF
      ENDIF

    END DO  ! ...end of the Layer DO loop

  END DO  ! ...end of Source Construction DO loop

  TotSourceConstructs = ConstrNum
  TotConstructs       = TotRegConstructs + TotSourceConstructs

  IF (TotConstructs > 0 .and. NoRegularMaterialsUsed) THEN
    !CALL ShowSevereError('This building has no thermal mass which can cause an unstable solution.')
    !CALL ShowContinueError('Use Material object for all opaque material definitions except very light insulation layers.')  !RS: Secret Search String
    WRITE(DebugFile,*) 'This building has no thermal mass which can cause an unstable solution.'
    WRITE(DebugFile,*) 'Use Material object for all opaque material definitions except very light insulation layers.'
  ENDIF

!-------------------------------------------------------------------------------
  ConstrNum = 0

  CurrentModuleObject='Construction:WindowDataFile'
  DO Loop = 1, TotWindow5Constructs  ! Loop through all Window5 constructions. These constructions come
                                     ! from the Window5 data file and can be referenced only by windows

      !Get the object names for each construction from the input processor
    CALL GetObjectItem(TRIM(CurrentModuleObject),Loop,ConstructAlphas,ConstructNumAlpha,DummyProps,DummyNumProp,IOSTAT,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

    ErrorInName=.false.
    IsBlank=.false.
    CALL VerifyName(ConstructAlphas(0),WConstructNames,ConstrNum,ErrorInName,IsBlank,  &
                     TRIM(CurrentModuleObject)//' Name')
    IF (ErrorInName .and. .not. IsBlank) THEN
      CALL ShowContinueError('...first instance will be used.')
      CYCLE
    ENDIF
    IF (IsBlank) THEN
      ErrorsFound=.true.
      CYCLE
    ENDIF

    ConstrNum=ConstrNum+1
    WConstructNames(ConstrNum)=ConstructAlphas(0)

          ! Obtain the data
    IF (DummyNumProp /= 0) THEN
      CALL ShowSevereError('Construction From Window5 Data File: there should be no numerical inputs for '//  &
                           TRIM(ConstructAlphas(0)))
      ErrorsFound = .TRUE.
      CYCLE
    END IF

    ! See if this construction is in the W5DataFile produced by the WINDOW 5 program;
    ! if so, ConstructionFound will be set to true and the Material objects
    ! associated with the construction will be created in subr. SearchWindow5DataFile.
    ! (If the matching construction on the Window5 data file has two glazing systems, a
    ! second construction and its associated materials will be created in subr.
    ! SearchWindow5DataFile and TotConstructs WILL BE INCREMENTED BY 1 in that routine.
    ! A FrameAndDivider object will also be created if window on data file has a
    ! frame or divider.)

    IF (ConstructAlphas(1) == ' ') THEN
      FullName=TRIM(CurrentWorkingFolder)//'Window5DataFile.dat'
    ELSE
      FullName=ConstructAlphas(1)
    ENDIF
    CALL DisplayString('Searching Window5 data file for Construction=' &
      //TRIM(ConstructAlphas(0)))

    CALL SearchWindow5DataFile(TRIM(FullName),ConstructAlphas(0),ConstructionFound,EOFonW5File,ErrorsFound)

    IF(EOFonW5File.OR..NOT.ConstructionFound) THEN
      CALL DisplayString('--Construction not found')
      ErrorsFound = .true.
      !CALL ShowSevereError('No match on WINDOW5 data file for Construction=' &
      ! //Trim(ConstructAlphas(0))//', or error in data file.')
      !CALL ShowContinueError('...Looking on file='//TRIM(FullName)) !RS: Secret Search String
      WRITE(DebugFile,*) 'No match on WINDOW5 data file for Construction=' &
        //TRIM(ConstructAlphas(0))//', or error in data file.'
      WRITE(DebugFile,*) '...Looking on file='//TRIM(FullName)
      CYCLE
    END IF

  END DO  ! ...end of Window5 Constructions DO loop

  DEALLOCATE(WConstructNames)

    ! set some (default) properties of the Construction Derived Type
  DO ConstrNum = 1, TotConstructs

    IF (NominalU(ConstrNum) /= 0.0) THEN
      NominalU(ConstrNum)=1.0/NominalU(ConstrNum)
    ELSE
      !CALL ShowSevereError('Nominal U is zero, for construction='//TRIM(Construct(ConstrNum)%Name)) !RS: Secret Search String
      WRITE(DebugFile,*) 'Nominal U is zero, for construction='//TRIM(Construct(ConstrNum)%Name)
      ErrorsFound=.true.
    ENDIF

    CALL CheckAndSetConstructionProperties(ConstrNum,ErrorsFound)

  END DO  ! End of ConstrNum DO loop

  RETURN

END SUBROUTINE GetConstructData

SUBROUTINE GetBuildingData(ErrorsFound)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Linda K. Lawrie
          !       DATE WRITTEN   November 1997
          !       MODIFIED       October 1998, FW; May 1999 FW; Oct 2004 LKL
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! This routine calls other routines to get the Zone, and Surface data
          !  from the input file.

          ! METHODOLOGY EMPLOYED:
          ! The GetObjectItem routines are employed to retrieve the data.

          ! REFERENCES:
          ! na


          ! USE STATEMENTS:
  USE SurfaceGeometry

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  LOGICAL, INTENT(INOUT) :: ErrorsFound ! If errors found in input

          ! SUBROUTINE PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:

    CALL GetZoneData(ErrorsFound)         ! Read Zone data from input file

    CALL SetupZoneGeometry(ErrorsFound)

  RETURN

END SUBROUTINE GetBuildingData

SUBROUTINE GetZoneData(ErrorsFound)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Linda K. Lawrie
          !       DATE WRITTEN   November 1997
          !       MODIFIED       PGE: Added ZONE LIST and ZONE GROUP objects, Nov 2003
          !                      RJH: Added init of DElight member of ZoneDaylight object, Jan 2004
          !                      JG: Added Part of Total Floor Area field March 2006
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! This subroutine gets the zone data for each zone in the input file.

          ! METHODOLOGY EMPLOYED:
          ! The GetObjectItem routines are employed to retrieve the data.

          ! REFERENCES:
          ! IDD Definition for Zone object


          ! USE STATEMENTS:
  USE DataDaylighting, ONLY: ZoneDaylight

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  LOGICAL, INTENT(INOUT) :: ErrorsFound ! If errors found in input

          ! SUBROUTINE PARAMETER DEFINITIONS:
  CHARACTER(len=*), PARAMETER :: BlankString=' '
  CHARACTER(len=*), PARAMETER :: RoutineName='GetZoneData: '
!  INTEGER, PARAMETER :: MaxZonesInList = 100 ! This is to allow DIMENSIONing below

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
!  CHARACTER(len=MaxNameLength), DIMENSION(MaxZonesInList + 1) :: Alphas
!  REAL(r64), DIMENSION(8)              :: Numbers
  INTEGER                         :: NumAlphas, NumNumbers
  INTEGER                         :: IOStatus
  INTEGER                         :: ZoneLoop
  INTEGER                         :: TMP
  INTEGER                         :: Loop
  INTEGER                         :: ListNum
  INTEGER                         :: ZoneNum
  CHARACTER(len=MaxNameLength)    :: ZoneName
  INTEGER                         :: GroupNum
  LOGICAL :: ErrorInName
  LOGICAL :: IsBlank
  
  INTEGER :: DebugFile       =150 !RS: Debugging file denotion, hopfully this works.
    
  OPEN(unit=DebugFile,file='Debug.txt')    !RS: Debugging

  cCurrentModuleObject='Zone'
  NumOfZones=GetNumObjectsFound(TRIM(cCurrentModuleObject))

  ALLOCATE(Zone(NumOfZones))

  ALLOCATE(ZoneDaylight(NumOfZones))

  ZoneLoop=0

  DO Loop=1,NumOfZones

    rNumericArgs=0.0       ! Zero out just in case
    CALL GetObjectItem(TRIM(cCurrentModuleObject),Loop,cAlphaArgs,NumAlphas,rNumericArgs,NumNumbers,IOStatus,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)
    TMP=INDEX(cAlphaArgs(1),CHAR(1))
    DO WHILE (TMP /= 0)
      cAlphaArgs(1)(TMP:TMP)=','
      TMP=INDEX(cAlphaArgs(1),CHAR(1))
    END DO
    TMP=INDEX(cAlphaArgs(1),CHAR(2))
    DO WHILE (TMP /= 0)
      cAlphaArgs(1)(TMP:TMP)='!'
      TMP=INDEX(cAlphaArgs(1),CHAR(2))
    END DO

    !    Make sure Zone Name is unique
    ErrorInName=.false.
    IsBlank=.false.
    CALL VerifyName(cAlphaArgs(1),Zone%Name,ZoneLoop,ErrorInName,IsBlank,TRIM(cCurrentModuleObject)//' Name')
    IF (ErrorInName) THEN
      ErrorsFound=.true.
      CYCLE
    ENDIF

    ZoneLoop=ZoneLoop+1
    Zone(ZoneLoop)%Name=cAlphaArgs(1)
    IF (NumNumbers >=1) &
      Zone(ZoneLoop)%RelNorth=rNumericArgs(1)
    IF (NumNumbers >=2) &
      Zone(ZoneLoop)%OriginX=rNumericArgs(2)
    IF (NumNumbers >=3) &
      Zone(ZoneLoop)%OriginY=rNumericArgs(3)
    IF (NumNumbers >=4) &
      Zone(ZoneLoop)%OriginZ=rNumericArgs(4)
    IF (NumNumbers >=5) &
      Zone(ZoneLoop)%OfType=rNumericArgs(5)
    Zone(ZoneLoop)%OfType=StandardZone
    IF (NumNumbers >=6) &
      Zone(ZoneLoop)%Multiplier=rNumericArgs(6)
    IF (NumNumbers >=7) &
      Zone(ZoneLoop)%CeilingHeight=rNumericArgs(7)
    IF (NumNumbers >=8) &
      Zone(ZoneLoop)%Volume=rNumericArgs(8)
    IF (NumNumbers >=9) &
      Zone(ZoneLoop)%UserEnteredFloorArea=rNumericArgs(9)

    IF (NumAlphas > 1 .and. .not. lAlphaFieldBlanks(2)) THEN
      SELECT CASE (cAlphaArgs(2))

        CASE ('SIMPLE')
          Zone(ZoneLoop)%InsideConvectionAlgo=ASHRAESimple

        CASE ('TARP', 'DETAILED')
         IF (cAlphaArgs(2) == 'DETAILED') THEN
           CALL ShowSevereError(RoutineName//TRIM(cCurrentModuleObject)//'="'//trim(Zone(ZoneLoop)%Name)//'".')
           CALL ShowContinueError('Deprecated value in '//TRIM(cAlphaFieldNames(2))//'="'//    &
              TRIM(cAlphaArgs(2))//'", defaulting to TARP.')
         ENDIF
          Zone(ZoneLoop)%InsideConvectionAlgo=ASHRAETARP

        CASE ('CEILINGDIFFUSER')
          Zone(ZoneLoop)%InsideConvectionAlgo=CeilingDiffuser

        CASE ('TROMBEWALL')
          Zone(ZoneLoop)%InsideConvectionAlgo=TrombeWall

        CASE ('ADAPTIVECONVECTIONALGORITHM ')
          Zone(ZoneLoop)%InsideConvectionAlgo=AdaptiveConvectionAlgorithm

        CASE DEFAULT
          CALL ShowWarningError(RoutineName//TRIM(cCurrentModuleObject)//'="'//trim(Zone(ZoneLoop)%Name)//'".')
          CALL ShowContinueError('Invalid value for '//TRIM(cAlphaFieldNames(2))//'="'//    &
              TRIM(cAlphaArgs(2))//'", defaulting to AdaptiveConvectionAlgorithm.')
          Zone(ZoneLoop)%InsideConvectionAlgo=AdaptiveConvectionAlgorithm

      END SELECT
    ELSE
      ! No zone specific algorithm specified, use default Inside Convection Algorithm
      Zone(ZoneLoop)%InsideConvectionAlgo=DefaultInsideConvectionAlgo

    ENDIF

    IF (NumAlphas > 2 .and. cAlphaArgs(3) /= BlankString) THEN
      SELECT CASE (cAlphaArgs(3))

        CASE ('SIMPLECOMBINED', 'SIMPLE')
          IF (cAlphaArgs(3) == 'SIMPLE') THEN
            CALL ShowSevereError(RoutineName//TRIM(cCurrentModuleObject)//'="'//trim(Zone(ZoneLoop)%Name)//'".')
            CALL ShowContinueError('Deprecated value in '//TRIM(cAlphaFieldNames(3))//'="'//    &
               TRIM(cAlphaArgs(3))//'", defaulting to SimpleCombined.')
          ENDIF
          Zone(ZoneLoop)%OutsideConvectionAlgo=ASHRAESimple

        CASE ('TARP', 'DETAILED', 'BLAST')
          IF (cAlphaArgs(3) == 'DETAILED') THEN
            CALL ShowSevereError(RoutineName//TRIM(cCurrentModuleObject)//'="'//trim(Zone(ZoneLoop)%Name)//'".')
            CALL ShowContinueError('Deprecated value in '//TRIM(cAlphaFieldNames(3))//'="'//    &
               TRIM(cAlphaArgs(3))//'", defaulting to TARP.')
          ENDIF
          IF (cAlphaArgs(3) == 'BLAST') THEN
            CALL ShowSevereError(RoutineName//TRIM(cCurrentModuleObject)//'="'//trim(Zone(ZoneLoop)%Name)//'".')
            CALL ShowContinueError('Deprecated value in '//TRIM(cAlphaFieldNames(3))//'="'//    &
               TRIM(cAlphaArgs(3))//'", defaulting to TARP.')
          ENDIF
          Zone(ZoneLoop)%OutsideConvectionAlgo=ASHRAETARP

        CASE ('MOWITT')
          Zone(ZoneLoop)%OutsideConvectionAlgo=MoWittHcOutside

        CASE ('DOE2','DOE-2')
          Zone(ZoneLoop)%OutsideConvectionAlgo=DOE2HcOutside

        CASE ('ADAPTIVECONVECTIONALGORITHM')
          Zone(ZoneLoop)%OutsideConvectionAlgo=AdaptiveConvectionAlgorithm

        CASE DEFAULT
          CALL ShowWarningError(RoutineName//TRIM(cCurrentModuleObject)//'="'//trim(Zone(ZoneLoop)%Name)//'".')
          CALL ShowContinueError('Invalid value for '//TRIM(cAlphaFieldNames(3))//'="'//    &
              TRIM(cAlphaArgs(3))//'", defaulting to AdaptiveConvectionAlgorithm.')
          Zone(ZoneLoop)%OutsideConvectionAlgo=AdaptiveConvectionAlgorithm

      END SELECT
    ELSE
      ! No zone specific algorithm specified, use default Outside Convection Algorithm
      Zone(ZoneLoop)%OutsideConvectionAlgo=DefaultOutsideConvectionAlgo

    ENDIF

    ! Process the input field:    Part of Total Floor Area
    !   The default value is YES and so only NO needs to be handled
    IF (NumAlphas > 3) THEN
      IF (SameString('No',cAlphaArgs(4))) THEN
        Zone(ZoneLoop)%isPartOfTotalArea = .FALSE.
      END IF
    END IF

    ! Zone outdoor environmental variables, used for zone infiltration/ventilation
    CALL SetupOutputVariable('Zone Outdoor Dry Bulb [C]',Zone(ZoneLoop)%OutDryBulbTemp, &
                               'Zone','Average',Zone(ZoneLoop)%Name)
    CALL SetupOutputVariable('Zone Outdoor Wet Bulb [C]',Zone(ZoneLoop)%OutWetBulbTemp, &
                               'Zone','Average',Zone(ZoneLoop)%Name)
    CALL SetupOutputVariable('Zone Outdoor Wind Speed [m/s]',Zone(ZoneLoop)%WindSpeed, &
                               'Zone','Average',Zone(ZoneLoop)%Name)
  END DO ! Loop

  DO Loop=1,NumOfZones
    ! Check to see if "nominally" controlled -- Zone Name appears in Zone Equip Configuration
    ! relies on zone name being the "name" of the Zone Controlled Equip Configuration
    IF (GetObjectItemNum('ZoneHVAC:EquipmentConnections',Zone(Loop)%Name) > 0) THEN
      Zone(Loop)%IsNominalControlled=.true.
    ELSE
      Zone(Loop)%IsNominalControlled=.false.
    ENDIF
  ENDDO

  ! Get ZONE LIST objects
  cCurrentModuleObject='ZoneList'
  NumOfZoneLists = GetNumObjectsFound(TRIM(cCurrentModuleObject))

  IF (NumOfZoneLists > 0) THEN

    ALLOCATE(ZoneList(NumOfZoneLists))

    DO ListNum = 1, NumOfZoneLists
      CALL GetObjectItem(TRIM(cCurrentModuleObject),ListNum,cAlphaArgs,NumAlphas,rNumericArgs,NumNumbers,IOStatus,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

      ! List name
      ErrorInName = .FALSE.
      IsBlank = .FALSE.
      CALL VerifyName(cAlphaArgs(1),ZoneList%Name,ListNum-1,ErrorInName,IsBlank,TRIM(cCurrentModuleObject)//' Name')

      IF (ErrorInName) THEN
        ErrorsFound = .TRUE.
      END IF

      ZoneList(ListNum)%Name = cAlphaArgs(1)

      IF (FindItemInList(ZoneList(ListNum)%Name,Zone%Name,NumOfZones) > 0) THEN
        CALL ShowWarningError(RoutineName//TRIM(cCurrentModuleObject)//'="'//TRIM(cAlphaArgs(1))//  &
           '":  is a duplicate of a zone name.')
        CALL ShowContinueError('This could be a problem in places where either a Zone Name or a Zone List can be used.')
      ENDIF

      ! List of zones
      ZoneList(ListNum)%NumOfZones = NumAlphas - 1

      IF (ZoneList(ListNum)%NumOfZones < 1) THEN
        CALL ShowSevereError(RoutineName//TRIM(cCurrentModuleObject)//'="'//TRIM(cAlphaArgs(1))//'":  No zones specified.')
        ErrorsFound = .TRUE.
      ELSE
        ALLOCATE(ZoneList(ListNum)%Zone(ZoneList(ListNum)%NumOfZones))
        ZoneList(ListNum)%Zone = 0

        DO ZoneNum = 1, ZoneList(ListNum)%NumOfZones
          ZoneName = cAlphaArgs(ZoneNum + 1)
          ZoneList(ListNum)%Zone(ZoneNum) = FindItemInList(ZoneName,Zone%Name,NumOfZones)
          IF (ZoneList(ListNum)%Zone(ZoneNum) == 0) THEN
            !CALL ShowSevereError(RoutineName//TRIM(cCurrentModuleObject)//'="'//TRIM(cAlphaArgs(1))//'":  '//  &
            !   TRIM(cAlphaFieldNames(ZoneNum+1))//' '//TRIM(ZoneName)//' not found.')   !RS: Secret String Search
            WRITE(DebugFile,*) RoutineName//TRIM(cCurrentModuleObject)//'="'//TRIM(cAlphaArgs(1))//'": '// &
                TRIM(cAlphaFieldNames(ZoneNum+1))//' '//TRIM(ZoneName)//' not found.'
            ErrorsFound = .TRUE.
          END IF

          ! Check for duplicate zones
          DO Loop = 1, ZoneNum - 1
            IF (ZoneList(ListNum)%Zone(ZoneNum) == ZoneList(ListNum)%Zone(Loop)) THEN
              !CALL ShowSevereError(RoutineName//TRIM(cCurrentModuleObject)//'="'//TRIM(cAlphaArgs(1))//'":  '//  &
              !   TRIM(cAlphaFieldNames(ZoneNum+1))//  &
              !   ' '//TRIM(ZoneName)//' appears more than once in list.')   !RS: Secret Search String
              WRITE(DebugFile,*) RoutineName//TRIM(cCurrentModuleObject)//'="'//TRIM(cAlphaArgs(1))//'": '// &
                TRIM(cAlphaFieldNames(ZoneNum+1))//' '//TRIM(ZoneName)//' appears more than once in list.'
              ErrorsFound = .TRUE.
            END IF
          END DO ! Loop
        END DO ! ZoneNum
      END IF

    END DO ! ListNum
  END IF

  ! Get ZONE GROUP objects
  cCurrentModuleObject='ZoneGroup'
  NumOfZoneGroups = GetNumObjectsFound(TRIM(cCurrentModuleObject))

  IF (NumOfZoneGroups > 0) THEN
    ALLOCATE(ZoneGroup(NumOfZoneGroups))

    DO GroupNum = 1, NumOfZoneGroups
      CALL GetObjectItem(TRIM(cCurrentModuleObject),GroupNum,cAlphaArgs,NumAlphas,rNumericArgs,NumNumbers,IOStatus,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

      ! Group name
      ErrorInName = .FALSE.
      IsBlank = .FALSE.
      CALL VerifyName(cAlphaArgs(1),ZoneGroup%Name,GroupNum-1,ErrorInName,IsBlank,TRIM(cCurrentModuleObject)//' Name')

      IF (ErrorInName) THEN
        ErrorsFound = .TRUE.
      END IF

      ZoneGroup(GroupNum)%Name = cAlphaArgs(1)

      ! Multiplier - checked already by IDD rules
      ZoneGroup(GroupNum)%Multiplier = rNumericArgs(1)

      ! Zone list
      ListNum = FindItemInList(cAlphaArgs(2),ZoneList%Name,NumOfZoneLists)
      ZoneGroup(GroupNum)%ZoneList = ListNum

      IF (ListNum == 0) THEN
        !CALL ShowSevereError(RoutineName//TRIM(cCurrentModuleObject)//'="'//TRIM(cAlphaArgs(1))//'":  '//  &
        !   TRIM(cAlphaFieldNames(2))//' named '//TRIM(cAlphaArgs(2))//' not found.')    !RS: Secret Search String
        WRITE(DebugFile,*) RoutineName//TRIM(cCurrentModuleObject)//'="'//TRIM(cAlphaArgs(1))//'": '// &
            TRIM(cAlphaFieldNames(2))//' named '//TRIM(cAlphaArgs(2))//' not found.'
        ErrorsFound = .TRUE.
      ELSE
        ! Check to make sure list is not in use by another ZONE GROUP
        DO Loop = 1, GroupNum - 1
          IF (ZoneGroup(GroupNum)%ZoneList == ZoneGroup(Loop)%ZoneList) THEN
            CALL ShowSevereError(RoutineName//TRIM(cCurrentModuleObject)//'="'//TRIM(cAlphaArgs(1))//'":  '//  &
              TRIM(cAlphaFieldNames(2))//' already used by '//TRIM(cCurrentModuleObject)//' named '// &
              TRIM(ZoneGroup(Loop)%Name)//'.')
            ErrorsFound = .TRUE.
          END IF
        END DO ! Loop

        ! Set group multiplier for each zone in the list
        DO Loop = 1, ZoneList(ListNum)%NumOfZones
          ZoneNum = ZoneList(ListNum)%Zone(Loop)

          IF (ZoneNum > 0 ) THEN
            ! Check to make sure group multiplier was not already set by another ZONE GROUP
            IF (Zone(ZoneNum)%ListGroup == 0) THEN
              Zone(ZoneNum)%ListMultiplier = ZoneGroup(GroupNum)%Multiplier
              Zone(ZoneNum)%ListGroup = ListNum
            ELSE
              CALL ShowSevereError(RoutineName//TRIM(cCurrentModuleObject)//'="'//TRIM(cAlphaArgs(1))//'":  Zone '//  &
                TRIM(Zone(ZoneNum)%Name)// &
                ' in ZoneList already exists in ZoneList of another ZoneGroup.')
              CALL ShowContinueError('Previous ZoneList='//TRIM(ZoneList(Zone(ZoneNum)%ListGroup)%Name))
              ErrorsFound = .TRUE.
            END IF
          END IF
        END DO ! Loop
      END IF

    END DO ! GroupNum
  END IF

  !allocate the array the holds the predefined report data
  ALLOCATE(ZonePreDefRep(NumOfZones))

  RETURN

END SUBROUTINE GetZoneData

! End of Get Input subroutines for the HB Module
!******************************************************************************




! Beginning Initialization Section of the Module
!******************************************************************************

SUBROUTINE InitHeatBalance  ! Heat Balance Initialization Manager

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Rick Strand
          !       DATE WRITTEN   April 1997
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! This subroutine is the main driver for initializations within the
          ! heat balance.

          ! METHODOLOGY EMPLOYED:
          ! Uses the status flags to trigger initialization events.  Some of the files
          !  have been moved to other heat balance managers.  More of these initializations
          !  will have to continue to be re-structured.

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE ConductionTransferFunctionCalc
  USE WindowManager
  USE SolarShading
  USE DaylightingDevices, ONLY: InitDaylightingDevices
!  USE DataRoomAirModel, ONLY: IsZoneDV,IsZoneCV,HVACMassFlow, ZoneDVMixedFlag

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
          ! na

          ! SUBROUTINE PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER :: StormWinNum     ! Number of StormWindow object
  INTEGER :: SurfNum         ! Surface number
  INTEGER :: Num             ! Loop counter
  LOGICAL,SAVE :: ChangeSet=.true.  ! Toggle for checking storm windows


          ! FLOW:

  IF (BeginSimFlag) THEN
    CAll AllocateHeatBalArrays ! Allocate the Module Arrays

    IF (SolutionAlgo == UseCTF .or. SolutionAlgo == UseEMPD) THEN
      CALL DisplayString('Initializing Response Factors')
      CALL InitConductionTransferFunctions ! Initialize the response factors
    END IF

    CALL DisplayString('Initializing Window Optical Properties')
    CALL InitGlassOpticalCalculations ! Initialize the window optical properties
    CALL InitDaylightingDevices ! Initialize any daylighting devices
    CALL DisplayString('Initializing Solar Calculations')
    CALL InitSolarCalculations ! Perform the shadowing calculations

  END IF


  IF (BeginEnvrnFlag) THEN

    MaxHeatLoadPrevDay = 0.0
    MaxCoolLoadPrevDay = 0.0
    MaxTempPrevDay = 0.0
    MinTempPrevDay = 0.0
    MaxHeatLoadZone=-9999.d0
    MaxCoolLoadZone=-9999.d0
    MaxTempZone=-9999.d0
    MinTempZone=1000.d0
    TempZone=-9999.d0
    LoadZone=-9999.d0
    TempZonePrevDay=1000.d0
    LoadZonePrevDay=-9999.d0
    TempZoneSecPrevDay=1000.d0
    TempZoneSecPrevDay=-9999.d0
    WarmupTempDiff=0.0
    WarmupLoadDiff=0.0
    TempZoneRpt=0.0
    LoadZoneRpt=0.0
    MaxLoadZoneRpt=0.0
    CountWarmupDayPoints=0

    DO Num=1,10
      SurfaceWindow%ThetaFace(Num) = 296.15d0
    ENDDO
    SurfaceWindow%EffInsSurfTemp = 23.d0

  END IF

  IF(TotStormWin > 0) THEN
    IF (BeginDayFlag) THEN
      CALL SetStormWindowControl
      ChangeSet=.false.
    ELSEIF (.not. ChangeSet) THEN
      StormWinChangeThisDay = .false.
      DO StormWinNum = 1,TotStormWin
        SurfNum = StormWindow(StormWinNum)%BaseWindowNum
        SurfaceWindow(SurfNum)%StormWinFlagPrevDay = SurfaceWindow(SurfNum)%StormWinFlag
      END DO
      ChangeSet=.true.
    END IF
  END IF

  IF (BeginDayFlag) THEN
    IF (.NOT. WarmupFlag) THEN
      IF (DayOfSim == 1) THEN
        MaxHeatLoadZone=-9999.d0
        MaxCoolLoadZone=-9999.d0
        MaxTempZone=-9999.d0
        MinTempZone=1000.d0
      END IF
    END IF
    CALL PerformSolarCalculations
  END IF

  RETURN

END SUBROUTINE InitHeatBalance

SUBROUTINE AllocateHeatBalArrays  ! Heat Balance Array Allocation

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Richard Liesen
          !       DATE WRITTEN   February 1998
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! This subroutine allocates the arrays to meet simulation requirements

          ! METHODOLOGY EMPLOYED:
          ! Uses the status flags to trigger variable allocation.

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
          ! na

          ! SUBROUTINE PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
          ! na

          ! FLOW:

! Use the total number of zones or surfaces to allocate variables to avoid a limit
    ! Allocate real Variables
    ! Following used for Calculations
    !  Allocate variables in DataHeatBalSys
  ALLOCATE(SumConvHTRadSys(NumOfZones))
  SumConvHTRadSys=0.0D0
  ALLOCATE(SumLatentHTRadSys(NumOfZones))
  SumLatentHTRadSys=0.0D0
  ALLOCATE(QHTRadSysToPerson(NumOfZones))
  QHTRadSysToPerson=0.0D0
  ALLOCATE(QHWBaseboardToPerson(NumOfZones))
  QHWBaseboardToPerson=0.0D0
  ALLOCATE(QSteamBaseboardToPerson(NumOfZones))
  QSteamBaseboardToPerson=0.0D0
  ALLOCATE(QElecBaseboardToPerson(NumOfZones))
  QElecBaseboardToPerson=0.0D0
  ALLOCATE(XMAT(NumOfZones))
  XMAT = 23.0d0
  ALLOCATE(XM2T(NumOfZones))
  XM2T = 23.0d0
  ALLOCATE(XM3T(NumOfZones))
  XM3T = 23.0d0
  ALLOCATE(XM4T(NumOfZones))
  XM4T = 23.0d0
  ALLOCATE(DSXMAT(NumOfZones))
  DSXMAT = 23.0D0
  ALLOCATE(DSXM2T(NumOfZones))
  DSXM2T = 23.0D0
  ALLOCATE(DSXM3T(NumOfZones))
  DSXM3T = 23.0D0
  ALLOCATE(DSXM4T(NumOfZones))
  DSXM4T = 23.0D0
  ALLOCATE(MCPI(NumOfZones))
  MCPI=0.0
  ALLOCATE(MCPTI(NumOfZones))
  MCPTI=0.0
  ALLOCATE(MCPV(NumOfZones))
  MCPV=0.0
  ALLOCATE(MCPTV(NumOfZones))
  MCPTV=0.0
  ALLOCATE(MCPM(NumOfZones))
  MCPM=0.0
  ALLOCATE(MCPTM(NumOfZones))
  MCPTM=0.0
  ALLOCATE(MixingMassFlowZone(NumOfZones))
  MixingMassFlowZone=0.0
  ALLOCATE(MixingMassFlowXHumRat(NumOfZones))
  MixingMassFlowXHumRat=0.0
  ALLOCATE(ZoneLatentGain(NumOfZones))
  ZoneLatentGain=0.0
  ALLOCATE(OAMFL(NumOfZones))
  OAMFL=0.0
  ALLOCATE(VAMFL(NumOfZones))
  VAMFL=0.0
  ALLOCATE(ZTAV(NumOfZones))
  ZTAV = 23.0d0
  ALLOCATE(ZTAVComf(NumOfZones))
  ZTAVComf = 23.0d0
  ALLOCATE(ZT(NumOfZones))
  ZT = 23.0d0
  ALLOCATE(TempTstatAir(NumOfZones))
  TempTstatAir = 23.0d0
  ALLOCATE(MAT(NumOfZones))
  MAT = 23.0d0
  ALLOCATE(ZoneTMX(NumOfZones))
  ZoneTMX = 23.0d0
  ALLOCATE(ZoneTM2(NumOfZones))
  ZoneTM2 = 23.0d0
! Allocate this zone air humidity ratio
  ALLOCATE(ZoneAirHumRatAvg(NumOfZones))
  ZoneAirHumRatAvg=0.01d0
  ALLOCATE(ZoneAirHumRatAvgComf(NumOfZones))
  ZoneAirHumRatAvgComf=0.01d0
  ALLOCATE(ZoneAirHumRat(NumOfZones))
  ZoneAirHumRat=0.01d0
  ALLOCATE(ZoneAirHumRatOld(NumOfZones))
  ZoneAirHumRatOld=0.01d0
  ALLOCATE(SumHmAW(NumOfZones))
  SumHmAW=0.0d0
  ALLOCATE(SumHmARa(NumOfZones))
  SumHmARa=0.0d0
  ALLOCATE(SumHmARaW(NumOfZones))
  SumHmARaW=0.0d0
  ALLOCATE(MCPTE(NumOfZones))
  MCPTE=0.0
  ALLOCATE(MCPE(NumOfZones))
  MCPE=0.0
  ALLOCATE(EAMFL(NumOfZones))
  EAMFL=0.0
  ALLOCATE(MCPTC(NumOfZones))
  MCPTC=0.0
  ALLOCATE(MCPC(NumOfZones))
  MCPC=0.0
  ALLOCATE(CTMFL(NumOfZones))
  CTMFL=0.0
  ALLOCATE(MDotCPOA(NumOfZones))
  MDotCPOA=0.0
  ALLOCATE(MDotOA(NumOfZones))
  MDotOA=0.0
  IF (Contaminant%CO2Simulation) Then
    OutdoorCO2 = GetCurrentScheduleValue(Contaminant%CO2OutdoorSchedPtr)
    ALLOCATE(ZoneAirCO2(NumOfZones))
    ZoneAirCO2=OutdoorCO2
    ALLOCATE(ZoneAirCO2Temp(NumOfZones))
    ZoneAirCO2Temp=OutdoorCO2
    ALLOCATE(ZoneAirCO2Avg(NumOfZones))
    ZoneAirCO2Avg=OutdoorCO2
  END IF
  IF (Contaminant%GenericContamSimulation) Then
    OutdoorGC = GetCurrentScheduleValue(Contaminant%GenericContamOutdoorSchedPtr)
    ALLOCATE(ZoneAirGC(NumOfZones))
    ZoneAirGC=OutdoorGC
    ALLOCATE(ZoneAirGCTemp(NumOfZones))
    ZoneAirGCTemp=OutdoorGC
    ALLOCATE(ZoneAirGCAvg(NumOfZones))
    ZoneAirGCAvg=OutdoorGC
  END IF
  ALLOCATE(MaxTempPrevDay(NumofZones))
           MaxTempPrevDay = 0.0
  ALLOCATE(MinTempPrevDay(NumofZones))
           MinTempPrevDay = 0.0
  ALLOCATE(MaxHeatLoadPrevDay(NumofZones))
           MaxHeatLoadPrevDay = 0.0
  ALLOCATE(MaxCoolLoadPrevDay(NumofZones))
           MaxCoolLoadPrevDay = 0.0
  ALLOCATE(MaxHeatLoadZone(NumofZones))
           MaxHeatLoadZone = -9999.d0
  ALLOCATE(MaxCoolLoadZone(NumofZones))
           MaxCoolLoadZone = -9999.d0
  ALLOCATE(MaxTempZone(NumofZones))
           MaxTempZone = -9999.d0
  ALLOCATE(MinTempZone(NumofZones))
           MinTempZone = 1000.d0
  ALLOCATE(TempZonePrevDay(NumofZones))
           TempZonePrevDay = 0.0
  ALLOCATE(LoadZonePrevDay(NumofZones))
           LoadZonePrevDay = 0.0
  ALLOCATE(TempZoneSecPrevDay(NumofZones))
           TempZoneSecPrevDay = 0.0
  ALLOCATE(LoadZoneSecPrevDay(NumofZones))
           LoadZoneSecPrevDay = 0.0
  ALLOCATE(WarmupTempDiff(NumofZones))
           WarmupTempDiff = 0.0
  ALLOCATE(WarmupLoadDiff(NumofZones))
           WarmupLoadDiff = 0.0
  ALLOCATE(TempZone(NumofZones))
           TempZone = 0.0
  ALLOCATE(LoadZone(NumofZones))
           LoadZone = 0.0
  ALLOCATE(TempZoneRpt(NumOfTimeStepInHour*24,NumofZones))
           TempZoneRpt=0.0
  ALLOCATE(LoadZoneRpt(NumOfTimeStepInHour*24,NumofZones))
           LoadZoneRpt=0.0
  ALLOCATE(MaxLoadZoneRpt(NumOfTimeStepInHour*24,NumofZones))
           MaxLoadZoneRpt=0.0
  ALLOCATE(WarmupConvergenceValues(NumOfZones))
  ALLOCATE(TempZoneRptStdDev(NumOfTimeStepInHour*24))
  ALLOCATE(LoadZoneRptStdDev(NumOfTimeStepInHour*24))

  CountWarmupDayPoints=0

RETURN

END SUBROUTINE AllocateHeatBalArrays


! End Initialization Section of the Module
!******************************************************************************


! Beginning of Record Keeping subroutines for the HB Module
! *****************************************************************************

SUBROUTINE RecKeepHeatBalance   ! Heat Balance Record Keeping Manager

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Rick Strand
          !       DATE WRITTEN   April 1997
          !       MODIFIED       June 2011, Daeho Kang for individual zone maximums & convergence outputs
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! This subroutine is the main driver for record keeping within the
          ! heat balance.

          ! METHODOLOGY EMPLOYED:
          ! Uses the status flags to trigger record keeping events.

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE General, ONLY: RoundSigDigits
  USE DataSystemVariables, ONLY: ReportDetailedWarmupConvergence

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
          ! na

          ! SUBROUTINE PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
!  CHARACTER(len=MaxNameLength) :: ZoneName
  INTEGER :: ZoneNum
  LOGICAL,SAVE :: FirstWarmupWrite=.true.

          ! FLOW:

  ! Always do the following record keeping (every time step):
  ! Record Maxs & Mins for individual zone
  DO ZoneNum = 1, NumOfZones
    IF (ZTAV(ZoneNum) > MaxTempZone(ZoneNum)) THEN
      MaxTempZone(ZoneNum) = ZTAV(ZoneNum)
    END IF
    IF (ZTAV(ZoneNum) < MinTempZone(ZoneNum)) THEN
      MinTempZone(ZoneNum) = ZTAV(ZoneNum)
    END IF
    IF (SNLoadHeatRate(ZoneNum) > MaxHeatLoadZone(ZoneNum)) THEN
      MaxHeatLoadZone(ZoneNum) = SNLoadHeatRate(ZoneNum)
    END IF
    IF (SNLoadCoolRate(ZoneNum) > MaxCoolLoadZone(ZoneNum)) THEN
      MaxCoolLoadZone(ZoneNum) = SNLoadCoolRate(ZoneNum)
    END IF

        ! Record temperature and load for individual zone
      TempZoneSecPrevDay(ZoneNum)=TempZonePrevDay(ZoneNum)
      LoadZoneSecPrevDay(ZoneNum)=LoadZonePrevDay(ZoneNum)
      TempZonePrevDay(ZoneNum)=TempZone(ZoneNum)
      LoadZonePrevDay(ZoneNum)=LoadZone(ZoneNum)
      TempZone(ZoneNum)=ZTAV(ZoneNum)
      LoadZone(ZoneNum)=Max(SNLoadHeatRate(ZoneNum),ABS(SNLoadCoolRate(ZoneNum)))

        ! Calculate differences in temperature and load for the last two warmup days
    IF (.NOT. WarmupFlag .AND. DayOfSim == 1 .AND. .NOT. DoingSizing) THEN
      WarmupTempDiff(ZoneNum)=ABS(TempZoneSecPrevDay(ZoneNum)-TempZonePrevDay(ZoneNum))
      WarmupLoadDiff(ZoneNum)=ABS(LoadZoneSecPrevDay(ZoneNum)-LoadZonePrevDay(ZoneNum))
      IF (ZoneNum == 1) CountWarmupDayPoints=CountWarmupDayPoints+1
      TempZoneRpt(CountWarmupDayPoints,ZoneNum)=WarmupTempDiff(ZoneNum)
      LoadZoneRpt(CountWarmupDayPoints,ZoneNum)=WarmupLoadDiff(ZoneNum)
      MaxLoadZoneRpt(CountWarmupDayPoints,ZoneNum)=LoadZone(ZoneNum)

      IF (ReportDetailedWarmupConvergence) THEN  ! only do this detailed thing when requested by user is on
          ! Write Warmup Convergence Information to the initialization output file
        IF (FirstWarmupWrite) THEN
          Write(OutputFileInits,732)
          FirstWarmupWrite=.false.
        ENDIF

        Write(OutputFileInits,731) TRIM(Zone(ZoneNum)%Name),trim(RoundSigDigits(TimeStep)),trim(RoundSigDigits(HourofDay)),  &
             trim(RoundSigDigits(WarmupTempDiff(ZoneNum),10)),trim(RoundSigDigits(WarmupLoadDiff(ZoneNum),10))

      ENDIF
    ENDIF

  END DO

  731 Format(' Warmup Convergence Information, ',A,',',A,',',A,',',A,',',A)
  732 Format('! <Warmup Convergence Information>,Zone Name,Time Step,Hour of Day,Warmup Temperature Difference {deltaC},'  &
         'Warmup Load Difference {W}')

  ! There is no hourly record keeping in the heat balance.

  ! There is no daily record keeping in the heat balance.

  ! There is no environment level record keeping in the heat balance.

  ! There is no simulation level record keeping in the heat balance.

  RETURN

END SUBROUTINE RecKeepHeatBalance


SUBROUTINE CheckWarmupConvergence

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Rick Strand
          !       DATE WRITTEN   April 1997
          !       MODIFIED       June 2011, Daeho Kang for individual zone comparison
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! This subroutine checks warmup convergence values.

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE General, ONLY: RoundSigDigits


  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
          ! na

          ! SUBROUTINE PARAMETER DEFINITIONS:
  REAL(r64), PARAMETER :: MinLoad = 100.d0     ! Minimum laods for convergence check
                                               ! To avoid big percentage difference in low load situations

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER :: ZoneNum
  LOGICAL,SAVE :: WarmupConvergenceWarning=.false.
  LOGICAL,SAVE :: SizingWarmupConvergenceWarning=.false.
  LOGICAL :: ConvergenceChecksFailed

          ! Convergence criteria for warmup days:
          ! Perform another warmup day unless both the % change in loads and
          ! absolute change in zone temp min & max are less than their criteria.

  ConvergenceChecksFailed=.false.

  IF (NumofZones <= 0) THEN ! if there are no zones, immediate convergence
      WarmupFlag=.FALSE.
  ELSE
    DO ZoneNum=1,NumofZones

      WarmupConvergenceValues(ZoneNum)%TestMaxTempValue=ABS(MaxTempPrevDay(ZoneNum)-MaxTempZone(ZoneNum))
      WarmupConvergenceValues(ZoneNum)%TestMinTempValue=ABS(MinTempPrevDay(ZoneNum)-MinTempZone(ZoneNum))
      IF (WarmupConvergenceValues(ZoneNum)%TestMaxTempValue <= TempConvergTol) THEN
        WarmupConvergenceValues(ZoneNum)%PassFlag(1)=2
      ELSE
        ConvergenceChecksFailed=.true.
        WarmupConvergenceValues(ZoneNum)%PassFlag(1)=1
      ENDIF

      IF (WarmupConvergenceValues(ZoneNum)%TestMinTempValue <= TempConvergTol) THEN
        WarmupConvergenceValues(ZoneNum)%PassFlag(2)=2
      ELSE
        ConvergenceChecksFailed=.true.
        WarmupConvergenceValues(ZoneNum)%PassFlag(2)=1
      ENDIF

      IF (MaxHeatLoadZone(ZoneNum) > 1.0d-4) THEN  ! make sure load big enough to divide
        MaxHeatLoadZone(ZoneNum) = ABS(Max(MaxHeatLoadZone(ZoneNum),MinLoad))
        WarmupConvergenceValues(ZoneNum)%TestMaxHeatLoadValue=  &
           ABS((MaxHeatLoadZone(ZoneNum)-MaxHeatLoadPrevDay(ZoneNum))/MaxHeatLoadZone(ZoneNum))
        IF (WarmupConvergenceValues(ZoneNum)%TestMaxHeatLoadValue <= LoadsConvergTol) THEN
          WarmupConvergenceValues(ZoneNum)%PassFlag(3)=2
        ELSE
          ConvergenceChecksFailed=.true.
          WarmupConvergenceValues(ZoneNum)%PassFlag(3)=1
        END IF
      ELSE
        WarmupConvergenceValues(ZoneNum)%PassFlag(3)=2
      END IF

      IF (MaxCoolLoadZone(ZoneNum) > 1.0d-4) THEN
        MaxCoolLoadZone(ZoneNum) = ABS(Max(MaxCoolLoadZone(ZoneNum),MinLoad))
        WarmupConvergenceValues(ZoneNum)%TestMaxCoolLoadValue=  &
           ABS((MaxCoolLoadZone(ZoneNum)-MaxCoolLoadPrevDay(ZoneNum))/MaxCoolLoadZone(ZoneNum))
        IF (WarmupConvergenceValues(ZoneNum)%TestMaxCoolLoadValue <= LoadsConvergTol) THEN
          WarmupConvergenceValues(ZoneNum)%PassFlag(4)=2
        ELSE
          ConvergenceChecksFailed=.true.
          WarmupConvergenceValues(ZoneNum)%PassFlag(4)=1
        END IF
      ELSE
        WarmupConvergenceValues(ZoneNum)%PassFlag(4)=2
      END IF

      IF (DayOfSim >= MaxNumberOfWarmupDays .and. WarmupFlag) THEN
            ! Check convergence for individual zone
        IF (SUM(WarmupConvergenceValues(ZoneNum)%PassFlag) /= 8) THEN ! pass=2 * 4 values for convergence
          CALL ShowSevereError('CheckWarmupConvergence: Loads Initialization, Zone="'//TRIM(Zone(ZoneNum)%Name)//  &
            '" did not converge after '//TRIM(RoundSigDigits(MaxNumberOfWarmupDays))//' warmup days.')
          IF (.not. WarmupConvergenceWarning .and. .not. DoingSizing) THEN
            CALL ShowContinueError('See Warmup Convergence Information in .eio file for details.')
            WarmupConvergenceWarning=.true.
          ELSEIF (.not. SizingWarmupConvergenceWarning .and. DoingSizing) THEN
            CALL ShowContinueError('Warmup Convergence failing during sizing.')
            SizingWarmupConvergenceWarning=.true.
          ENDIF
          IF (RunPeriodEnvironment) THEN
            CALL ShowContinueError('...Environment(RunPeriod)="'//TRIM(EnvironmentName)//'"')
          ELSE
            CALL ShowContinueError('...Environment(SizingPeriod)="'//TRIM(EnvironmentName)//'"')
          ENDIF

          CALL ShowContinueError('..Max Temp Comparison = '//  &
                       TRIM(RoundSigDigits(WarmupConvergenceValues(ZoneNum)%TestMaxTempValue,2))//  &
                       ' vs Temperature Convergence Tolerance='//TRIM(RoundSigDigits(TempConvergTol,2))//  &
                       ' - '//PassFail(WarmupConvergenceValues(ZoneNum)%PassFlag(1))//' Convergence')
          CALL ShowContinueError('..Min Temp Comparison = '//  &
                       TRIM(RoundSigDigits(WarmupConvergenceValues(ZoneNum)%TestMinTempValue,2))//  &
                       ' vs Temperature Convergence Tolerance='//TRIM(RoundSigDigits(TempConvergTol,2))//  &
                       ' - '//PassFail(WarmupConvergenceValues(ZoneNum)%PassFlag(2))//' Convergence')
          CALL ShowContinueError('..Max Heat Load Comparison = '//    &
                TRIM(RoundSigDigits(WarmupConvergenceValues(ZoneNum)%TestMaxHeatLoadValue,4))//  &
                ' vs Loads Convergence Tolerance='//TRIM(RoundSigDigits(LoadsConvergTol,2))//  &
                ' - '//PassFail(WarmupConvergenceValues(ZoneNum)%PassFlag(3))//' Convergence')
            CALL ShowContinueError('..Max Cool Load Comparison = '// &
                TRIM(RoundSigDigits(WarmupConvergenceValues(ZoneNum)%TestMaxCoolLoadValue,4))//  &
               ' vs Loads Convergence Tolerance='//TRIM(RoundSigDigits(LoadsConvergTol,2))//  &
               ' - '//PassFail(WarmupConvergenceValues(ZoneNum)%PassFlag(4))//' Convergence')
        END IF

      END IF

          ! Transfer current daily max and min loads and temperatures to the
          ! variables containing the last day's values
      MaxHeatLoadPrevDay(ZoneNum)=MaxHeatLoadZone(ZoneNum)
      MaxCoolLoadPrevDay(ZoneNum)=MaxCoolLoadZone(ZoneNum)
      MaxTempPrevDay(ZoneNum)=MaxTempZone(ZoneNum)
      MinTempPrevDay(ZoneNum)=MinTempZone(ZoneNum)

      MaxHeatLoadZone(ZoneNum)=-9999.d0
      MaxCoolLoadZone(ZoneNum)=-9999.d0
      MaxTempZone(ZoneNum)=-9999.d0
      MinTempZone(ZoneNum)=1000.d0

    END DO

          ! Limit the number of warmup days, regardless of the number of zones
          ! in the building, to some arbitrary value based on common sense and
          ! experience with the (I)BLAST program.  If too many warmup days were
          ! required, notify the program user.

    IF ((DayOfSim >= MaxNumberOfWarmupDays) .and. WarmupFlag .and. ConvergenceChecksFailed) THEN
      IF (MaxNumberOfWarmupDays < DefaultMaxNumberOfWarmupDays) THEN
        CALL ShowSevereError('CheckWarmupConvergence: User supplied maximum warmup days='//  &
          TRIM(RoundSigDigits(MaxNumberOfWarmupDays))//' is insufficient.')
        CALL ShowContinueError('Suggest setting maximum number of warmup days to at least '//  &
           trim(RoundSigDigits(DefaultMaxNumberOfWarmupDays))//'.')
      ENDIF
    ENDIF

          ! Set warmup flag to true depending on value of ConvergenceChecksFailed (true=fail)
          ! and minimum number of warmup days
    IF (.not. ConvergenceChecksFailed .and. DayOfSim >= MinNumberOfWarmupDays) THEN
      WarmupFlag=.false.
    ELSEIF (.not. ConvergenceChecksFailed .and. DayOfSim < MinNumberOfWarmupDays) THEN
      WarmupFlag=.true.
    END IF

          ! If max warmup days reached and still warmupflag, then go to non-warmup state.
          ! prior messages will have been displayed
    IF ((DayOfSim >= MaxNumberOfWarmupDays) .and. WarmupFlag) THEN
      WarmupFlag=.false.
    ENDIF

  END IF

  RETURN

END SUBROUTINE CheckWarmupConvergence

SUBROUTINE ReportWarmupConvergence

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Linda Lawrie
          !       DATE WRITTEN   October 2011
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! na

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE General, ONLY: RoundSigDigits

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
          ! na

          ! SUBROUTINE PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER :: ZoneNum
  LOGICAL,SAVE :: FirstWarmupWrite=.true.
  REAL(r64) :: AverageZoneTemp
  REAL(r64) :: AverageZoneLoad
  REAL(r64) :: StdDevZoneTemp
  REAL(r64) :: StdDevZoneLoad
  CHARACTER(len=15) :: EnvHeader
  INTEGER :: Num  ! loop control


    IF (.not. WarmupFlag) THEN   ! Report out average/std dev
        ! Write Warmup Convervence Information to the initialization output file
      IF (FirstWarmupWrite .and. NumOfZones > 0) THEN
        Write(OutputFileInits,730)
        FirstWarmupWrite=.false.
      ENDIF

      TempZoneRptStdDev=0.0d0
      LoadZoneRptStdDev=0.0d0

      IF (RunPeriodEnvironment) THEN
        EnvHeader='RunPeriod:'
      ELSE
        EnvHeader='SizingPeriod:'
      ENDIF

      DO ZoneNum=1,NumofZones
        AverageZoneTemp=SUM(TempZoneRpt(1:CountWarmupDayPoints,ZoneNum))/REAL(CountWarmupDayPoints,r64)
        DO Num=1,CountWarmupDayPoints
          IF (MaxLoadZoneRpt(Num,ZoneNum) > 1.d-4) THEN
            LoadZoneRpt(Num,ZoneNum)=LoadZoneRpt(Num,ZoneNum)/MaxLoadZoneRpt(Num,ZoneNum)
          ELSE
            LoadZoneRpt(Num,ZoneNum)=0.0
          ENDIF
        ENDDO
        AverageZoneLoad=SUM(LoadZoneRpt(1:CountWarmupDayPoints,ZoneNum))/REAL(CountWarmupDayPoints,r64)
        StdDevZoneTemp=0.0d0
        StdDevZoneLoad=0.0d0
        DO Num=1,CountWarmupDayPoints
          TempZoneRptStdDev(Num)=(TempZoneRpt(Num,ZoneNum)-AverageZoneTemp)**2
          LoadZoneRptStdDev(Num)=(LoadZoneRpt(Num,ZoneNum)-AverageZoneLoad)**2
        ENDDO
        StdDevZoneTemp=SQRT(SUM(TempZoneRptStdDev(1:CountWarmupDayPoints))/REAL(CountWarmupDayPoints,r64))
        StdDevZoneLoad=SQRT(SUM(LoadZoneRptStdDev(1:CountWarmupDayPoints))/REAL(CountWarmupDayPoints,r64))

        Write(OutputFileInits,731) TRIM(Zone(ZoneNum)%Name),  &
                                   trim(EnvHeader)//' '//trim(EnvironmentName),  &
                                   trim(RoundSigDigits(AverageZoneTemp,10)),  &
                                   trim(RoundSigDigits(StdDevZoneTemp,10)),   &
                                   trim(PassFail(WarmupConvergenceValues(ZoneNum)%PassFlag(1))),  &
                                   trim(PassFail(WarmupConvergenceValues(ZoneNum)%PassFlag(2))),  &
                                   trim(RoundSigDigits(AverageZoneLoad,10)),  &
                                   trim(RoundSigDigits(StdDevZoneLoad,10)),  &
                                   trim(PassFail(WarmupConvergenceValues(ZoneNum)%PassFlag(3))),  &
                                   trim(PassFail(WarmupConvergenceValues(ZoneNum)%PassFlag(4)))
      ENDDO

    END IF


  RETURN


  730 Format('! <Warmup Convergence Information>,Zone Name,Environment Type/Name,'  &
         'Average Warmup Temperature Difference {deltaC},'  &
         'Std Dev Warmup Temperature Difference {deltaC},Max Temperature Pass/Fail Convergence,'  &
         'Min Temperature Pass/Fail Convergence,Average Warmup Load Difference {W},Std Dev Warmup Load Difference {W},',  &
         'Heating Load Pass/Fail Convergence,Cooling Load Pass/Fail Convergence')
  731 Format(' Warmup Convergence Information',10(',',A))

END SUBROUTINE ReportWarmupConvergence

!        End of Record Keeping subroutines for the HB Module
! *****************************************************************************


! Beginning of Reporting subroutines for the HB Module
! *****************************************************************************
SUBROUTINE ReportHeatBalance  ! Heat Balance Reporting Manager

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Rick Strand
          !       DATE WRITTEN   July 1997
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! This subroutine is the main driver for reporting within the heat
          ! balance.

          ! METHODOLOGY EMPLOYED:
          ! Uses the status flags to trigger record keeping events.

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE OutputReportTabular, ONLY: UpdateTabularReports
  USE ScheduleManager, ONLY: ReportScheduleValues
  USE NodeInputManager, ONLY: CalcMoreNodeInfo
  USE EconomicTariff, ONLY: UpdateUtilityBills        !added for computing annual utility costs
  USE DataSystemVariables, ONLY: ReportDuringWarmup, UpdateDataDuringWarmupExternalInterface ! added for FMI
  USE DataReportingFlags

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
          ! na

          ! SUBROUTINE PARAMETER DEFINITIONS:
  CHARACTER(len=*), PARAMETER :: EndOfHeaderFormat = "('End of Data Dictionary')"    ! End of data dictionary marker
  CHARACTER(len=*), PARAMETER :: EnvironmentStampFormat = "(a,',',a,3(',',f7.2),',',f7.2)" ! Format descriptor for environ stamp

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
!  LOGICAL, SAVE :: PrintEnvrnStamp=.false.

          ! FLOW:

          ! Time step level reporting:

  CALL ReportScheduleValues
  IF (.NOT. WarmupFlag .AND. DoOutputReporting) THEN
    CALL CalcMoreNodeInfo
    CALL UpdateDataandReport(ZoneTSReporting)
    CALL UpdateTabularReports(ZoneTSReporting)
    CALL UpdateUtilityBills
  ELSEIF (.not. KickOffSimulation .and. DoOutputReporting .and. ReportDuringWarmup) THEN
    IF (BeginDayFlag .and. .not. PrintEnvrnStampWarmupPrinted) THEN
      PrintEnvrnStampWarmup=.true.
      PrintEnvrnStampWarmupPrinted=.true.
    ENDIF
    IF (.not. BeginDayFlag) PrintEnvrnStampWarmupPrinted=.false.
    IF (PrintEnvrnStampWarmup) THEN
      IF (PrintEndDataDictionary .AND. DoOutputReporting) THEN
        WRITE (OutputFileStandard,EndOfHeaderFormat)
        WRITE (OutputFileMeters,EndOfHeaderFormat)
        PrintEndDataDictionary = .FALSE.
      ENDIF
      IF (DoOutputReporting) THEN
        WRITE (OutputFileStandard,EnvironmentStampFormat) '1', &
                  'Warmup {'//trim(cWarmupDay)//'} '//Trim(EnvironmentName),Latitude,Longitude,TimeZoneNumber,Elevation
        WRITE (OutputFileMeters,EnvironmentStampFormat) '1', &
                  'Warmup {'//trim(cWarmupDay)//'} '//Trim(EnvironmentName),Latitude,Longitude,TimeZoneNumber,Elevation
        PrintEnvrnStampWarmup=.FALSE.
      END IF
    END IF
    CALL CalcMoreNodeInfo
    CALL UpdateDataandReport(ZoneTSReporting)
  ELSEIF (UpdateDataDuringWarmupExternalInterface) THEN ! added for FMI
      CALL UpdateDataandReport(ZoneTSReporting)
  END IF
          ! There is no hourly reporting in the heat balance.

          ! There is no daily reporting in the heat balance.

          ! There is no simulation level record keeping in the heat balance.

  RETURN

END SUBROUTINE ReportHeatBalance

!        End of Reporting subroutines for the HB Module

SUBROUTINE GetFrameAndDividerData(ErrorsFound)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Fred Winkelmann
          !       DATE WRITTEN   May 2000
          !       MODIFIED       April 2002 (FCW): get window reveal data
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! Gets input data for window frame and/or divider and/or window
          ! inside/outside reveal.

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor, ONLY: GetNumObjectsFound, GetObjectItem, FindItemInList, VerifyName

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  LOGICAL, INTENT(INOUT) :: ErrorsFound  ! set to true if errors found in input

          ! SUBROUTINE PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:

  INTEGER :: IOStat               ! IO Status when calling get input subroutine
  CHARACTER(len=MaxNameLength),DIMENSION(2) &
          :: FrameDividerNames    ! Frame/Divider Alpha names
  INTEGER :: FrameDividerNum      ! Counter to keep track of the frame/divider number
  INTEGER :: FrameDividerNumAlpha ! Number of frame/divider alpha names being passed
  INTEGER :: FrameDividerNumProp  ! Number of frame/divider properties being passed
  REAL(r64), DIMENSION(23) :: FrameDividerProps !Temporary array to transfer frame/divider properties
  INTEGER :: Loop
  LOGICAL :: ErrorInName
  LOGICAL :: IsBlank

  CurrentModuleObject='WindowProperty:FrameAndDivider'
  TotFrameDivider=GetNumObjectsFound(TRIM(CurrentModuleObject))
  ALLOCATE (FrameDivider(TotFrameDivider))
  IF(TotFrameDivider == 0) RETURN

  FrameDividerNum=0

  DO Loop=1,TotFrameDivider

    !Call Input Get routine to retrieve frame/divider data
    CALL GetObjectItem(TRIM(CurrentModuleObject),Loop,FrameDividerNames,FrameDividerNumAlpha, &
                       FrameDividerProps,FrameDividerNumProp,IOSTAT,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

    ErrorInName=.false.
    IsBlank=.false.
    CALL VerifyName(FrameDividerNames(1),FrameDivider%Name,FrameDividerNum,ErrorInName,IsBlank, &
           TRIM(CurrentModuleObject)//' Name')
    IF (ErrorInName) THEN
      ErrorsFound=.true.
      CYCLE
    ENDIF

    !Load the frame/divider derived type from the input data.
    FrameDividerNum=FrameDividerNum+1
    FrameDivider(FrameDividerNum)%Name = FrameDividerNames(1)
    FrameDivider(FrameDividerNum)%FrameWidth = FrameDividerProps(1)
    FrameDivider(FrameDividerNum)%FrameProjectionOut = FrameDividerProps(2)
    FrameDivider(FrameDividerNum)%FrameProjectionIn = FrameDividerProps(3)
    IF(FrameDivider(FrameDividerNum)%FrameWidth == 0.0) THEN
      FrameDivider(FrameDividerNum)%FrameProjectionOut = 0.0
      FrameDivider(FrameDividerNum)%FrameProjectionIn = 0.0
    END IF
    FrameDivider(FrameDividerNum)%FrameConductance = FrameDividerProps(4)
    FrameDivider(FrameDividerNum)%FrEdgeToCenterGlCondRatio = FrameDividerProps(5)
    FrameDivider(FrameDividerNum)%FrameSolAbsorp = FrameDividerProps(6)
    FrameDivider(FrameDividerNum)%FrameVisAbsorp = FrameDividerProps(7)
    FrameDivider(FrameDividerNum)%FrameEmis = FrameDividerProps(8)
    IF (SameString(FrameDividerNames(2),'DividedLite')) THEN
      FrameDivider(FrameDividerNum)%DividerType = DividedLite
    ELSEIF (SameString(FrameDividerNames(2),'Suspended')) THEN
      FrameDivider(FrameDividerNum)%DividerType = Suspended
    ELSE
      CALL ShowWarningError(TRIM(CurrentModuleObject)//'="'//trim(FrameDividerNames(1))//  &
         '", Invalid '//TRIM(cAlphaFieldNames(2)))
      CALL ShowContinueError('Entered="'//trim(FrameDividerNames(2))//  &
         '", must be DividedLite or Suspended.  Will be set to DividedLite.')
      FrameDivider(FrameDividerNum)%DividerType = DividedLite
    ENDIF
    FrameDivider(FrameDividerNum)%DividerWidth = FrameDividerProps(9)
    FrameDivider(FrameDividerNum)%HorDividers = FrameDividerProps(10)
    FrameDivider(FrameDividerNum)%VertDividers = FrameDividerProps(11)
    FrameDivider(FrameDividerNum)%DividerProjectionOut = FrameDividerProps(12)
    FrameDivider(FrameDividerNum)%DividerProjectionIn = FrameDividerProps(13)
    IF(FrameDivider(FrameDividerNum)%DividerWidth == 0.0 .OR. &
       FrameDivider(FrameDividerNum)%DividerType == Suspended) THEN
      FrameDivider(FrameDividerNum)%DividerProjectionOut = 0.0
      FrameDivider(FrameDividerNum)%DividerProjectionIn = 0.0
    END IF
    FrameDivider(FrameDividerNum)%DividerConductance = FrameDividerProps(14)
    FrameDivider(FrameDividerNum)%DivEdgeToCenterGlCondRatio = FrameDividerProps(15)
    FrameDivider(FrameDividerNum)%DividerSolAbsorp = FrameDividerProps(16)
    FrameDivider(FrameDividerNum)%DividerVisAbsorp = FrameDividerProps(17)
    FrameDivider(FrameDividerNum)%DividerEmis = FrameDividerProps(18)
    FrameDivider(FrameDividerNum)%OutsideRevealSolAbs = FrameDividerProps(19)
    FrameDivider(FrameDividerNum)%InsideSillDepth = FrameDividerProps(20)
    FrameDivider(FrameDividerNum)%InsideSillSolAbs = FrameDividerProps(21)
    FrameDivider(FrameDividerNum)%InsideReveal = FrameDividerProps(22)
    FrameDivider(FrameDividerNum)%InsideRevealSolAbs = FrameDividerProps(23)

    IF (FrameDivider(FrameDividerNum)%DividerWidth > 0.0 .and.   &
        (FrameDivider(FrameDividerNum)%HorDividers == 0 .and. FrameDivider(FrameDividerNum)%VertDividers == 0)) THEN
      CALL ShowWarningError(TRIM(CurrentModuleObject)//': In FrameAndDivider '//TRIM(FrameDivider(FrameDividerNum)%Name)// &
        ' '//TRIM(cNumericFieldNames(9))//' > 0 ')
      CALL ShowContinueError('...but '//TRIM(cNumericFieldNames(10))//' = 0 and '//  &
        TRIM(cNumericFieldNames(11))//' = 0.')
      CALL ShowContinueError('...'//TRIM(cNumericFieldNames(9))//' set to 0.')
      FrameDivider(FrameDividerNum)%DividerWidth=0.0
    ENDIF
    ! Prevent InsideSillDepth < InsideReveal
    IF(FrameDivider(FrameDividerNum)%InsideSillDepth < FrameDivider(FrameDividerNum)%InsideReveal) THEN
      CALL ShowWarningError(TRIM(CurrentModuleObject)//': In FrameAndDivider '//TRIM(FrameDivider(FrameDividerNum)%Name)// &
        ' '//TRIM(cNumericFieldNames(20))//' is less than '//TRIM(cNumericFieldNames(22))//'; it will be set to '// &
        TRIM(cNumericFieldNames(22))//'.')
      FrameDivider(FrameDividerNum)%InsideSillDepth = FrameDivider(FrameDividerNum)%InsideReveal
    END IF

!    ! Warn if InsideSillDepth OR InsideReveal > 0.2meters to warn of inaccuracies
!    IF(FrameDivider(FrameDividerNum)%InsideSillDepth > 0.2d0) THEN
!      CALL ShowWarningError(TRIM(CurrentModuleObject)//': In FrameAndDivider '//TRIM(FrameDivider(FrameDividerNum)%Name)// &
!        ' '//TRIM(cNumericFieldNames(20))//' is greater than 0.2 meters, which could cause inaccuracies in zone cooling energy.')
!    END IF
!    IF(FrameDivider(FrameDividerNum)%InsideReveal > 0.2d0) THEN
!      CALL ShowWarningError(TRIM(CurrentModuleObject)//': In FrameAndDivider '//TRIM(FrameDivider(FrameDividerNum)%Name)// &
!        ' '//TRIM(cNumericFieldNames(22))//' is greater than 0.2 meters, which could cause inaccuracies in zone cooling energy.')
!    END IF

  END DO

  RETURN

END SUBROUTINE GetFrameAndDividerData

SUBROUTINE SearchWindow5DataFile(DesiredFileName,DesiredConstructionName,ConstructionFound,EOFonFile,ErrorsFound)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Fred Winkelmann
          !       DATE WRITTEN   August 2001
          !       MODIFIED       June 2002, FW: do all reallocation here for constructions found on
          !                        data file; 1 new construction of entry has one glazing system;
          !                        2 new constructions if entry has two glazing systems.
          !                      Nov 2002, FW: skip read of mullion data line if one glazing system;
          !                        add error messages for bad data; increase length of input line
          !                        from 132 to 200 to handle case where Window5 puts in extra blanks
          !                        in gas data line.
          !                      Feb 2007, LKL: Add more checks on Window5DataFile
          !                      Jan 2008, LKL: Change Edge/Cond ratio check.
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! Searches the WINDOW5 data file for a window with the name "DesiredConstructionName,"
          ! which is the name of an idf Construction input using CONSTRUCTION FROM WINDOW5 DATA FILE.
          ! (The WINDOW5 data file contains data for one or more complete windows --
          ! glazing, frame, mullion, and divider.
          ! WINDOW5 writes the data file for export to EnergyPlus so that an annual energy
          ! analysis can be done on exactly the same window without having to re-input into
          ! EnergyPlus.)

          ! If a match is found, a Construction is created and the Material objects associated with
          ! the Construction are created. If there is an associated frame or
          ! divider in the Window5 data file for this Construction, a FrameAndDivider object will
          ! also be created.

          ! If the window on the data file has two glazing systems, a second Construction (and its
          ! associated materials) corresponding to the second glazing system is created.

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor, ONLY: MakeUPPERcase
  USE DataStringGlobals
  USE General, ONLY: POLYF,TrimSigDigits ! POLYF       ! Polynomial in cosine of angle of incidence
  USE DataSystemVariables, ONLY: iASCII_CR, iUnicode_end,GoodIOStatValue,TempFullFileName,CheckForActualFileName

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*),INTENT(IN) :: DesiredFileName ! File name that contains the Window5 constructions.
  CHARACTER(len=*),INTENT(IN) :: DesiredConstructionName ! Name that will be searched for in the Window5 data file
  LOGICAL, INTENT(OUT)        :: ConstructionFound       ! True if DesiredConstructionName is in the Window5 data file
  LOGICAL, INTENT(INOUT)      :: ErrorsFound             ! True if there is a problem with the entry requested from the data file
  LOGICAL, INTENT(OUT)        :: EOFonFile               ! True if EOF during file read

          ! SUBROUTINE PARAMETER DEFINITIONS:
  CHARACTER(len=*),PARAMETER :: NumName(5)=(/'1','2','3','4','5'/)

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER, SAVE      :: W5DataFileNum
  INTEGER            :: FileLineCount       ! counter for number of lines read (used in some error messages)
  CHARACTER(len=200) :: DataLine(100)       ! Array of data lines
  CHARACTER(len=200) :: NextLine            ! Line of data
  CHARACTER(len=MaxNameLength) :: WindowNameInW5DataFile, W5Name
  CHARACTER(len=10)  :: GasName(3)          ! Gas name from data file
  CHARACTER(len=20)  :: LayerName           ! Layer name from data file
  CHARACTER(len=10)  :: MullionOrientation  ! Horizontal, vertical or none
  INTEGER            :: LineNum
  INTEGER            :: ReadStat            ! File read status
  INTEGER            :: NGlass(2)           ! Number of glass layers in glazing system
  INTEGER            :: NumGases(2,4)       ! Number of gases in each gap of a glazing system
  INTEGER            :: MaterNumSysGlass(2,5) ! Material numbers for glazing system / glass combinations
  INTEGER            :: MaterNumSysGap(2,4) ! Material numbers for glazing system / gap combinations
  INTEGER            :: TotMaterialsPrev    ! Number of materials before adding ones from W5DataFile
  INTEGER            :: TotFrameDividerPrev ! Number of FrameAndDivider objects before adding ones from W5DataFile
  INTEGER            :: NGaps(2)            ! Number of gaps in window construction
  INTEGER            :: NGlSys              ! Number of glazing systems (normally 1, but 2 for mullioned window
                                            !  with two different glazing systems
  INTEGER            :: loop                ! DO loop counter
  INTEGER            :: ILine               ! Line counter
  INTEGER            :: ConstrNum           ! Construction number
  INTEGER            :: IGlass              ! Glass layer counter
  INTEGER            :: IGap                ! Gap counter
  INTEGER            :: IGas                ! Gas counter
!  INTEGER            :: ICoeff              ! Gas property coefficient counter
  INTEGER            :: IGlSys              ! Glazing system counter
  INTEGER            :: MaterNum,MatNum     ! Material number
  INTEGER            :: FrDivNum            ! FrameDivider number
  LOGICAL            :: exists              ! True if Window5 data file exists
  REAL(r64)          :: WinHeight(2),WinWidth(2)  ! Height, width for glazing system (m)
  REAL(r64)          :: UValCenter(2)       ! Center of glass U-value (W/m2-K) for glazing system
  REAL(r64)          :: SCCenter(2)         ! Center of glass shading coefficient for glazing system
  REAL(r64)          :: SHGCCenter(2)       ! Center of glass solar heat gain coefficient for glazing system
  REAL(r64)          :: TVisCenter(2)       ! Center of glass visible transmittance for glazing system
  REAL(r64)          :: Tsol(11)            ! Solar transmittance vs incidence angle; diffuse trans.
  REAL(r64)          :: AbsSol(5,11)        ! Solar absorptance vs inc. angle in each glass layer
  REAL(r64)          :: Rfsol(11)           ! Front solar reflectance vs inc. angle
  REAL(r64)          :: Rbsol(11)           ! Back solar reflectance vs inc. angle
  REAL(r64)          :: Tvis(11)            ! Visible transmittance vs inc. angle
  REAL(r64)          :: Rfvis(11)           ! Front visible reflectance vs inc. angle
  REAL(r64)          :: Rbvis(11)           ! Back visible reflectance vs inc. angle
  REAL(r64)          :: CosPhiIndepVar(10)  ! Cosine of incidence angle from 0 to 90 deg in 10 deg increments
  INTEGER            :: IPhi                ! Incidence angle counter
  REAL(r64)          :: Phi                 ! Incidence angle (deg)
  REAL(r64)          :: CosPhi              ! Cosine of incidence angle
  REAL(r64)          :: tsolFit(10)         ! Fitted solar transmittance vs incidence angle
  REAL(r64)          :: tvisFit(10)         ! Fitted visible transmittance vs incidence angle
  REAL(r64)          :: rfsolFit(10)        ! Fitted solar front reflectance vs incidence angle
  REAL(r64)          :: solabsFit(10,5)     ! Fitted solar absorptance vs incidence angle for each glass layer
  INTEGER,EXTERNAL   :: GetNewUnitNumber
  CHARACTER(len=20)  :: DividerType(2)      ! Divider type: DividedLite or Suspended
  REAL(r64)          :: FrameWidth
  REAL(r64)          :: MullionWidth
  REAL(r64)          :: FrameProjectionOut
  REAL(r64)          :: FrameProjectionIn
  REAL(r64)          :: FrameConductance
  REAL(r64)          :: FrEdgeToCenterGlCondRatio
  REAL(r64)          :: FrameSolAbsorp
  REAL(r64)          :: FrameVisAbsorp
  REAL(r64)          :: FrameEmis
  INTEGER            :: HorDividers(2)      ! For divider: number horizontal for each glazing system
  INTEGER            :: VertDividers(2)     ! For divider: number vertical for each glazing system
  REAL(r64)          :: DividerWidth(2)
  REAL(r64)          :: DividerProjectionOut(2)
  REAL(r64)          :: DividerProjectionIn(2)
  REAL(r64)          :: DividerConductance(2)
  REAL(r64)          :: DivEdgeToCenterGlCondRatio(2)
  REAL(r64)          :: DividerSolAbsorp(2)
  REAL(r64)          :: DividerVisAbsorp(2)
  REAL(r64)          :: DividerEmis(2)
  INTEGER :: endcol
  LOGICAL :: StripCR
  TYPE (FrameDividerProperties), ALLOCATABLE, DIMENSION(:) :: FrameDividerSave
  
  INTEGER :: DebugFile       =150 !RS: Debugging file denotion, hopefully this works.
    
  OPEN(unit=DebugFile,file='Debug.txt')    !RS: Debugging

                                            ! In the following four gas-related data sets, the first
                                            !  index is gas type (1=air, 2=Argon, 3=Krypton, 4=Xenon)
                                            !  and the second index gives a,b,c in the expression
                                            !  property value = a + bT(K) + cT(K)**2, where T is mean
                                            !  gap temperature in deg K.

ConstructionFound = .FALSE.
!ErrorsFound = .FALSE.
EOFonFile = .FALSE.

CALL CheckForActualFileName(DesiredFileName,exists,TempFullFileName)
!INQUIRE(FILE=TRIM(DesiredFileName), EXIST=exists)
IF(.NOT.exists) THEN
  !CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: '//   &
  !      'Could not locate Window5 Data File, expecting it as file name='//TRIM(DesiredFileName))
  !CALL ShowContinueError('Certain run environments require a full path to be included with the file name in the input field.')
  !CALL ShowContinueError('Try again with putting full path and file name in the field.')
  !CALL ShowFatalError('Program terminates due to these conditions.')    !RS: Secret Search String
  WRITE(DebugFile,*) 'HeatBalanceManager: SearchWindow5DataFile: Could not locate Window5 Data File,'// &
    ' expecting it as file name='//TRIM(DesiredFileName)
  WRITE(DebugFile,*) 'Certain run environments require a full path to be included with the file name in the input field.'
  WRITE(DebugFile,*) 'Try again with putting full path and file name in the field.'
  WRITE(DebugFile,*) 'Program terminates due to these conditions. (Wishful thinking!)'
ENDIF

W5DataFileNum = GetNewUnitNumber()
OPEN(UNIT=W5DataFileNum, FILE=TempFullFileName, Action='read', Err=999)
StripCR=.false.
READ(Unit=W5DataFileNum, FMT=fmtA) NextLine
endcol=LEN_TRIM(NextLine)
IF (endcol > 0) THEN
  IF (ICHAR(NextLine(endcol:endcol)) == iASCII_CR) THEN
    StripCR=.true.
  ENDIF
  IF (ICHAR(NextLine(endcol:endcol)) == iUnicode_end) THEN
    CALL ShowSevereError('SearchWindow5DataFile: For "'//TRIM(DesiredConstructionName)//'" in '//TRIM(DesiredFileName)//  &
     ' file, appears to be a Unicode file.')
    CALL ShowContinueError('...This file cannot be read by this program. Please save as PC or Unix file and try again')
    CALL ShowFatalError('Program terminates due to previous condition.')
  ENDIF
ENDIF

REWIND(W5DataFileNum)
FileLineCount=0

READ(W5DataFileNum,fmtA,IOSTAT=ReadStat) NextLine
IF(ReadStat < GoodIOStatValue) GO TO 1000
IF (StripCR) THEN
  endcol=LEN_TRIM(NextLine)
  IF (endcol > 0) NextLine(endcol:endcol)=Blank
ENDIF
FileLineCount=FileLineCount+1
IF(MakeUPPERCase(NextLine(1:7)) /= 'WINDOW5') THEN
  CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Error in Data File='//TRIM(DesiredFileName))
  CALL ShowFatalError('Error reading Window5 Data File: first word of window entry is "' &
         //TRIM(NextLine(1:7))//'", should be Window5.')
END IF

10 DO LineNum = 2,5
      READ(W5DataFileNum,'(A)',IOSTAT=ReadStat) DataLine(LineNum)
      IF(ReadStat < GoodIOStatValue) GO TO 1000
      IF (StripCR) THEN
        endcol=LEN_TRIM(DataLine(LineNum))
        IF (endcol > 0) DataLine(LineNum)(endcol:endcol)=Blank
      ENDIF
      FileLineCount=FileLineCount+1
  END DO

  ! Get window name and check for match
  READ(Dataline(4)(20:),fmtA) W5Name
  WindowNameInW5DataFile = MakeUPPERcase(W5Name)
  IF(TRIM(DesiredConstructionName) /= TRIM(WindowNameInW5DataFile)) THEN
    ! Doesn't match; read through file until next window entry is found
 20 READ(W5DataFileNum,fmtA,IOSTAT=ReadStat) NextLine
    IF(ReadStat < GoodIOStatValue) GO TO 1000
    IF (StripCR) THEN
      endcol=LEN_TRIM(NextLine)
      IF (endcol > 0) NextLine(endcol:endcol)=Blank
    ENDIF
    FileLineCount=FileLineCount+1
    IF (MakeUPPERCase(NextLine(1:7)) /= 'WINDOW5') GOTO 20
    ! Beginning of next window entry found
    GO TO 10
  ELSE
    ! Match found
    ConstructionFound = .TRUE.

    ! Create Material:WindowGlass, Material:WindowGas, Construction
    ! and WindowFrameAndDividerObjects for this window

    READ(W5DataFileNum,fmtA,IOSTAT=ReadStat) NextLine
    IF(ReadStat < GoodIOStatValue) GO TO 1000
    IF (StripCR) THEN
      endcol=LEN_TRIM(NextLine)
      IF (endcol > 0) NextLine(endcol:endcol)=Blank
    ENDIF
    FileLineCount=FileLineCount+1
    READ(NextLine(20:),*) NGlSys
      IF(NGlSys <= 0 .OR. NGlSys > 2) THEN
        CALL ShowFatalError('Construction='//TRIM(DesiredConstructionName)// &
          ' from the Window5 data file cannot be used: it has '&
          //TRIM(TrimSigDigits(NGlSys))//' glazing systems; only 1 or 2 are allowed.')
      END IF
    READ(W5DataFileNum,fmtA,IOSTAT=ReadStat) NextLine
    IF(ReadStat < GoodIOStatValue) GO TO 1000
    IF (StripCR) THEN
      endcol=LEN_TRIM(NextLine)
      IF (endcol > 0) NextLine(endcol:endcol)=Blank
    ENDIF
    FileLineCount=FileLineCount+1
    DO IGlSys = 1, NGlSys
      READ(W5DataFileNum,fmtA,IOSTAT=ReadStat) NextLine
      IF(ReadStat < GoodIOStatValue) GO TO 1000
      IF (StripCR) THEN
        endcol=LEN_TRIM(NextLine)
        IF (endcol > 0) NextLine(endcol:endcol)=Blank
      ENDIF
      FileLineCount=FileLineCount+1
      READ(NextLine(20:),*,IOSTAT=ReadStat)  WinHeight(IGlSys),WinWidth(IGlSys),NGlass(IGlSys),UvalCenter(IGlSys), &
                              SCCenter(IGlSys),SHGCCenter(IGlSys),TvisCenter(IGlSys)
        IF (ReadStat /= 0) THEN
          CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Error in Read of glazing system values.'//  &
              ' For glazing system='//TRIM(TrimSigDigits(IGlSys)))
          CALL ShowContinueError('Line (~'//TRIM(TrimSigDigits(FileLineCount))//') in error (first 100 characters)='//  &
                 TRIM(NextLine(1:100)))
          ErrorsFound=.true.
        ENDIF
        IF(WinHeight(IGlSys) == 0.0 .OR. WinWidth(IGlSys) == 0.0) THEN
          CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Construction='//TRIM(DesiredConstructionName)// &
            ' from the Window5 data file cannot be used:'// &
            ' it has window height or width = 0 for glazing system '//TRIM(TrimSigDigits(IGlSys)))
          ErrorsFound = .TRUE.
        END IF
        IF(NGlass(IGlSys) <= 0 .OR. NGlass(IGlSys) > 4) THEN
          CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Construction='//TRIM(DesiredConstructionName)// &
            ' from the Window5 data file cannot be used:'// &
            ' it has 0 or more than 4 glass layers in glazing system '//TRIM(TrimSigDigits(IGlSys)))
          ErrorsFound = .TRUE.
        END IF
        IF( UvalCenter(IGlSys) <= 0.0) THEN
          CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Construction='//TRIM(DesiredConstructionName)// &
            ' from the Window5 data file cannot be used:'// &
            ' it has Center-of-Glass U-value <= 0 in glazing system '//TRIM(TrimSigDigits(IGlSys)))
          ErrorsFound = .TRUE.
        END IF
        IF(SCCenter(IGlSys) <= 0.0) THEN
          CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Construction='//TRIM(DesiredConstructionName)// &
            ' from the Window5 data file cannot be used:'// &
            ' it has Shading Coefficient <= 0 in glazing system '//TRIM(TrimSigDigits(IGlSys)))
          ErrorsFound = .TRUE.
        END IF
        IF(SHGCCenter(IGlSys) <= 0.0) THEN
          CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Construction='//TRIM(DesiredConstructionName)// &
            ' from the Window5 data file cannot be used:'// &
            ' it has SHGC <= 0 in glazing system '//TRIM(TrimSigDigits(IGlSys)))
          ErrorsFound = .TRUE.
        END IF
      WinHeight(IGlSys) = 0.001d0*WinHeight(IGlSys)
      WinWidth(IGlSys)  = 0.001d0*WinWidth(IGlSys)
    END DO
    DO LineNum = 1,11
      READ(W5DataFileNum,fmtA,IOSTAT=ReadStat) DataLine(LineNum)
      IF(ReadStat == -1) GO TO 1000
      IF (StripCR) THEN
        endcol=LEN_TRIM(DataLine(LineNum))
        IF (endcol > 0) DataLine(LineNum)(endcol:endcol)=Blank
      ENDIF
    END DO

    ! Mullion width and orientation
    MullionWidth = 0.0
    MullionOrientation = 'Vertical'
    IF(NGlSys == 2) THEN
      READ(Dataline(10)(20:),*,IOSTAT=ReadStat) MullionWidth
        IF (ReadStat /= 0) THEN
          CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Error in Read of Mullion Width.')
          CALL ShowContinueError('Line (~'//TRIM(TrimSigDigits(FileLineCount+10))//') in error (first 100 characters)='//  &
                 TRIM(DataLine(10)(1:100)))
          ErrorsFound=.true.
        ENDIF
      MullionWidth = 0.001*MullionWidth
      READ(Dataline(10)(89:),*,IOSTAT=ReadStat) MullionOrientation
        IF (ReadStat /= 0) THEN
          CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Error in Read of Mullion Orientation.')
          CALL ShowContinueError('Line (~'//TRIM(TrimSigDigits(FileLineCount+10))//') in error (first 100 characters)='//  &
                 TRIM(DataLine(10)(1:100)))
          ErrorsFound=.true.
        ENDIF
    END IF

    ! Frame data; if there are two glazing systems, the frame is assumed to be
    ! the same for both.
    FrameWidth=0.0
    FrameProjectionOut=0.0
    FrameProjectionIn=0.0
    FrameConductance=0.0
    FrEdgeToCenterGlCondRatio=0.0
    FrameSolAbsorp=0.0
    FrameVisAbsorp=0.0
    FrameEmis=0.0
    READ(DataLine(11)(20:),*,IOStat=ReadStat) FrameWidth,FrameProjectionOut,FrameProjectionIn,FrameConductance, &
      FrEdgeToCenterGlCondRatio,FrameSolAbsorp,FrameVisAbsorp,FrameEmis
      IF (ReadStat /= 0) THEN
        CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Error in Read of frame data values.')
        CALL ShowContinueError('Line (~'//TRIM(TrimSigDigits(FileLineCount+11))//') in error (first 100 characters)='//  &
               TRIM(DataLine(11)(1:100)))
        ErrorsFound=.true.
      ENDIF
      IF(FrameWidth > 0.0) THEN
        IF(FrameConductance <= 0.0) THEN
            CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Construction='//TRIM(DesiredConstructionName)// &
            ' from the Window5 data file cannot be used: it has Frame Conductance <= 0.0')
          ErrorsFound = .TRUE.
        END IF
! Relax this check for Window5 data: 1/28/2008.
!        IF(FrEdgeToCenterGlCondRatio < 1.0) THEN
!            CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Construction='//TRIM(DesiredConstructionName)// &
!            ' from the Window5 data file cannot be used: it has Frame Edge-of-Glass Conduction Ratio < 1.0')
!          ErrorsFound = .TRUE.
!        END IF
        IF(FrameSolAbsorp < 0.0 .OR. FrameSolAbsorp > 1.0) THEN
            CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Construction='//TRIM(DesiredConstructionName)// &
            ' from the Window5 data file cannot be used: it has Frame Solar Absorptance < 0.0 or > 1.0')
          ErrorsFound = .TRUE.
        END IF
        IF(FrameEmis <= 0.0 .OR. FrameEmis >= 1.0) THEN
            CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Construction='//TRIM(DesiredConstructionName)// &
            ' from the Window5 data file cannot be used: it has Frame Emissivity <= 0.0 or >= 1.0')
          ErrorsFound = .TRUE.
        END IF
      END IF
      FrameWidth         = 0.001d0*FrameWidth
      FrameProjectionOut = 0.001d0*FrameProjectionOut
      FrameProjectionIn  = 0.001d0*FrameProjectionIn
    FileLineCount=FileLineCount+11

    READ(W5DataFileNum,fmtA,IOSTAT=ReadStat) NextLine
    IF(ReadStat < GoodIOStatValue) GO TO 1000
    IF (StripCR) THEN
      endcol=LEN_TRIM(NextLine)
      IF (endcol > 0) NextLine(endcol:endcol)=Blank
    ENDIF
    FileLineCount=FileLineCount+1

    ! Divider data for each glazing system
    DO IGlSys = 1,NGlSys
      READ(W5DataFileNum,fmtA,IOSTAT=ReadStat) NextLine
      IF(ReadStat < GoodIOStatValue) GO TO 1000
      FileLineCount=FileLineCount+1
      READ(NextLine(20:),*,IOSTAT=ReadStat) DividerWidth(IGlSys),DividerProjectionOut(IGlSys),DividerProjectionIn(IGlSys), &
        DividerConductance(IGlSys),DivEdgeToCenterGlCondRatio(IGlSys),DividerSolAbsorp(IGlSys),DividerVisAbsorp(IGlSys), &
        DividerEmis(IGlSys),DividerType(IGlSys),HorDividers(IGlSys),VertDividers(IGlSys)
      IF (ReadStat /= 0) THEN
        CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Error in Read of divider data values. '//  &
           'For Glazing System='//TRIM(TrimSigDigits(IGLSys)))
        CALL ShowContinueError('Line (~'//TRIM(TrimSigDigits(FileLineCount+11))//') in error (first 100 characters)='//  &
               TRIM(NextLine(1:100)))
        ErrorsFound=.true.
      ENDIF
      DividerType(IGlSys) = MakeUpperCase(DividerType(IGlSys))
      IF(DividerWidth(IGlSys) > 0.0) THEN
        IF(HorDividers(IGlSys) == 0 .AND. VertDividers(IGlSys) == 0) THEN
          CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Construction='//TRIM(DesiredConstructionName)// &
            ' from the Window5 data file cannot be used:')
          CALL ShowContinueError('glazing system '//TRIM(TrimSigDigits(IGLSys))// &
            ' has a divider but number of horizontal and vertical divider elements = 0')
          ErrorsFound = .TRUE.
        END IF
        IF(DividerConductance(IGlSys) <= 0.0) THEN
          CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Construction='//TRIM(DesiredConstructionName)// &
            ' from the Window5 data file cannot be used:')
          CALL ShowContinueError('glazing system '//TRIM(TrimSigDigits(IGLSys))//' has Divider Conductance <= 0.0')
          ErrorsFound = .TRUE.
        END IF
        IF(DivEdgeToCenterGlCondRatio(IGlSys) < 1.0) THEN
          CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Construction='//TRIM(DesiredConstructionName)// &
            ' from the Window5 data file cannot be used:')
          CALL ShowContinueError('glazing system '//TRIM(TrimSigDigits(IGLSys))// &
            ' has Divider Edge-Of-Glass Conduction Ratio < 1.0')
          ErrorsFound = .TRUE.
        END IF
        IF(DividerSolAbsorp(IGlSys) < 0.0 .OR. DividerSolAbsorp(IGlSys) > 1.0) THEN
          CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Construction='//TRIM(DesiredConstructionName)// &
            ' from the Window5 data file cannot be used:')
          CALL ShowContinueError('glazing system '//TRIM(TrimSigDigits(IGLSys))// &
            ' has Divider Solar Absorptance < 0.0 or > 1.0')
          ErrorsFound = .TRUE.
        END IF
        IF(DividerEmis(IGlSys) <= 0.0 .OR. DividerEmis(IGlSys) >= 1.0) THEN
          CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Construction='//TRIM(DesiredConstructionName)// &
            ' from the Window5 data file cannot be used:')
          CALL ShowContinueError('glazing system '//TRIM(TrimSigDigits(IGLSys))// &
            ' has Divider Emissivity <= 0.0 or >= 1.0')
          ErrorsFound = .TRUE.
        END IF
        IF(DividerType(IGlSys) /= 'DIVIDEDLITE' .AND. DividerType(IGlSys) /= 'SUSPENDED') THEN
          CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Construction='//TRIM(DesiredConstructionName)// &
            ' from the Window5 data file cannot be used:')
          CALL ShowContinueError('glazing system '//TRIM(TrimSigDigits(IGLSys))// &
            ' has Divider Type = '//TRIM(DividerType(IGlSys))// &
            '; it should be DIVIDEDLITE or SUSPENDED.')
          ErrorsFound = .TRUE.
        END IF
      END IF
      DividerWidth(IGlSys)         = 0.001d0*DividerWidth(IGlSys)
      IF(DividerType(IGlSys) == 'DIVIDEDLITE') THEN
        DividerProjectionOut(IGlSys) = 0.001d0*DividerProjectionOut(IGlSys)
        DividerProjectionIn(IGlSys)  = 0.001d0*DividerProjectionIn(IGlSys)
      ELSE
        DividerProjectionOut(IGlSys) = 0.0
        DividerProjectionIn(IGlSys)  = 0.0
      END IF
    END DO

    IF(ErrorsFound)   &
       CALL ShowFatalError('HeatBalanceManager: SearchWindow5DataFile: Construction='//TRIM(DesiredConstructionName)// &
            ' from the Window5 data file cannot be used because of above errors')

    TotMaterialsPrev = TotMaterials
    DO IGlSys = 1,NGlSys
      NGaps(IGlSys) = NGlass(IGlSys)-1
      TotMaterials = TotMaterials + NGlass(IGlSys) + NGaps(IGlSys)
    END DO

    ! Create Material objects

    ! reallocate Material type

    ALLOCATE(MaterialSave(TotMaterialsPrev))
    ALLOCATE(NominalRSave(TotMaterialsPrev))
    NominalRSave=NominalR
    DO loop = 1,TotMaterialsPrev
      MaterialSave(loop) = Material(loop)
!      NominalRSave(loop) = NominalR(loop)
    END DO
    DEALLOCATE(Material)
    DEALLOCATE(NominalR)
    ALLOCATE(Material(TotMaterials))
    ALLOCATE(NominalR(TotMaterials))
    NominalR=0.0
    NominalR(1:TotMaterialsPrev)=NominalRSave
    DO loop = 1,TotMaterialsPrev
      Material(loop) = MaterialSave(loop)
!      NominalR(loop) = NominalRSave(loop)
    END DO
    DEALLOCATE(MaterialSave)
    DEALLOCATE(NominalRSave)

    ! Initialize new materials
    DO loop = TotMaterialsPrev+1,TotMaterials
      Material(loop)%Name=' '
      Material(loop)%Group=-1
      Material(loop)%Roughness=0
      Material(loop)%Conductivity=0.0
      Material(loop)%Density=0.0
      Material(loop)%IsoMoistCap=0.0
      Material(loop)%Porosity=0.0
      Material(loop)%Resistance=0.0
      Material(loop)%SpecHeat=0.0
      Material(loop)%ThermGradCoef=0.0
      Material(loop)%Thickness=0.0
      Material(loop)%VaporDiffus=0.0
      Material(loop)%AbsorpSolar=0.0
      Material(loop)%AbsorpThermal=0.0
      Material(loop)%AbsorpVisible=0.0
      Material(loop)%ReflectShade=0.0
      Material(loop)%Trans=0.0
      Material(loop)%ReflectShadeVis=0.0
      Material(loop)%TransVis=0.0
      Material(loop)%GlassTransDirtFactor=1.0
      Material(loop)%SolarDiffusing=.false.
      Material(loop)%AbsorpThermalBack=0.0
      Material(loop)%AbsorpThermalFront=0.0
      Material(loop)%ReflectSolBeamBack=0.0
      Material(loop)%ReflectSolBeamFront=0.0
      Material(loop)%ReflectSolDiffBack=0.0
      Material(loop)%ReflectSolDiffFront=0.0
      Material(loop)%ReflectVisBeamBack=0.0
      Material(loop)%ReflectVisBeamFront=0.0
      Material(loop)%ReflectVisDiffBack=0.0
      Material(loop)%ReflectVisDiffFront=0.0
      Material(loop)%TransSolBeam=0.0
      Material(loop)%TransThermal=0.0
      Material(loop)%TransVisBeam=0.0
      Material(loop)%GlassSpectralDataPtr=0
      Material(loop)%NumberOfGasesInMixture=0
      Material(loop)%GasCon=0.0
      Material(loop)%GasVis=0.0
      Material(loop)%GasCp=0.0
      Material(loop)%GasType=0
      Material(loop)%GasWght=0.0
      Material(loop)%GasFract=0.0
      Material(loop)%WinShadeToGlassDist=0.0
      Material(loop)%WinShadeTopOpeningMult=0.0
      Material(loop)%WinShadeBottomOpeningMult=0.0
      Material(loop)%WinShadeLeftOpeningMult=0.0
      Material(loop)%WinShadeRightOpeningMult=0.0
      Material(loop)%WinShadeAirFlowPermeability=0.0
      Material(loop)%BlindDataPtr=0
      Material(loop)%EMPDVALUE=0.0
      Material(loop)%MoistACoeff=0.0
      Material(loop)%MoistBCoeff=0.0
      Material(loop)%MoistCCoeff=0.0
      Material(loop)%MoistDCoeff=0.0
    END DO

    ! Glass objects
    READ(W5DataFileNum,fmtA,IOSTAT=ReadStat) NextLine
    IF(ReadStat < GoodIOStatValue) GO TO 1000
    IF (StripCR) THEN
      endcol=LEN_TRIM(NextLine)
      IF (endcol > 0) NextLine(endcol:endcol)=Blank
    ENDIF
    FileLineCount=FileLineCount+1
    MaterNum = TotMaterialsPrev
    DO IGlSys = 1,NGlSys
      DO IGlass = 1,NGlass(IGlSys)
        MaterNum = MaterNum + 1
        MaterNumSysGlass(IGlSys,IGlass) = MaterNum
        Material(MaterNum)%Group = WindowGlass
        READ(W5DataFileNum,'(A)',IOSTAT=ReadStat) NextLine
        FileLineCount=FileLineCount+1
        READ(NextLine(26:),*) &
          Material(MaterNum)%Thickness,                   &
          Material(MaterNum)%Conductivity,                &
          Material(MaterNum)%Trans,                       &
          Material(MaterNum)%ReflectSolBeamFront,         &
          Material(MaterNum)%ReflectSolBeamBack,          &
          Material(MaterNum)%TransVis,                    &
          Material(MaterNum)%ReflectVisBeamFront,         &
          Material(MaterNum)%ReflectVisBeamBack,          &
          Material(MaterNum)%TransThermal,                &
          Material(MaterNum)%AbsorpThermalFront,          &
          Material(MaterNum)%AbsorpThermalBack,           &
          LayerName
        Material(MaterNum)%Thickness = 0.001*Material(MaterNum)%Thickness
        IF (Material(MaterNum)%Thickness <= 0.0) THEN
        ENDIF
        IF(NGlSys == 1) THEN
          Material(MaterNum)%Name = &
            'W5:'//TRIM(DesiredConstructionName)//':GLASS'//NumName(IGlass)
        ELSE
          Material(MaterNum)%Name = &
            'W5:'//TRIM(DesiredConstructionName)//':'//NumName(IGlSys)//':GLASS'//NumName(IGlass)
        END IF
        Material(MaterNum)%Roughness = VerySmooth
        Material(MaterNum)%AbsorpThermal = Material(MaterNum)%AbsorpThermalBack
        IF (Material(MaterNum)%Thickness <= 0.0) THEN
          CALL ShowSevereError('SearchWindow5DataFile: Material="'//trim(Material(MaterNum)%Name)//  &
            '" has thickness of 0.0.  Will be set to thickness = .001 but inaccuracies may result.')
          CALL ShowContinueError('Line being read='//trim(NextLine))
          CALL ShowContinueError('Thickness field starts at column 26='//trim(NextLine(26:)))
          Material(MaterNum)%Thickness=.001d0
        ENDIF
      END DO
    END DO

    ! Gap objects
    READ(W5DataFileNum,fmtA,IOSTAT=ReadStat) NextLine
    IF(ReadStat < GoodIOStatValue) GO TO 1000
    IF (StripCR) THEN
      endcol=LEN_TRIM(NextLine)
      IF (endcol > 0) NextLine(endcol:endcol)=Blank
    ENDIF
    FileLineCount=FileLineCount+1
    DO IGlSys = 1,NGlSys
      DO IGap = 1,NGaps(IGlSys)
        MaterNum = MaterNum + 1
        MaterNumSysGap(IGlSys,IGap) = MaterNum
        READ(W5DataFileNum,fmtA,IOSTAT=ReadStat) NextLine
        FileLineCount=FileLineCount+1
        READ(NextLine(24:),*) Material(MaterNum)%Thickness,NumGases(IGlSys,IGap)
        IF(NGlSys == 1) THEN
          Material(MaterNum)%Name = &
            'W5:'//TRIM(DesiredConstructionName)//':GAP'//NumName(IGap)
        ELSE
          Material(MaterNum)%Name = &
            'W5:'//TRIM(DesiredConstructionName)//':'//NumName(IGlSys)//':GAP'//NumName(IGap)
        END IF
        Material(MaterNum)%Thickness = 0.001d0*Material(MaterNum)%Thickness
        Material(MaterNum)%RoughNess = MediumRough  ! Unused
      END DO
    END DO

    READ(W5DataFileNum,fmtA,IOSTAT=ReadStat) NextLine
    IF(ReadStat < GoodIOStatValue) GO TO 1000
    IF (StripCR) THEN
      endcol=LEN_TRIM(NextLine)
      IF (endcol > 0) NextLine(endcol:endcol)=Blank
    ENDIF
    FileLineCount=FileLineCount+1
    DO IGlSys = 1,NGlSys
      DO IGap = 1,NGaps(IGlSys)
        MaterNum = MaterNumSysGap(IGlSys,IGap)
        Material(MaterNum)%NumberOfGasesInMixture = NumGases(IGlSys,IGap)
        Material(MaterNum)%Group = WindowGas
        IF(NumGases(IGlSys,IGap) > 1) Material(MaterNum)%Group = WindowGasMixture
        DO IGas = 1,NumGases(IGlSys,IGap)
          READ(W5DataFileNum,'(A)',IOSTAT=ReadStat) NextLine
          FileLineCount=FileLineCount+1
          READ(NextLine(20:),*) &
            GasName(IGas),Material(MaterNum)%GasFract(IGas),Material(MaterNum)%GasWght(IGas), &
            Material(MaterNum)%GasCon(IGas,:),Material(MaterNum)%GasVis(IGas,:),Material(MaterNum)%GasCp(IGas,:)
            ! Nominal resistance of gap at room temperature (based on first gas in mixture)
            NominalR(MaterNum)=Material(MaterNum)%Thickness/(Material(MaterNum)%GasCon(1,1) + &
              Material(MaterNum)%GasCon(1,2)*300.0d0 + Material(MaterNum)%GasCon(1,3)*90000.d0)
        END DO
      END DO
    END DO

    ! Construction objects

    ! reallocate Construct types
    ALLOCATE(ConstructSave(TotConstructs))
    ALLOCATE(NominalUSave(TotConstructs))
    DO loop = 1,TotConstructs
      ConstructSave(loop) = Construct(loop)
      NominalUSave(loop) = NominalU(loop)
    END DO
    DEALLOCATE(Construct)
    DEALLOCATE(NominalU)
    TotConstructs = TotConstructs + NGlSys
    ALLOCATE(Construct(TotConstructs))
    ALLOCATE(NominalU(TotConstructs))
    DO loop = 1,TotConstructs-NGlSys
      Construct(loop) = ConstructSave(loop)
      NominalU(loop) = NominalUSave(loop)
    END DO
    DEALLOCATE(ConstructSave)
    DEALLOCATE(NominalUSave)

    READ(W5DataFileNum,fmtA,IOSTAT=ReadStat) NextLine
    IF(ReadStat < GoodIOStatValue) GO TO 1000
    IF (StripCR) THEN
      endcol=LEN_TRIM(NextLine)
      IF (endcol > 0) NextLine(endcol:endcol)=Blank
    ENDIF
    FileLineCount=FileLineCount+1

    DO IGlSys = 1,NGlSys
      ConstrNum = TotConstructs - NGlSys + IGlSys
      IF(IGlSys == 1) THEN
        Construct(ConstrNum)%Name = TRIM(DesiredConstructionName)
      ELSE
        Construct(ConstrNum)%Name = TRIM(DesiredConstructionName)//':2'
      END IF
      DO loop = 1,MaxLayersInConstruct
        Construct(ConstrNum)%LayerPoint(loop) = 0
      END DO
      Construct(ConstrNum)%InsideAbsorpSolar   = 0.0
      Construct(ConstrNum)%OutsideAbsorpSolar  = 0.0
      Construct(ConstrNum)%DayltPropPtr        = 0
      Construct(ConstrNum)%CTFCross            = 0.0D0
      Construct(ConstrNum)%CTFFlux             = 0.0D0
      Construct(ConstrNum)%CTFInside           = 0.0D0
      Construct(ConstrNum)%CTFOutside          = 0.0D0
      Construct(ConstrNum)%CTFSourceIn         = 0.0D0
      Construct(ConstrNum)%CTFSourceOut        = 0.0D0
      Construct(ConstrNum)%CTFTimeStep         = 0.0D0
      Construct(ConstrNum)%CTFTSourceOut       = 0.0D0
      Construct(ConstrNum)%CTFTSourceIn        = 0.0D0
      Construct(ConstrNum)%CTFTSourceQ         = 0.0D0
      Construct(ConstrNum)%CTFTUserOut         = 0.0D0
      Construct(ConstrNum)%CTFTUserIn          = 0.0D0
      Construct(ConstrNum)%CTFTUserSource      = 0.0D0
      Construct(ConstrNum)%NumHistories        = 0
      Construct(ConstrNum)%NumCTFTerms         = 0
      Construct(ConstrNum)%UValue              = 0.0
      Construct(ConstrNum)%SourceSinkPresent   = .FALSE.
      Construct(ConstrNum)%SolutionDimensions  = 0
      Construct(ConstrNum)%SourceAfterLayer    = 0
      Construct(ConstrNum)%TempAfterLayer      = 0
      Construct(ConstrNum)%ThicknessPerpend    = 0.0
      Construct(ConstrNum)%AbsDiff             = 0.0
      Construct(ConstrNum)%AbsDiffBack         = 0.0
      Construct(ConstrNum)%AbsDiffShade        = 0.0
      Construct(ConstrNum)%AbsDiffBackShade    = 0.0
      Construct(ConstrNum)%ShadeAbsorpThermal  = 0.0
      Construct(ConstrNum)%AbsBeamCoef         = 0.0
      Construct(ConstrNum)%AbsBeamBackCoef     = 0.0
      Construct(ConstrNum)%AbsBeamShadeCoef    = 0.0
      Construct(ConstrNum)%AbsDiffIn           = 0.0
      Construct(ConstrNum)%AbsDiffOut          = 0.0
      Construct(ConstrNum)%TransDiff           = 0.0
      Construct(ConstrNum)%TransDiffVis        = 0.0
      Construct(ConstrNum)%ReflectSolDiffBack  = 0.0
      Construct(ConstrNum)%ReflectSolDiffFront = 0.0
      Construct(ConstrNum)%ReflectVisDiffBack  = 0.0
      Construct(ConstrNum)%ReflectVisDiffFront = 0.0
      Construct(ConstrNum)%TransSolBeamCoef    = 0.0
      Construct(ConstrNum)%TransVisBeamCoef    = 0.0
      Construct(ConstrNum)%ReflSolBeamFrontCoef= 0.0
      Construct(ConstrNum)%ReflSolBeamBackCoef = 0.0
      Construct(ConstrNum)%W5FrameDivider      = 0
      Construct(ConstrNum)%TotLayers = NGlass(IGlSys) + NGaps(IGlSys)
      Construct(ConstrNum)%TotGlassLayers = NGlass(IGlSys)
      Construct(ConstrNum)%TotSolidLayers = NGlass(IGlSys)

      DO IGlass = 1,NGlass(IGlSys)
        Construct(ConstrNum)%LayerPoint(2*IGlass-1) = MaterNumSysGlass(IGlSys,IGlass)
        IF(IGlass < NGlass(IGlSys)) Construct(ConstrNum)%LayerPoint(2*IGlass) = MaterNumSysGap(IGlSys,IGlass)
      END DO

      Construct(ConstrNum)%OutsideRoughness    = VerySmooth
      Construct(ConstrNum)%InsideAbsorpThermal = Material(TotMaterialsPrev+NGlass(IGlSys))%AbsorpThermalBack
      Construct(ConstrNum)%OutsideAbsorpThermal= Material(TotMaterialsPrev+1)%AbsorpThermalFront
      Construct(ConstrNum)%TypeIsWindow        = .TRUE.
      Construct(ConstrNum)%FromWindow5DataFile = .TRUE.
      Construct(ConstrNum)%W5FileGlazingSysHeight   = WinHeight(IGlSys)
      Construct(ConstrNum)%W5FileGlazingSysWidth    = WinWidth(IGlSys)
      IF (SameString(MullionOrientation,'Vertical')) THEN
        Construct(ConstrNum)%W5FileMullionOrientation = Vertical
      ELSEIF (SameString(MullionOrientation,'Horizontal')) THEN
        Construct(ConstrNum)%W5FileMullionOrientation = Horizontal
      ELSE
      ENDIF
      Construct(ConstrNum)%W5FileMullionWidth       = MullionWidth

      ! Fill Construct with system transmission, reflection and absorption properties

      DO IPhi = 1,10
        CosPhiIndepVar(IPhi) = COS((IPhi-1)*10.*DegToRadians)
      END DO

      READ(W5DataFileNum,fmtA,IOSTAT=ReadStat) NextLine
      IF(ReadStat < GoodIOStatValue) GO TO 1000
      IF (StripCR) THEN
        endcol=LEN_TRIM(NextLine)
        IF (endcol > 0) NextLine(endcol:endcol)=Blank
      ENDIF
      FileLineCount=FileLineCount+1
      IF(IGlSys == 1) THEN
        READ(W5DataFileNum,fmtA,IOSTAT=ReadStat) NextLine
        IF(ReadStat < GoodIOStatValue) GO TO 1000
        IF (StripCR) THEN
          endcol=LEN_TRIM(NextLine)
          IF (endcol > 0) NextLine(endcol:endcol)=Blank
        ENDIF
        FileLineCount=FileLineCount+1
      ENDIF
      READ(W5DataFileNum,fmtA,IOSTAT=ReadStat) NextLine
      IF(ReadStat < GoodIOStatValue) GO TO 1000
      IF (StripCR) THEN
        endcol=LEN_TRIM(NextLine)
        IF (endcol > 0) NextLine(endcol:endcol)=Blank
      ENDIF
      FileLineCount=FileLineCount+1
      READ(NextLine(6:),*,IOSTAT=ReadStat) Tsol
      IF (ReadStat /= 0) THEN
        CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Error in Read of TSol values.')
          CALL ShowContinueError('Line (~'//TRIM(TrimSigDigits(FileLineCount))//') in error (first 100 characters)='//  &
                 TRIM(NextLine(1:100)))
        ErrorsFound=.true.
      ELSEIF (ANY(Tsol < 0.0) .or. ANY(TSol > 1.0)) THEN
        CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Error in Read of TSol values. (out of range [0,1])')
          CALL ShowContinueError('Line (~'//TRIM(TrimSigDigits(FileLineCount))//') in error (first 100 characters)='//  &
                 TRIM(NextLine(1:100)))
        ErrorsFound=.true.
      ENDIF
      DO IGlass = 1,NGlass(IGlSys)
        READ(W5DataFileNum,fmtA,IOSTAT=ReadStat) NextLine
        IF (StripCR) THEN
          endcol=LEN_TRIM(NextLine)
          IF (endcol > 0) NextLine(endcol:endcol)=Blank
        ENDIF
        FileLineCount=FileLineCount+1
        READ(NextLine(6:),*,IOSTAT=ReadStat) AbsSol(IGlass,:)
        IF (ReadStat /= 0) THEN
          CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Error in Read of AbsSol values. For Glass='//  &
                                  TRIM(TrimSigDigits(IGlass)))
          CALL ShowContinueError('Line (~'//TRIM(TrimSigDigits(FileLineCount))//') in error (first 100 characters)='//  &
                 TRIM(NextLine(1:100)))
          ErrorsFound=.true.
        ELSEIF (ANY(AbsSol(IGlass,:) < 0.0) .or. ANY(AbsSol(IGlass,:) > 1.0)) THEN
          CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Error in Read of AbsSol values. '//  &
                    '(out of range [0,1]) For Glass='//TRIM(TrimSigDigits(IGlass)))
          CALL ShowContinueError('Line (~'//TRIM(TrimSigDigits(FileLineCount))//') in error (first 100 characters)='//  &
                 TRIM(NextLine(1:100)))
          ErrorsFound=.true.
        ENDIF
      END DO
      DO ILine = 1,5
        READ(W5DataFileNum,fmtA,IOSTAT=ReadStat) DataLine(ILine)
        IF (StripCR) THEN
          endcol=LEN_TRIM(DataLine(ILine))
          IF (endcol > 0) DataLine(ILine)(endcol:endcol)=Blank
        ENDIF
      END DO
      READ(DataLine(1)(6:),*,IOSTAT=ReadStat)  Rfsol
      IF (ReadStat /= 0) THEN
        CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Error in Read of RfSol values.')
        CALL ShowContinueError('Line (~'//TRIM(TrimSigDigits(FileLineCount+1))//') in error (first 100 characters)='//  &
                              TRIM(DataLine(1)(1:100)))
        ErrorsFound=.true.
      ELSEIF (ANY(Rfsol < 0.0) .or. ANY(RfSol > 1.0)) THEN
        CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Error in Read of RfSol values. (out of range [0,1])')
        CALL ShowContinueError('Line (~'//TRIM(TrimSigDigits(FileLineCount+1))//') in error (first 100 characters)='//  &
                              TRIM(DataLine(1)(1:100)))
        ErrorsFound=.true.
      ENDIF
      READ(DataLine(2)(6:),*,IOSTAT=ReadStat)  Rbsol
      IF (ReadStat /= 0) THEN
        CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Error in Read of RbSol values.')
        CALL ShowContinueError('Line (~'//TRIM(TrimSigDigits(FileLineCount+2))//') in error (first 100 characters)='//  &
                              TRIM(DataLine(2)(1:100)))
        ErrorsFound=.true.
      ELSEIF (ANY(Rbsol < 0.0) .or. ANY(RbSol > 1.0)) THEN
        CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Error in Read of RbSol values. (out of range [0,1])')
        CALL ShowContinueError('Line (~'//TRIM(TrimSigDigits(FileLineCount+2))//') in error (first 100 characters)='//  &
                              TRIM(DataLine(2)(1:100)))
        ErrorsFound=.true.
      ENDIF
      READ(DataLine(3)(6:),*,IOSTAT=ReadStat)  Tvis
      IF (ReadStat /= 0) THEN
        CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Error in Read of Tvis values.')
        CALL ShowContinueError('Line (~'//TRIM(TrimSigDigits(FileLineCount+3))//') in error (first 100 characters)='//  &
                              TRIM(DataLine(3)(1:100)))
        ErrorsFound=.true.
      ELSEIF (ANY(Tvis < 0.0) .or. ANY(Tvis > 1.0)) THEN
        CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Error in Read of Tvis values. (out of range [0,1])')
        CALL ShowContinueError('Line (~'//TRIM(TrimSigDigits(FileLineCount+3))//') in error (first 100 characters)='//  &
                              TRIM(DataLine(3)(1:100)))
        ErrorsFound=.true.
      ENDIF
      READ(DataLine(4)(6:),*,IOSTAT=ReadStat)  Rfvis
      IF (ReadStat /= 0) THEN
        CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Error in Read of Rfvis values.')
        CALL ShowContinueError('Line (~'//TRIM(TrimSigDigits(FileLineCount+4))//') in error (first 100 characters)='//  &
                              TRIM(DataLine(4)(1:100)))
        ErrorsFound=.true.
      ELSEIF (ANY(Rfvis < 0.0) .or. ANY(Rfvis > 1.0)) THEN
        CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Error in Read of Rfvis values. (out of range [0,1])')
        CALL ShowContinueError('Line (~'//TRIM(TrimSigDigits(FileLineCount+4))//') in error (first 100 characters)='//  &
                              TRIM(DataLine(4)(1:100)))
        ErrorsFound=.true.
      ENDIF
      READ(DataLine(5)(6:),*,IOSTAT=ReadStat)  Rbvis
      IF (ReadStat /= 0) THEN
        CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Error in Read of Rbvis values.')
        CALL ShowContinueError('Line (~'//TRIM(TrimSigDigits(FileLineCount+5))//') in error (first 100 characters)='//  &
                              TRIM(DataLine(5)(1:100)))
        ErrorsFound=.true.
      ELSEIF (ANY(Rbvis < 0.0) .or. ANY(Rbvis > 1.0)) THEN
        CALL ShowSevereError('HeatBalanceManager: SearchWindow5DataFile: Error in Read of Rbvis values. (out of range [0,1])')
        CALL ShowContinueError('Line (~'//TRIM(TrimSigDigits(FileLineCount+5))//') in error (first 100 characters)='//  &
                              TRIM(DataLine(5)(1:100)))
        ErrorsFound=.true.
      ENDIF
      FileLineCount=FileLineCount+5

      IF(ErrorsFound)   &
        CALL ShowFatalError('HeatBalanceManager: SearchWindow5DataFile: Construction='//TRIM(DesiredConstructionName)// &
             ' from the Window5 data file cannot be used because of above errors')

      ! Hemis
      Construct(ConstrNum)%TransDiff = Tsol(11)
      Construct(ConstrNum)%TransDiffVis = Tvis(11)
      Construct(ConstrNum)%ReflectSolDiffFront = Rfsol(11)
      Construct(ConstrNum)%ReflectSolDiffBack = Rbsol(11)
      Construct(ConstrNum)%ReflectVisDiffFront = Rfvis(11)
      Construct(ConstrNum)%ReflectVisDiffBack = Rbvis(11)

      CALL W5LsqFit(CosPhiIndepVar,Tsol,6,1,10,Construct(ConstrNum)%TransSolBeamCoef)
      CALL W5LsqFit(CosPhiIndepVar,Tvis,6,1,10,Construct(ConstrNum)%TransVisBeamCoef)
      CALL W5LsqFit(CosPhiIndepVar,Rfsol,6,1,10,Construct(ConstrNum)%ReflSolBeamFrontCoef)
      DO IGlass = 1,NGlass(IGlSys)
        CALL W5LsqFit(CosPhiIndepVar,AbsSol(IGlass,:),6,1,10,Construct(ConstrNum)%AbsBeamCoef(IGlass,:))
      END DO

      ! For comparing fitted vs. input distribution in incidence angle
      DO IPhi = 1,10
        Phi = REAL(IPhi-1,r64)*10.d0
        CosPhi = COS(Phi*DegToRadians)
        if (abs(CosPhi) < .0001d0) CosPhi=0.0
        tsolFit(IPhi) = POLYF(CosPhi,Construct(ConstrNum)%TransSolBeamCoef(1:6))
        tvisFit(IPhi) = POLYF(CosPhi,Construct(ConstrNum)%TransVisBeamCoef(1:6))
        rfsolFit(IPhi)= POLYF(CosPhi,Construct(ConstrNum)%ReflSolBeamFrontCoef(1:6))
        DO IGlass = 1,NGlass(IGlSys)
          solabsFit(IPhi,IGlass) = POLYF(CosPhi,Construct(ConstrNum)%AbsBeamCoef(IGlass,1:6))
        END DO
      END DO
      ! end

      ! NominalU of this construction (actually the total resistance of all of its layers; gas layer
      ! conductivity here ignores convective efffects in gap.)
      NominalU(ConstrNum) = 0.
      DO loop = 1,NGlass(IGlSys)+NGaps(IGlSys)
        MatNum = Construct(ConstrNum)%LayerPoint(loop)
        IF(Material(MatNum)%Group == WindowGlass) THEN
          NominalU(ConstrNum) = NominalU(ConstrNum) + Material(MatNum)%Thickness/Material(MatNum)%Conductivity
        ELSE IF(Material(MatNum)%Group == WindowGas .OR. Material(MatNum)%Group == WindowGasMixture) THEN
            ! If mixture, use conductivity of first gas in mixture
            NominalU(ConstrNum) = NominalU(ConstrNum) + Material(MatNum)%Thickness / &
                (Material(MatNum)%GasCon(1,1) + Material(MatNum)%GasCon(1,2)*300.0d0 + &
                 Material(MatNum)%GasCon(1,3)*90000.0d0)
        END IF
      END DO

    END DO  ! End of loop over glazing systems

    ! WindowFrameAndDivider objects

    TotFrameDividerPrev = TotFrameDivider
    DO IGlSys = 1,NGlSys
      IF(FrameWidth > 0.0 .OR. DividerWidth(IGlSys) > 0.0) THEN
        TotFrameDivider = TotFrameDivider + 1
        Construct(TotConstructs-NGlSys+IGlSys)%W5FrameDivider = TotFrameDivider
      END IF
    END DO

    IF(TotFrameDivider > TotFrameDividerPrev) THEN
      ALLOCATE(FrameDividerSave(TotFrameDividerPrev))
      DO loop = 1,TotFrameDividerPrev
        FrameDividerSave(loop) = FrameDivider(loop)
      END DO
      DEALLOCATE(FrameDivider)
      ALLOCATE(FrameDivider(TotFrameDivider))
      DO loop = 1,TotFrameDividerPrev
        FrameDivider(loop) = FrameDividerSave(loop)
      END DO
      DEALLOCATE(FrameDividerSave)
    END IF

    DO IGlSys = 1,NGlSys
      IF(FrameWidth > 0.0 .OR. DividerWidth(IGlSys) > 0.0) THEN
        FrDivNum = Construct(TotConstructs-NGlSys+IGlSys)%W5FrameDivider
        FrameDivider(FrDivNum)%FrameWidth                = FrameWidth
        FrameDivider(FrDivNum)%FrameProjectionOut        = FrameProjectionOut
        FrameDivider(FrDivNum)%FrameProjectionIn         = FrameProjectionIn
        FrameDivider(FrDivNum)%FrameConductance          = FrameConductance
        FrameDivider(FrDivNum)%FrEdgeToCenterGlCondRatio = FrEdgeToCenterGlCondRatio
        FrameDivider(FrDivNum)%FrameSolAbsorp            = FrameSolAbsorp
        FrameDivider(FrDivNum)%FrameVisAbsorp            = FrameVisAbsorp
        FrameDivider(FrDivNum)%FrameEmis                 = FrameEmis
        FrameDivider(FrDivNum)%FrameEdgeWidth            = 0.06355d0  ! 2.5 in
        IF (SameString(MullionOrientation,'Vertical')) THEN
          FrameDivider(FrDivNum)%MullionOrientation      = Vertical
        ELSEIF (SameString(MullionOrientation,'Horizontal')) THEN
          FrameDivider(FrDivNum)%MullionOrientation      = Horizontal
        ENDIF
        IF (SameString(DividerType(IGlSys),'DividedLite')) THEN
          FrameDivider(FrDivNum)%DividerType             = DividedLite
        ELSEIF (SameString(DividerType(IGlSys),'Suspended')) THEN
          FrameDivider(FrDivNum)%DividerType             = Suspended
        ENDIF
        FrameDivider(FrDivNum)%DividerWidth              = DividerWidth(IGlSys)
        FrameDivider(FrDivNum)%HorDividers               = HorDividers(IGlSys)
        FrameDivider(FrDivNum)%VertDividers              = VertDividers(IGlSys)
        FrameDivider(FrDivNum)%DividerProjectionOut      = DividerProjectionOut(IGlSys)
        FrameDivider(FrDivNum)%DividerProjectionIn       = DividerProjectionIn(IGlSys)
        FrameDivider(FrDivNum)%DividerConductance        = DividerConductance(IGlSys)
        FrameDivider(FrDivNum)%DivEdgeToCenterGlCondRatio = DivEdgeToCenterGlCondRatio(IGlSys)
        FrameDivider(FrDivNum)%DividerSolAbsorp          = DividerSolAbsorp(IGlSys)
        FrameDivider(FrDivNum)%DividerVisAbsorp          = DividerVisAbsorp(IGlSys)
        FrameDivider(FrDivNum)%DividerEmis               = DividerEmis(IGlSys)
        FrameDivider(FrDivNum)%DividerEdgeWidth          = 0.06355d0 ! 2.5 in
        IF(NGlSys == 1) THEN
          FrameDivider(FrDivNum)%Name = 'W5:'//TRIM(DesiredConstructionName)
        ELSE
          FrameDivider(FrDivNum)%Name = 'W5:'//TRIM(DesiredConstructionName)//':'//NumName(IGlSys)
        END IF
      END IF
    END DO

   IF(FrameWidth > 0.0 .AND. DividerWidth(1) > 0.0) THEN
     CALL DisplayString('--Construction and associated frame and divider found')
   ELSE IF(FrameWidth > 0.0) THEN
     CALL DisplayString('--Construction and associated frame found')
   ELSE IF(DividerWidth(1) > 0.0) THEN
     CALL DisplayString('--Construction and associated divider found')
   ELSE
     CALL DisplayString('--Construction without frame or divider found')
   END IF

 END IF

 CLOSE (W5DataFileNum)
 RETURN

 !999   CALL ShowFatalError('HeatBalanceManager: SearchWindow5DataFile: '//   &
 !       'Could not open Window5 Data File, expecting it as file name='//TRIM(DesiredFileName))  !RS: Secret Search String
999       WRITE(DebugFile,*) 'HeatBalanceManager: Search Window5DataFile: '// &
        'Could not open Window5 Data File, expecting it as file name='//TRIM(DesiredFileName)
 RETURN

 1000 EOFonFile = .TRUE.
      CLOSE (W5DataFileNum)
 RETURN

 END SUBROUTINE SearchWindow5DataFile

SUBROUTINE SetStormWindowControl

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Fred Winkelmann
          !       DATE WRITTEN   Jan 2004
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! Sets the storm window flag for each window, which is:
          !  -1: if storm window is not applicable (this will always be the value for interior
          !        windows since storm windows can only be applied to exterior windows
          !   0: if the window has a storm window but it is off
          !   1: if the window has a storm window and it is on

          ! A "storm window" is a single layer of exterior glass separated from the main window by air gap.
          ! Whether the storm window is in place is determined by the following values, which
          ! which are specified in the Storm Window object for the window:
          !  -Month that Storm Window Is Put On
          !  -Day of Month that Storm Window Is Put On
          !  -Month that Storm Window Is Taken Off
          !  -Day of Month that Storm Window Is Taken Off

          ! REFERENCES:na
          ! USE STATEMENTS:
  USE General, ONLY: BetweenDates

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE PARAMETER DEFINITIONS:na

          ! INTERFACE BLOCK SPECIFICATIONS:na

          ! DERIVED TYPE DEFINITIONS:na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:

  INTEGER            :: SurfNum         ! Surface number
  INTEGER            :: StormWinNum     ! Number of storm window object
  INTEGER            :: StormWinFlag    ! Storm window flag; this routine sets the following values:
                                        !   0: if the storm window is off this time step
                                        !   1: if the storm window is on this time step
  INTEGER            :: DateOff         ! Date Off for calculation

StormWinChangeThisDay = .false.

DO StormWinNum = 1,TotStormWin
  SurfNum = StormWindow(StormWinNum)%BaseWindowNum
  SurfaceWindow(SurfNum)%StormWinFlagPrevDay = SurfaceWindow(SurfNum)%StormWinFlag
  DateOff=StormWindow(StormWinNum)%DateOff-1
  ! Note: Dateon = Dateoff is not allowed and will have produced an error in getinput.
  IF (DateOff == 0) DateOff=366
  IF (BetweenDates(DayOfYear_Schedule,StormWindow(StormWinNum)%DateOn,DateOff)) THEN
    StormWinFlag=1
  ELSE
    StormWinFlag=0
  ENDIF
  SurfaceWindow(SurfNum)%StormWinFlag = StormWinFlag
  IF(BeginSimFlag) SurfaceWindow(SurfNum)%StormWinFlagPrevDay = StormWinFlag
  IF(SurfaceWindow(SurfNum)%StormWinFlag /= SurfaceWindow(SurfNum)%StormWinFlagPrevDay) StormWinChangeThisDay = .true.
END DO

RETURN
END SUBROUTINE SetStormWindowControl

SUBROUTINE CreateFCfactorConstructions(ConstrNum,ErrorsFound)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Tianzhen Hong
          !       DATE WRITTEN   July 2009
          !       MODIFIED
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! This subroutine goes through each construction defined with Ffactor or Cfactor method,
          ! and creates a construction (concrete + insulation) used in the heat transfer calculation.
          ! This subroutine only gets called once in the GetConstructionData subroutine

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE DataStringGlobals
  USE General, ONLY: RoundSigDigits

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  INTEGER, INTENT(INOUT) :: ConstrNum   ! Counter for Constructions
  LOGICAL, INTENT(INOUT) :: ErrorsFound ! If errors found in input

          ! SUBROUTINE PARAMETER DEFINITIONS:
   !Thermal resistance of the inside air film, m2.K/W. Average of 0.14 (heat flow up) and 0.11 (heat flow down)
  REAL(r64),PARAMETER :: Rfilm_in = 0.125d0
    !Thermal resistance of the outside air film used in calculating the Ffactor, m2.K/W. 0.17/5.678
  REAL(r64),PARAMETER  :: Rfilm_out = 0.03d0

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER :: ConstructNumAlpha ! Number of construction alpha names being passed
  INTEGER :: DummyNumProp      ! dummy variable for properties being passed
  INTEGER :: IOStat            ! IO Status when calling get input subroutine
  CHARACTER(len=MaxNameLength),DIMENSION(1) :: ConstructAlphas   ! Construction Alpha names defined
  REAL(r64), DIMENSION(4) :: DummyProps !Temporary array to transfer construction properties
  LOGICAL :: ErrorInName
  LOGICAL :: IsBlank
  INTEGER :: Loop

  INTEGER :: TotFfactorConstructs  ! Number of slabs-on-grade or underground floor constructions defined with F factors
  INTEGER :: TotCfactorConstructs  ! Number of underground wall constructions defined with C factors

  REAL(r64) :: Ffactor             !Ffactor in W/m-K, applies to deltaT of outside - indoor air temperature
  REAL(r64) :: Cfactor             !Cfactor in W/m2-K, does not include soil or air films
  REAL(r64) :: Area                !floor area in m2
  REAL(r64) :: PerimeterExposed    !perimeter exposed in m
  REAL(r64) :: Height              !Height of the underground wall in m

  REAL(r64) :: Reff                !Effective thermal resistance, m2.K/W
  REAL(r64) :: Rcon                !Concrete layer thermal resistance, m2.K/W
  REAL(r64) :: Rfic                !Thermal resistance of the fictitious material, m2.K/W
  INTEGER   :: MaterNum            !Material index
  REAL(r64) :: Rsoilequ            !Effective R-value of soil for underground walls
  INTEGER   :: iFCConcreteLayer    !Layer pointer to the materials array

  ! First get the concrete layer
  iFCConcreteLayer = FindIteminList('~FC_Concrete',Material%Name,TotMaterials)
  Rcon = Material(iFCConcreteLayer)%Resistance

  ! Count number of constructions defined with Ffactor or Cfactor method
  TotFfactorConstructs = GetNumObjectsFound('Construction:FfactorGroundFloor')
  TotCfactorConstructs = GetNumObjectsFound('Construction:CfactorUndergroundWall')

  ! First create ground floor constructions defined with F factor method if any
  CurrentModuleObject='Construction:FfactorGroundFloor'

  ! Loop through all constructs defined with Ffactor method
  DO Loop = 1, TotFfactorConstructs

      !Get the object names for each construction from the input processor
    CALL GetObjectItem(TRIM(CurrentModuleObject),Loop,ConstructAlphas,ConstructNumAlpha,DummyProps,DummyNumProp,IOSTAT,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

    ErrorInName=.false.
    IsBlank=.false.
    CALL VerifyName(ConstructAlphas(1),Construct%Name,ConstrNum,ErrorInName,IsBlank,TRIM(CurrentModuleObject)//' Name')
    IF (ErrorInName) THEN
      ErrorsFound=.true.
      CYCLE
    ENDIF

    ConstrNum = ConstrNum+1

    Construct(ConstrNum)%Name = ConstructAlphas(1)
    Construct(ConstrNum)%TypeIsFfactorFloor = .true.

    Ffactor = DummyProps(1)
    Area = DummyProps(2)
    PerimeterExposed = DummyProps(3)

    Construct(ConstrNum)%Area = Area
    Construct(ConstrNum)%PerimeterExposed = PerimeterExposed
    Construct(ConstrNum)%Ffactor = Ffactor

    IF (Ffactor <= 0 ) THEN
      CALL ShowSevereError(CurrentModuleObject // '="' //Trim(ConstructAlphas(1))// &
      '" has '//trim(cNumericFieldNames(1))//' <= 0.0, must be > 0.0.')
      CALL ShowContinueError('Entered value=['//trim(RoundSigDigits(Cfactor,2))//']')
      ErrorsFound=.true.
    ENDIF

    IF (Area <= 0 ) THEN
      CALL ShowSevereError(CurrentModuleObject // '="' //Trim(ConstructAlphas(1))// &
      '" has '//trim(cNumericFieldNames(2))//' <= 0.0, must be > 0.0.')
      CALL ShowContinueError('Entered value=['//trim(RoundSigDigits(Area,2))//']')
      ErrorsFound=.true.
    ENDIF

    IF (PerimeterExposed < 0 ) THEN
      CALL ShowSevereError(CurrentModuleObject // '="' //Trim(ConstructAlphas(1))// &
      '" has '//trim(cNumericFieldNames(3))//' <= 0.0, must be > 0.0.')
      CALL ShowContinueError('Entered value=['//trim(RoundSigDigits(PerimeterExposed,2))//']')
      ErrorsFound=.true.
    ENDIF

    ! The construction has two layers which have been created in GetMaterialData
    Construct(ConstrNum)%TotLayers = 2

    ! The concrete is the inside layer
    Construct(ConstrNum)%LayerPoint(2) = iFCConcreteLayer

    ! The fictitious insulation is the outside layer
    MaterNum = FindIteminList('~FC_Insulation_' // RoundSigDigits(Loop),Material%Name,TotMaterials)
    Construct(ConstrNum)%LayerPoint(1) = MaterNum

    ! Calculate the thermal resistance of the fictitious insulation layer
    IF (PerimeterExposed > 0.0) THEN
      Reff = Area / (PerimeterExposed * Ffactor)
    ELSE  ! PerimeterExposed = 0 for underground floor, assume R-1000 (IP)
      Reff = 177
    ENDIF

    Rfic = Reff - Rfilm_in - Rfilm_out - Rcon
    IF (Rfic <=0 ) THEN
      CALL ShowSevereError(CurrentModuleObject // '="' //Trim(ConstructAlphas(1))// &
      '" has calculated R value <= 0.0, must be > 0.0.')
      CALL ShowContinueError('Calculated value=['//trim(RoundSigDigits(Rfic,2))//'] Check definition.')
      ErrorsFound=.true.
    ENDIF

    Material(MaterNum)%Resistance = Rfic
    NominalR(MaterNum) = Rfic
    NominalU(ConstrNum) = NominalU(ConstrNum) + Reff
  END DO

  ! Then create underground wall constructions defined with C factor method if any
  CurrentModuleObject = 'Construction:CfactorUndergroundWall'

  DO Loop = 1, TotCfactorConstructs ! Loop through all constructs defined with Ffactor method

      !Get the object names for each construction from the input processor
    CALL GetObjectItem(TRIM(CurrentModuleObject),Loop,ConstructAlphas,ConstructNumAlpha,DummyProps,DummyNumProp,IOSTAT,  &
                   AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
                   AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

    ErrorInName=.false.
    IsBlank=.false.
    CALL VerifyName(ConstructAlphas(1),Construct%Name,ConstrNum,ErrorInName,IsBlank,TRIM(CurrentModuleObject)//' Name')
    IF (ErrorInName) THEN
      ErrorsFound=.true.
      CYCLE
    ENDIF

    ConstrNum = ConstrNum+1

    Construct(ConstrNum)%Name = ConstructAlphas(1)
    Construct(ConstrNum)%TypeIsCfactorWall = .true.

    Cfactor = DummyProps(1)
    Height = DummyProps(2)

    Construct(ConstrNum)%Height = Height
    Construct(ConstrNum)%Cfactor = Cfactor

    IF (Cfactor <= 0 ) THEN
      CALL ShowSevereError(CurrentModuleObject // ' ' //Trim(ConstructAlphas(1))// &
      ' has '//trim(cNumericFieldNames(1))//' <= 0.0, must be > 0.0.')
      CALL ShowContinueError('Entered value=['//trim(RoundSigDigits(Cfactor,2))//']')
      ErrorsFound=.true.
    ENDIF

    IF (Height <= 0 ) THEN
      CALL ShowSevereError(CurrentModuleObject // ' ' //Trim(ConstructAlphas(1))// &
      ' has '//trim(cNumericFieldNames(2))//' <= 0.0, must be > 0.0.')
      CALL ShowContinueError('Entered value=['//trim(RoundSigDigits(Height,2))//']')
      ErrorsFound=.true.
    ENDIF

    ! The construction has two layers which have been created in GetMaterialData
    Construct(ConstrNum)%TotLayers = 2

    ! The concrete is the inside layer
    Construct(ConstrNum)%LayerPoint(2) = iFCConcreteLayer

    ! The fictitious insulation is the outside layer
    MaterNum = FindIteminList('~FC_Insulation_' // RoundSigDigits(Loop+TotFfactorConstructs,0),Material%Name,TotMaterials)
    Construct(ConstrNum)%LayerPoint(1) = MaterNum

    IF (Height <= 0.305d0) THEN      ! 1 ft
      Rsoilequ = 4.88d0
    ELSEIF (Height >= 3.048d0) THEN  ! 10 ft
      Rsoilequ = 34.64d0
    ELSE  ! regression from ASHRAE 90.1-2007 TABLE C6.10.1 Effective R-Value of Soil, R2 = 0.9972
      Rsoilequ = 2.592d0 + 10.736d0 * Height
    ENDIF

    Reff = 1.0d0/Cfactor + Rsoilequ    ! Cfactor does not include air films

    Rfic = Reff - Rcon
    IF (Rfic <=0 ) THEN
      CALL ShowSevereError(CurrentModuleObject // '="' //Trim(ConstructAlphas(1))// &
      '" has calculated R value <= 0.0, must be > 0.0.')
      CALL ShowContinueError('Calculated value=['//trim(RoundSigDigits(Rfic,2))//'] Check definition.')
      ErrorsFound = .true.
    ENDIF

    Material(MaterNum)%Resistance = Rfic
    NominalR(MaterNum) = Rfic
    NominalU(ConstrNum) = NominalU(ConstrNum) + Reff
  END DO

  RETURN

END SUBROUTINE CreateFCfactorConstructions


SUBROUTINE CreateTCConstructions(ErrorsFound)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Tianzhen Hong
          !       DATE WRITTEN   January 2009
          !       MODIFIED
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! This subroutine goes through each TC master construction and creates a complete series
          ! of the slave thermochromic constructions.
          ! This subroutine only gets called once in the GetHeatBalanceInput subroutine
          !  after materials, constructions and building geometry data are read.

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:

  USE DataStringGlobals
  USE General, ONLY: RoundSigDigits

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  LOGICAL, INTENT(INOUT) :: ErrorsFound ! If errors found in input

          ! SUBROUTINE PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:

  INTEGER :: Loop
  INTEGER :: iTC = 0
  INTEGER :: iMat = 0
  INTEGER :: NumNewConst = 0
  INTEGER :: iTCG = 0

  TYPE (ConstructionData),  ALLOCATABLE, DIMENSION(:) :: ConstructExtra

  ! create a complete series of slave TC constructions for each TCGlazings object
  !   with the master construction as the basis
  ALLOCATE (ConstructExtra(500))   ! temporary assume a maximum of 500 extra TC constructions

  NumNewConst = 0
  DO Loop = 1, TotConstructs
    IF (Construct(Loop)%TCFlag == 1) THEN
      iTCG = Material(Construct(Loop)%TCLayer)%TCParent
      iMat = TCGlazings(iTCG)%NumGlzMat
      DO iTC = 1, iMat
        ! create a new TC construction
        IF (NumNewConst >=500 ) THEN
          CALL ShowFatalError('Errors found in creating thermochromic window constructions: '//  &
             'There are more than 500 thermochromic glazings.')
          ErrorsFound=.true.
        ENDIF

        NumNewConst = NumNewConst +1
        ConstructExtra(NumNewConst) = Construct(Loop)  ! copy data
        ConstructExtra(NumNewConst)%Name = TRIM(Construct(Loop)%Name) // '_TC_'//RoundSigDigits(TCGlazings(iTCG)%SpecTemp(iTC),0)
        ConstructExtra(NumNewConst)%TCLayer = TCGlazings(iTCG)%LayerPoint(iTC)
        ConstructExtra(NumNewConst)%LayerPoint(Construct(Loop)%TCLayerID) = ConstructExtra(NumNewConst)%TCLayer
        ConstructExtra(NumNewConst)%TCFlag = 1
        ConstructExtra(NumNewConst)%TCMasterConst = Loop
        ConstructExtra(NumNewConst)%TCLayerID = Construct(Loop)%TCLayerID
        ConstructExtra(NumNewConst)%TCGlassID = Construct(Loop)%TCGlassID
        ConstructExtra(NumNewConst)%TypeIsWindow = .True.
      ENDDO
    ENDIF
  ENDDO

  ! Increase Construct() and copy the extra constructions
  ALLOCATE(ConstructSave(TotConstructs))
  ALLOCATE(NominalUSave(TotConstructs))
  DO Loop = 1,TotConstructs
    ConstructSave(Loop) = Construct(Loop)
    NominalUSave(Loop) = NominalU(Loop)
  END DO
  DEALLOCATE(Construct)
  DEALLOCATE(NominalU)
  TotConstructs = TotConstructs + NumNewConst
  ALLOCATE(Construct(TotConstructs))
  ALLOCATE(NominalU(TotConstructs))
  DO Loop = 1, TotConstructs - NumNewConst
    Construct(Loop) = ConstructSave(Loop)
    NominalU(Loop) = NominalUSave(Loop)
  END DO
  DEALLOCATE(ConstructSave)
  DEALLOCATE(NominalUSave)

  DO Loop =1, NumNewConst
    Construct(TotConstructs - NumNewConst + Loop) = ConstructExtra(Loop)
  END DO

  DEALLOCATE(ConstructExtra)

  RETURN

END SUBROUTINE CreateTCConstructions

SUBROUTINE SetupSimpleWindowGlazingSystem(MaterNum)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         B. Griffith
          !       DATE WRITTEN   January 2009
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! Convert simple window performance indices into all the properties needed to
          ! describe a single, equivalent glass layer

          ! METHODOLOGY EMPLOYED:
          ! The simple window indices are converted to a single materal layer using a "block model"
          !

          ! REFERENCES:
          ! draft paper by Arasteh, Kohler, and Griffith

          ! USE STATEMENTS:
          ! na

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  INTEGER  :: MaterNum

          ! SUBROUTINE PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  REAL(r64) :: Riw = 0.0D0 ! thermal resistance of interior film coefficient under winter conditions (m2-K/W)
  REAL(r64) :: Row = 0.0D0  ! theraml resistance of exterior film coefficient under winter conditions (m2-K/W)
  REAL(r64) :: Rlw = 0.0D0  ! thermal resistance of block model layer (m2-K/W)
  REAL(r64) :: Ris = 0.0D0  ! thermal resistance of interior film coefficient under summer conditions (m2-K/W)
  REAL(r64) :: Ros = 0.0D0  ! theraml resistance of exterior film coefficient under summer conditions (m2-K/W)
  REAL(r64) :: InflowFraction = 0.0d0  ! inward flowing fraction for SHGC, intermediate value non dimensional
  REAL(r64) :: SolarAbsorb = 0.0d0  ! solar aborptance
  LOGICAL   :: ErrorsFound = .false.
  REAL(r64) :: TsolLowSide = 0.0d0   ! intermediate solar transmission for interpolating
  REAL(r64) :: TsolHiSide = 0.0d0   ! intermediate solar transmission for interpolating
  REAL(r64) :: DeltaSHGCandTsol = 0.0d0  ! intermediate difference
  REAL(r64) :: RLowSide   = 0.0d0
  REAL(r64) :: RHiSide    = 0.0d0


  ! first fill out defaults
  Material(MaterNum)%GlassSpectralDataPtr = 0
  Material(MaterNum)%SolarDiffusing       = .FALSE.
  Material(MaterNum)%Roughness            = VerySmooth
  Material(MaterNum)%TransThermal         = 0.0d0
  Material(MaterNum)%AbsorpThermalBack    = 0.84d0
  Material(MaterNum)%AbsorpThermalFront   = 0.84d0
  Material(MaterNum)%AbsorpThermal        = Material(MaterNum)%AbsorpThermalBack

  ! step 1. Determine U-factor without film coefficients
  ! Simple window model has its own correlation for film coefficients (m2-K/W) under Winter conditions as function of U-factor
  IF ( Material(MaterNum)%SimpleWindowUfactor < 5.85d0 ) THEN
    Riw = 1.0d0 / ( 0.359073d0 * Log( Material(MaterNum)%SimpleWindowUfactor ) + 6.949915d0)
  Else
    Riw = 1.0d0 / (1.788041d0 *  Material(MaterNum)%SimpleWindowUfactor - 2.886625d0 )
  Endif
  Row = 1.0d0 / (0.025342d0 * Material(MaterNum)%SimpleWindowUfactor + 29.163853d0 )

  ! determine 1/U without film coefficients
  Rlw = (1.0d0/Material(MaterNum)%SimpleWindowUfactor) - Riw - Row
  IF (Rlw <= 0.0d0) THEN ! U factor of film coefficients is better than user input.
    Rlw = MAX(Rlw, 0.001d0)
    CALL ShowWarningError('WindowMaterial:SimpleGlazingSystem: ' //Trim(Material(MaterNum)%Name)// &
          ' has U-factor higher than that provided by surface film resistances, Check value of U-factor')
  ENDIF

  ! Step 2. determine layer thickness.

  IF ( (1.0d0 / Rlw) > 7.0d0 ) THEN
    Material(MaterNum)%Thickness = 0.002d0
  ELSE
    Material(MaterNum)%Thickness = 0.05914d0 - (0.00714d0 / Rlw )
  ENDIF

  ! Step 3. determine effective conductivity

  Material(MaterNum)%Conductivity = Material(MaterNum)%Thickness / Rlw
  IF (Material(MaterNum)%Conductivity > 0.0) THEN
    NominalR(MaterNum) = Rlw
  ELSE
    ErrorsFound = .true.
    CALL ShowSevereError('WindowMaterial:SimpleGlazingSystem: ' //Trim(Material(MaterNum)%Name)// &
          ' has Conductivity <= 0.0, must be >0.0, Check value of U-factor')
  END IF

  !step 4. determine solar transmission (revised to 10-1-2009 version from LBNL.)

  IF (Material(MaterNum)%SimpleWindowUfactor > 4.5d0) THEN

    IF (Material(MaterNum)%SimpleWindowSHGC < 0.7206d0 ) THEN

      Material(MaterNum)%Trans = 0.939998d0  * Material(MaterNum)%SimpleWindowSHGC**2 &
                                 + 0.20332d0 * Material(MaterNum)%SimpleWindowSHGC
    ELSE ! >= 0.7206

      Material(MaterNum)%Trans = 1.30415d0 * Material(MaterNum)%SimpleWindowSHGC - 0.30515d0

    ENDIF

  ELSEIF (Material(MaterNum)%SimpleWindowUfactor < 3.4d0) THEN

    IF (Material(MaterNum)%SimpleWindowSHGC <= 0.15d0) THEN
      Material(MaterNum)%Trans = 0.41040d0 * Material(MaterNum)%SimpleWindowSHGC
    ELSE ! > 0.15
      Material(MaterNum)%Trans =  0.085775d0*(Material(MaterNum)%SimpleWindowSHGC**2) &
                                  + 0.963954d0*Material(MaterNum)%SimpleWindowSHGC - 0.084958d0
    ENDIF
  ELSE ! interpolate. 3.4 <= Ufactor <= 4.5

    IF (Material(MaterNum)%SimpleWindowSHGC < 0.7206d0 ) THEN
      TsolHiSide = 0.939998d0    * Material(MaterNum)%SimpleWindowSHGC**2 &
                                 + 0.20332d0 * Material(MaterNum)%SimpleWindowSHGC
    ELSE ! >= 0.7206
      TsolHiSide = 1.30415d0 * Material(MaterNum)%SimpleWindowSHGC - 0.30515d0
    ENDIF

    IF (Material(MaterNum)%SimpleWindowSHGC <= 0.15d0) THEN
      TsolLowSide = 0.41040d0 * Material(MaterNum)%SimpleWindowSHGC
    ELSE ! > 0.15
      TsolLowSide =  0.085775d0*(Material(MaterNum)%SimpleWindowSHGC**2) &
                                  + 0.963954d0*Material(MaterNum)%SimpleWindowSHGC - 0.084958d0
    ENDIF

    Material(MaterNum)%Trans = ((Material(MaterNum)%SimpleWindowUfactor - 3.4d0) &
                                  / (4.5d0 - 3.4d0) )   &
                                  * (TsolHiSide - TsolLowSide)  + TsolLowSide

  ENDIF
  If (Material(MaterNum)%Trans < 0.0d0) Material(MaterNum)%Trans = 0.0d0

  !step 5.  determine solar reflectances

  DeltaSHGCandTsol = Material(MaterNum)%SimpleWindowSHGC - Material(MaterNum)%Trans

  IF (Material(MaterNum)%SimpleWindowUfactor > 4.5d0) THEN

   Ris = 1.0d0 / (29.436546d0*DeltaSHGCandTsol**3.0d0 - 21.943415d0*DeltaSHGCandTsol**2 &
         + 9.945872d0*DeltaSHGCandTsol + 7.426151d0 )
   Ros = 1.0d0 / (2.225824d0*DeltaSHGCandTsol + 20.577080d0 )
  ELSEIF (Material(MaterNum)%SimpleWindowUfactor < 3.4d0) THEN

   Ris = 1.0d0 / (199.8208128d0*DeltaSHGCandTsol**3.0d0 - 90.639733d0*DeltaSHGCandTsol**2 &
         + 19.737055d0*DeltaSHGCandTsol + 6.766575d0 )
   Ros = 1.0d0 / (5.763355d0*DeltaSHGCandTsol + 20.541528d0 )
  ELSE ! interpolate. 3.4 <= Ufactor <= 4.5
   !inside first
   RLowSide = 1.0d0 / (199.8208128d0*DeltaSHGCandTsol**3.0d0 - 90.639733d0*DeltaSHGCandTsol**2 &
         + 19.737055d0*DeltaSHGCandTsol + 6.766575d0 )
   RHiSide  = 1.0d0 / (29.436546d0*DeltaSHGCandTsol**3 - 21.943415d0*DeltaSHGCandTsol**2 &
         + 9.945872d0*DeltaSHGCandTsol + 7.426151d0 )
   Ris = ((Material(MaterNum)%SimpleWindowUfactor - 3.4d0) &
                                  / (4.5d0 - 3.4d0) )   &
                                  * (RLowSide - RHiSide)  + RLowSide
   ! then outside
   RLowSide = 1.0d0 / (5.763355d0*DeltaSHGCandTsol + 20.541528d0 )
   RHiSide  = 1.0d0 / (2.225824d0*DeltaSHGCandTsol + 20.577080d0 )
   Ros = ((Material(MaterNum)%SimpleWindowUfactor - 3.4d0) &
                                  / (4.5d0 - 3.4d0) )   &
                                  * (RLowSide - RHiSide)  + RLowSide

  ENDIF

  InflowFraction = (Ros + 0.5d0*Rlw)/(Ros + Rlw + Ris)

  SolarAbsorb = (Material(MaterNum)%SimpleWindowSHGC - Material(MaterNum)%Trans) / InflowFraction
  Material(MaterNum)%ReflectSolBeamBack  = 1.0 - Material(MaterNum)%Trans - SolarAbsorb
  Material(MaterNum)%ReflectSolBeamFront = Material(MaterNum)%ReflectSolBeamBack

  !step 6. determine visible properties.
  IF (Material(MaterNum)%SimpleWindowVTinputByUser) THEN
    Material(MaterNum)%TransVis = Material(MaterNum)%SimpleWindowVisTran
    Material(MaterNum)%ReflectVisBeamBack  = - 0.7409d0 * Material(MaterNum)%TransVis**3  &
                                             + 1.6531d0 * Material(MaterNum)%TransVis**2  &
                                             - 1.2299d0 * Material(MaterNum)%TransVis + 0.4545d0
    IF (Material(MaterNum)%TransVis + Material(MaterNum)%ReflectVisBeamBack >= 1.0d0) THEN
      Material(MaterNum)%ReflectVisBeamBack = 0.999d0 - Material(MaterNum)%TransVis
    ENDIF

    Material(MaterNum)%ReflectVisBeamFront = - 0.0622d0 * Material(MaterNum)%TransVis**3  &
                                             + 0.4277d0 * Material(MaterNum)%TransVis**2  &
                                             - 0.4169d0 * Material(MaterNum)%TransVis + 0.2399d0
    IF (Material(MaterNum)%TransVis + Material(MaterNum)%ReflectVisBeamFront >= 1.0d0) THEN
      Material(MaterNum)%ReflectVisBeamFront = 0.999d0 - Material(MaterNum)%TransVis
    ENDIF
  ELSE
    Material(MaterNum)%TransVis = Material(MaterNum)%Trans
    Material(MaterNum)%ReflectVisBeamBack  = Material(MaterNum)%ReflectSolBeamBack
    Material(MaterNum)%ReflectVisBeamFront = Material(MaterNum)%ReflectSolBeamFront
  ENDIF

  !step 7. The dependence on incident angle is in subroutine TransAndReflAtPhi

  !step 8.  Hemispherical terms are averaged using standard method


  IF (ErrorsFound) THEN
    CALL ShowFatalError('Program halted because of input problem(s) in WindowMaterial:SimpleGlazingSystem')
  ENDIF

  RETURN

END SUBROUTINE SetupSimpleWindowGlazingSystem


! *****************************************************************************

!     NOTICE
!
!     Copyright � 1996-2012 The Board of Trustees of the University of Illinois
!     and The Regents of the University of California through Ernest Orlando Lawrence
!     Berkeley National Laboratory.  All rights reserved.
!
!     Portions of the EnergyPlus software package have been developed and copyrighted
!     by other individuals, companies and institutions.  These portions have been
!     incorporated into the EnergyPlus software package under license.   For a complete
!     list of contributors, see "Notice" located in EnergyPlus.f90.
!
!     NOTICE: The U.S. Government is granted for itself and others acting on its
!     behalf a paid-up, nonexclusive, irrevocable, worldwide license in this data to
!     reproduce, prepare derivative works, and perform publicly and display publicly.
!     Beginning five (5) years after permission to assert copyright is granted,
!     subject to two possible five year renewals, the U.S. Government is granted for
!     itself and others acting on its behalf a paid-up, non-exclusive, irrevocable
!     worldwide license in this data to reproduce, prepare derivative works,
!     distribute copies to the public, perform publicly and display publicly, and to
!     permit others to do so.
!
!     TRADEMARKS: EnergyPlus is a trademark of the US Department of Energy.
!

END MODULE HeatBalanceManager
