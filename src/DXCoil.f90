MODULE DXCoils

  ! Module containing the DX coil simulation routines

  ! MODULE INFORMATION:
  !       AUTHOR         Fred Buhl
  !       DATE WRITTEN   May 2000
  !       MODIFIED       Don Shirey, Aug/Sept 2000, Feb/Oct 2001, Sept 2003, Jan 2004
  !                      Feb 2005 M. J. Witte, GARD Analytics, Inc.
  !                        Add new coil type COIL:DX:MultiMode:CoolingEmpirical:
  !                        Work supported by ASHRAE research project 1254-RP
  !                      Aug 2006 B Griffith, NREL
  !                        Added water system interactions for new water manager,
  !                      Feb 2010 B Nigusse, FSEC
  !                        Added Standard Rating for Coil:Cooling:DX:SingleSpeed
  !                      April 2010 Chandan Sharma, FSEC
  !                        Added basin heater routines for Coil:Cooling:DX:SingleSpeed,
  !                        Coil:Cooling:DX:TwoSpeed, Coil:Cooling:DX:MultiSpeed
  !                        and Coil:Cooling:DX:TwoStageWithHumidityControlMode
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS MODULE:
  ! To encapsulate the data and algorithms required to simulate DX cooling coils in
  ! EnergyPlus. Module currently models air-cooled or evap-cooled direct expansion systems
  ! (split or packaged). Air-side performance is modeled to determine coil discharge
  ! air conditions. The module also determines the DX unit's electrical energy usage.
  ! Neither the air-side performance nor the electrical energy usage includes the effect
  ! of supply air fan heat/energy usage. The supply air fan is modeled by other modules.

  ! METHODOLOGY EMPLOYED:
  !

  ! REFERENCES:


  ! OTHER NOTES:
  !

  ! USE STATEMENTS:
  ! Use statements for data only modules
USE DataPrecisionGlobals
USE DataLoopNode
USE DataGlobals
USE DataHVACGlobals
USE Psychrometrics
Use DataEnvironment, ONLY: StdBaroPress, EnvironmentName, CurMnDy, OutDryBulbTemp, OutHumRat, OutBaroPress, OutWetBulbTemp
USE DataHeatBalance, ONLY: HeatReclaimDXCoil
USE DataInterfaces

  ! Use statements for access to subroutines in other modules
USE ScheduleManager

IMPLICIT NONE         ! Enforce explicit typing of all variables

PRIVATE ! Everything private unless explicitly made public

  !MODULE PARAMETER DEFINITIONS
! Defrost strategy (heat pump only)
INTEGER, PARAMETER :: ReverseCycle     = 1 ! uses reverse cycle defrost strategy
INTEGER, PARAMETER :: Resistive        = 2 ! uses electric resistance heater for defrost
! Defrost control  (heat pump only)
INTEGER, PARAMETER :: Timed            = 1 ! defrost cycle is timed
INTEGER, PARAMETER :: OnDemand         = 2 ! defrost cycle occurs only when required
! Compressor operation
INTEGER, PARAMETER :: On               = 1 ! normal compressor operation
INTEGER, PARAMETER :: Off              = 0 ! signal DXCoil that compressor shouldn't run
! Condenser Type (using same numbering scheme as for chillers)
INTEGER, PARAMETER :: AirCooled        = 1 ! Air-cooled condenser
!INTEGER, PARAMETER :: WaterCooled      = 2 ! Water-cooled condenser (not yet implemented)
INTEGER, PARAMETER :: EvapCooled       = 3 ! Evaporatively-cooled condenser
INTEGER, PARAMETER :: WaterHeater      = 4 ! Condenser heats water (e.g., in water heater tank)

REAL(r64), PARAMETER ::    RatedInletAirTemp   = 26.6667d0   ! 26.6667C or 80F
REAL(r64), PARAMETER ::    RatedInletWetbulbTemp = 19.44d0   ! 19.44 or 67F
REAL(r64), PARAMETER ::    RatedInletAirHumRat = 0.01125d0   ! Humidity ratio corresponding to 80F dry bulb/67F wet bulb
REAL(r64), PARAMETER ::    RatedOutdoorAirTemp = 35.0d0      ! 35 C or 95F
REAL(r64), PARAMETER ::    RatedInletAirTempHeat = 21.11d0   ! 21.11C or 70F
REAL(r64), PARAMETER ::    RatedOutdoorAirTempHeat = 8.33d0  ! 8.33 C or 47F
REAL(r64), PARAMETER ::    RatedInletWetbulbTempHeat = 15.55d0 ! 15.55 or 60F

! Airflow per total capacity range
REAL(r64), PARAMETER ::    MaxRatedVolFlowPerRatedTotCap = 0.00006041d0 ! m3/s per watt = 450 cfm/ton
REAL(r64), PARAMETER ::    MinRatedVolFlowPerRatedTotCap = 0.00004027d0 ! m3/s per watt = 300 cfm/ton
REAL(r64), PARAMETER ::    MaxHeatVolFlowPerRatedTotCap  = 0.00008056d0 ! m3/s per watt = 600 cfm/ton
REAL(r64), PARAMETER ::    MaxCoolVolFlowPerRatedTotCap  = 0.00006713d0 ! m3/s per watt = 500 cfm/ton
REAL(r64), PARAMETER ::    MinOperVolFlowPerRatedTotCap  = 0.00002684d0 ! m3/s per watt = 200 cfm/ton

! Curve Types
INTEGER, PARAMETER :: Linear      = 1
INTEGER, PARAMETER :: Bilinear    = 2
INTEGER, PARAMETER :: Quadratic   = 3
INTEGER, PARAMETER :: Biquadratic = 4
INTEGER, PARAMETER :: Cubic       = 5

! Multimode DX Coil
INTEGER, PARAMETER :: MaxCapacityStages = 2  ! Maximum number of capacity stages supported
INTEGER, PARAMETER :: MaxDehumidModes   = 1  ! Maximum number of enhanced dehumidification modes supported
INTEGER, PARAMETER :: MaxModes=MaxCapacityStages*(MaxDehumidModes+1)  ! Maximum number of performance modes

!Water Systems
INTEGER, PARAMETER :: CondensateDiscarded = 1001 ! default mode where water is "lost"
INTEGER, PARAMETER :: CondensateToTank    = 1002 ! collect coil condensate from air and store in water storage tank

INTEGER, PARAMETER :: WaterSupplyFromMains = 101
INTEGER, PARAMETER :: WaterSupplyFromTank  = 102

INTEGER, PARAMETER :: NumValidFuelTypes=7
CHARACTER(len=*), PARAMETER, DIMENSION(NumValidFuelTypes) :: cValidFuelTypes=    &
                 (/'Electricity',  &
                   'NaturalGas ',  &
                   'PropaneGas ',  &
                   'Diesel     ',  &
                   'Gasoline   ',  &
                   'FuelOil#1  ',  &
                   'FuelOil#2  '/)

! Fuel Types
INTEGER, PARAMETER :: FuelTypeElectricity = 1     ! Fuel type for electricity
INTEGER, PARAMETER :: FuelTypeNaturalGas  = 2     ! Fuel type for natural gas
INTEGER, PARAMETER :: FuelTypePropaneGas  = 3     ! Fuel type for propane gas
INTEGER, PARAMETER :: FuelTypeDiesel      = 4     ! Fuel type for diesel
INTEGER, PARAMETER :: FuelTypeGasoline    = 5     ! Fuel type for gasoline
INTEGER, PARAMETER :: FuelTypeFuelOil1    = 6     ! Fuel type for fuel oil #1
INTEGER, PARAMETER :: FuelTypeFuelOil2    = 7     ! Fuel type for fuel oil #2

  ! DERIVED TYPE DEFINITIONS
TYPE, PUBLIC :: DXCoilData
!          Some variables in this type are arrays (dimension=MaxModes) to support coil type
!          COIL:DX:MultiMode:CoolingEmpirical.  Other coil types only use the first element.
  CHARACTER(len=MaxNameLength) :: Name           =' '    ! Name of the DX Coil
  CHARACTER(len=MaxNameLength) :: DXCoilType     =' '    ! type of coil
  INTEGER                      :: DXCoilType_Num = 0     ! Integer equivalent to DXCoilType
  CHARACTER(len=MaxNameLength) :: Schedule       =' '    ! WaterCoil Operation Schedule
  INTEGER :: SchedPtr               = 0 ! Pointer to the correct schedule
!          RatedCoolCap, RatedSHR and RatedCOP do not include the thermal or electrical
!          effects due to the supply air fan
  REAL(r64) :: RatedTotCap(MaxModes)  =0.0 ! Gross total cooling capacity at rated conditions [watts]
  LOGICAL   :: RatedTotCapEMSOverrideOn(MaxModes) = .FALSE.  !if true, then EMS is calling to override rated total capacity
  REAL(r64) :: RatedTotCapEMSOverrideValue(MaxModes) = 0.0  ! value to use for EMS override

  REAL(r64) :: RatedSHR(MaxModes)     =0.0 ! Sensible heat ratio (sens cap/total cap) at rated conditions
  LOGICAL   :: RatedSHREMSOverrideOn(MaxModes)     =.false. ! if true, then EMS is calling to override Sensible heat ratio
  REAL(r64) :: RatedSHREMSOverrideValue(MaxModes)     =0.0 ! value to use for EMS override forSensible heat ratio
  REAL(r64) :: RatedCOP(MaxModes)     =0.0 ! Coefficient of performance at rated conditions
  REAL(r64) :: RatedAirVolFlowRate(MaxModes) =0.0  ! Air volume flow rate through coil at rated conditions [m3/s]
                                                   ! This is adjusted for bypassed air if any (see BypassedFlowFrac)
  LOGICAL   :: RatedAirVolFlowRateEMSOverrideON(MaxModes) =.false.  ! if true, then EMS is calling to override Air volume flow rate
  REAL(r64) :: RatedAirVolFlowRateEMSOverrideValue(MaxModes) =0.0  ! value to use for EMS override Air volume flow rate
  REAL(r64) :: FanPowerPerEvapAirFlowRate(MaxModes)=0.0  ! Fan Power Per Air volume flow rate through the
                                                         ! Evaporator coil at rated conditions [W/(m3/s)]
  REAL(r64) :: RatedAirMassFlowRate(MaxModes) =0.0 ! Air mass flow rate through coil at rated conditions [kg/s]
                                                 ! This is adjusted for bypassed air if any (see BypassedFlowFrac)
  REAL(r64) :: BypassedFlowFrac(MaxModes) =0.0     ! Fraction of air flow bypassed around coil
  REAL(r64) :: RatedCBF(MaxModes)     =0.0 ! rated coil bypass factor, determined using RatedTotCap and RatedSHR
  INTEGER :: AirInNode              = 0  ! Air inlet node number
  INTEGER :: AirOutNode             = 0  ! Air outlet node number
  INTEGER :: CCapFTemp(MaxModes)    = 0  ! index of total cooling capacity modifier curve
                                         ! (function of entering wetbulb, outside drybulb)
  INTEGER :: CCapFTempErrorIndex    = 0  ! Used for warning messages when output of CCapFTemp is negative
  INTEGER :: TotCapTempModFacCurveType(MaxModes) = 0 !type of curve for CCapFTemp (cubic,quadratic,bi-quadratic)
  INTEGER :: CCapFFlow(MaxModes)    = 0  ! index of total cooling capacity modifier curve
                                         ! (function of actual supply air flow vs rated air flow)
  INTEGER :: CCapFFlowErrorIndex    = 0  ! Used for warning messages when output of CCapFFlow is negative
  INTEGER :: EIRFTemp(MaxModes)     = 0  ! index of energy input ratio modifier curve
                                         ! (function of entering wetbulb, outside drybulb)
  INTEGER :: EIRFTempErrorIndex     = 0  ! Used for warning messages when output of EIRFTemp is negative
  INTEGER :: EIRTempModFacCurveType(MaxModes) = 0  !type of curve for EIRFTemp (cubic,quadratic,bi-quadratic)
  INTEGER :: EIRFFlow(MaxModes)     = 0  ! index of energy input ratio modifier curve
                                         ! (function of actual supply air flow vs rated air flow)
  INTEGER :: EIRFFlowErrorIndex     = 0  ! Used for warning messages when output of EIRFFlow is negative
  INTEGER :: PLFFPLR(MaxModes)      = 0  ! index of part-load factor vs part-load ratio curve
  LOGICAL :: ReportCoolingCoilCrankcasePower =.true. ! logical determines if the cooling coil crankcase heater power is reported
  REAL(r64) :: CrankcaseHeaterCapacity =0.0 ! total crankcase heater capacity [W]
  REAL(r64) :: CrankcaseHeaterPower    =0.0 ! report variable for average crankcase heater power [W]
  REAL(r64) :: MaxOATCrankcaseHeater   =0.0 ! maximum OAT for crankcase heater operation [C]
  REAL(r64) :: CrankcaseHeaterConsumption  = 0.0 ! report variable for total crankcase heater energy consumption [J]
  REAL(r64) :: BasinHeaterPowerFTempDiff   = 0.0 ! Basin heater capacity per degree C below setpoint (W/C)
  REAL(r64) :: BasinHeaterSetPointTemp     = 0.0 ! setpoint temperature for basin heater operation (C)
  INTEGER :: CompanionUpstreamDXCoil = 0  ! index number of the DX coil that is "upstream" of this DX coil. Currently used for
                                          ! UnitarySystem:HeatPump:AirToAir for proper calculation of crankcase heater energy
                                          ! consumption
  LOGICAL :: FindCompanionUpStreamCoil = .TRUE. ! Flag to get the companion coil in Init.
  INTEGER :: CondenserInletNodeNum(MaxModes) = 0  ! Node number of outdoor condenser(s) (actually an evaporator for heating coils)
  INTEGER :: LowOutletTempIndex     =0   ! used for low outlet temperature warnings
  REAL(r64) :: FullLoadOutAirTempLast =0.0 ! used for low outlet temperature warnings
  REAL(r64) :: FullLoadInletAirTempLast =0.0 ! used for low outlet temperature warnings
  LOGICAL :: PrintLowOutTempMessage= .FALSE.  ! used to print warning message for low outlet air dry-bulb conditions
  CHARACTER(len=300) :: LowOutTempBuffer1=' ' ! holds warning message until next iteration (only prints 1 message/iteration)
  CHARACTER(len=300) :: LowOutTempBuffer2=' ' ! holds warning message until next iteration (only prints 1 message/iteration)
  INTEGER :: HeatingCoilPLFCurvePTR  = 0      ! PLF curve index to gas or electric heating coil (used in latent degradation model)
  INTEGER :: BasinHeaterSchedulePtr  = 0      ! Pointer to basin heater schedule

! start of multi-speed compressor variables
  REAL(r64) :: RatedTotCap2           =0.0 ! Gross total cooling capacity at rated conditions, low speed [watts]
                                         ! Note: For HPWHs, RatedTotCap2   = Water Heating Capacity for Coil:DX:HPWH and
                                         !                  RatedTotCap(1) = Air Cooling Coil Capacity for Coil:DX:HPWH
  REAL(r64) :: RatedSHR2              =0.0 ! Sensible heat ratio (sens cap/total cap) at rated conditions, low speed
  REAL(r64) :: RatedCOP2              =0.0 ! Coefficient of performance at rated conditions, low speed
  REAL(r64) :: RatedAirVolFlowRate2   =0.0 ! Air volume flow rate through unit at rated conditions, low speed [m3/s]
  REAL(r64) :: RatedAirMassFlowRate2  =0.0 ! Air mass flow rate through unit at rated conditions, low speed [kg/s]
  REAL(r64) :: RatedCBF2              =0.0 ! rated coil bypass factor (low speed), determined using RatedTotCap2 and RatedSHR2
  INTEGER :: CCapFTemp2             = 0  ! index of total cooling capacity modifier curve (low speed)
  INTEGER :: EIRFTemp2              = 0  ! index of energy input ratio modifier curve (low speed)
                                         ! (function of entering wetbulb, outside drybulb)
  REAL(r64) :: RatedEIR2              =0.0 ! rated energy input ratio (low speed, inverse of COP2)
! end of multi-speed compressor variables

  REAL(r64) :: RatedEIR(MaxModes)     =0.0 ! rated energy input ratio (inverse of COP)
  REAL(r64) :: InletAirMassFlowRate   =0.0
  REAL(r64) :: InletAirMassFlowRateMax=0.0
  REAL(r64) :: InletAirTemp           =0.0
  REAL(r64) :: InletAirHumRat         =0.0
  REAL(r64) :: InletAirEnthalpy       =0.0
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!  REAL(r64) :: InletAirPressure       =0.0
  REAL(r64) :: OutletAirTemp          =0.0
  REAL(r64) :: OutletAirHumRat        =0.0
  REAL(r64) :: OutletAirEnthalpy      =0.0
  REAL(r64) :: PartLoadRatio          =0.0 ! Ratio of actual sensible cooling load to steady-state sensible cooling capacity
  REAL(r64) :: TotalCoolingEnergy     =0.0
  REAL(r64) :: SensCoolingEnergy      =0.0
  REAL(r64) :: LatCoolingEnergy       =0.0
  REAL(r64) :: TotalCoolingEnergyRate =0.0
  REAL(r64) :: SensCoolingEnergyRate  =0.0
  REAL(r64) :: LatCoolingEnergyRate   =0.0
  REAL(r64) :: ElecCoolingConsumption =0.0
  REAL(r64) :: ElecCoolingPower       =0.0
  REAL(r64) :: CoolingCoilRuntimeFraction =0.0 ! Run time fraction of the DX cooling unit

! start of variables used in heat pump heating coils only
  REAL(r64) :: TotalHeatingEnergy     =0.0
  REAL(r64) :: TotalHeatingEnergyRate =0.0
  REAL(r64) :: ElecHeatingConsumption =0.0
  REAL(r64) :: ElecHeatingPower       =0.0
  REAL(r64) :: HeatingCoilRuntimeFraction =0.0 ! Run time fraction of the DX heating unit
  INTEGER :: DefrostStrategy        = 0   ! defrost strategy; 1=reverse-cycle, 2=resistive
  INTEGER :: DefrostControl         = 0   ! defrost control; 1=timed, 2=on-demand
  INTEGER :: EIRFPLR                = 0   ! index of energy input ratio vs part-load ratio curve
  INTEGER :: DefrostEIRFT           = 0   ! index of defrost mode total cooling capacity for reverse cycle heat pump
  INTEGER :: RegionNum              = 0   ! Region number for calculating HSPF of single speed DX heating coil
  REAL(r64) :: MinOATCompressor       =0.0  ! Minimum OAT for heat pump compressor operation
  REAL(r64) :: OATempCompressorOn   = 0.0  ! The outdoor temperature when the compressor is automatically turned back on,
                                           ! if applicable, following automatic shut off. This field is used only for
                                           ! HSPF calculation.
  REAL(r64) :: MaxOATCompressor       =0.0  ! Maximum OAT for VRF heat pump compressor operation
  REAL(r64) :: MaxOATDefrost          =0.0  ! Maximum OAT for defrost operation
  REAL(r64) :: DefrostTime            =0.0  ! Defrost time period in hours
  REAL(r64) :: DefrostCapacity        =0.0  ! Resistive defrost to nominal capacity (at 21.11C/8.33C) ratio
  REAL(r64) :: HPCompressorRuntime    =0.0  ! keep track of compressor runtime
  REAL(r64) :: HPCompressorRuntimeLast=0.0  ! keep track of last time step compressor runtime (if simulation downshifts)
  REAL(r64) :: TimeLeftToDefrost      =0.0  ! keep track of time left to defrost heat pump
  REAL(r64) :: DefrostPower           =0.0  ! power used during defrost
  REAL(r64) :: DefrostConsumption     =0.0  ! energy used during defrost
  INTEGER   :: HeatingPerformanceOATType=DryBulbIndicator  ! Heating performance curve OAT type (1-wetbulb, 2-drybulb)
  LOGICAL   :: HPCoilIsInCoilSystemHeatingDX = .FALSE.
  LOGICAL   :: OATempCompressorOnOffBlank = .FALSE.
! end of variables used in heat pump heating coils only

! start of variables for DX cooling coil latent degradation model
  REAL(r64) :: Twet_Rated(MaxModes)                 =0.0 ! Nominal time for condensate to begin leaving the coil's
                                                       ! condensate drain line (sec)
  REAL(r64) :: Gamma_Rated(MaxModes)                =0.0 ! Initial moisture evaporation rate divided by steady-state
                                                       ! AC latent capacity (dimensionless)
  REAL(r64) :: MaxONOFFCyclesperHour(MaxModes)      =0.0 ! Maximum ON/OFF cycles per hour for the compressor (cycles/hour)
  REAL(r64) :: LatentCapacityTimeConstant(MaxModes) =0.0 ! Time constant for latent capacity to reach steady state
                                                       ! after startup (sec)
! end of variables for DX cooling coil latent degradation model

  INTEGER :: CondenserType(MaxModes) = AirCooled ! Type of condenser for DX cooling coil: AIR COOLED or EVAP COOLED

! start of variables for DX cooling coil evaporative condenser option
  LOGICAL :: ReportEvapCondVars =.false. ! true if any performance mode includes an evap condenser
  REAL(r64) :: EvapCondEffect(MaxModes) =0.0  ! effectiveness of the evaporatively cooled condenser
                                            ! [high speed for multi-speed unit] (-)
  REAL(r64) :: CondInletTemp  =0.0            ! Evap condenser inlet temperature [C], report variable
  REAL(r64) :: EvapCondAirFlow(MaxModes)=0.0  ! Air flow rate through the evap condenser at high speed,
                                                 ! for water use calcs [m3/s]
  REAL(r64) :: EvapCondPumpElecNomPower(MaxModes)=0.0  ! Nominal power input to the evap condenser water circulation pump
                                                     ! at high speed [W]
  REAL(r64) :: EvapCondPumpElecPower =0.0    ! Average power consumed by the evap condenser water circulation pump over
                                           ! the time step [W]
  REAL(r64) :: EvapCondPumpElecConsumption =0.0 ! Electric energy consumed by the evap condenser water circulation pump [J]
  REAL(r64) :: EvapWaterConsumpRate =0.0 ! Evap condenser water consumption rate [m3/s]
  REAL(r64) :: EvapWaterConsump =0.0 ! Evap condenser water consumption [m3]
  REAL(r64) :: EvapCondAirFlow2 =0.0 ! Air flow rate through the evap condenser at low speed, for water use calcs [m3/s]
  REAL(r64) :: EvapCondEffect2  =0.0 ! effectiveness of the evaporatively cooled condenser at low speed (-)
  REAL(r64) :: EvapCondPumpElecNomPower2 = 0.0  ! Nominal power input to the evap condenser water circulation pump at low speed [W]
  REAL(r64) :: BasinHeaterPower          = 0.0  ! Basin heater power (W)
  REAL(r64) :: BasinHeaterConsumption    = 0.0  ! Basin heater energy consumption (J)
! end of variables for DX cooling coil evaporative condenser option

! start of variables for Multimode DX cooling coil
  INTEGER :: NumCapacityStages=1 ! number of capacity stages, up to MaxCapacityStages for Multimode DX coil,
                                 ! always 1 for other coils
  INTEGER :: NumDehumidModes=0   ! number of enhanced dehumidification modes, up to MaxDehumidModes for Multimode DX coil,
                                 ! always 0 for other coils)
  CHARACTER(len=MaxNameLength) :: CoilPerformanceType(MaxModes)    =' '  ! Coil Performance object type
  INTEGER                      :: CoilPerformanceType_Num(MaxModes)= 0   ! Coil Performance object type number
  CHARACTER(len=MaxNameLength) :: CoilPerformanceName(MaxModes)    =' '  ! Coil Performance object names
  REAL(r64) :: CoolingCoilStg2RuntimeFrac =0.0 ! Run time fraction of stage 2
  INTEGER :: DehumidificationMode               =0   ! Dehumidification mode for multimode coil,
                                                     ! 0=normal, 1+=enhanced dehumidification mode
! end of variables for Multimode DX cooling coil

! start of variables for heat pump water heater DX coil
  INTEGER :: WaterInNode                 = 0       ! Condenser water inlet node number for HPWH DX coil
  INTEGER :: WaterOutNode                = 0       ! Condenser water outlet node number for HPWH DX coil
  INTEGER :: HCOPFTemp                   = 0       ! COP as a function of temperature curve index
  INTEGER :: HCOPFTempErrorIndex         = 0       ! Used for warning messages when output of HCOPFTemp is negative
  INTEGER :: HCOPFTempCurveType          = 0       ! COP as a function of temperature curve type
  INTEGER :: HCOPFAirFlow                = 0       ! COP as a function of air flow rate ratio curve index
  INTEGER :: HCOPFAirFlowErrorIndex      = 0       ! Used for warning messages when output of HCOPFAirFlow is negative
  INTEGER :: HCOPFWaterFlow              = 0       ! COP as a function of water flow rate ratio curve index
  INTEGER :: HCOPFWaterFlowErrorIndex    = 0       ! Used for warning messages when output of HCOPFWaterFlow is negative
  INTEGER :: HCapFTemp                   = 0       ! Heating capacity as a function of temperature curve index
  INTEGER :: HCapFTempErrorIndex         = 0       ! Used for warning messages when output of HCapFTemp is negative
  INTEGER :: HCapFTempCurveType          = 0       ! Heating capacity as a function of temperature curve type
  INTEGER :: HCapFAirFlow                = 0       ! Heating capacity as a function of air flow rate ratio curve index
  INTEGER :: HCapFAirFlowErrorIndex      = 0       ! Used for warning messages when output of HCapFAirFlow is negative
  INTEGER :: HCapFWaterFlow              = 0       ! Heating capacity as a function of water flow rate ratio curve index
  INTEGER :: HCapFWaterFlowErrorIndex    = 0       ! Used for warning messages when output of HCapFWaterFlow is negative
  INTEGER :: InletAirTemperatureType     = 0       ! Specifies to use either air wet-bulb or dry-bulb temp for curve objects
  REAL(r64) :: RatedInletDBTemp            = 0.0     ! Rated inlet air dry-bulb temperature [C]
  REAL(r64) :: RatedInletWBTemp            = 0.0     ! Rated inlet air wet-bulb temperature [C]
  REAL(r64) :: RatedInletWaterTemp         = 0.0     ! Rated condenser water inlet temperature [C]
!  REAL(r64) :: CondenserInletWaterTemp     = 0.0     ! Actual inlet water temperature to condenser of the HPWH DX coil [C]
  REAL(r64) :: HPWHCondPumpElecNomPower    = 0.0     ! Nominal power input to the condenser water circulation pump [W]
  REAL(r64) :: HPWHCondPumpFracToWater     = 0.0     ! Nominal power fraction to water for the condenser water circulation pump
  REAL(r64) :: RatedHPWHCondWaterFlow      = 0.0     ! Rated water flow rate through the condenser of the HPWH DX coil [m3/s]
  REAL(r64) :: ElecWaterHeatingPower       = 0.0     ! Total electric power consumed by compressor and condenser pump [W]
  REAL(r64) :: ElecWaterHeatingConsumption = 0.0     ! Total electric consumption by compressor and condenser pump [J]
  LOGICAL :: FanPowerIncludedInCOP       = .TRUE.  ! Indicates that fan heat is included in heating capacity and COP
  LOGICAL :: CondPumpHeatInCapacity      = .FALSE. ! Indicates that condenser pump heat is included in heating capacity
  LOGICAL :: CondPumpPowerInCOP          = .FALSE. ! Indicates that condenser pump power is included in heating COP
  LOGICAL :: AirVolFlowAutoSized         = .FALSE. ! Used to report autosizing info for the HPWH DX coil
  LOGICAL :: WaterVolFlowAutoSized       = .FALSE. ! Used to report autosizing info for the HPWH DX coil
! end of variables for heat pump water heater DX coil

! Error tracking
  REAL(r64) :: LowTempLast=0.0 ! low ambient temp entering condenser when warning message occurred
  REAL(r64) :: HighTempLast=0.0 ! high ambient temp entering condenser when warning message occurred
  INTEGER :: ErrIndex1=0     ! index/pointer to recurring error structure for Air volume flow rate per watt of
                             ! rated total cooling capacity error
  INTEGER :: ErrIndex2=0     ! index/pointer to recurring error structure for PLF curve values must be >= 0.7. error
  INTEGER :: ErrIndex3=0     ! index/pointer to recurring error structure for DX cooling coil runtime fraction > 1.0 warning
  INTEGER :: ErrIndex4=0     ! index/pointer to recurring error structure for DX heating coil runtime fraction > 1.0 warning
  INTEGER :: LowAmbErrIndex=0   ! index/pointer to recurring error structure for low ambient temp entering condenser
  INTEGER :: HighAmbErrIndex=0  ! index/pointer to recurring error structure for high ambient temp entering condenser
  INTEGER :: PLFErrIndex=0   ! index/pointer to recurring error structure for PLF <> 1 at speed 1 for a multiple speed coil
  INTEGER :: PLRErrIndex=0   ! index/pointer to recurring error structure for PLR < .7
  LOGICAL :: PrintLowAmbMessage= .FALSE.  ! used to print warning message for low ambient conditions
  CHARACTER(len=300) :: LowAmbBuffer1=' ' ! holds warning message until next iteration (only prints 1 message/iteration)
  CHARACTER(len=300) :: LowAmbBuffer2=' ' ! holds warning message until next iteration (only prints 1 message/iteration)
  LOGICAL :: PrintHighAmbMessage= .FALSE.  ! used to print warning message for high ambient conditions
  CHARACTER(len=300) :: HighAmbBuffer1=' ' ! holds warning message until next iteration (only prints 1 message/iteration)
  CHARACTER(len=300) :: HighAmbBuffer2=' ' ! holds warning message until next iteration (only prints 1 message/iteration)

  !begin variables for Water System interactions
  INTEGER ::EvapWaterSupplyMode                   = WaterSupplyFromMains !  where does water come from
  CHARACTER(len=MaxNameLength) :: EvapWaterSupplyName = ' ' ! name of water source e.g. water storage tank
  INTEGER ::EvapWaterSupTankID                    = 0 !
  INTEGER ::EvapWaterTankDemandARRID              = 0 !
  INTEGER ::CondensateCollectMode                 = CondensateDiscarded !  where does water come from
  CHARACTER(len=MaxNameLength) :: CondensateCollectName = ' ' ! name of water source e.g. water storage tank
  INTEGER ::CondensateTankID                      = 0 !
  INTEGER ::CondensateTankSupplyARRID             = 0 !

  REAL(r64)   :: CondensateVdot = 0.0 ! rate of water condensation from air stream [m3/s]
  REAL(r64)   :: CondensateVol  = 0.0 ! amount of water condensed from air stream [m3]

  !end variables for water system interactions

  ! used to print low ambient warning message for DOE2 coil only after time step has incremented
  REAL(r64) :: CurrentEndTimeLast  = 0.0 ! end time of time step for last simulation time step
  REAL(r64) :: TimeStepSysLast     = 0.0 ! last system time step (used to check for downshifting)
  ! for multispeed DX coil type
  INTEGER :: FuelType       =0   ! Fuel type
  INTEGER :: NumOfSpeeds    =0   ! Number of speeds
  LOGICAL :: PLRImpact      =.FALSE.   ! Part load fraction applied to Speed Number > 1
  LOGICAL :: LatentImpact   =.FALSE.   ! Latent degradation applied to Speed Number > 1
  INTEGER, DIMENSION(:), ALLOCATABLE :: MSErrIndex ! index flag for num speeds/recurring messages
  REAL(r64), DIMENSION(:), ALLOCATABLE :: MSRatedTotCap ! Rated cooling capacity for MS heat pump [W]
  REAL(r64), DIMENSION(:), ALLOCATABLE :: MSRatedSHR    ! Rated SHR for MS heat pump [dimensionless]
  REAL(r64), DIMENSION(:), ALLOCATABLE :: MSRatedCOP    ! Rated COP for MS heat pump [dimensionless]
  REAL(r64), DIMENSION(:), ALLOCATABLE :: MSRatedAirVolFlowRate  ! Air volume flow rate through unit at rated conditions [m3/s]
  REAL(r64), DIMENSION(:), ALLOCATABLE :: MSRatedAirMassFlowRate ! Air mass flow rate through unit at rated conditions [m3/s]
  REAL(r64), DIMENSION(:), ALLOCATABLE :: MSRatedCBF    ! rated coil bypass factor
  INTEGER, DIMENSION(:), ALLOCATABLE :: MSCCapFTemp   ! index of total cooling capacity modifier curve
  INTEGER, DIMENSION(:), ALLOCATABLE :: MSCCapFFlow   ! index of total cooling capacity modifier curve
  INTEGER, DIMENSION(:), ALLOCATABLE :: MSEIRFTemp    ! index of energy input ratio modifier curve as a function of temperature
  INTEGER, DIMENSION(:), ALLOCATABLE :: MSEIRFFlow    ! index of energy input ratio modifier curve as a function of flow fraction
  INTEGER, DIMENSION(:), ALLOCATABLE :: MSPLFFPLR     ! index of part load factor as a function of part load ratio
  INTEGER, DIMENSION(:), ALLOCATABLE :: MSWasteHeat   ! index of waste heat as a function of temperature
  REAL(r64), DIMENSION(:), ALLOCATABLE :: MSWasteHeatFrac  ! Waste heat fraction
  REAL(r64), DIMENSION(:), ALLOCATABLE :: MSEvapCondEffect ! effectiveness of the evaporatively cooled condenser
  REAL(r64), DIMENSION(:), ALLOCATABLE :: MSEvapCondAirFlow ! Air flow rate through the evap condenser for water use calcs [m3/s]
  REAL(r64), DIMENSION(:), ALLOCATABLE :: MSEvapCondPumpElecNomPower ! Nominal power input to the evap condenser
                                                                     ! water circulation pump
  INTEGER, DIMENSION(:), ALLOCATABLE :: MSTotCapTempModFacCurveType ! type of curve for CCapFTemp (cubic,quadratic,bi-quadratic)
  INTEGER, DIMENSION(:), ALLOCATABLE :: MSEIRTempModFacCurveType !type of curve for EIRFTemp (cubic,quadratic,bi-quadratic)

  REAL(r64), DIMENSION(:), ALLOCATABLE :: MSTwet_Rated      ! Nominal time for condensate to begin leaving the coil's
                                                   ! condensate drain line (sec)
  REAL(r64), DIMENSION(:), ALLOCATABLE :: MSGamma_Rated     ! Initial moisture evaporation rate divided by steady-state
                                                   ! AC latent capacity (dimensionless)
  REAL(r64), DIMENSION(:), ALLOCATABLE :: MSMaxONOFFCyclesperHour ! Maximum ON/OFF cycles per hour for the compressor (cycles/hour)
  REAL(r64), DIMENSION(:), ALLOCATABLE :: MSLatentCapacityTimeConstant ! Time constant for latent capacity to reach steady state
  REAL(r64) :: FuelUsed            ! Energy used, in addition to electricity [W]
  REAL(r64) :: FuelConsumed        ! Energy consumed, in addition to electricity [J]
  ! End of multispeed DX coil input

  ! VRF system variables used for sizing
  LOGICAL :: CoolingCoilPresent = .TRUE.   ! FALSE if coil not present
  LOGICAL :: HeatingCoilPresent = .TRUE.   ! FALSE if coil not present
  
  REAL(r64) :: CndFanPwr !RS: Debugging: Trying to pass through power from HPSim (9/3/14)
  REAL(R64) :: CompPwr   !RS: Debugging: Trying to pass through power from HPSim (9/3/14)
  
END TYPE DXCoilData

  ! MODULE VARIABLE DECLARATIONS:
TYPE (DXCoilData) , PUBLIC, ALLOCATABLE, DIMENSION(:) :: DXCoil
REAL(r64), PUBLIC, ALLOCATABLE, DIMENSION(:) :: DXCoilOutletTemp    ! DX coil outlet dry bulb temperature [C]
REAL(r64), PUBLIC, ALLOCATABLE, DIMENSION(:) :: DXCoilOutletHumRat  ! DX coil outlet humidity ratio [kgWater/kgDryAir]
REAL(r64), PUBLIC, ALLOCATABLE, DIMENSION(:) :: DXCoilPartLoadRatio ! DX coil part-load ratio
INTEGER,   PUBLIC, ALLOCATABLE, DIMENSION(:) :: DXCoilFanOpMode       ! supply air fan operating mode
REAL(r64), PUBLIC, ALLOCATABLE, DIMENSION(:) :: DXCoilFullLoadOutAirTemp ! DX coil full load outlet dry bulb temperature [C]
REAL(r64), PUBLIC, ALLOCATABLE, DIMENSION(:) :: DXCoilFullLoadOutAirHumRat ! DX coil full load outlet humidity ratio [kgWater/kgDryAir]
REAL(r64), PUBLIC, ALLOCATABLE, DIMENSION(:) :: DXCoilTotalCooling ! DX cooling coil total cooling output [W]
REAL(r64), PUBLIC, ALLOCATABLE, DIMENSION(:) :: DXCoilTotalHeating ! DX heating coil total heating output [W]
REAL(r64), PUBLIC, ALLOCATABLE, DIMENSION(:) :: DXCoilCoolInletAirWBTemp ! DX cooling coil inlet air wet-bulb temp [C]
REAL(r64), PUBLIC, ALLOCATABLE, DIMENSION(:) :: DXCoilHeatInletAirDBTemp ! DX heating coil inlet air dry-bulb temp [C]
REAL(r64), PUBLIC, ALLOCATABLE, DIMENSION(:) :: DXCoilHeatInletAirWBTemp ! DX heating coil inlet air wet-bulb temp [C]
INTEGER, PUBLIC :: CurDXCoilNum=0

INTEGER, PUBLIC :: NumDXCoils             = 0      ! Total number of DX coils
REAL(r64), PUBLIC    :: HPWHHeatingCapacity    = 0.0    ! Used by Heat Pump:Water Heater object as total water heating capacity [W]
REAL(r64), PUBLIC    :: HPWHHeatingCOP         = 0.0    ! Used by Heat Pump:Water Heater object as water heating COP [W/W]
LOGICAL         :: GetCoilsInputFlag      = .TRUE. ! First time, input is "gotten"
INTEGER, PRIVATE :: NumVRFHeatingCoils     = 0     ! number of VRF heat pump heating coils
INTEGER, PRIVATE :: NumVRFCoolingCoils     = 0     ! number of VRF heat pump cooling coils
INTEGER, PRIVATE :: NumDXHeatingCoils     = 0      ! number of DX heat pump heating coils
INTEGER, PRIVATE :: NumDoe2DXCoils        = 0      ! number of doe2 DX  coils
INTEGER, PRIVATE :: NumDXHeatPumpWaterHeaterCoils = 0 ! number of DX  water heater coils
INTEGER, PRIVATE :: NumDXMulSpeedCoils    = 0      ! number of DX coils with multi-speed compressor
INTEGER, PRIVATE :: NumDXMulModeCoils     = 0      ! number of DX coils with multi-mode performance

INTEGER, PRIVATE :: NumDXMulSpeedCoolCoils= 0      ! number of multispeed DX cooling coils
INTEGER, PRIVATE :: NumDXMulSpeedHeatCoils= 0      ! number of multispeed DX heating coils
!INTEGER, PRIVATE :: NumHPSimCoils= 0    !RS: Implement
LOGICAL, ALLOCATABLE, DIMENSION(:) :: CheckEquipName

INTEGER, PUBLIC :: DXCoilHPSimNum = 0  !RS: Implementation: Making another variable for HPSim implementation
INTEGER, PUBLIC :: DXCoilHPSimNum2 = 0  !RS: Implementation: Making another variable for HPSim implementation (10/31/14)
REAL, PUBLIC :: WeathOutDryBulb !RS: Debugging: Getting the outside air properties for the HPSim Condenser (11/3/14)
REAL, PUBLIC :: WeathOutHumRat !RS: Debugging: Getting the outside air properties for the HPSim Condenser (11/3/14)
REAL, PUBLIC :: WeathOutBaroPress !RS: Debugging: Getting the outside air properties for the HPSim Condenser (11/3/14)

  ! SUBROUTINE SPECIFICATIONS FOR MODULE

          ! Driver/Manager Routines
PUBLIC  SimDXCoil
PUBLIC  SimDXCoilMultiSpeed
PUBLIC  SimDXCoilMultiMode

          ! Get Input routines for module
PRIVATE GetDXCoils

          ! Initialization routines for module
PRIVATE InitDXCoil
PRIVATE SizeDXCoil
PRIVATE CalcBasinHeaterPowerForMultiModeDXCoil

          ! Update routines to check convergence and update nodes
PUBLIC  CalcHPWHDXCoil
PUBLIC  CalcDoe2DXCoil
PUBLIC  CalcHPSimDXCoil !RS: Implementation: Trying to create a new subroutine to call HPSim
PUBLIC  CalcVRFCoolingCoil
PUBLIC  CalcDXHeatingCoil
PUBLIC  CalcMultiSpeedDXCoil
PUBLIC  CalcMultiSpeedDXCoilCooling
PUBLIC  CalcMultiSpeedDXCoilHeating
PRIVATE UpdateDXCoil
PRIVATE ReportDXCoil

          ! Common routines
PRIVATE CalcTotCapSHR
PRIVATE CalcCBF
PRIVATE AdjustCBF

          ! External function calls
PUBLIC  GetDXCoilIndex
PUBLIC  GetCoilCapacity
PUBLIC  GetCoilCapacityByIndexType
PUBLIC  GetCoilTypeNum
PUBLIC  GetMinOATCompressor
PUBLIC  GetCoilInletNode
PUBLIC  GetCoilOutletNode
PUBLIC  GetCoilCondenserInletNode
PUBLIC  GetDXCoilBypassedFlowFrac
PRIVATE GetHPCoolingCoilIndex
PUBLIC  GetDXCoilNumberOfSpeeds
PUBLIC  SetDXCoolingCoilData
PUBLIC  SetCoilSystemHeatingDXFlag

CONTAINS

SUBROUTINE SimDXCoil(CompName,CompOp,FirstHVACIteration,PartLoadRatio,CompIndex,FanOpMode,OnOffAFR, CoilCoolingHeatingPLRRatio, &
                     MaxCap, CompCyclingRatio)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Fred Buhl
          !       DATE WRITTEN   May 2000
          !       MODIFIED       Don Shirey, Sept 2000, October 2001, June 2005
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! Manages the simulation of a single speed on/off DX coil.

          ! METHODOLOGY EMPLOYED:
          ! NA

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor, ONLY: FindItemInList
  USE General,        ONLY: TrimSigDigits

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT (IN)           :: CompName            ! name of the fan coil unit
  INTEGER         , INTENT (IN)           :: CompOp              ! compressor operation; 1=on, 0=off
  LOGICAL         , INTENT(IN)            :: FirstHVACIteration  ! True when first HVAC iteration
  REAL(r64)       , INTENT (IN), OPTIONAL :: PartLoadRatio       ! part load ratio (for single speed cycling unit)
  INTEGER         , INTENT (INOUT)        :: CompIndex
  INTEGER         , INTENT (IN)           :: FanOpMode  ! allows parent object to control fan mode
  REAL(r64)       , INTENT (IN), OPTIONAL :: OnOffAFR   ! ratio of compressor on airflow to compressor off airflow
  REAL(r64)       , INTENT (IN), OPTIONAL :: CoilCoolingHeatingPLRRatio ! used for cycling fan RH control
  REAL(r64)       , INTENT (IN), OPTIONAL :: MaxCap ! maximum cooling capacity of VRF terminal units
  REAL(r64)       , INTENT (IN), OPTIONAL :: CompCyclingRatio ! cycling ratio of VRF condenser connected to this TU

          ! SUBROUTINE PARAMETER DEFINITIONS:
  CHARACTER(len=*), PARAMETER :: Blank = ' '

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER      :: DXCoilNum        ! index of fan coil unit being simulated
  REAL(r64) :: AirFlowRatio        ! ratio of compressor on airflow to compressor off airflow
  REAL(r64) :: CompCycRatio        ! compressor cycling ratio of VRF condenser

          ! FLOW

! First time SimDXCoil is called, get the input for all the DX coils (condensing units)
IF (GetCoilsInputFlag) THEN
  CALL GetDXCoils
  GetCoilsInputFlag = .FALSE. ! Set GetInputFlag false so you don't get coil inputs again
END IF

IF (CompIndex == 0) THEN
  DXCoilNum = FindItemInList(CompName,DXCoil%Name,NumDXCoils)
  IF (DXCoilNum == 0) THEN
    CALL ShowFatalError('DX Coil not found='//TRIM(CompName))
  ENDIF
  CompIndex=DXCoilNum
ELSE
  DXCoilNum=CompIndex
  IF (DXCoilNum > NumDXCoils .or. DXCoilNum < 1) THEN
    CALL ShowFatalError('SimDXCoil: Invalid CompIndex passed='//  &
                        TRIM(TrimSigDigits(DXCoilNum))// &
                        ', Number of DX Coils='//TRIM(TrimSigDigits(NumDXCoils))//  &
                        ', Coil name='//TRIM(CompName))
  ENDIF
  IF (CheckEquipName(DXCoilNum)) THEN
    IF (CompName /= Blank .AND. CompName /= DXCoil(DXCoilNum)%Name) THEN
      CALL ShowFatalError('SimDXCoil: Invalid CompIndex passed='//  &
                          TRIM(TrimSigDigits(DXCoilNum))// &
                          ', Coil name='//TRIM(CompName)//', stored Coil Name for that index='//  &
                          TRIM(DXCoil(DXCoilNum)%Name))
    ENDIF
    CheckEquipName(DXCoilNum)=.false.
  ENDIF
ENDIF

IF (PRESENT(OnOffAFR)) THEN
  AirFlowRatio = OnOffAFR
ELSE
  AirFlowRatio = 1.0d0
END IF

IF(PRESENT(CompCyclingRatio))THEN
  CompCycRatio = CompCyclingRatio
ELSE
  CompCycRatio = 1.d0
END IF

CurDXCoilNum = DXCoilNum

IF (HPSimFlag .EQ. 1) THEN  !RS: Implementation: Debugging: Setting up a case when HPSim is called
    !IF (WarmUpFlag) THEN  !RS: Debugging: Trying to get it to run as a curve-fit model during warmup (8/11/14)
    !    DXCoil(DXCoilNum)%DXCoilType_Num=CoilDX_CoolingSingleSpeed  !RS: Make sure it's same CASE as in IDF (8/11/14)
    !ELSE
    !    DXCoil(DXCoilNum)%DXCoilType_Num=CoilDX_HPSim
    !END IF
END IF

! Initialize the DX coil unit
CALL InitDXCoil(DXCoilNum)

! Select the correct unit type
SELECT CASE(DXCoil(DXCoilNum)%DXCoilType_Num)

CASE (CoilDX_CoolingSingleSpeed)

    IF (WARMUPFLAG) THEN   
    IF (PRESENT(CoilCoolingHeatingPLRRatio)) THEN
      CALL CalcDoe2DXCoil(DXCoilNum,CompOp,FirstHVACIteration,PartLoadRatio, FanOpMode, &
                          OnOffAirFlowRatio=AirFlowRatio, CoolingHeatingPLR=CoilCoolingHeatingPLRRatio)
    ELSE
      CALL CalcDoe2DXCoil(DXCoilNum,CompOp,FirstHVACIteration,PartLoadRatio, FanOpMode,OnOffAirFlowRatio=AirFlowRatio)
    END IF
    ELSE    !RS: Debugging: Trying to handle case when the coil isn't called because it's not needed (12/17/14)
        !IF (CurrentTime .LT. 12 .OR. CurrentTime .GT. 13) THEN    !RS: Debugging: Hardcoding so only one hour of HPSim is being run (3/10/18)    
            IF (PRESENT(CoilCoolingHeatingPLRRatio)) THEN
                CALL CalcDoe2DXCoil(DXCoilNum,CompOp,FirstHVACIteration,PartLoadRatio, FanOpMode, &
                          OnOffAirFlowRatio=AirFlowRatio, CoolingHeatingPLR=CoilCoolingHeatingPLRRatio)
            ELSE
                DXCoil(DXCoilNum)%OutletAirEnthalpy=DXCoil(DXCoilNum)%InletAirEnthalpy
                DXCoil(DXCoilNum)%OutletAirTemp=DXCoil(DXCoilNum)%InletAirTemp
                DXCoil(DXCoilNum)%OutletAirHumRat=DXCoil(DXCoilNum)%InletAirHumRat
                DXCoil(DXCoilNum)%TotalCoolingEnergyRate=0  !RS: Debugging: The coil is "off" so it's 0 (8/6/16)
                DXCoil(DXCoilNum)%SensCoolingEnergyRate=0   !RS: Debugging: The coil is "off" so it's 0 (8/6/16)
                DXCoil(DXCoilNum)%TotalCoolingEnergy=0  !RS: Debugging: The coil is "off" so it's 0 (8/6/16)
                DXCoil(DXCoilNum)%ElecCoolingPower=0    !RS: Debugging: The coil is "off" so it's 0 (8/6/16)
                DXCoil(DXCoilNum)%LatCoolingEnergyRate=0    !RS: Debugging: The coil is "off" so it's 0 (8/6/16)
            END IF
        !END IF    !RS: Debugging: Hardcoding so only one hour of HPSim is being run (3/10/18)
    END IF

  CASE (CoilDX_HeatingEmpirical)

    CALL CalcDXHeatingCoil(DXCoilNum,PartLoadRatio, FanOpMode,OnOffAirFlowRatio=AirFlowRatio)

  CASE (CoilDX_HeatPumpWaterHeater)

!   call the HPWHDXCoil routine to calculate water side performance set up the DX coil info for air-side calcs
    CALL CalcHPWHDXCoil(DXCoilNum,PartLoadRatio)
!    CALL CalcDOE2DXCoil(DXCoilNum, CompOp, FirstHVACIteration,PartLoadRatio), perform air-side calculations
    CALL CalcDOE2DXCoil(DXCoilNum, 1, FirstHVACIteration,PartLoadRatio, FanOpMode)

  CASE (CoilVRF_Cooling)

    CALL CalcVRFCoolingCoil(DXCoilNum, 1, FirstHVACIteration,PartLoadRatio, FanOpMode, CompCycRatio, MaxCoolCap=MaxCap)

  CASE (CoilVRF_Heating)

    CALL CalcDXHeatingCoil(DXCoilNum,PartLoadRatio, FanOpMode, MaxHeatCap=MaxCap)
    
  CASE (CoilDX_HPSim)   !RS: Implementation: Debugging: Setting up a case when HPSim is called
      
      DXCoilHPSimNum=DXCoilNum
      CALL CalcHPSimDXCoil(1, FirstHVACIteration,PartLoadRatio, FanOpMode)

  CASE DEFAULT
    CALL ShowSevereError('Error detected in DX Coil='//TRIM(CompName))
    CALL ShowContinueError('Invalid DX Coil Type='//TRIM(DXCoil(DXCoilNum)%DXCoilType))
    CALL ShowFatalError('Preceding condition causes termination.')


END SELECT

! Update the unit outlet nodes
CALL UpdateDXCoil(DXCoilNum)

! Report the result of the simulation
CALL ReportDXCoil(DXCoilNum)


RETURN
END SUBROUTINE SimDXCoil

SUBROUTINE SimDXCoilMultiSpeed(CompName,SpeedRatio,CycRatio,CompIndex,SpeedNum,FanOpMode,CompOp)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Fred Buhl
          !       DATE WRITTEN   September 2002
          !       MODIFIED       Lixing Gu, Sep. 2007
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! Manages the simulation of a multi speed DX coil.

          ! METHODOLOGY EMPLOYED:
          ! NA

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor, ONLY: FindItemInList
  USE General,        ONLY: TrimSigDigits

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT (IN) :: CompName            ! name of the fan coil unit
  REAL(r64)       , INTENT (IN) :: SpeedRatio          ! = (CompressorSpeed - CompressorSpeedMin) /
                                                       !   (CompressorSpeedMax - CompressorSpeedMin)
                                                       ! for variable speed or 2 speed compressors
  REAL(r64)       , INTENT (IN) :: CycRatio            ! cycling part load ratio for variable speed
                                                       ! or 2 speed compressors
  INTEGER, INTENT(INOUT)        :: CompIndex
  INTEGER, INTENT(IN), OPTIONAL :: SpeedNum            ! Speed number for multispeed cooling coil onlyn
  INTEGER, INTENT(IN), OPTIONAL :: FanOpMode           ! Fan operation mode
  INTEGER, INTENT(IN), OPTIONAL :: CompOp              ! Compressor on/off; 1=on, 0=off

          ! SUBROUTINE PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER      :: DXCoilNum              ! index of fan coil unit being simulated

          ! FLOW

! First time SimDXCoil is called, get the input for all the DX coils (condensing units)
IF (GetCoilsInputFlag) THEN
  CALL GetDXCoils
  GetCoilsInputFlag = .FALSE. ! Set GetInputFlag false so you don't get coil inputs again
END IF

!  find correct DX Coil

IF (CompIndex == 0) THEN
  DXCoilNum = FindItemInList(CompName,DXCoil%Name,NumDXCoils)
  IF (DXCoilNum == 0) THEN
    CALL ShowFatalError('DX Coil not found='//TRIM(CompName))
  ENDIF
  CompIndex=DXCoilNum
ELSE
  DXCoilNum=CompIndex
  IF (DXCoilNum > NumDXCoils .or. DXCoilNum < 1) THEN
    CALL ShowFatalError('SimDXCoilMultiSpeed: Invalid CompIndex passed='//  &
                        TRIM(TrimSigDigits(DXCoilNum))// &
                        ', Number of DX Coils='//TRIM(TrimSigDigits(NumDXCoils))//  &
                        ', Coil name='//TRIM(CompName))
  ENDIF
  IF (CheckEquipName(DXCoilNum)) THEN
    IF (CompName /= DXCoil(DXCoilNum)%Name) THEN
      CALL ShowFatalError('SimDXCoilMultiSpeed: Invalid CompIndex passed='//  &
                          TRIM(TrimSigDigits(DXCoilNum))// &
                          ', Coil name='//TRIM(CompName)//', stored Coil Name for that index='//  &
                          TRIM(DXCoil(DXCoilNum)%Name))
    ENDIF
    CheckEquipName(DXCoilNum)=.false.
  ENDIF
ENDIF

CurDXCoilNum = DXCoilNum

! Initialize the DX coil unit
CALL InitDXCoil(DXCoilNum)

! Select the correct unit type
SELECT CASE(DXCoil(DXCoilNum)%DXCoilType_Num)

  CASE (CoilDX_CoolingTwoSpeed)

    CALL CalcMultiSpeedDXCoil(DXCoilNum,SpeedRatio,CycRatio)

  CASE (CoilDX_MultiSpeedCooling)
    If (PRESENT(SpeedNum)) CALL CalcMultiSpeedDXCoilCooling(DXCoilNum,SpeedRatio,CycRatio,SpeedNum,FanOpMode,CompOp)

  CASE (CoilDX_MultiSpeedHeating)
    If (PRESENT(SpeedNum)) CALL CalcMultiSpeedDXCoilHeating(DXCoilNum,SpeedRatio,CycRatio,SpeedNum,FanOpMode)

  CASE DEFAULT
    CALL ShowSevereError('Error detected in DX Coil='//TRIM(CompName))
    CALL ShowContinueError('Invalid DX Coil Type='//TRIM(DXCoil(DXCoilNum)%DXCoilType))
    CALL ShowFatalError('Preceding condition causes termination.')


END SELECT

! Update the unit outlet nodes
CALL UpdateDXCoil(DXCoilNum)

! Report the result of the simulation
CALL ReportDXCoil(DXCoilNum)


RETURN

END SUBROUTINE SimDXCoilMultiSpeed

SUBROUTINE SimDXCoilMultiMode(CompName,CompOp,FirstHVACIteration,PartLoadRatio,DehumidMode,CompIndex,FanOpMode)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         M. J. Witte (based on SimDXCoilMultiSpeed by Fred Buhl)
          !       DATE WRITTEN   February 2005
          !       MODIFIED       April 2010, Chandan sharma, added basin heater
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! Manages the simulation of a DX coil with multiple performance modes, such as
          ! multiple stages, or sub-cool reheat for humidity control.

          ! METHODOLOGY EMPLOYED:
          ! NA

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor, ONLY: FindItemInList
  USE General,        ONLY: TrimSigDigits

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT (IN) :: CompName            ! name of the fan coil unit
  INTEGER,          INTENT (IN) :: CompOp              ! compressor operation; 1=on, 0=off !unused1208
  LOGICAL,          INTENT (IN) :: FirstHVACIteration  ! true if first hvac iteration
  REAL(r64),        INTENT (IN) :: PartLoadRatio       ! part load ratio
  INTEGER,          INTENT (IN) :: DehumidMode         ! dehumidification mode (0=normal, 1=enhanced)
  INTEGER,        INTENT(INOUT) :: CompIndex
  INTEGER,          INTENT (IN) :: FanOpMode           ! allows parent object to control fan mode

          ! SUBROUTINE PARAMETER DEFINITIONS:
  CHARACTER(len=*), PARAMETER :: RoutineName='SimDXCoilMultiMode'

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER :: DXCoilNum                 ! index of coil being simulated
  INTEGER :: PerfMode                  ! Performance mode for MultiMode DX coil; Always 1 for other coil types
                                       ! 1-2=normal mode: 1=stage 1 only, 2=stage 1&2
                                       ! 3-4=enhanced dehumidification mode: 3=stage 1 only, 4=stage 1&2
  REAL(r64) :: AirMassFlow               ! Dry air mass flow rate through coil [kg/s]

  REAL(r64) :: S1OutletAirTemp           ! Stage 1   Outlet air dry bulb temp [C]
  REAL(r64) :: S1OutletAirHumRat         ! Stage 1   Outlet air humidity ratio [kgWater/kgDryAir]
  REAL(r64) :: S1OutletAirEnthalpy       ! Stage 1   Outlet air enthalpy
  REAL(r64) :: S1PLR                     ! Stage 1   Ratio of actual sensible cooling load to
                                         !           steady-state sensible cooling capacity
  REAL(r64) :: S1TotalCoolingEnergyRate  ! Stage 1   Total cooling rate [W]
  REAL(r64) :: S1SensCoolingEnergyRate   ! Stage 1   Sensible cooling rate [W]
  REAL(r64) :: S1LatCoolingEnergyRate    ! Stage 1   Latent cooling rate [W]
  REAL(r64) :: S1ElecCoolingPower        ! Stage 1   Electric power input [W]
  REAL(r64) :: S1RuntimeFraction =0.0    ! Stage 1   Run time fraction (overlaps with stage1&2 run time)
  REAL(r64) :: S1EvapCondPumpElecPower   ! Stage 1   Evaporative condenser pump electric power input [W]
  REAL(r64) :: S1EvapWaterConsumpRate    ! Stage 1   Evap condenser water consumption rate [m3/s]
  REAL(r64) :: S1CrankcaseHeaterPower    ! Stage 1   Report variable for average crankcase heater power [W]
  REAL(r64) :: S1FFullLoadOutAirTemp     ! Stage 1   Full load outlet temperature [C]
  REAL(r64) :: S1FullLoadOutAirHumRat    ! Stage 1   Full load outlet humidity ratio [kgWater/kgDryAir]

  REAL(r64) :: S12OutletAirTemp          ! Stage 1&2 Outlet air dry bulb temp [C]
  REAL(r64) :: S12OutletAirHumRat        ! Stage 1&2 Outlet air humidity ratio [kgWater/kgDryAir]
  REAL(r64) :: S12OutletAirEnthalpy      ! Stage 1&2 Outlet air enthalpy
!  REAL(r64) :: S12PLR                  ! Stage 1&2 Ratio of actual sensible cooling load to
!                                       !           steady-state sensible cooling capacity
  REAL(r64) :: S12TotalCoolingEnergyRate ! Stage 1&2 Total cooling rate [W]
  REAL(r64) :: S12SensCoolingEnergyRate  ! Stage 1&2 Sensible cooling rate [W]
  REAL(r64) :: S12LatCoolingEnergyRate   ! Stage 1&2 Latent cooling rate [W]
  REAL(r64) :: S12ElecCoolingPower       ! Stage 1&2 Electric power input [W]
  REAL(r64) :: S12ElecCoolFullLoadPower  ! Stage 1&2 Electric power input at full load (PLR=1) [W]
  REAL(r64) :: S12RuntimeFraction =0.0   ! Stage 1&2 Run time fraction (overlaps with stage1 run time)
  REAL(r64) :: S12EvapCondPumpElecPower  ! Stage 1&2 Evaporative condenser pump electric power input [W]
  REAL(r64) :: S12EvapWaterConsumpRate   ! Stage 1&2 Evap condenser water consumption rate [m3/s]
  REAL(r64) :: S12CrankcaseHeaterPower   ! Stage 1&2 Report variable for average crankcase heater power [W]
  REAL(r64) :: S2PLR                     ! Stage 2   Ratio of actual sensible cooling load to
                                         !           steady-state sensible cooling capacity
  REAL(r64) :: MinAirHumRat = 0.0        ! minimum of the inlet air humidity ratio and the outlet air humidity ratio
  REAL(r64) :: TSat                      ! calculation to avoid calling psych routines twice

          ! FLOW

! First time SimDXCoil is called, get the input for all the DX coils (condensing units)
IF (GetCoilsInputFlag) THEN
  CALL GetDXCoils
  GetCoilsInputFlag = .FALSE. ! Set GetInputFlag false so you don't get coil inputs again
END IF

!  find correct DX Coil
IF (CompIndex == 0) THEN
  DXCoilNum = FindItemInList(CompName,DXCoil%Name,NumDXCoils)
  IF (DXCoilNum == 0) THEN
    CALL ShowFatalError('DX Coil not found='//TRIM(CompName))
  ENDIF
  CompIndex=DXCoilNum
ELSE
  DXCoilNum=CompIndex
  IF (DXCoilNum > NumDXCoils .or. DXCoilNum < 1) THEN
    CALL ShowFatalError('SimDXCoilMultiMode: Invalid CompIndex passed='//  &
                        TRIM(TrimSigDigits(DXCoilNum))// &
                        ', Number of DX Coils='//TRIM(TrimSigDigits(NumDXCoils))//  &
                        ', Coil name='//TRIM(CompName))
  ENDIF
  IF (CheckEquipName(DXCoilNum)) THEN
    IF ((CompName /='') .and. (CompName /= DXCoil(DXCoilNum)%Name)) THEN
      CALL ShowFatalError('SimDXCoilMultiMode: Invalid CompIndex passed='//  &
                          TRIM(TrimSigDigits(DXCoilNum))// &
                          ', Coil name='//TRIM(CompName)//', stored Coil Name for that index='//  &
                          TRIM(DXCoil(DXCoilNum)%Name))
    ENDIF
    CheckEquipName(DXCoilNum)=.false.
  ENDIF
ENDIF

CurDXCoilNum = DXCoilNum

! Initialize the DX coil unit
CALL InitDXCoil(DXCoilNum)

! Select the correct unit type
SELECT CASE(DXCoil(DXCoilNum)%DXCoilType_Num)

  CASE (CoilDX_CoolingTwoStageWHumControl)

    ! Initialize local variables
    S1RuntimeFraction            = 0.0
    S1OutletAirEnthalpy          = DXCoil(DXCoilNum)%InletAirEnthalpy
    S1OutletAirHumRat            = DXCoil(DXCoilNum)%InletAirHumRat
    S1OutletAirTemp              = DXCoil(DXCoilNum)%InletAirTemp
    S1ElecCoolingPower           = 0.0
    S1TotalCoolingEnergyRate     = 0.0
    S1SensCoolingEnergyRate      = 0.0
    S1LatCoolingEnergyRate       = 0.0
    S1CrankcaseHeaterPower       = 0.0
    S1EvapWaterConsumpRate       = 0.0
    S1EvapCondPumpElecPower      = 0.0

    S12RuntimeFraction           = 0.0
    S12OutletAirEnthalpy         = DXCoil(DXCoilNum)%InletAirEnthalpy
    S12OutletAirHumRat           = DXCoil(DXCoilNum)%InletAirHumRat
    S12OutletAirTemp             = DXCoil(DXCoilNum)%InletAirTemp
    S12ElecCoolingPower          = 0.0
    S12TotalCoolingEnergyRate    = 0.0
    S12SensCoolingEnergyRate     = 0.0
    S12LatCoolingEnergyRate      = 0.0
    S12CrankcaseHeaterPower      = 0.0
    S12EvapWaterConsumpRate      = 0.0
    S12EvapCondPumpElecPower     = 0.0

    DXCoil(DXCoilNum)%DehumidificationMode = DehumidMode
    IF (DehumidMode .GT. DXCoil(DXCoilNum)%NumDehumidModes) THEN
      CALL ShowFatalError(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//&
        '" - Requested enhanced dehumidification mode not available.')
    END IF

    ! If a single-stage coil OR If part load is zero,
    ! run stage 1 at zero part load to set leaving conditions
    IF ((DXCoil(DXCoilNum)%NumCapacityStages .EQ. 1) .OR. (PartLoadRatio .LE. 0.0)) THEN
      ! Run stage 1 at its part load
      PerfMode = DehumidMode*2 + 1
      CALL CalcDoe2DXCoil(DXCoilNum,On,FirstHVACIteration,PartLoadRatio,FanOpMode,PerfMode)
      S1PLR = PartLoadRatio
      S2PLR = 0.0
    ELSE
      ! If a two-stage coil
      ! Run stage 1 at full load
      PerfMode = DehumidMode*2 + 1
      CALL CalcDoe2DXCoil(DXCoilNum,On,FirstHVACIteration,1.0d0, FanOpMode,PerfMode)
      S1SensCoolingEnergyRate = DXCoil(DXCoilNum)%SensCoolingEnergyRate
      IF (S1SensCoolingEnergyRate > 0.0) THEN
        S1PLR = PartLoadRatio
      ELSE
        S1PLR = 0.0
      END IF
      ! Run stage 1+2 at full load
      IF (DXCoil(DXCoilNum)%NumCapacityStages .GE. 2) THEN
        PerfMode = DehumidMode*2 + 2
        CALL CalcDoe2DXCoil(DXCoilNum,On,FirstHVACIteration,1.0d0,FanOpMode,PerfMode)
        S12SensCoolingEnergyRate = DXCoil(DXCoilNum)%SensCoolingEnergyRate
        S12ElecCoolFullLoadPower = DXCoil(DXCoilNum)%ElecCoolingPower
      END IF

      ! Determine run-time fractions for each stage based on sensible capacities
      !   Relationships:
      !     Stage 1   PLR1=   Load/Cap1
      !     Stage1+2  PLR12=  Load/Cap12
      !     Stage 2   PLR2=   (Load-Cap1)/(Cap2)
      !     PLR = Load/(Cap1+Cap2)
      !     Load= PLR*(Cap1+Cap2)
      !     PLR1= Min(1,(PLR*(Cap1+Cap2)/Cap1))
      !     PLR2= Min(1,((PLR*(Cap1+Cap2)-Cap1)/Cap2))

      IF (S1SensCoolingEnergyRate > 0.0) THEN
        S1PLR = PartLoadRatio*S12SensCoolingEnergyRate/S1SensCoolingEnergyRate
      ELSE
        S1PLR = 0.0
      END IF
      S1PLR = MIN(1.0d0,S1PLR)
      S1PLR = MAX(0.0d0,S1PLR)
      IF ((S12SensCoolingEnergyRate-S1SensCoolingEnergyRate) > 0.0) THEN
        S2PLR = (PartLoadRatio*S12SensCoolingEnergyRate-S1SensCoolingEnergyRate)/ &
                (S12SensCoolingEnergyRate-S1SensCoolingEnergyRate)
      ELSE
        S2PLR = 0.0
      END IF
      S2PLR = MIN(1.0d0,S2PLR)
      S2PLR = MAX(0.0d0,S2PLR)

      ! Run stage 1 at its part load
      PerfMode = DehumidMode*2 + 1
      CALL CalcDoe2DXCoil(DXCoilNum,On,FirstHVACIteration,S1PLR,FanOpMode,PerfMode)
    ENDIF
    ! For stage-1 only operation, all outputs are set by CalcDoe2DXCoil.
    ! No further adjustments are necessary.

    ! Run stage 2 if needed and available
    IF ((S2PLR > 0.0) .AND. (DXCoil(DXCoilNum)%NumCapacityStages .GE. 2)) THEN
      ! Store stage 1 outputs
      S1RuntimeFraction            = DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction
      S1OutletAirEnthalpy          = DXCoil(DXCoilNum)%OutletAirEnthalpy
      S1OutletAirHumRat            = DXCoil(DXCoilNum)%OutletAirHumRat
      S1OutletAirTemp              = DXCoil(DXCoilNum)%OutletAirTemp
      S1ElecCoolingPower           = DXCoil(DXCoilNum)%ElecCoolingPower
      S1TotalCoolingEnergyRate     = DXCoil(DXCoilNum)%TotalCoolingEnergyRate
      S1SensCoolingEnergyRate      = DXCoil(DXCoilNum)%SensCoolingEnergyRate
      S1LatCoolingEnergyRate       = DXCoil(DXCoilNum)%LatCoolingEnergyRate
      S1CrankcaseHeaterPower       = DXCoil(DXCoilNum)%CrankcaseHeaterPower
      S1EvapWaterConsumpRate       = DXCoil(DXCoilNum)%EvapWaterConsumpRate
      S1EvapCondPumpElecPower      = DXCoil(DXCoilNum)%EvapCondPumpElecPower

      ! Save first stage full load outlet conditions to pass to heat recovery
      S1FFullLoadOutAirTemp        = DXCoilFullLoadOutAirTemp(DXCoilNum)
      S1FullLoadOutAirHumRat       = DXCoilFullLoadOutAirHumRat(DXCoilNum)

      ! Run stage 1+2 at its part load
      PerfMode = DehumidMode*2 + 2
      CALL CalcDoe2DXCoil(DXCoilNum,On,FirstHVACIteration,S2PLR,FanOpMode,PerfMode)
      S12RuntimeFraction            = DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction
      S12OutletAirEnthalpy          = DXCoil(DXCoilNum)%OutletAirEnthalpy
      S12OutletAirHumRat            = DXCoil(DXCoilNum)%OutletAirHumRat
      S12OutletAirTemp              = DXCoil(DXCoilNum)%OutletAirTemp
      S12ElecCoolingPower           = DXCoil(DXCoilNum)%ElecCoolingPower
      S12TotalCoolingEnergyRate     = DXCoil(DXCoilNum)%TotalCoolingEnergyRate
      S12SensCoolingEnergyRate      = DXCoil(DXCoilNum)%SensCoolingEnergyRate
      S12LatCoolingEnergyRate       = DXCoil(DXCoilNum)%LatCoolingEnergyRate
      S12CrankcaseHeaterPower       = DXCoil(DXCoilNum)%CrankcaseHeaterPower
      S12EvapWaterConsumpRate       = DXCoil(DXCoilNum)%EvapWaterConsumpRate
      S12EvapCondPumpElecPower      = DXCoil(DXCoilNum)%EvapCondPumpElecPower

      ! Determine combined performance
      DXCoil(DXCoilNum)%OutletAirEnthalpy = (1.d0-S12RuntimeFraction)*S1OutletAirEnthalpy &
                                             +S12RuntimeFraction*S12OutletAirEnthalpy
      DXCoil(DXCoilNum)%OutletAirHumRat = (1.d0-S12RuntimeFraction)*S1OutletAirHumRat &
                                            +S12RuntimeFraction*S12OutletAirHumRat
      DXCoil(DXCoilNum)%OutletAirTemp =   &
         PsyTdbFnHW(DXCoil(DXCoilNum)%OutletAirEnthalpy,DXCoil(DXCoilNum)%OutletAirHumRat,RoutineName)
      ! Check for saturation error and modify temperature at constant enthalpy
      IF (DXCoil(DXCoilNum)%CondenserInletNodeNum(PerfMode) /= 0) THEN
        TSat=PsyTsatFnHPb(DXCoil(DXCoilNum)%OutletAirEnthalpy,  &
           Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(PerfMode))%Press,RoutineName)
        IF(DXCoil(DXCoilNum)%OutletAirTemp .LT. TSat) THEN
          DXCoil(DXCoilNum)%OutletAirTemp = TSat
        ENDIF
        DXCoil(DXCoilNum)%OutletAirHumRat  = PsyWFnTdbH(DXCoil(DXCoilNum)%OutletAirTemp,  &
           DXCoil(DXCoilNum)%OutletAirEnthalpy,RoutineName)
      ELSE
        TSat=PsyTsatFnHPb(DXCoil(DXCoilNum)%OutletAirEnthalpy,OutBaroPress,RoutineName)
        IF(DXCoil(DXCoilNum)%OutletAirTemp .LT. TSat) THEN
          DXCoil(DXCoilNum)%OutletAirTemp = TSat
        ENDIF
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!      IF(DXCoil(DXCoilNum)%OutletAirTemp .LT. PsyTsatFnHPb(DXCoil(DXCoilNum)%OutletAirEnthalpy, &
!                 Node(DXCoil(DXCoilNum)%AirInNode)%Press)) THEN
!        DXCoil(DXCoilNum)%OutletAirTemp = PsyTsatFnHPb(DXCoil(DXCoilNum)%OutletAirEnthalpy, &
!                 Node(DXCoil(DXCoilNum)%AirInNode)%Press)
        DXCoil(DXCoilNum)%OutletAirHumRat  = PsyWFnTdbH(DXCoil(DXCoilNum)%OutletAirTemp,  &
           DXCoil(DXCoilNum)%OutletAirEnthalpy,RoutineName)
      END IF

!      DXCoil(DXCoilNum)%ElecCoolingPower = (1-S12RuntimeFraction)*S1ElecCoolingPower &
!                                             +S12RuntimeFraction*S12ElecCoolingPower
      !  S12ElecCoolingPower overstates S1 portion of power, because it is also adjust by S12PLR
      !  So, must make an adjustment for S12ElecCoolingPower/S12ElecCoolFullLoadPower
      !  when subtracting off S1ElecCoolingPower
      IF(S12ElecCoolFullLoadPower .GT. 0.0)THEN
        DXCoil(DXCoilNum)%ElecCoolingPower = S1RuntimeFraction*S1ElecCoolingPower &
         +S12RuntimeFraction*(S12ElecCoolingPower-S1ElecCoolingPower*S12ElecCoolingPower/S12ElecCoolFullLoadPower)
      ELSE
        DXCoil(DXCoilNum)%ElecCoolingPower = 0.0
      END IF

      DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction = S1RuntimeFraction

      AirMassFlow = DXCoil(DXCoilNum)%InletAirMassFlowRate
      DXCoil(DXCoilNum)%TotalCoolingEnergyRate = AirMassFlow * (DXCoil(DXCoilNum)%InletAirEnthalpy - &
                                              DXCoil(DXCoilNum)%OutletAirEnthalpy)
      MinAirHumRat = MIN(DXCoil(DXCoilNum)%InletAirHumRat,DXCoil(DXCoilNum)%OutletAirHumRat)
      DXCoil(DXCoilNum)%SensCoolingEnergyRate = AirMassFlow * &
              (PsyHFnTdbW(DXCoil(DXCoilNum)%InletAirTemp,MinAirHumRat,RoutineName) - &
               PsyHFnTdbW(DXCoil(DXCoilNum)%OutletAirTemp,MinAirHumRat,RoutineName))
      !  Don't let sensible capacity be greater than total capacity
      IF (DXCoil(DXCoilNum)%SensCoolingEnergyRate .GT. DXCoil(DXCoilNum)%TotalCoolingEnergyRate) THEN
        DXCoil(DXCoilNum)%SensCoolingEnergyRate = DXCoil(DXCoilNum)%TotalCoolingEnergyRate
      END IF

      DXCoil(DXCoilNum)%LatCoolingEnergyRate = DXCoil(DXCoilNum)%TotalCoolingEnergyRate - DXCoil(DXCoilNum)%SensCoolingEnergyRate

      DXCoil(DXCoilNum)%EvapWaterConsumpRate = (1.d0-S12RuntimeFraction)*S1EvapWaterConsumpRate &
                                                 +S12RuntimeFraction*S12EvapWaterConsumpRate
      DXCoil(DXCoilNum)%EvapCondPumpElecPower = (1.d0-S12RuntimeFraction)*S1EvapCondPumpElecPower &
                                                  +S12RuntimeFraction*S12EvapCondPumpElecPower


      ! Stage 1 runtime sets the crankcase heater power
      DXCoil(DXCoilNum)%CrankcaseHeaterPower = S1CrankcaseHeaterPower

      DXCoilOutletTemp(DXCoilNum)   = DXCoil(DXCoilNum)%OutletAirTemp
      DXCoilOutletHumRat(DXCoilNum) = DXCoil(DXCoilNum)%OutletAirHumRat

!     calculate average full load outlet conditions for second stage operation
      DXCoilFullLoadOutAirTemp(DXCoilNum) = (1.0d0 - S2PLR)*S1FFullLoadOutAirTemp + S2PLR*DXCoilFullLoadOutAirTemp(DXCoilNum)
      DXCoilFullLoadOutAirHumRat(DXCoilNum) = (1.0d0 - S2PLR)*S1FullLoadOutAirHumRat + S2PLR*DXCoilFullLoadOutAirHumRat(DXCoilNum)

    ENDIF ! End if stage 2 is operating

!   set the part load ratio and heat reclaim capacity for use by desuperheater heating coils
    DXCoil(DXCoilNum)%PartLoadRatio            = S1PLR
    DXCoilPartLoadRatio(DXCoilNum)             = S1PLR

!   Calculation for heat reclaim needs to be corrected to use compressor power (not including condenser fan power)
    HeatReclaimDXCoil(DXCoilNum)%AvailCapacity = DXCoil(DXCoilNum)%TotalCoolingEnergyRate + DXCoil(DXCoilNum)%ElecCoolingPower

    DXCoil(DXCoilNum)%CoolingCoilStg2RuntimeFrac = S12RuntimeFraction

!   Calculate basin heater power
    CALL CalcBasinHeaterPowerForMultiModeDXCoil(DXCoilNum, DehumidMode)

  CASE DEFAULT
    CALL ShowSevereError('Error detected in DX Coil='//TRIM(CompName))
    CALL ShowContinueError('Invalid DX Coil Type='//TRIM(DXCoil(DXCoilNum)%DXCoilType))
    CALL ShowFatalError('Preceding condition causes termination.')


END SELECT

! Update the unit outlet nodes
CALL UpdateDXCoil(DXCoilNum)

! Report the result of the simulation
CALL ReportDXCoil(DXCoilNum)


RETURN

END SUBROUTINE SimDXCoilMultiMode


SUBROUTINE GetDXCoils

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Fred Buhl
          !       DATE WRITTEN   May 2000
          !       MODIFIED       Don Shirey, Aug/Sept 2000, Feb/Oct 2001, Sept 2003, Jan/July 2004
          !                      Feb 2005 M. J. Witte, GARD Analytics, Inc.
          !                        Add new coil type COIL:DX:MultiMode:CoolingEmpirical:
          !                      May 2005, Rich Raustad, FSEC, Added COIL:DX:HeatPumpWaterHeater
          !                      June 2007 L. Gu, FSEC
          !                      Added new coil type COIL:DX:MULTISPEED:COOLING and COIL:DX:MULTISPEED:HEATING
          !                      April 2010, Chandan Sharma, FSEC, added basin heater inputs
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! Obtains input data for DX coils and stores it in DX coil data structure

          ! METHODOLOGY EMPLOYED:
          ! Uses "Get" routines to read in data.

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor,        ONLY: GetNumObjectsFound, GetObjectItem, GetObjectItemNum, VerifyName, SameString,GetObjectDefMaxArgs
  USE CurveManager,          ONLY: GetCurveIndex, GetCurveType, CurveValue, SetCurveOutputMinMaxValues
  USE BranchNodeConnections, ONLY: TestCompSet
  USE NodeInputManager,      ONLY: GetOnlySingleNode
  USE DataSizing,            ONLY: AutoSize
  USE General,               ONLY: TrimSigDigits
  USE WaterManager,          ONLY: SetupTankDemandComponent, SetupTankSupplyComponent
  USE OutAirNodeManager,     ONLY: CheckOutAirNodeNumber
  USE ScheduleManager,       ONLY: GetScheduleIndex
  USE DataGlobals,           ONLY: AnyEnergyManagementSystemInModel, emsCallFromComponentGetInput
  USE EMSManager,            ONLY: ManageEMS

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
          ! na

          ! SUBROUTINE PARAMETER DEFINITIONS:
  CHARACTER(len=*), PARAMETER :: Blank = ' '
  CHARACTER(len=*), PARAMETER :: RoutineName='GetDXCoils: ' ! include trailing blank space

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
INTEGER :: DXCoilIndex          ! loop index
INTEGER :: DXCoilNum            ! current DX coil number
INTEGER :: NumAlphas            ! Number of alphas in input
INTEGER :: NumNumbers           ! Number of numeric items in input
CHARACTER(Len=MaxNameLength), ALLOCATABLE, DIMENSION(:) :: Alphas2         ! Alpha input items for object
REAL(r64), ALLOCATABLE, DIMENSION(:) :: Numbers2           ! Numeric input items for object
CHARACTER(Len=MaxNameLength), ALLOCATABLE, DIMENSION(:) :: cAlphaFields2   ! Alpha field names
CHARACTER(Len=MaxNameLength), ALLOCATABLE, DIMENSION(:) :: cNumericFields2 ! Numeric field names
LOGICAL, ALLOCATABLE, DIMENSION(:)   :: lAlphaBlanks2      ! Logical array, alpha field input BLANK = .true.
LOGICAL, ALLOCATABLE, DIMENSION(:)   :: lNumericBlanks2    ! Logical array, numeric field input BLANK = .true.
INTEGER :: NumAlphas2           ! Number of alphas in input for performance object
INTEGER :: NumNumbers2          ! Number of numeric items in input for performance object
INTEGER :: IOStatus             ! Input status returned from GetObjectItem
LOGICAL :: IsNotOK              ! Flag to verify name
LOGICAL :: IsBlank              ! Flag for blank name
LOGICAL :: ErrorsFound=.false.  ! Set to true if errors in input, fatal at end of routine
INTEGER :: DXHPWaterHeaterCoilNum ! Loop index for 1,NumDXHeatPumpWaterHeaterCoils
INTEGER :: CapacityStageNum     ! Loop index for 1,Number of capacity stages
INTEGER :: DehumidModeNum       ! Loop index for 1,Number of enhanced dehumidification modes
INTEGER :: PerfModeNum          ! Performance mode index
INTEGER :: PerfObjectNum        ! Item number for performance object
INTEGER :: AlphaIndex           ! Index for current alpha field
CHARACTER(len=MaxNameLength) :: CurrentModuleObject     ! Object type for getting and error messages
CHARACTER(len=MaxNameLength) :: PerfObjectType   ! Performance object type for getting and error messages
CHARACTER(len=MaxNameLength) :: PerfObjectName   ! Performance object name for getting and error messages
REAL(r64) :: InletAirTemp       ! Used to pass proper inlet air temp to HPWH DX coil performance curves
REAL(r64) :: InletWaterTemp     ! Used to pass proper inlet water temp to HPWH DX coil performance curves
REAL(r64) :: HeatCapFTemp       ! Used to verify HPWH DX coil heating capacity (function of temp) performance curve
REAL(r64) :: HeatCOPFTemp       ! Used to verify HPWH DX coil heating COP (function of temp) performance curve
REAL(r64) :: HeatCapFAirFlow    ! Used to verify HPWH DX coil heating capacity (function of air flow) performance curve
REAL(r64) :: HeatCOPFAirFlow    ! Used to verify HPWH DX coil heating COP (function of air flow) performance curve
REAL(r64) :: HeatCapFWaterFlow  ! Used to verify HPWH DX coil heating capacity (function of water flow) performance curve
REAL(r64) :: HeatCOPFWaterFlow  ! Used to verify HPWH DX coil heating COP (function of water flow) performance curve
INTEGER   :: I                  ! Index of speeds
REAL(r64) :: CurveVal           ! Used to verify modifier curves equal 1 at rated conditions
CHARACTER(Len=MaxNameLength), ALLOCATABLE, DIMENSION(:) :: Alphas         ! Alpha input items for object
CHARACTER(Len=MaxNameLength), ALLOCATABLE, DIMENSION(:) :: cAlphaFields   ! Alpha field names
CHARACTER(Len=MaxNameLength), ALLOCATABLE, DIMENSION(:) :: cNumericFields ! Numeric field names
REAL(r64), ALLOCATABLE, DIMENSION(:) :: Numbers           ! Numeric input items for object
LOGICAL, ALLOCATABLE, DIMENSION(:)   :: lAlphaBlanks      ! Logical array, alpha field input BLANK = .true.
LOGICAL, ALLOCATABLE, DIMENSION(:)   :: lNumericBlanks    ! Logical array, numeric field input BLANK = .true.
INTEGER   :: MaxNumbers=0          ! Maximum number of numeric input fields
INTEGER   :: MaxAlphas=0        ! Maximum number of alpha input fields
INTEGER   :: TotalArgs=0        ! Total number of alpha and numeric arguments (max) for a
                                !   certain object in the input file
REAL(r64) :: MinCurveVal        ! used for testing PLF curve output
REAL(r64) :: MinCurvePLR        ! used for testing PLF curve output
REAL(r64) :: MaxCurveVal        ! used for testing PLF curve output
REAL(r64) :: MaxCurvePLR        ! used for testing PLF curve output
REAL(r64) :: CurveInput         ! index used for testing PLF curve output

! find number of each type of DX coil and calculate the total number
NumDoe2DXCoils = GetNumObjectsFound('Coil:Cooling:DX:SingleSpeed')
NumDXHeatingCoils = GetNumObjectsFound('Coil:Heating:DX:SingleSpeed')
NumDXMulSpeedCoils = GetNumObjectsFound('Coil:Cooling:DX:TwoSpeed')
NumDXMulModeCoils = GetNumObjectsFound('Coil:Cooling:DX:TwoStageWithHumidityControlMode')
NumDXHeatPumpWaterHeaterCoils = GetNumObjectsFound('Coil:WaterHeating:AirToWaterHeatPump')
NumDXMulSpeedCoolCoils = GetNumObjectsFound('Coil:Cooling:DX:MultiSpeed')
NumDXMulSpeedHeatCoils = GetNumObjectsFound('Coil:Heating:DX:MultiSpeed')
NumVRFCoolingCoils = GetNumObjectsFound(cAllCoilTypes(CoilVRF_Cooling))
NumVRFHeatingCoils = GetNumObjectsFound(cAllCoilTypes(CoilVRF_Heating))
DXCoilHPSimNum2=GetNumObjectsFound('ZoneHVAC:HPSim')
DXCoilHPSimNum=GetNumObjectsFound('Coil:HPSim') !RS: Debugging: Moving to an air loop instead of packaged (10/31/14)

NumDXCoils = NumDoe2DXCoils + NumDXHeatingCoils + NumDXMulSpeedCoils + NumDXMulModeCoils + NumDXHeatPumpWaterHeaterCoils &
             + NumDXMulSpeedCoolCoils + NumDXMulSpeedHeatCoils + NumVRFCoolingCoils + NumVRFHeatingCoils + DXCoilHPSimNum &
             + DXCoilHPSimNum2

! Determine max number of alpha and numeric arguments for all objects being read, in order to allocate local arrays
CALL GetObjectDefMaxArgs('Coil:Cooling:DX:SingleSpeed',TotalArgs,NumAlphas,NumNumbers)
MaxNumbers=NumNumbers
MaxAlphas=NumAlphas
CALL GetObjectDefMaxArgs('Coil:Heating:DX:SingleSpeed',TotalArgs,NumAlphas,NumNumbers)
MaxNumbers=MAX(MaxNumbers,NumNumbers)
MaxAlphas=MAX(MaxAlphas,NumAlphas)
CALL GetObjectDefMaxArgs('Coil:Cooling:DX:TwoSpeed',TotalArgs,NumAlphas,NumNumbers)
MaxNumbers=MAX(MaxNumbers,NumNumbers)
MaxAlphas=MAX(MaxAlphas,NumAlphas)
CALL GetObjectDefMaxArgs('Coil:Cooling:DX:TwoStageWithHumidityControlMode',TotalArgs,NumAlphas,NumNumbers)
MaxNumbers=MAX(MaxNumbers,NumNumbers)
MaxAlphas=MAX(MaxAlphas,NumAlphas)
CALL GetObjectDefMaxArgs('Coil:WaterHeating:AirToWaterHeatPump',TotalArgs,NumAlphas,NumNumbers)
MaxNumbers=MAX(MaxNumbers,NumNumbers)
MaxAlphas=MAX(MaxAlphas,NumAlphas)
CALL GetObjectDefMaxArgs('Coil:Cooling:DX:MultiSpeed',TotalArgs,NumAlphas,NumNumbers)
MaxNumbers=MAX(MaxNumbers,NumNumbers)
MaxAlphas=MAX(MaxAlphas,NumAlphas)
CALL GetObjectDefMaxArgs('Coil:Heating:DX:MultiSpeed',TotalArgs,NumAlphas,NumNumbers)
MaxNumbers=MAX(MaxNumbers,NumNumbers)
MaxAlphas=MAX(MaxAlphas,NumAlphas)
CALL GetObjectDefMaxArgs(cAllCoilTypes(CoilVRF_Cooling),TotalArgs,NumAlphas,NumNumbers)
MaxNumbers=MAX(MaxNumbers,NumNumbers)
MaxAlphas=MAX(MaxAlphas,NumAlphas)
CALL GetObjectDefMaxArgs(cAllCoilTypes(CoilVRF_Heating),TotalArgs,NumAlphas,NumNumbers)
MaxNumbers=MAX(MaxNumbers,NumNumbers)
MaxAlphas=MAX(MaxAlphas,NumAlphas)
CALL GetObjectDefMaxArgs('CoilPerformance:DX:Cooling',TotalArgs,NumAlphas,NumNumbers)
MaxNumbers=MAX(MaxNumbers,NumNumbers)
MaxAlphas=MAX(MaxAlphas,NumAlphas)


ALLOCATE(Alphas(MaxAlphas))
Alphas=' '
ALLOCATE(cAlphaFields(MaxAlphas))
cAlphaFields=' '
ALLOCATE(cNumericFields(MaxNumbers))
cNumericFields=' '
ALLOCATE(Numbers(MaxNumbers))
Numbers=0.0d0
ALLOCATE(lAlphaBlanks(MaxAlphas))
lAlphaBlanks=.TRUE.
ALLOCATE(lNumericBlanks(MaxNumbers))
lNumericBlanks=.TRUE.

ALLOCATE(Alphas2(MaxAlphas))
Alphas2=' '
ALLOCATE(cAlphaFields2(MaxAlphas))
cAlphaFields2=' '
ALLOCATE(cNumericFields2(MaxNumbers))
cNumericFields2=' '
ALLOCATE(Numbers2(MaxNumbers))
Numbers2=0.0d0
ALLOCATE(lAlphaBlanks2(MaxAlphas))
lAlphaBlanks2=.TRUE.
ALLOCATE(lNumericBlanks2(MaxNumbers))
lNumericBlanks2=.TRUE.

! allocate the data structure

! Derived types
ALLOCATE(DXCoil(NumDXCoils))
ALLOCATE(HeatReclaimDXCoil(NumDXCoils))
ALLOCATE(CheckEquipName(NumDXCoils))
CheckEquipName=.true.


! Module level variable arrays
ALLOCATE(DXCoilOutletTemp(NumDXCoils))
ALLOCATE(DXCoilOutletHumRat(NumDXCoils))
ALLOCATE(DXCoilPartLoadRatio(NumDXCoils))
ALLOCATE(DXCoilFanOpMode(NumDXCoils))
ALLOCATE(DXCoilFullLoadOutAirTemp(NumDXCoils))
ALLOCATE(DXCoilFullLoadOutAirHumRat(NumDXCoils))
ALLOCATE(DXCoilTotalCooling(NumDXCoils))
ALLOCATE(DXCoilTotalHeating(NumDXCoils))
ALLOCATE(DXCoilCoolInletAirWBTemp(NumDXCoils))
ALLOCATE(DXCoilHeatInletAirDBTemp(NumDXCoils))
ALLOCATE(DXCoilHeatInletAirWBTemp(NumDXCoils))
!
! initialize the module level arrays
!
DXCoilOutletTemp           = 0.0
DXCoilOutletHumRat         = 0.0
DXCoilPartLoadRatio        = 0.0
DXCoilFanOpMode            = 0
DXCoilFullLoadOutAirTemp   = 0.0
DXCoilFullLoadOutAirHumRat = 0.0

! initialize the coil counter
DXCoilNum = 0

! Loop over the Doe2 DX Coils and get & load the data
CurrentModuleObject='Coil:Cooling:DX:SingleSpeed'
DO DXCoilIndex = 1, NumDoe2DXCoils

  CALL GetObjectItem(TRIM(CurrentModuleObject),DXCoilIndex,Alphas,NumAlphas,Numbers,NumNumbers,IOStatus, &
                     NumBlank=lNumericBlanks,AlphaBlank=lAlphaBlanks, &
                     AlphaFieldNames=cAlphaFields,NumericFieldNames=cNumericFields)

  DXCoilNum = DXCoilNum+1
  IsNotOK=.FALSE.
  IsBlank=.FALSE.
  CALL VerifyName(Alphas(1),DXCoil%Name,DXCoilNum-1,IsNotOK,IsBlank,TRIM(CurrentModuleObject)//' Name')
  IF (IsNotOK) THEN
    ErrorsFound=.true.
    IF (IsBlank) Alphas(1)='xxxxx'
  ENDIF
  DXCoil(DXCoilNum)%Name = Alphas(1)
! Initialize DataHeatBalance heat reclaim variable name for use by heat reclaim coils
  HeatReclaimDXCoil(DXCoilNum)%Name = DXCoil(DXCoilNum)%Name
  HeatReclaimDXCoil(DXCoilNum)%SourceType = TRIM(CurrentModuleObject)
  DXCoil(DXCoilNum)%DXCoilType = TRIM(CurrentModuleObject)
  DXCoil(DXCoilNum)%DXCoilType_Num = CoilDX_CoolingSingleSpeed
  DXCoil(DXCoilNum)%Schedule = Alphas(2)
  DXCoil(DXCoilNum)%SchedPtr = GetScheduleIndex(Alphas(2))  ! convert schedule name to pointer
  IF (DXCoil(DXCoilNum)%SchedPtr .EQ. 0) THEN
    IF (lAlphaBlanks(2)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(2))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...'//TRIM(cAlphaFields(2))//'="'//TRIM(Alphas(2))//'".')
    END IF
    ErrorsFound=.TRUE.
  END IF
  DXCoil(DXCoilNum)%RatedTotCap(1) = Numbers(1)
  DXCoil(DXCoilNum)%RatedSHR(1)    = Numbers(2)
  DXCoil(DXCoilNum)%RatedCOP(1)    = Numbers(3)
  IF (DXCoil(DXCoilNum)%RatedCOP(1) .LE. 0.0d0) THEN
     CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
     CALL ShowContinueError('...'//TRIM(cNumericFields(3))//' must be > 0.0,'//  &
        ' entered value=['//trim(TrimSigDigits(Numbers(3),2))//'].')
     ErrorsFound = .TRUE.
  END IF

  DXCoil(DXCoilNum)%RatedAirVolFlowRate(1) = Numbers(4)
  DXCoil(DXCoilNum)%FanPowerPerEvapAirFlowRate(1) = Numbers(5)

  DXCoil(DXCoilNum)%AirInNode = &
               GetOnlySingleNode(Alphas(3),ErrorsFound,TRIM(CurrentModuleObject),Alphas(1), &
               NodeType_Air,NodeConnectionType_Inlet,1,ObjectIsNotParent)

  DXCoil(DXCoilNum)%AirOutNode = &
               GetOnlySingleNode(Alphas(4),ErrorsFound,TRIM(CurrentModuleObject),Alphas(1), &
               NodeType_Air,NodeConnectionType_Outlet,1,ObjectIsNotParent)

  CALL TestCompSet(TRIM(CurrentModuleObject),Alphas(1),Alphas(3),Alphas(4),'Air Nodes')

  DXCoil(DXCoilNum)%CCapFTemp(1) = GetCurveIndex(Alphas(5)) ! convert curve name to number
  IF (DXCoil(DXCoilNum)%CCapFTemp(1) .EQ. 0) THEN
    IF (lAlphaBlanks(5)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(5))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(5))//'="'//TRIM(Alphas(5))//'".')
    END IF
    ErrorsFound = .TRUE.
  ELSE
    ! Verify Curve Object, only legal type is BiQuadratic
    SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%CCapFTemp(1)))

    CASE('BIQUADRATIC')
      DXCoil(DXCoilNum)%TotCapTempModFacCurveType(1)=Biquadratic

    CASE DEFAULT
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(5))//' type for this object = '// &
                           TRIM(GetCurveType(DXCoil(DXCoilNum)%CCapFTemp(1))))
      CALL ShowContinueError('Curve type must be Biquadratic.')
      ErrorsFound=.TRUE.
    END SELECT
  END IF

  DXCoil(DXCoilNum)%CCapFFlow(1) = GetCurveIndex(Alphas(6)) ! convert curve name to number
  IF (DXCoil(DXCoilNum)%CCapFFlow(1) .EQ. 0) THEN
    IF (lAlphaBlanks(6)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(6))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(6))//'="'//TRIM(Alphas(6))//'".')
    END IF
    ErrorsFound = .TRUE.
  ELSE
    ! Verify Curve Object, only legal type is Quadratic
    SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%CCapFFlow(1)))
    CASE('QUADRATIC')

    CASE('CUBIC')

    CASE DEFAULT
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(6))//' type for this object = '// &
                           TRIM(GetCurveType(DXCoil(DXCoilNum)%CCapFFlow(1))))
      CALL ShowContinueError('Curve type must be Quadratic or Cubic.')
      ErrorsFound=.TRUE.
    END SELECT
  END IF

  DXCoil(DXCoilNum)%EIRFTemp(1) = GetCurveIndex(Alphas(7)) ! convert curve name to number
  IF (DXCoil(DXCoilNum)%EIRFTemp(1) .EQ. 0) THEN
    IF (lAlphaBlanks(7)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(7))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(7))//'="'//TRIM(Alphas(7))//'".')
    END IF
    ErrorsFound = .TRUE.
  ELSE
    ! Verify Curve Object, only legal type is Biquadratic
    SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%EIRFTemp(1)))

    CASE('BIQUADRATIC')
      DXCoil(DXCoilNum)%EIRTempModFacCurveType(1)=Biquadratic

    CASE DEFAULT
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(7))//' type for this object = '// &
                           TRIM(GetCurveType(DXCoil(DXCoilNum)%EIRFTemp(1))))
      CALL ShowContinueError('Curve type must be Biquadratic.')
      ErrorsFound=.TRUE.
    END SELECT
  END IF

  DXCoil(DXCoilNum)%EIRFFlow(1) = GetCurveIndex(Alphas(8)) ! convert curve name to number
  IF (DXCoil(DXCoilNum)%EIRFFlow(1) .EQ. 0) THEN
    IF (lAlphaBlanks(8)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(8))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(8))//'="'//TRIM(Alphas(8))//'".')
    END IF
    ErrorsFound = .TRUE.
  ELSE
    ! Verify Curve Object, only legal type is Quadratic
    SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%EIRFFlow(1)))

    CASE('QUADRATIC')

    CASE('CUBIC')

    CASE DEFAULT
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(8))//' type for this object = '// &
                           TRIM(GetCurveType(DXCoil(DXCoilNum)%EIRFFlow(1))))
      CALL ShowContinueError('Curve type must be Quadratic or Cubic.')
      ErrorsFound=.TRUE.
    END SELECT
  END IF

  DXCoil(DXCoilNum)%PLFFPLR(1) = GetCurveIndex(Alphas(9)) ! convert curve name to number
  IF (DXCoil(DXCoilNum)%PLFFPLR(1) .EQ. 0) THEN
    IF (lAlphaBlanks(9)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(9))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(9))//'="'//TRIM(Alphas(9))//'".')
    END IF
    ErrorsFound = .TRUE.
  ELSE
    ! Verify Curve Object, only legal types are Quadratic or Cubic
    SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%PLFFPLR(1)))

    CASE('QUADRATIC')

    CASE('CUBIC')

    CASE DEFAULT
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(9))//' type for this object = '// &
                            TRIM(GetCurveType(DXCoil(DXCoilNum)%PLFFPLR(1))))
      CALL ShowContinueError('Curve type must be Quadratic or Cubic.')
      ErrorsFound=.TRUE.
    END SELECT

    IF(.NOT. ErrorsFound)THEN
!     Test PLF curve minimum and maximum. Cap if less than 0.7 or greater than 1.0.
      MinCurveVal = 999.0d0
      MaxCurveVal = -999.0d0
      CurveInput = 0.0d0
      DO WHILE (CurveInput <= 1.0d0)
        CurveVal = CurveValue(DXCoil(DXCoilNum)%PLFFPLR(1),CurveInput)
        IF(CurveVal .LT. MinCurveVal)THEN
          MinCurveVal = CurveVal
          MinCurvePLR = CurveInput
        END IF
        IF(CurveVal .GT. MaxCurveVal)THEN
          MaxCurveVal = CurveVal
          MaxCurvePLR = CurveInput
        END IF
        CurveInput=CurveInput+0.01d0
      END DO
      IF(MinCurveVal .LT. 0.7d0)THEN
        CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...'//TRIM(cAlphaFields(9))//'="'//TRIM(Alphas(9))//'" has out of range values.')
        CALL ShowContinueError('...Curve minimum must be >= 0.7, '// &
                        'curve min at PLR = '//TRIM(TrimSigDigits(MinCurvePLR,2))//' is '//TRIM(TrimSigDigits(MinCurveVal,3)))
        CALL ShowContinueError('...Setting curve minimum to 0.7 and simulation continues.')
        CALL SetCurveOutputMinMaxValues(DXCoil(DXCoilNum)%PLFFPLR(1),ErrorsFound,CurveMin=0.7d0)
      END IF

      IF(MaxCurveVal .GT. 1.0d0)THEN
        CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...'//TRIM(cAlphaFields(9))//' = '//TRIM(Alphas(9))//' has out of range value.')
        CALL ShowContinueError('...Curve maximum must be <= 1.0, '// &
                        'curve max at PLR = '//TRIM(TrimSigDigits(MaxCurvePLR,2))//' is '//TRIM(TrimSigDigits(MaxCurveVal,3)))
        CALL ShowContinueError('...Setting curve maximum to 1.0 and simulation continues.')
        CALL SetCurveOutputMinMaxValues(DXCoil(DXCoilNum)%PLFFPLR(1),ErrorsFound,CurveMax=1.0d0)
      END IF

    END IF
  END IF

  DXCoil(DXCoilNum)%Twet_Rated(1)            = Numbers(6)
  DXCoil(DXCoilNum)%Gamma_Rated(1)           = Numbers(7)
  DXCoil(DXCoilNum)%MaxONOFFCyclesperHour(1) = Numbers(8)
  DXCoil(DXCoilNum)%LatentCapacityTimeConstant(1) = Numbers(9)

  ! Numbers (6) through (9) must all be greater than zero to use the latent capacity degradation model
  IF ((Numbers(6) .GT. 0.0d0 .OR. Numbers (7) .GT. 0.0d0 .OR. Numbers (8) .GT. 0.0d0 .OR. Numbers (9) .GT. 0.0d0) &
    .AND. (Numbers(6) .LE. 0.0d0 .OR. Numbers (7) .LE. 0.0d0 .OR. Numbers (8) .LE. 0.0d0 .OR. Numbers (9) .LE. 0.0d0)) THEN
       CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'":')
       CALL ShowContinueError('...At least one of the four input parameters for the latent capacity degradation model')
       CALL ShowContinueError('...is set to zero. Therefore, the latent degradation model will not be used for this simulation.')
  END IF

! outdoor condenser node
  IF (lAlphaBlanks(10)) THEN
    DXCoil(DXCoilNum)%CondenserInletNodeNum(1) = 0
  ELSE
    DXCoil(DXCoilNum)%CondenserInletNodeNum(1) = &
       GetOnlySingleNode(Alphas(10),ErrorsFound,TRIM(CurrentModuleObject),DXCoil(DXCoilNum)%Name, &
                         NodeType_Air,NodeConnectionType_OutsideAirReference,1,ObjectIsNotParent)

    IF (.not. CheckOutAirNodeNumber(DXCoil(DXCoilNum)%CondenserInletNodeNum(1))) THEN
      CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", may be invalid')
      CALL ShowContinueError(TRIM(cAlphaFields(10))//'="'//TRIM(Alphas(10))// &
                               '", node does not appear in an OutdoorAir:NodeList or as an OutdoorAir:Node.')
      CALL ShowContinueError('This node needs to be included in an air system or the coil model will not be valid' &
                             //', and the simulation continues')
    END IF
  ENDIF

  IF ((SameString(Alphas(11),'AirCooled')) .OR. lAlphaBlanks(11)) THEN
    DXCoil(DXCoilNum)%CondenserType(1) = AirCooled
  ELSEIF (SameString(Alphas(11),'EvaporativelyCooled')) THEN
    DXCoil(DXCoilNum)%CondenserType(1) = EvapCooled
    DXCoil(DXCoilNum)%ReportEvapCondVars = .TRUE.
  ELSE
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//TRIM(cAlphaFields(11))//'="'//TRIM(Alphas(11))//'":')
    CALL ShowContinueError('...must be AirCooled or EvaporativelyCooled.')
    ErrorsFound = .TRUE.
  END IF

  DXCoil(DXCoilNum)%EvapCondEffect(1) = Numbers(10)
  IF (DXCoil(DXCoilNum)%EvapCondEffect(1) .LT. 0.0d0 .OR. DXCoil(DXCoilNum)%EvapCondEffect(1) .GT. 1.0d0) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//trim(cNumericFields(10))//' cannot be < 0.0 or > 1.0.')
    CALL ShowContinueError('...entered value=['//trim(TrimSigDigits(Numbers(10),2))//'].')
    ErrorsFound = .TRUE.
  END IF

  DXCoil(DXCoilNum)%EvapCondAirFlow(1) = Numbers(11)
  IF (DXCoil(DXCoilNum)%EvapCondAirFlow(1) .LT. 0.0d0 .AND. DXCoil(DXCoilNum)%EvapCondAirFlow(1) /= AutoSize) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//trim(cNumericFields(11))//' cannot be < 0.0.')
    CALL ShowContinueError('...entered value=['//trim(TrimSigDigits(Numbers(11),2))//'].')
    ErrorsFound = .TRUE.
  END IF

  DXCoil(DXCoilNum)%EvapCondPumpElecNomPower(1) = Numbers(12)
  IF (DXCoil(DXCoilNum)%EvapCondPumpElecNomPower(1) .LT. 0.0d0 .AND. DXCoil(DXCoilNum)%EvapCondAirFlow(1) /= AutoSize) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//trim(cNumericFields(12))//' cannot be < 0.0.')
    CALL ShowContinueError('...entered value=['//trim(TrimSigDigits(Numbers(12),2))//'].')
    ErrorsFound = .TRUE.
  END IF

  !Set crankcase heater capacity
  DXCoil(DXCoilNum)%CrankcaseHeaterCapacity = Numbers(13)
  IF (DXCoil(DXCoilNum)%CrankcaseHeaterCapacity .LT. 0.0d0) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//trim(cNumericFields(13))//' cannot be < 0.0.')
    CALL ShowContinueError('...entered value=['//trim(TrimSigDigits(Numbers(13),2))//'].')
    ErrorsFound = .TRUE.
  END IF

  !Set crankcase heater cutout temperature
  DXCoil(DXCoilNum)%MaxOATCrankcaseHeater = Numbers(14)

  IF (DXCoil(DXCoilNum)%RatedCOP(1) .GT. 0.0d0) THEN
     DXCoil(DXCoilNum)%RatedEIR(1) = 1.d0 / DXCoil(DXCoilNum)%RatedCOP(1)
  END IF

  ! Get Water System tank connections
  !  A12, \field Name of Water Storage Tank for Supply
  DXCoil(DXCoilNum)%EvapWaterSupplyName = Alphas(12)
  IF (lAlphaBlanks(12)) THEN
    DXCoil(DXCoilNum)%EvapWaterSupplyMode = WaterSupplyFromMains
  ELSE
    DXCoil(DXCoilNum)%EvapWaterSupplyMode = WaterSupplyFromTank
    CALL SetupTankDemandComponent(DXCoil(DXCoilNum)%Name,TRIM(CurrentModuleObject), &
                 DXCoil(DXCoilNum)%EvapWaterSupplyName, ErrorsFound, DXCoil(DXCoilNum)%EvapWaterSupTankID, &
                 DXCoil(DXCoilNum)%EvapWaterTankDemandARRID )
  ENDIF

  !A13; \field Name of Water Storage Tank for Condensate Collection
  DXCoil(DXCoilNum)%CondensateCollectName = Alphas(13)
  IF (lAlphaBlanks(13)) THEN
    DXCoil(DXCoilNum)%CondensateCollectMode = CondensateDiscarded
  ELSE
    DXCoil(DxCoilNum)%CondensateCollectMode = CondensateToTank
    CALL SetupTankSupplyComponent(DXCoil(DXCoilNum)%Name,TRIM(CurrentModuleObject), &
                 DXCoil(DXCoilNum)%CondensateCollectName, ErrorsFound, DXCoil(DXCoilNum)%CondensateTankID, &
                 DXCoil(DXCoilNum)%CondensateTankSupplyARRID )
  END IF

  !   Basin heater power as a function of temperature must be greater than or equal to 0
  DXCoil(DxCoilNum)%BasinHeaterPowerFTempDiff = Numbers(15)
  IF(Numbers(15) .LT. 0.0d0) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//trim(cNumericFields(15))//' must be >= 0.0.')
    CALL ShowContinueError('...entered value=['//trim(TrimSigDigits(Numbers(15),2))//'].')
    ErrorsFound = .TRUE.
  END IF

  DXCoil(DxCoilNum)%BasinHeaterSetPointTemp = Numbers(16)
  IF(DXCoil(DxCoilNum)%BasinHeaterPowerFTempDiff .GT. 0.0d0) THEN
    IF(NumNumbers .LT. 16) THEN
      DXCoil(DxCoilNum)%BasinHeaterSetPointTemp = 2.0d0
    ENDIF
    IF(DXCoil(DxCoilNum)%BasinHeaterSetPointTemp < 2.0d0) THEN
      CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", freeze possible')
      CALL ShowContinueError('...'//trim(cNumericFields(16))//' is < 2 {C}. Freezing could occur.')
      CALL ShowContinueError('...entered value=['//trim(TrimSigDigits(Numbers(16),2))//'].')
    END IF
  END IF

  IF(.NOT. lAlphaBlanks(14))THEN
    DXCoil(DxCoilNum)%BasinHeaterSchedulePtr   = GetScheduleIndex(Alphas(14))
    IF(DXCoil(DxCoilNum)%BasinHeaterSchedulePtr .EQ. 0)THEN
      CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(14))//'="'//TRIM(Alphas(14))//'".')
      CALL ShowContinueError('Basin heater will be available to operate throughout the simulation.')
    END IF
  END IF

END DO ! end of the Doe2 DX coil loop

IF (ErrorsFound) THEN
  CALL ShowFatalError(RoutineName//'Errors found in getting '//TRIM(CurrentModuleObject)//' input. '//&
                      'Preceding condition(s) causes termination.')
END IF

! Loop over the Multimode DX Coils and get & load the data
CurrentModuleObject='Coil:Cooling:DX:TwoStageWithHumidityControlMode'
DO DXCoilIndex = 1, NumDXMulModeCoils

  CALL GetObjectItem(TRIM(CurrentModuleObject),DXCoilIndex,Alphas,NumAlphas,Numbers,NumNumbers,IOStatus, &
                     NumBlank=lNumericBlanks,AlphaBlank=lAlphaBlanks, &
                     AlphaFieldNames=cAlphaFields,NumericFieldNames=cNumericFields)

  DXCoilNum = DXCoilNum+1
  IsNotOK=.FALSE.
  IsBlank=.FALSE.
  CALL VerifyName(Alphas(1),DXCoil%Name,DXCoilNum-1,IsNotOK,IsBlank,TRIM(CurrentModuleObject)//' Name')
  IF (IsNotOK) THEN
    ErrorsFound=.true.
    IF (IsBlank) Alphas(1)='xxxxx'
  ENDIF
  DXCoil(DXCoilNum)%Name = Alphas(1)
! Initialize DataHeatBalance heat reclaim variable name for use by heat reclaim coils
  HeatReclaimDXCoil(DXCoilNum)%Name = DXCoil(DXCoilNum)%Name
  HeatReclaimDXCoil(DXCoilNum)%SourceType = TRIM(CurrentModuleObject)
  DXCoil(DXCoilNum)%DXCoilType = TRIM(CurrentModuleObject)
  DXCoil(DXCoilNum)%DXCoilType_Num = CoilDX_CoolingTwoStageWHumControl
  DXCoil(DXCoilNum)%Schedule = Alphas(2)
  DXCoil(DXCoilNum)%SchedPtr = GetScheduleIndex(Alphas(2))  ! convert schedule name to pointer
  IF (DXCoil(DXCoilNum)%SchedPtr .EQ. 0) THEN
    IF (lAlphaBlanks(2)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(2))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...'//TRIM(cAlphaFields(2))//'="'//TRIM(Alphas(2))//'".')
    END IF
    ErrorsFound=.TRUE.
  END IF

  DXCoil(DXCoilNum)%AirInNode = &
               GetOnlySingleNode(Alphas(3),ErrorsFound,TRIM(CurrentModuleObject),Alphas(1), &
               NodeType_Air,NodeConnectionType_Inlet,1,ObjectIsNotParent)

  DXCoil(DXCoilNum)%AirOutNode = &
               GetOnlySingleNode(Alphas(4),ErrorsFound,TRIM(CurrentModuleObject),Alphas(1), &
               NodeType_Air,NodeConnectionType_Outlet,1,ObjectIsNotParent)

  CALL TestCompSet(TRIM(CurrentModuleObject),Alphas(1),Alphas(3),Alphas(4),'Air Nodes')

  !Set crankcase heater capacity
  DXCoil(DXCoilNum)%CrankcaseHeaterCapacity = Numbers(1)
  IF (DXCoil(DXCoilNum)%CrankcaseHeaterCapacity .LT. 0.0d0) THEN
     CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
     CALL ShowContinueError('...'//TRIM(cNumericFields(1))//' must be >= 0.0,'//  &
        ' entered value=['//trim(TrimSigDigits(Numbers(1),2))//'].')
    ErrorsFound = .TRUE.
  END IF

  !Set crankcase heater cutout temperature
  DXCoil(DXCoilNum)%MaxOATCrankcaseHeater = Numbers(2)

  !  Number of capacity stages
  DXCoil(DXCoilNum)%NumCapacityStages = Numbers(3)
  !  Check if requested number of capacity stages exceeds limits
  IF ((DXCoil(DXCoilNum)%NumCapacityStages .GT. MaxCapacityStages) .OR. (DXCoil(DXCoilNum)%NumCapacityStages .LT. 1)) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...illegal '//TRIM(cNumericFields(3))//' = '//TRIM(TrimSigDigits(DXCoil(DXCoilNum)%NumCapacityStages)))
    CALL ShowContinueError('...Valid range is 1 to '//TRIM(TrimSigDigits(MaxCapacityStages)))
    ErrorsFound=.TRUE.
  END IF

  !  Number of enhanced dehumidification modes
  DXCoil(DXCoilNum)%NumDehumidModes = Numbers(4)
  !  Check if requested number of enhanced dehumidification modes exceeds limits
  IF ((DXCoil(DXCoilNum)%NumDehumidModes .GT. MaxDehumidModes) .OR. (DXCoil(DXCoilNum)%NumDehumidModes .LT. 0)) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...illegal '//TRIM(cNumericFields(4))//' = '//TRIM(TrimSigDigits(DXCoil(DXCoilNum)%NumDehumidModes)))
    CALL ShowContinueError('...Valid range is 0 to '//TRIM(TrimSigDigits(MaxDehumidModes)))
    ErrorsFound=.TRUE.
  END IF

  !  Set starting alpha index for coil performance inputs
  AlphaIndex = 5

  !  Loop through capacity stages and dehumidification modes
  DO DehumidModeNum = 0, DXCoil(DXCoilNum)%NumDehumidModes
    DO CapacityStageNum = 1, DXCoil(DXCoilNum)%NumCapacityStages
      !  Check if sufficient number of fields entered
      IF ((AlphaIndex+1) .GT. NumAlphas) THEN
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...not enough remaining fields for specified Number of Operating Modes.')
        CALL ShowContinueError('...Need additional Coil Performance Object Type and Coil Performance Object Name fields.')
        ErrorsFound=.TRUE.
      ELSE
        PerfObjectType = TRIM(Alphas(AlphaIndex))
        PerfObjectName = Alphas(AlphaIndex+1)
        PerfModeNum    = DehumidModeNum*2 + CapacityStageNum
        DXCoil(DXCoilNum)%CoilPerformanceType(PerfModeNum) = PerfObjectType
        IF (SameString(PerfObjectType,'CoilPerformance:DX:Cooling')) THEN
          DXCoil(DXCoilNum)%CoilPerformanceType_Num(PerfModeNum) = CoilPerfDX_CoolByPassEmpirical
        ELSE
          CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
          CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(AlphaIndex))//'="'//TRIM(PerfObjectType)//'".')
          CALL ShowContinueError('Must be "CoilPerformance:DX:Cooling".')
          ErrorsFound=.TRUE.
        END IF
        DXCoil(DXCoilNum)%CoilPerformanceName(PerfModeNum) = PerfObjectName
        ! Get for CoilPerformance object
        PerfObjectNum = GetObjectItemNum(PerfObjectType,PerfObjectName)
        IF (PerfObjectNum > 0) THEN

          CALL GetObjectItem(PerfObjectType,PerfObjectNum, &
                             Alphas2,NumAlphas2,Numbers2,NumNumbers2,IOStatus, &
                             NumBlank=lNumericBlanks2,AlphaBlank=lAlphaBlanks2, &
                             AlphaFieldNames=cAlphaFields2,NumericFieldNames=cNumericFields2)

          DXCoil(DXCoilNum)%RatedTotCap(PerfModeNum) = Numbers2(1)
          DXCoil(DXCoilNum)%RatedSHR(PerfModeNum)    = Numbers2(2)
          DXCoil(DXCoilNum)%RatedCOP(PerfModeNum)    = Numbers2(3)
               ! Rated flow is immediately adjusted for bypass fraction if not autosized
          DXCoil(DXCoilNum)%BypassedFlowFrac(PerfModeNum)    = Numbers2(5)
          DXCoil(DXCoilNum)%RatedAirVolFlowRate(PerfModeNum) = Numbers2(4)
          IF (DXCoil(DXCoilNum)%RatedAirVolFlowRate(PerfModeNum) /= Autosize) THEN
            DXCoil(DXCoilNum)%RatedAirVolFlowRate(PerfModeNum)= &
               DXCoil(DXCoilNum)%RatedAirVolFlowRate(PerfModeNum) * (1.d0-DXCoil(DXCoilNum)%BypassedFlowFrac(PerfModeNum))
          ENDIF

          DXCoil(DXCoilNum)%CCapFTemp(PerfModeNum) = GetCurveIndex(Alphas2(2)) ! convert curve name to number
          IF (DXCoil(DXCoilNum)%CCapFTemp(PerfModeNum) .EQ. 0) THEN
            IF (lAlphaBlanks2(2)) THEN
              CALL ShowSevereError(RoutineName//TRIM(PerfObjectType)//'="'//TRIM(PerfObjectName)//'", invalid')
              CALL ShowContinueError('...required '//trim(cAlphaFields2(2))//' is blank.')
            ELSE
              CALL ShowSevereError(RoutineName//TRIM(PerfObjectType)//'="'//TRIM(PerfObjectName)//'", invalid')
              CALL ShowContinueError('...not found '//TRIM(cAlphaFields2(2))//'="'//TRIM(Alphas2(2))//'".')
            END IF
            ErrorsFound = .TRUE.
          ELSE
            ! Verify Curve Object, only legal type is BiQuadratic
            SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%CCapFTemp(PerfModeNum)))

            CASE('BIQUADRATIC')
              DXCoil(DXCoilNum)%TotCapTempModFacCurveType(PerfModeNum)=Biquadratic

            CASE DEFAULT
              CALL ShowSevereError(RoutineName//TRIM(PerfObjectType)//'="'//TRIM(PerfObjectName)//'", invalid')
              CALL ShowContinueError('...illegal '//TRIM(cAlphaFields2(2))//' type for this object = '// &
                                   TRIM(GetCurveType(DXCoil(DXCoilNum)%CCapFTemp(PerfModeNum))))
              CALL ShowContinueError('Curve type must be BiQuadratic.')
              ErrorsFound=.TRUE.
            END SELECT
          END IF

          DXCoil(DXCoilNum)%CCapFFlow(PerfModeNum) = GetCurveIndex(Alphas2(3)) ! convert curve name to number
          IF (DXCoil(DXCoilNum)%CCapFFlow(PerfModeNum) .EQ. 0) THEN
            IF (lAlphaBlanks2(3)) THEN
              CALL ShowSevereError(RoutineName//TRIM(PerfObjectType)//'="'//TRIM(PerfObjectName)//'", invalid')
              CALL ShowContinueError('...required '//trim(cAlphaFields2(3))//' is blank.')
            ELSE
              CALL ShowSevereError(RoutineName//TRIM(PerfObjectType)//'="'//TRIM(PerfObjectName)//'", invalid')
              CALL ShowContinueError('...not found '//TRIM(cAlphaFields2(3))//'="'//TRIM(Alphas2(3))//'".')
            END IF
            ErrorsFound = .TRUE.
          ELSE
            ! Verify Curve Object, only legal type is Quadratic
            SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%CCapFFlow(PerfModeNum)))

            CASE('QUADRATIC')

            CASE('CUBIC')

            CASE DEFAULT
              CALL ShowSevereError(RoutineName//TRIM(PerfObjectType)//'="'//TRIM(PerfObjectName)//'", invalid')
              CALL ShowContinueError('...illegal '//TRIM(cAlphaFields2(3))//' type for this object = '// &
                                   TRIM(GetCurveType(DXCoil(DXCoilNum)%CCapFFlow(PerfModeNum))))
              CALL ShowContinueError('Curve type must be Quadratic or Cubic.')
              ErrorsFound=.TRUE.
            END SELECT
          END IF

          DXCoil(DXCoilNum)%EIRFTemp(PerfModeNum) = GetCurveIndex(Alphas2(4)) ! convert curve name to number
          IF (DXCoil(DXCoilNum)%EIRFTemp(PerfModeNum) .EQ. 0) THEN
            IF (lAlphaBlanks2(4)) THEN
              CALL ShowSevereError(RoutineName//TRIM(PerfObjectType)//'="'//TRIM(PerfObjectName)//'", invalid')
              CALL ShowContinueError('...required '//trim(cAlphaFields2(4))//' is blank.')
            ELSE
              CALL ShowSevereError(RoutineName//TRIM(PerfObjectType)//'="'//TRIM(PerfObjectName)//'", invalid')
              CALL ShowContinueError('...not found '//TRIM(cAlphaFields2(4))//'="'//TRIM(Alphas2(4))//'".')
            END IF
            ErrorsFound = .TRUE.
          ELSE
            ! Verify Curve Object, only legal type is Biquadratic
            SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%EIRFTemp(PerfModeNum)))

            CASE('BIQUADRATIC')
              DXCoil(DXCoilNum)%EIRTempModFacCurveType(PerfModeNum)=Biquadratic

            CASE DEFAULT
              CALL ShowSevereError(RoutineName//TRIM(PerfObjectType)//'="'//TRIM(PerfObjectName)//'", invalid')
              CALL ShowContinueError('...illegal '//TRIM(cAlphaFields2(4))//' type for this object = '// &
                                   TRIM(GetCurveType(DXCoil(DXCoilNum)%EIRFTemp(PerfModeNum))))
              CALL ShowContinueError('Curve type must be BiQuadratic.')
              ErrorsFound=.TRUE.
            END SELECT
          END IF

          DXCoil(DXCoilNum)%EIRFFlow(PerfModeNum) = GetCurveIndex(Alphas2(5)) ! convert curve name to number
          IF (DXCoil(DXCoilNum)%EIRFFlow(PerfModeNum) .EQ. 0) THEN
            IF (lAlphaBlanks2(5)) THEN
              CALL ShowSevereError(RoutineName//TRIM(PerfObjectType)//'="'//TRIM(PerfObjectName)//'", invalid')
              CALL ShowContinueError('...required '//trim(cAlphaFields2(5))//' is blank.')
            ELSE
              CALL ShowSevereError(RoutineName//TRIM(PerfObjectType)//'="'//TRIM(PerfObjectName)//'", invalid')
              CALL ShowContinueError('...not found '//TRIM(cAlphaFields2(5))//'="'//TRIM(Alphas2(5))//'".')
            END IF
            ErrorsFound = .TRUE.
          ELSE
            ! Verify Curve Object, only legal type is Quadratic
            SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%EIRFFlow(PerfModeNum)))

            CASE('QUADRATIC')

            CASE('CUBIC')

            CASE DEFAULT
              CALL ShowSevereError(RoutineName//TRIM(PerfObjectType)//'="'//TRIM(PerfObjectName)//'", invalid')
              CALL ShowContinueError('...illegal '//TRIM(cAlphaFields2(5))//' type for this object = '// &
                                   TRIM(GetCurveType(DXCoil(DXCoilNum)%EIRFFlow(PerfModeNum))))
              CALL ShowContinueError('Curve type must be Quadratic or Cubic.')
              ErrorsFound=.TRUE.
            END SELECT
          END IF

          DXCoil(DXCoilNum)%PLFFPLR(PerfModeNum) = GetCurveIndex(Alphas2(6)) ! convert curve name to number
          IF (DXCoil(DXCoilNum)%PLFFPLR(PerfModeNum) .EQ. 0) THEN
            IF (lAlphaBlanks2(6)) THEN
              CALL ShowSevereError(RoutineName//TRIM(PerfObjectType)//'="'//TRIM(PerfObjectName)//'", invalid')
              CALL ShowContinueError('...required '//trim(cAlphaFields2(6))//' is blank.')
            ELSE
              CALL ShowSevereError(RoutineName//TRIM(PerfObjectType)//'="'//TRIM(PerfObjectName)//'", invalid')
              CALL ShowContinueError('...not found '//TRIM(cAlphaFields2(6))//'="'//TRIM(Alphas2(6))//'".')
            END IF
            ErrorsFound = .TRUE.
          ELSE
            ! Verify Curve Object, only legal types are Quadratic or Cubic
            SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%PLFFPLR(PerfModeNum)))

            CASE('QUADRATIC')

            CASE('CUBIC')

            CASE DEFAULT
              CALL ShowSevereError(RoutineName//TRIM(PerfObjectType)//'="'//TRIM(PerfObjectName)//'", invalid')
              CALL ShowContinueError('...illegal '//TRIM(cAlphaFields2(6))//' type for this object = '// &
                                   TRIM(GetCurveType(DXCoil(DXCoilNum)%PLFFPLR(PerfModeNum))))
              CALL ShowContinueError('Curve type must be Quadratic or Cubic.')
              ErrorsFound=.TRUE.
            END SELECT

            IF(.NOT. ErrorsFound)THEN
!             Test PLF curve minimum and maximum. Cap if less than 0.7 or greater than 1.0.
              MinCurveVal = 999.0d0
              MaxCurveVal = -999.0d0
              CurveInput = 0.0d0
              DO WHILE (CurveInput <= 1.0d0)
                CurveVal = CurveValue(DXCoil(DXCoilNum)%PLFFPLR(PerfModeNum),CurveInput)
                IF(CurveVal .LT. MinCurveVal)THEN
                  MinCurveVal = CurveVal
                  MinCurvePLR = CurveInput
                END IF
                IF(CurveVal .GT. MaxCurveVal)THEN
                  MaxCurveVal = CurveVal
                  MaxCurvePLR = CurveInput
                END IF
                CurveInput=CurveInput+0.01d0
              END DO
              IF(MinCurveVal .LT. 0.7d0)THEN
                CALL ShowWarningError(RoutineName//TRIM(PerfObjectType)//'="'//TRIM(PerfObjectName)//'", invalid')
                CALL ShowContinueError('...'//TRIM(cAlphaFields2(6))//' = '//TRIM(Alphas2(6))//' has out of range value.')
                CALL ShowContinueError('...Curve minimum must be >= 0.7, '// &
                        'curve min at PLR = '//TRIM(TrimSigDigits(MinCurvePLR,2))//' is '//TRIM(TrimSigDigits(MinCurveVal,3)))
                CALL ShowContinueError('...Setting curve minimum to 0.7 and simulation continues.')
                CALL SetCurveOutputMinMaxValues(DXCoil(DXCoilNum)%PLFFPLR(PerfModeNum),ErrorsFound,CurveMin=0.7d0)
              END IF

              IF(MaxCurveVal .GT. 1.0d0)THEN
                CALL ShowWarningError(RoutineName//TRIM(PerfObjectType)//'="'//TRIM(PerfObjectName)//'", invalid')
                CALL ShowContinueError('...'//TRIM(cAlphaFields2(6))//' = '//TRIM(Alphas2(6))//' has out of range value.')
                CALL ShowContinueError('...Curve maximum must be <= 1.0, '// &
                        'curve max at PLR = '//TRIM(TrimSigDigits(MaxCurvePLR,2))//' is '//TRIM(TrimSigDigits(MaxCurveVal,3)))
                CALL ShowContinueError('...Setting curve maximum to 1.0 and simulation continues.')
                CALL SetCurveOutputMinMaxValues(DXCoil(DXCoilNum)%PLFFPLR(PerfModeNum),ErrorsFound,CurveMax=1.0d0)
              END IF

            END IF
          END IF

          DXCoil(DXCoilNum)%Twet_Rated(PerfModeNum)  = Numbers2(6)
          DXCoil(DXCoilNum)%Gamma_Rated(PerfModeNum) = Numbers2(7)
          DXCoil(DXCoilNum)%MaxONOFFCyclesperHour(PerfModeNum)      = Numbers2(8)
          DXCoil(DXCoilNum)%LatentCapacityTimeConstant(PerfModeNum) = Numbers2(9)
          ! Numbers2 (6) through (9) must all be greater than zero to use the latent capacity degradation model
          IF ((Numbers2(6) .GT. 0.0 .OR. Numbers2 (7) .GT. 0.0 .OR. Numbers2 (8) .GT. 0.0 .OR. Numbers2 (9) .GT. 0.0) &
              .AND. (Numbers2(6) .LE. 0.0 .OR. Numbers2 (7) .LE. 0.0 .OR. Numbers2 (8) .LE. 0.0 .OR. Numbers2 (9) .LE. 0.0)) THEN
            CALL ShowWarningError(RoutineName//TRIM(PerfObjectType)//'="'//TRIM(PerfObjectName)//'":')
            CALL ShowContinueError('...At least one of the four input parameters for the latent capacity degradation model')
            CALL ShowContinueError('...is set to zero. Therefore, the latent degradation model will not be used '//  &
               'for this simulation.')
          END IF

          ! outdoor condenser node
          IF (lAlphaBlanks2(7)) THEN
            DXCoil(DXCoilNum)%CondenserInletNodeNum(PerfModeNum) = 0
          ELSE
            DXCoil(DXCoilNum)%CondenserInletNodeNum(PerfModeNum) = &
               GetOnlySingleNode(Alphas2(7),ErrorsFound,PerfObjectType,PerfObjectName, &
                                 NodeType_Air,NodeConnectionType_OutsideAirReference,1,ObjectIsNotParent)
            IF (.not. CheckOutAirNodeNumber(DXCoil(DXCoilNum)%CondenserInletNodeNum(PerfModeNum))) THEN
              CALL ShowWarningError(RoutineName//TRIM(PerfObjectType)//'="'//TRIM(PerfObjectName)//'":')
              CALL ShowContinueError('may not be valid '//TRIM(cAlphaFields2(7))//'="'//TRIM(Alphas2(7))//'".')
              CALL ShowContinueError('node does not appear in an OutdoorAir:NodeList or as an OutdoorAir:Node.')
              CALL ShowContinueError('This node needs to be included in an air system or the coil model will not be valid' &
                             //', and the simulation continues')
            END IF
          ENDIF
          IF ((SameString(Alphas2(8),'AirCooled')) .OR. lAlphaBlanks2(8)) THEN
            DXCoil(DXCoilNum)%CondenserType(PerfModeNum) = AirCooled
          ELSEIF (SameString(Alphas2(8),'EvaporativelyCooled')) THEN
            DXCoil(DXCoilNum)%CondenserType(PerfModeNum) = EvapCooled
            DXCoil(DXCoilNum)%ReportEvapCondVars = .true.
          ELSE
            CALL ShowSevereError(RoutineName//TRIM(PerfObjectType)//'="'//TRIM(PerfObjectName)//'", invalid')
            CALL ShowContinueError('...'//TRIM(cAlphaFields2(8))//'="'//TRIM(Alphas2(8))//'":')
            CALL ShowContinueError('...must be AirCooled or EvaporativelyCooled.')
            ErrorsFound = .TRUE.
          END IF

          DXCoil(DXCoilNum)%EvapCondEffect(PerfModeNum) = Numbers2(10)
          IF (DXCoil(DXCoilNum)%EvapCondEffect(PerfModeNum) .LT. 0.0d0 .OR.   &
              DXCoil(DXCoilNum)%EvapCondEffect(PerfModeNum) .GT. 1.0d0) THEN
            CALL ShowSevereError(RoutineName//trim(PerfObjectType)//'="'//trim(PerfObjectName)//'", invalid')
            CALL ShowContinueError('...'//trim(cNumericFields2(10))//' cannot be < 0.0 or > 1.0.')
            CALL ShowContinueError('...entered value=['//trim(TrimSigDigits(Numbers2(10),2))//'].')
            ErrorsFound = .TRUE.
          END IF

          DXCoil(DXCoilNum)%EvapCondAirFlow(PerfModeNum) = Numbers2(11)
          IF (DXCoil(DXCoilNum)%EvapCondAirFlow(PerfModeNum) .LT. 0.0d0 .AND.  &
              DXCoil(DXCoilNum)%EvapCondAirFlow(PerfModeNum) /= AutoSize) THEN
            CALL ShowSevereError(RoutineName//trim(PerfObjectType)//'="'//trim(PerfObjectName)//'", invalid')
            CALL ShowContinueError('...'//trim(cNumericFields2(11))//' cannot be < 0.0.')
            CALL ShowContinueError('...entered value=['//trim(TrimSigDigits(Numbers2(11),2))//'].')
            ErrorsFound = .TRUE.
          END IF

          DXCoil(DXCoilNum)%EvapCondPumpElecNomPower(PerfModeNum) = Numbers2(12)
          IF (DXCoil(DXCoilNum)%EvapCondPumpElecNomPower(PerfModeNum) .LT. 0.0d0 .AND.  &
              DXCoil(DXCoilNum)%EvapCondAirFlow(PerfModeNum) /= AutoSize) THEN
            CALL ShowSevereError(RoutineName//trim(PerfObjectType)//'="'//trim(PerfObjectName)//'", invalid')
            CALL ShowContinueError('...'//trim(cNumericFields2(12))//' cannot be less than zero.')
            CALL ShowContinueError('...entered value=['//trim(TrimSigDigits(Numbers2(12),2))//'].')
            ErrorsFound = .TRUE.
          END IF

          DXCoil(DXCoilNum)%RatedEIR(PerfModeNum) = 1.0d0 / DXCoil(DXCoilNum)%RatedCOP(PerfModeNum)
        ELSE ! invalid performance object
          CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
          CALL ShowContinueError('... not found '//TRIM(PerfObjectType)//'="'//TRIM(PerfObjectName)//'".')
          ErrorsFound=.TRUE.
        END IF ! end of valid performance object check
        AlphaIndex = AlphaIndex + 2
      END IF ! end of sufficient number of fields entered check
    END DO ! End of multimode DX capacity stages loop
    ! Warn if inputs entered for unused capacity stages
    DO CapacityStageNum = (DXCoil(DXCoilNum)%NumCapacityStages+1), MaxCapacityStages
      IF ((AlphaIndex .LE. NumAlphas) .AND. &
         ((Alphas(AlphaIndex) .NE. Blank) .OR. (Alphas(AlphaIndex+1) .NE. Blank))) THEN
        CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'":')
        CALL ShowContinueError('...Capacity Stage '//TRIM(TrimSigDigits(CapacityStageNum))//' not active. Therefore,'//  &
            trim(cAlphaFields(AlphaIndex)))
        CALL ShowContinueError('... and '//trim(cAlphaFields(AlphaIndex+1))//' fields will be ignored.')
      END IF
      AlphaIndex = AlphaIndex + 2
    END DO ! End of unused capacity stages loop
  END DO ! End of multimode DX dehumidification modes loop

!  ! Warn if excess fields entered
!  IF (NumAlphas .GE. AlphaIndex .and. ANY(Alphas(AlphaIndex:) /= Blank)) THEN
!    CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'":')
!    CALL ShowContinueError('...too many remaining fields for specified Capacity Stages and Dehumidification Modes.')
!    CALL ShowContinueError('...Excess Coil Performance Object Type and Coil Performance Object Name fields will be ignored.')
!  ENDIF

  ! Get Water System tank connections
  !  A13, \field Name of Water Storage Tank for Supply
  DXCoil(DXCoilNum)%EvapWaterSupplyName = Alphas(13)
  IF (lAlphaBlanks(13)) THEN
    DXCoil(DXCoilNum)%EvapWaterSupplyMode = WaterSupplyFromMains
  ELSE
    DXCoil(DXCoilNum)%EvapWaterSupplyMode = WaterSupplyFromTank
    CALL SetupTankDemandComponent(DXCoil(DXCoilNum)%Name, TRIM(CurrentModuleObject), &
                 DXCoil(DXCoilNum)%EvapWaterSupplyName, ErrorsFound, DXCoil(DXCoilNum)%EvapWaterSupTankID, &
                 DXCoil(DXCoilNum)%EvapWaterTankDemandARRID )
  ENDIF

  !A14; \field Name of Water Storage Tank for Condensate Collection
  DXCoil(DXCoilNum)%CondensateCollectName = Alphas(14)
  IF (lAlphaBlanks(14)) THEN
    DXCoil(DXCoilNum)%CondensateCollectMode = CondensateDiscarded
  ELSE
    DXCoil(DxCoilNum)%CondensateCollectMode = CondensateToTank
    CALL SetupTankSupplyComponent(DXCoil(DXCoilNum)%Name, TRIM(CurrentModuleObject), &
                 DXCoil(DXCoilNum)%CondensateCollectName, ErrorsFound, DXCoil(DXCoilNum)%CondensateTankID, &
                 DXCoil(DXCoilNum)%CondensateTankSupplyARRID )
  ENDIF

  !Basin heater power as a function of temperature must be greater than or equal to 0
  DXCoil(DxCoilNum)%BasinHeaterPowerFTempDiff = Numbers(5)
  IF(Numbers(5) .LT. 0.0d0) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//trim(cNumericFields(5))//' must be >= 0.')
    CALL ShowContinueError('...entered value=['//trim(TrimSigDigits(Numbers(5),2))//'].')
    ErrorsFound = .TRUE.
  END IF

  DXCoil(DxCoilNum)%BasinHeaterSetPointTemp = Numbers(6)
  IF(DXCoil(DxCoilNum)%BasinHeaterPowerFTempDiff .GT. 0.0d0) THEN
    IF(NumNumbers .LT. 6) THEN
      DXCoil(DxCoilNum)%BasinHeaterSetPointTemp = 2.0d0
    ENDIF
    IF(DXCoil(DxCoilNum)%BasinHeaterSetPointTemp < 2.0d0) THEN
      CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", freeze possible')
      CALL ShowContinueError('...'//trim(cNumericFields(6))//' is < 2 {C}. Freezing could occur.')
      CALL ShowContinueError('...entered value=['//trim(TrimSigDigits(Numbers(6),2))//'].')
    END IF
  END IF

  IF(.NOT. lAlphaBlanks(15))THEN
    DXCoil(DxCoilNum)%BasinHeaterSchedulePtr   = GetScheduleIndex(Alphas(15))
    IF(DXCoil(DxCoilNum)%BasinHeaterSchedulePtr .EQ. 0)THEN
      CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(15))//'="'//TRIM(Alphas(15))//'".')
      CALL ShowContinueError('Basin heater will be available to operate throughout the simulation.')
    END IF
  END IF

END DO ! end of the Multimode DX coil loop

IF (ErrorsFound) THEN
  CALL ShowFatalError(RoutineName//'Errors found in getting '//TRIM(CurrentModuleObject)//' input.  '//&
                      'Preceding condition(s) causes termination.')
END IF

!************* Read Heat Pump (DX Heating Coil) Input **********
CurrentModuleObject='Coil:Heating:DX:SingleSpeed'
DO DXCoilIndex = 1,NumDXHeatingCoils

  DXCoilNum = DXCoilNum+1

  CALL GetObjectItem(TRIM(CurrentModuleObject),DXCoilIndex,Alphas,NumAlphas,Numbers,NumNumbers,IOStatus, &
                     NumBlank=lNumericBlanks,AlphaBlank=lAlphaBlanks, &
                     AlphaFieldNames=cAlphaFields,NumericFieldNames=cNumericFields)
  IsNotOK=.FALSE.
  IsBlank=.FALSE.
  CALL VerifyName(Alphas(1),DXCoil%Name,DXCoilNum-1,IsNotOK,IsBlank,TRIM(CurrentModuleObject)//' Name')
  IF (IsNotOK) THEN
    ErrorsFound=.true.
    IF (IsBlank) Alphas(1)='xxxxx'
  ENDIF
  DXCoil(DXCoilNum)%Name       = Alphas(1)
  DXCoil(DXCoilNum)%DXCoilType = TRIM(CurrentModuleObject)
  DXCoil(DXCoilNum)%DXCoilType_Num = CoilDX_HeatingEmpirical
  DXCoil(DXCoilNum)%Schedule = Alphas(2)
  DXCoil(DXCoilNum)%SchedPtr = GetScheduleIndex(Alphas(2))  ! convert schedule name to pointer
  IF (DXCoil(DXCoilNum)%SchedPtr .EQ. 0) THEN
    IF (lAlphaBlanks(2)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(2))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...'//TRIM(cAlphaFields(2))//'="'//TRIM(Alphas(2))//'".')
    END IF
    ErrorsFound=.TRUE.
  END IF

  DXCoil(DXCoilNum)%AirInNode = &
               GetOnlySingleNode(Alphas(3),ErrorsFound,TRIM(CurrentModuleObject),Alphas(1), &
               NodeType_Air,NodeConnectionType_Inlet,1,ObjectIsNotParent)

  DXCoil(DXCoilNum)%AirOutNode = &
               GetOnlySingleNode(Alphas(4),ErrorsFound,TRIM(CurrentModuleObject),Alphas(1), &
               NodeType_Air,NodeConnectionType_Outlet,1,ObjectIsNotParent)

  CALL TestCompSet(TRIM(CurrentModuleObject),Alphas(1),Alphas(3),Alphas(4),'Air Nodes')

  DXCoil(DXCoilNum)%CCapFTemp(1) = GetCurveIndex(Alphas(5)) ! convert curve name to number
  IF (DXCoil(DXCoilNum)%CCapFTemp(1) .EQ. 0) THEN
    IF (lAlphaBlanks(5)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(5))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(5))//'="'//TRIM(Alphas(5))//'".')
    END IF
    ErrorsFound = .TRUE.
  ELSE
    ! only legal types are Quadratic, Biquadratic and Cubic
    SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%CCapFTemp(1)))

    CASE('QUADRATIC')
      DXCoil(DXCoilNum)%TotCapTempModFacCurveType(1)=Quadratic

    CASE('BIQUADRATIC')
      DXCoil(DXCoilNum)%TotCapTempModFacCurveType(1)=Biquadratic

    CASE('CUBIC')
      DXCoil(DXCoilNum)%TotCapTempModFacCurveType(1)=Cubic

    CASE DEFAULT
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(5))//' type for this object = '// &
                           TRIM(GetCurveType(DXCoil(DXCoilNum)%CCapFTemp(1))))
      CALL ShowContinueError('Curve type must be Biquadratic, Quadratic or Cubic.')
      ErrorsFound=.TRUE.
    END SELECT
  END IF

  DXCoil(DXCoilNum)%CCapFFlow(1) = GetCurveIndex(Alphas(6)) ! convert curve name to number
  IF (DXCoil(DXCoilNum)%CCapFFlow(1) .EQ. 0) THEN
    IF (lAlphaBlanks(6)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(6))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(6))//'="'//TRIM(Alphas(6))//'".')
    END IF
    ErrorsFound = .TRUE.
  ELSE
    ! Verify Curve Object, only legal type is Quadratic
    SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%CCapFFlow(1)))

    CASE('QUADRATIC')

    CASE('CUBIC')

    CASE DEFAULT
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(6))//' type for this object = '// &
                           TRIM(GetCurveType(DXCoil(DXCoilNum)%CCapFFlow(1))))
      CALL ShowContinueError('Curve type must be Quadratic or Cubic.')
      ErrorsFound=.TRUE.
    END SELECT
  END IF

  DXCoil(DXCoilNum)%EIRFTemp(1) = GetCurveIndex(Alphas(7)) ! convert curve name to number
  IF (DXCoil(DXCoilNum)%EIRFTemp(1) .EQ. 0) THEN
    IF (lAlphaBlanks(7)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(7))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(7))//'="'//TRIM(Alphas(7))//'".')
    END IF
    ErrorsFound = .TRUE.
  ELSE
    ! only legal types are Quadratic, Biquadratic and Cubic
    SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%EIRFTemp(1)))

    CASE('QUADRATIC')
      DXCoil(DXCoilNum)%EIRTempModFacCurveType(1)=Quadratic

    CASE('BIQUADRATIC')
      DXCoil(DXCoilNum)%EIRTempModFacCurveType(1)=Biquadratic

    CASE('CUBIC')
      DXCoil(DXCoilNum)%EIRTempModFacCurveType(1)=Cubic

    CASE DEFAULT
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(7))//' type for this object = '// &
                           TRIM(GetCurveType(DXCoil(DXCoilNum)%EIRFTemp(1))))
      CALL ShowContinueError('Curve type must be Biquadratic, Quadratic or Cubic.')
      ErrorsFound=.TRUE.
    END SELECT
  END IF

  DXCoil(DXCoilNum)%EIRFFlow(1) = GetCurveIndex(Alphas(8)) ! convert curve name to number
  IF (DXCoil(DXCoilNum)%EIRFFlow(1) .EQ. 0) THEN
    IF (lAlphaBlanks(8)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(8))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(8))//'="'//TRIM(Alphas(8))//'".')
    END IF
    ErrorsFound = .TRUE.
  ELSE
    ! Verify Curve Object, only legal type is Quadratic or Cubic
    SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%EIRFFlow(1)))

    CASE('QUADRATIC')

    CASE('CUBIC')

    CASE DEFAULT
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(8))//' type for this object = '// &
                           TRIM(GetCurveType(DXCoil(DXCoilNum)%EIRFFlow(1))))
      CALL ShowContinueError('Curve type must be Quadratic or Cubic.')
      ErrorsFound=.TRUE.
    END SELECT
  END IF

  DXCoil(DXCoilNum)%PLFFPLR(1) = GetCurveIndex(Alphas(9)) ! convert curve name to number
  IF (DXCoil(DXCoilNum)%PLFFPLR(1) .EQ. 0) THEN
    IF (lAlphaBlanks(9)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(9))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(9))//'="'//TRIM(Alphas(9))//'".')
    END IF
    ErrorsFound = .TRUE.
  ELSE
    ! Verify Curve Object, only legal types are Quadratic or Cubic
    SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%PLFFPLR(1)))

    CASE('QUADRATIC')

    CASE('CUBIC')

    CASE DEFAULT
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(9))//' type for this object = '// &
                           TRIM(GetCurveType(DXCoil(DXCoilNum)%PLFFPLR(1))))
      CALL ShowContinueError('Curve type must be Quadratic or Cubic.')
      ErrorsFound=.TRUE.
    END SELECT

    IF(.NOT. ErrorsFound)THEN
!     Test PLF curve minimum and maximum. Cap if less than 0.7 or greater than 1.0.
      MinCurveVal = 999.0d0
      MaxCurveVal = -999.0d0
      CurveInput = 0.0d0
      DO WHILE (CurveInput <= 1.0d0)
        CurveVal = CurveValue(DXCoil(DXCoilNum)%PLFFPLR(1),CurveInput)
        IF(CurveVal .LT. MinCurveVal)THEN
          MinCurveVal = CurveVal
          MinCurvePLR = CurveInput
        END IF
        IF(CurveVal .GT. MaxCurveVal)THEN
          MaxCurveVal = CurveVal
          MaxCurvePLR = CurveInput
        END IF
        CurveInput=CurveInput+0.01d0
      END DO
      IF(MinCurveVal .LT. 0.7d0)THEN
        CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...'//TRIM(cAlphaFields(9))//' = '//TRIM(Alphas(9))//' has out of range value.')
        CALL ShowContinueError('...Curve minimum must be >= 0.7, '// &
                        'curve min at PLR = '//TRIM(TrimSigDigits(MinCurvePLR,2))//' is '//TRIM(TrimSigDigits(MinCurveVal,3)))
        CALL ShowContinueError('...Setting curve minimum to 0.7 and simulation continues.')
        CALL SetCurveOutputMinMaxValues(DXCoil(DXCoilNum)%PLFFPLR(1),ErrorsFound,CurveMin=0.7d0)
      END IF

      IF(MaxCurveVal .GT. 1.0d0)THEN
        CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...'//TRIM(cAlphaFields(9))//' = '//TRIM(Alphas(9))//' has out of range value.')
        CALL ShowContinueError('...Curve maximum must be <= 1.0, '// &
                        'curve max at PLR = '//TRIM(TrimSigDigits(MaxCurvePLR,2))//' is '//TRIM(TrimSigDigits(MaxCurveVal,3)))
        CALL ShowContinueError('...Setting curve maximum to 1.0 and simulation continues.')
        CALL SetCurveOutputMinMaxValues(DXCoil(DXCoilNum)%PLFFPLR(1),ErrorsFound,CurveMax=1.0d0)
      END IF

    END IF

  END IF

! Only required for reverse cycle heat pumps
  DXCoil(DXCoilNum)%DefrostEIRFT = GetCurveIndex(Alphas(10)) ! convert curve name to number
  IF (SameString(Alphas(11),'ReverseCycle')) THEN
    IF (DXCoil(DXCoilNum)%DefrostEIRFT .EQ. 0) THEN
      IF (lAlphaBlanks(10)) THEN
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
        CALL ShowContinueError('...required '//trim(cAlphaFields(10))//' is blank.')
        CALL ShowContinueError('...field is required because '//trim(cAlphaFields(11))//' is "ReverseCycle".')
      ELSE
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...not found '//TRIM(cAlphaFields(10))//'="'//TRIM(Alphas(10))//'".')
      END IF
      ErrorsFound = .TRUE.
    ELSE
      ! Verify Curve Object, only legal type is BiQuadratic
      SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%DefrostEIRFT))

        CASE('BIQUADRATIC')

        CASE DEFAULT
          CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
          CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(10))//' type for this object = '// &
                             TRIM(GetCurveType(DXCoil(DXCoilNum)%DefrostEIRFT)))
          CALL ShowContinueError('Curve type must be BiQuadratic.')
          ErrorsFound=.TRUE.
      END SELECT
    END IF
  END IF

  IF (SameString(Alphas(11),'ReverseCycle'))  DXCoil(DXCoilNum)%DefrostStrategy = ReverseCycle
  IF (SameString(Alphas(11),'Resistive')) DXCoil(DXCoilNum)%DefrostStrategy = Resistive
  IF (DXCoil(DXCoilNum)%DefrostStrategy .EQ.0) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(11))//'="'//TRIM(Alphas(11))//'".')
    CALL ShowContinueError('...valid values for this field are ReverseCycle or Resistive.')
    ErrorsFound = .TRUE.
  END IF

  IF (SameString(Alphas(12),'Timed'))     DXCoil(DXCoilNum)%DefrostControl = Timed
  IF (SameString(Alphas(12),'OnDemand')) DXCoil(DXCoilNum)%DefrostControl = OnDemand
  IF (DXCoil(DXCoilNum)%DefrostControl .EQ.0) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(12))//'="'//TRIM(Alphas(12))//'".')
    CALL ShowContinueError('...valid values for this field are Timed or OnDemand.')
    ErrorsFound = .TRUE.
  END IF

  DXCoil(DXCoilNum)%RatedSHR(1)    = 1.0d0
  DXCoil(DXCoilNum)%RatedTotCap(1) = Numbers(1)
  DXCoil(DXCoilNum)%RatedCOP(1)    = Numbers(2)
  DXCoil(DXCoilNum)%RatedAirVolFlowRate(1) = Numbers(3)
  DXCoil(DXCoilNum)%FanPowerPerEvapAirFlowRate(1) = Numbers(4)

  !Set minimum OAT for heat pump compressor operation
  DXCoil(DXCoilNum)%MinOATCompressor = Numbers(5)

  DXCoil(DXCoilNum)%OATempCompressorOn = Numbers(6)

  IF ( lNumericBlanks(6) .OR. lNumericBlanks(5)) THEN
    DXCoil(DXCoilNum)%OATempCompressorOnOffBlank = .TRUE.
  ELSE
    DXCoil(DXCoilNum)%OATempCompressorOnOffBlank = .FALSE.
  ENDIF

  IF (DXCoil(DXCoilNum)%OATempCompressorOn .LT. DXCoil(DXCoilNum)%MinOATCompressor) &
  DXCoil(DXCoilNum)%OATempCompressorOn = DXCoil(DXCoilNum)%MinOATCompressor

  !Set maximum outdoor temp for defrost to occur
  DXCoil(DXCoilNum)%MaxOATDefrost = Numbers(7)

  !Set crankcase heater capacity
  DXCoil(DXCoilNum)%CrankcaseHeaterCapacity = Numbers(8)
  IF (DXCoil(DXCoilNum)%CrankcaseHeaterCapacity .LT. 0.0d0) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//TRIM(cNumericFields(8))//' cannot be < 0.0.')
    CALL ShowContinueError('...entered value=['//trim(TrimSigDigits(Numbers(8),2))//'].')
    ErrorsFound = .TRUE.
  END IF

  !Set crankcase heater cutout temperature
  DXCoil(DXCoilNum)%MaxOATCrankcaseHeater = Numbers(9)

  !Set defrost time period
  DXCoil(DXCoilNum)%DefrostTime = Numbers(10)
  IF(DXCoil(DXCoilNum)%DefrostTime .EQ. 0.0d0 .AND. DXCoil(DXCoilNum)%DefrostControl .EQ. 1) THEN
    CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", ')
    CALL ShowContinueError('...'//TRIM(cNumericFields(10))//' = 0.0 for defrost control = TIMED.')
  END IF

  !Set defrost capacity (for resistive defrost)
  DXCoil(DXCoilNum)%DefrostCapacity = Numbers(11)
  IF(DXCoil(DXCoilNum)%DefrostCapacity .EQ. 0.0d0 .AND. DXCoil(DXCoilNum)%DefrostStrategy .EQ. 2) THEN
    CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", ')
    CALL ShowContinueError('...'//TRIM(cNumericFields(11))//' = 0.0 for defrost strategy = RESISTIVE.')
  END IF

  !Set Region number for calculating HSPF
  DXCoil(DXCoilNum)%RegionNum = Numbers(12)

  IF ( lNumericBlanks(12)) THEN
    DXCoil(DXCoilNum)%RegionNum = 4
  ENDIF

  DXCoil(DXCoilNum)%RatedEIR(1) = 1.d0 / DXCoil(DXCoilNum)%RatedCOP(1)

  !A13 is optional evaporator node name
  IF ( lAlphaBlanks(13) ) THEN
    DXCoil(DXCoilNum)%CondenserInletNodeNum(1) = 0
  ELSE
    DXCoil(DXCoilNum)%CondenserInletNodeNum(1) = &
       GetOnlySingleNode(Alphas(13),ErrorsFound,TRIM(CurrentModuleObject),DXCoil(DXCoilNum)%Name, &
                         NodeType_Air,NodeConnectionType_OutsideAirReference,1,ObjectIsNotParent)
     ! warn if not an outdoor node, but allow
    IF (.not. CheckOutAirNodeNumber(DXCoil(DXCoilNum)%CondenserInletNodeNum(1))) THEN
      CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", may be invalid')
      CALL ShowContinueError(TRIM(cAlphaFields(13))//'="'//TRIM(Alphas(13))//&
                               '", node does not appear in an OutdoorAir:NodeList or as an OutdoorAir:Node.')
      CALL ShowContinueError('This node needs to be included in an air system or the coil model will not be valid' &
                             //', and the simulation continues')
    END IF
  ENDIF

END DO ! end of the DX heating coil loop

IF (ErrorsFound) THEN
  CALL ShowFatalError(RoutineName//'Errors found in getting '//TRIM(CurrentModuleObject)//' input. '//&
                      'Preceding condition(s) causes termination.')
END IF

CurrentModuleObject='Coil:Cooling:DX:TwoSpeed'
DO DXCoilIndex = 1,NumDXMulSpeedCoils

  DXCoilNum = DXCoilNum+1

  CALL GetObjectItem(TRIM(CurrentModuleObject),DXCoilIndex,Alphas,NumAlphas,Numbers,NumNumbers,IOStatus, &
                     NumBlank=lNumericBlanks,AlphaBlank=lAlphaBlanks, &
                     AlphaFieldNames=cAlphaFields,NumericFieldNames=cNumericFields)

  IsNotOK=.FALSE.
  IsBlank=.FALSE.
  CALL VerifyName(Alphas(1),DXCoil%Name,DXCoilNum-1,IsNotOK,IsBlank,TRIM(CurrentModuleObject)//' Name')
  IF (IsNotOK) THEN
    ErrorsFound=.true.
    IF (IsBlank) Alphas(1)='xxxxx'
  ENDIF
  DXCoil(DXCoilNum)%Name = Alphas(1)
! Initialize DataHeatBalance heat reclaim variable name for use by heat reclaim coils
  HeatReclaimDXCoil(DXCoilNum)%Name = DXCoil(DXCoilNum)%Name
  HeatReclaimDXCoil(DXCoilNum)%SourceType = TRIM(CurrentModuleObject)
  DXCoil(DXCoilNum)%DXCoilType = TRIM(CurrentModuleObject)
  DXCoil(DXCoilNum)%DXCoilType_Num = CoilDX_CoolingTwoSpeed
  DXCoil(DXCoilNum)%Schedule = Alphas(2)
  DXCoil(DXCoilNum)%SchedPtr = GetScheduleIndex(Alphas(2))  ! convert schedule name to pointer
  IF (DXCoil(DXCoilNum)%SchedPtr .EQ. 0) THEN
    IF (lAlphaBlanks(2)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(2))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...'//TRIM(cAlphaFields(2))//'="'//TRIM(Alphas(2))//'".')
    END IF
    ErrorsFound=.TRUE.
  END IF
  DXCoil(DXCoilNum)%RatedTotCap(1) = Numbers(1)
  DXCoil(DXCoilNum)%RatedSHR(1)    = Numbers(2)
  DXCoil(DXCoilNum)%RatedCOP(1)    = Numbers(3)
  DXCoil(DXCoilNum)%RatedAirVolFlowRate(1) = Numbers(4)

  DXCoil(DXCoilNum)%AirInNode = &
               GetOnlySingleNode(Alphas(3),ErrorsFound,TRIM(CurrentModuleObject),Alphas(1), &
               NodeType_Air,NodeConnectionType_Inlet,1,ObjectIsNotParent)

  DXCoil(DXCoilNum)%AirOutNode = &
               GetOnlySingleNode(Alphas(4),ErrorsFound,TRIM(CurrentModuleObject),Alphas(1), &
               NodeType_Air,NodeConnectionType_Outlet,1,ObjectIsNotParent)

  CALL TestCompSet(TRIM(CurrentModuleObject),Alphas(1),Alphas(3),Alphas(4),'Air Nodes')

  DXCoil(DXCoilNum)%CCapFTemp(1) = GetCurveIndex(Alphas(5)) ! convert curve name to number
  IF (DXCoil(DXCoilNum)%CCapFTemp(1) .EQ. 0) THEN
    IF (lAlphaBlanks(5)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(5))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(5))//'="'//TRIM(Alphas(5))//'".')
    END IF
    ErrorsFound = .TRUE.
  ELSE
    ! Verify Curve Object, only legal type is BiQuadratic
    SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%CCapFTemp(1)))

    CASE('BIQUADRATIC')
      DXCoil(DXCoilNum)%TotCapTempModFacCurveType(1)=Biquadratic

    CASE DEFAULT
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(5))//' type for this object = '// &
                           TRIM(GetCurveType(DXCoil(DXCoilNum)%CCapFTemp(1))))
      CALL ShowContinueError('Curve type must be BiQuadratic.')
      ErrorsFound=.true.
    END SELECT
  END IF

  DXCoil(DXCoilNum)%CCapFFlow(1) = GetCurveIndex(Alphas(6)) ! convert curve name to number
  IF (DXCoil(DXCoilNum)%CCapFFlow(1) .EQ. 0) THEN
    IF (lAlphaBlanks(6)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(6))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(6))//'="'//TRIM(Alphas(6))//'".')
    END IF
    ErrorsFound = .TRUE.
  ELSE
    ! Verify Curve Object, only legal type is Quadratic
    SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%CCapFFlow(1)))

    CASE('QUADRATIC')

    CASE('CUBIC')

    CASE DEFAULT
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(6))//' type for this object = '// &
                            TRIM(GetCurveType(DXCoil(DXCoilNum)%CCapFFlow(1))))
      CALL ShowContinueError('Curve type must be Quadratic or Cubic.')
      ErrorsFound=.TRUE.
    END SELECT
  END IF

  DXCoil(DXCoilNum)%EIRFTemp(1) = GetCurveIndex(Alphas(7)) ! convert curve name to number
  IF (DXCoil(DXCoilNum)%EIRFTemp(1) .EQ. 0) THEN
    IF (lAlphaBlanks(7)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(7))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(7))//'="'//TRIM(Alphas(7))//'".')
    END IF
    ErrorsFound = .TRUE.
  ELSE
    ! Verify Curve Object, only legal type is Biquadratic
    SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%EIRFTemp(1)))

    CASE('BIQUADRATIC')
      DXCoil(DXCoilNum)%EIRTempModFacCurveType(1)=Biquadratic

    CASE DEFAULT
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(7))//' type for this object = '// &
                           TRIM(GetCurveType(DXCoil(DXCoilNum)%EIRFTemp(1))))
      CALL ShowContinueError('Curve type must be BiQuadratic.')
      ErrorsFound=.TRUE.
    END SELECT
  END IF

  DXCoil(DXCoilNum)%EIRFFlow(1) = GetCurveIndex(Alphas(8)) ! convert curve name to number
  IF (DXCoil(DXCoilNum)%EIRFFlow(1) .EQ. 0) THEN
    IF (lAlphaBlanks(8)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(8))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(8))//'="'//TRIM(Alphas(8))//'".')
    END IF
    ErrorsFound = .TRUE.
  ELSE
    ! Verify Curve Object, only legal type is Quadratic
    SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%EIRFFlow(1)))

    CASE('QUADRATIC')

    CASE('CUBIC')

    CASE DEFAULT
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(8))//' type for this object = '// &
                           TRIM(GetCurveType(DXCoil(DXCoilNum)%EIRFFlow(1))))
      CALL ShowContinueError('Curve type must be Quadratic or Cubic.')
      ErrorsFound=.TRUE.
    END SELECT
  END IF

  DXCoil(DXCoilNum)%PLFFPLR(1) = GetCurveIndex(Alphas(9)) ! convert curve name to number
  IF (DXCoil(DXCoilNum)%PLFFPLR(1) .EQ. 0) THEN
    IF (lAlphaBlanks(9)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(9))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(9))//'="'//TRIM(Alphas(9))//'".')
    END IF
    ErrorsFound = .TRUE.
  ELSE
    ! Verify Curve Object, only legal types are Quadratic or Cubic
    SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%PLFFPLR(1)))

    CASE('QUADRATIC')

    CASE('CUBIC')

    CASE DEFAULT
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(9))//' type for this object = '// &
                           TRIM(GetCurveType(DXCoil(DXCoilNum)%PLFFPLR(1))))
      CALL ShowContinueError('Curve type must be Quadratic or Cubic.')
      ErrorsFound=.TRUE.
    END SELECT

    IF(.NOT. ErrorsFound)THEN
!     Test PLF curve minimum and maximum. Cap if less than 0.7 or greater than 1.0.
      MinCurveVal = 999.0d0
      MaxCurveVal = -999.0d0
      CurveInput = 0.0d0
      DO WHILE (CurveInput <= 1.0d0)
        CurveVal = CurveValue(DXCoil(DXCoilNum)%PLFFPLR(1),CurveInput)
        IF(CurveVal .LT. MinCurveVal)THEN
          MinCurveVal = CurveVal
          MinCurvePLR = CurveInput
        END IF
        IF(CurveVal .GT. MaxCurveVal)THEN
          MaxCurveVal = CurveVal
          MaxCurvePLR = CurveInput
        END IF
        CurveInput=CurveInput+0.01d0
      END DO
      IF(MinCurveVal .LT. 0.7d0)THEN
        CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...'//TRIM(cAlphaFields(9))//' = '//TRIM(Alphas(9))//' has out of range value.')
        CALL ShowContinueError('...Curve minimum must be >= 0.7, '// &
                        'curve min at PLR = '//TRIM(TrimSigDigits(MinCurvePLR,2))//' is '//TRIM(TrimSigDigits(MinCurveVal,3)))
        CALL ShowContinueError('...Setting curve minimum to 0.7 and simulation continues.')
        CALL SetCurveOutputMinMaxValues(DXCoil(DXCoilNum)%PLFFPLR(1),ErrorsFound,CurveMin=0.7d0)
      END IF

      IF(MaxCurveVal .GT. 1.0d0)THEN
        CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...'//TRIM(cAlphaFields(9))//' = '//TRIM(Alphas(9))//' has out of range value.')
        CALL ShowContinueError('...Curve maximum must be <= 1.0, '// &
                        'curve max at PLR = '//TRIM(TrimSigDigits(MaxCurvePLR,2))//' is '//TRIM(TrimSigDigits(MaxCurveVal,3)))
        CALL ShowContinueError('...Setting curve maximum to 1.0 and simulation continues.')
        CALL SetCurveOutputMinMaxValues(DXCoil(DXCoilNum)%PLFFPLR(1),ErrorsFound,CurveMax=1.0d0)
      END IF

    END IF

  END IF

  DXCoil(DXCoilNum)%RatedEIR(1) = 1.d0 / DXCoil(DXCoilNum)%RatedCOP(1)

  DXCoil(DXCoilNum)%RatedTotCap2 = Numbers(5)
  DXCoil(DXCoilNum)%RatedSHR2    = Numbers(6)
  DXCoil(DXCoilNum)%RatedCOP2    = Numbers(7)
  DXCoil(DXCoilNum)%RatedAirVolFlowRate2 = Numbers(8)

  DXCoil(DXCoilNum)%CCapFTemp2 = GetCurveIndex(Alphas(10)) ! convert curve name to number
  IF (DXCoil(DXCoilNum)%CCapFTemp2 .EQ. 0) THEN
    IF (lAlphaBlanks(10)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(10))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(10))//'="'//TRIM(Alphas(10))//'".')
    END IF
    ErrorsFound = .TRUE.
  ELSE
    ! Verify Curve Object, only legal type is BiQuadratic
    SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%CCapFTemp2))

    CASE('BIQUADRATIC')
      DXCoil(DXCoilNum)%TotCapTempModFacCurveType(2)=Biquadratic

    CASE DEFAULT
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(10))//' type for this object = '// &
                           TRIM(GetCurveType(DXCoil(DXCoilNum)%CCapFTemp2)))
      CALL ShowContinueError('Curve type must be BiQuadratic.')
      ErrorsFound=.TRUE.
    END SELECT
  END IF

  DXCoil(DXCoilNum)%EIRFTemp2 = GetCurveIndex(Alphas(11)) ! convert curve name to number
  IF (DXCoil(DXCoilNum)%EIRFTemp2 .EQ. 0) THEN
    IF (lAlphaBlanks(11)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(11))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(11))//'="'//TRIM(Alphas(11))//'".')
    END IF
    ErrorsFound = .TRUE.
  ELSE
    ! Verify Curve Object, only legal type is Biquadratic
    SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%EIRFTemp2))

    CASE('BIQUADRATIC')
      DXCoil(DXCoilNum)%EIRTempModFacCurveType(2)=Biquadratic

    CASE DEFAULT
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(11))//' type for this object = '// &
                           TRIM(GetCurveType(DXCoil(DXCoilNum)%EIRFTemp2)))
      CALL ShowContinueError('Curve type must be BiQuadratic.')
      ErrorsFound=.TRUE.
    END SELECT
  END IF

! outdoor condenser node
  IF (lAlphaBlanks(12)) THEN
    DXCoil(DXCoilNum)%CondenserInletNodeNum(1) = 0
  ELSE
    DXCoil(DXCoilNum)%CondenserInletNodeNum(1) = &
       GetOnlySingleNode(Alphas(12),ErrorsFound,TRIM(CurrentModuleObject),DXCoil(DXCoilNum)%Name, &
                         NodeType_Air,NodeConnectionType_OutsideAirReference,1,ObjectIsNotParent)
    IF (.not. CheckOutAirNodeNumber(DXCoil(DXCoilNum)%CondenserInletNodeNum(1))) THEN
      CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", may be invalid')
      CALL ShowContinueError(TRIM(cAlphaFields(12))//'="'//TRIM(Alphas(12))//&
                               '", node does not appear in an OutdoorAir:NodeList or as an OutdoorAir:Node.')
      CALL ShowContinueError('This node needs to be included in an air system or the coil model will not be valid' &
                             //', and the simulation continues')
    END IF
  ENDIF

  IF ((SameString(Alphas(13),'AirCooled')) .OR. lAlphaBlanks(13))THEN
    DXCoil(DXCoilNum)%CondenserType(1) = AirCooled
  ELSEIF (SameString(Alphas(13),'EvaporativelyCooled')) THEN
    DXCoil(DXCoilNum)%CondenserType(1) = EvapCooled
    DXCoil(DXCoilNum)%ReportEvapCondVars = .TRUE.
  ELSE
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//TRIM(cAlphaFields(13))//'="'//TRIM(Alphas(13))//'":')
    CALL ShowContinueError('...must be AirCooled or EvaporativelyCooled.')
    ErrorsFound = .TRUE.
  END IF

  DXCoil(DXCoilNum)%EvapCondEffect(1) = Numbers(9)
  IF (DXCoil(DXCoilNum)%EvapCondEffect(1) .LT. 0.0d0 .OR. DXCoil(DXCoilNum)%EvapCondEffect(1) .GT. 1.0d0) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//trim(cNumericFields(9))//' cannot be < 0.0 or > 1.0.')
    CALL ShowContinueError('...entered value=['//trim(TrimSigDigits(Numbers(9),2))//'].')
    ErrorsFound = .TRUE.
  END IF

  DXCoil(DXCoilNum)%EvapCondAirFlow(1) = Numbers(10)
  IF (DXCoil(DXCoilNum)%EvapCondAirFlow(1) .LT. 0.0d0 .AND. DXCoil(DXCoilNum)%EvapCondAirFlow(1) /= AutoSize) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//trim(cNumericFields(10))//' cannot be < 0.0.')
    CALL ShowContinueError('...entered value=['//trim(TrimSigDigits(Numbers(10),2))//'].')
    ErrorsFound = .TRUE.
  END IF

  DXCoil(DXCoilNum)%EvapCondPumpElecNomPower(1) = Numbers(11)
  IF (DXCoil(DXCoilNum)%EvapCondPumpElecNomPower(1) .LT. 0.0d0 .AND. &
      DXCoil(DXCoilNum)%EvapCondPumpElecNomPower(1) /= AutoSize) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//trim(cNumericFields(11))//' cannot be < 0.0.')
    CALL ShowContinueError('...entered value=['//trim(TrimSigDigits(Numbers(11),2))//'].')
    ErrorsFound = .TRUE.
  END IF

  DXCoil(DXCoilNum)%EvapCondEffect2 = Numbers(12)
  IF (DXCoil(DXCoilNum)%EvapCondEffect2 .LT. 0.0d0 .OR. DXCoil(DXCoilNum)%EvapCondEffect2 .GT. 1.0d0) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//trim(cNumericFields(12))//' cannot be cannot be < 0.0 or > 1.0.')
    CALL ShowContinueError('...entered value=['//trim(TrimSigDigits(Numbers(12),2))//'].')
    ErrorsFound = .TRUE.
  END IF

  DXCoil(DXCoilNum)%EvapCondAirFlow2 = Numbers(13)
  IF (DXCoil(DXCoilNum)%EvapCondAirFlow2 .LT. 0.0d0 .AND. DXCoil(DXCoilNum)%EvapCondAirFlow2 /= AutoSize) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//trim(cNumericFields(13))//' cannot be < 0.0.')
    CALL ShowContinueError('...entered value=['//trim(TrimSigDigits(Numbers(13),2))//'].')
    ErrorsFound = .TRUE.
  END IF

  DXCoil(DXCoilNum)%EvapCondPumpElecNomPower2 = Numbers(14)
  IF (DXCoil(DXCoilNum)%EvapCondPumpElecNomPower2 .LT. 0.0d0 .AND. DXCoil(DXCoilNum)%EvapCondPumpElecNomPower2 /= AutoSize) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//trim(cNumericFields(14))//' cannot be < 0.0.')
    CALL ShowContinueError('...entered value=['//trim(TrimSigDigits(Numbers(14),2))//'].')
    ErrorsFound = .TRUE.
  END IF

  DXCoil(DXCoilNum)%RatedEIR2 = 1.d0 / DXCoil(DXCoilNum)%RatedCOP2

  ! Get Water System tank connections
  !  A14, \field Name of Water Storage Tank for Supply
  DXCoil(DXCoilNum)%EvapWaterSupplyName = Alphas(14)
  IF (lAlphaBlanks(14)) THEN
    DXCoil(DXCoilNum)%EvapWaterSupplyMode = WaterSupplyFromMains
  ELSE
    DXCoil(DXCoilNum)%EvapWaterSupplyMode = WaterSupplyFromTank
    CALL SetupTankDemandComponent(DXCoil(DXCoilNum)%Name, TRIM(CurrentModuleObject), &
                 DXCoil(DXCoilNum)%EvapWaterSupplyName, ErrorsFound, DXCoil(DXCoilNum)%EvapWaterSupTankID, &
                 DXCoil(DXCoilNum)%EvapWaterTankDemandARRID )
  ENDIF

  !A15; \field Name of Water Storage Tank for Condensate Collection
  DXCoil(DXCoilNum)%CondensateCollectName = Alphas(15)
  IF (lAlphaBlanks(15)) THEN
    DXCoil(DXCoilNum)%CondensateCollectMode = CondensateDiscarded
  ELSE
    DXCoil(DxCoilNum)%CondensateCollectMode = CondensateToTank
    CALL SetupTankSupplyComponent(DXCoil(DXCoilNum)%Name, TRIM(CurrentModuleObject), &
                 DXCoil(DXCoilNum)%CondensateCollectName, ErrorsFound, DXCoil(DXCoilNum)%CondensateTankID, &
                 DXCoil(DXCoilNum)%CondensateTankSupplyARRID )
  ENDIF

  ! Basin heater power as a function of temperature must be greater than or equal to 0
  DXCoil(DxCoilNum)%BasinHeaterPowerFTempDiff = Numbers(15)
  IF(Numbers(15) .LT. 0.0d0) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//trim(cNumericFields(15))//' must be >= 0.0.')
    CALL ShowContinueError('...entered value=['//trim(TrimSigDigits(Numbers(15),2))//'].')
    ErrorsFound = .TRUE.
  END IF

  DXCoil(DxCoilNum)%BasinHeaterSetPointTemp = Numbers(16)
  IF(DXCoil(DxCoilNum)%BasinHeaterPowerFTempDiff .GT. 0.0d0) THEN
    IF(NumNumbers .LT. 16) THEN
      DXCoil(DxCoilNum)%BasinHeaterSetPointTemp = 2.0d0
    ENDIF
    IF(DXCoil(DxCoilNum)%BasinHeaterSetPointTemp < 2.0d0) THEN
      CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", freeze possible')
      CALL ShowContinueError('...'//trim(cNumericFields(16))//' is < 2 {C}. Freezing could occur.')
      CALL ShowContinueError('...entered value=['//trim(TrimSigDigits(Numbers(16),2))//'].')
    END IF
  END IF

  IF(.NOT. lAlphaBlanks(16))THEN
    DXCoil(DxCoilNum)%BasinHeaterSchedulePtr   = GetScheduleIndex(Alphas(16))
    IF(DXCoil(DxCoilNum)%BasinHeaterSchedulePtr .EQ. 0)THEN
      CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(16))//'="'//TRIM(Alphas(16))//'".')
      CALL ShowContinueError('Basin heater will be available to operate throughout the simulation.')
    END IF
  END IF

END DO

IF (ErrorsFound) THEN
  CALL ShowFatalError(RoutineName//'Errors found in getting '//TRIM(CurrentModuleObject)//' input.  '//&
                      'Preceding condition(s) causes termination.')
END IF

! Loop over the DX Water Heater Coils and get & load the data
CurrentModuleObject='Coil:WaterHeating:AirToWaterHeatPump'
DO DXHPWaterHeaterCoilNum = 1, NumDXHeatPumpWaterHeaterCoils

  CALL GetObjectItem(TRIM(CurrentModuleObject),DXHPWaterHeaterCoilNum,Alphas,NumAlphas,Numbers,NumNumbers,IOStatus, &
                     NumBlank=lNumericBlanks,AlphaBlank=lAlphaBlanks, &
                     AlphaFieldNames=cAlphaFields,NumericFieldNames=cNumericFields)

  DXCoilNum = DXCoilNum+1
  IsNotOK=.FALSE.
  IsBlank=.FALSE.
  CALL VerifyName(Alphas(1),DXCoil%Name,DXCoilNum-1,IsNotOK,IsBlank,TRIM(CurrentModuleObject)//' Name')
  IF (IsNotOK) THEN
    ErrorsFound=.true.
    IF (IsBlank) Alphas(1)='xxxxx'
  ENDIF
  DXCoil(DXCoilNum)%Name = Alphas(1)
  DXCoil(DXCoilNum)%DXCoilType = TRIM(CurrentModuleObject)
  DXCoil(DXCoilNum)%DXCoilType_Num = CoilDX_HeatPumpWaterHeater
  DXCoil(DXCoilNum)%SchedPtr = 0  ! heat pump water heater DX coil has no schedule

! Store the HPWH DX coil heating capacity in RatedTotCap2. After backing off pump and fan heat,
! move to RatedTotCap() for use by DX coil
  DXCoil(DXCoilNum)%RatedTotCap2                = Numbers(1)
  IF(DXCoil(DXCoilNum)%RatedTotCap2 .LE. 0.0d0) THEN
     CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
     CALL ShowContinueError('...'//TRIM(cNumericFields(1))//' must be > 0.0,'//  &
        ' entered value=['//trim(TrimSigDigits(Numbers(1),2))//'].')
    ErrorsFound=.TRUE.
  END IF

  DXCoil(DXCoilNum)%RatedCOP(1)                 = Numbers(2)
  IF(DXCoil(DXCoilNum)%RatedCOP(1) .LE. 0.0d0) THEN
     CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
     CALL ShowContinueError('...'//TRIM(cNumericFields(2))//' must be > 0.0,'//  &
        ' entered value=['//trim(TrimSigDigits(Numbers(2),2))//'].')
    ErrorsFound=.TRUE.
  END IF

  DXCoil(DXCoilNum)%RatedSHR(1) = Numbers(3)
  IF(DXCoil(DXCoilNum)%RatedSHR(1) .LE. 0.0d0 .OR. DXCoil(DXCoilNum)%RatedSHR(1) .GT. 1.0d0) THEN
     CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
     CALL ShowContinueError('...'//TRIM(cNumericFields(3))//' must be > 0 and <= 1. '//  &
                ' entered value=['//trim(TrimSigDigits(Numbers(3),3))//'].')

    ErrorsFound=.TRUE.
  END IF

  DXCoil(DXCoilNum)%RatedInletDBTemp         = Numbers(4)
  IF(DXCoil(DXCoilNum)%RatedInletDBTemp .LE. 5.0d0) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//TRIM(cNumericFields(4))//' must be > 5 {C}. '//  &
                ' entered value=['//trim(TrimSigDigits(Numbers(4),1))//'].')
    ErrorsFound=.TRUE.
  END IF

  DXCoil(DXCoilNum)%RatedInletWBTemp         = Numbers(5)
  IF(DXCoil(DXCoilNum)%RatedInletWBTemp .LE. 5.0d0) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//TRIM(cNumericFields(5))//' must be > 5 {C}. '//  &
                ' entered value=['//trim(TrimSigDigits(Numbers(5),1))//'].')
    ErrorsFound=.TRUE.
  END IF

  DXCoil(DXCoilNum)%RatedInletWaterTemp      = Numbers(6)
  IF(DXCoil(DXCoilNum)%RatedInletWaterTemp .LE. 25.0d0) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//TRIM(cNumericFields(6))//' must be > 25 {C}. '//  &
                ' entered value=['//trim(TrimSigDigits(Numbers(6),1))//'].')
    ErrorsFound=.TRUE.
  END IF

  DXCoil(DXCoilNum)%RatedAirVolFlowRate(1)    = Numbers(7)
  IF(DXCoil(DXCoilNum)%RatedAirVolFlowRate(1) /= AutoCalculate) THEN
    IF(DXCoil(DXCoilNum)%RatedAirVolFlowRate(1) .LE. 0.0d0) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//TRIM(cNumericFields(7))//' must be > 0.0. '//  &
                ' entered value=['//trim(TrimSigDigits(Numbers(7),3))//'].')
      ErrorsFound=.TRUE.
    END IF
  END IF

  DXCoil(DXCoilNum)%RatedHPWHCondWaterFlow    = Numbers(8)
! move to init
  IF(DXCoil(DXCoilNum)%RatedHPWHCondWaterFlow /= AutoCalculate) THEN
    IF(DXCoil(DXCoilNum)%RatedHPWHCondWaterFlow .LE. 0.0) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...'//TRIM(cNumericFields(8))//' must be > 0.0 '//  &
                ' entered value=['//trim(TrimSigDigits(Numbers(8),3))//'].')
      ErrorsFound=.TRUE.
    END IF
!   check the range of flow rate to be >= 1 gpm/ton and <= 5 gpm/ton
    IF(DXCoil(DXCoilNum)%RatedHPWHCondWaterFlow/DXCoil(DXCoilNum)%RatedTotCap2 .LT. 1.79405d-8 .OR. &
       DXCoil(DXCoilNum)%RatedHPWHCondWaterFlow/DXCoil(DXCoilNum)%RatedTotCap2 .GT. 8.97024d-8)THEN
      CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", outside range')
      CALL ShowContinueError('...'//TRIM(cNumericFields(8))//' per watt of '//TRIM(cNumericFields(1))//  &
         ' is outside the recommended range of >= 1.79405E-8 m3/s/W (0.083 gpm/MBH) and <= 8.97024E-8 m3/s/W (0.417 gpm/MBH).')
      CALL ShowContinueError('...Entered Flow rate per watt = ['//  &
             TRIM(TrimSigDigits((DXCoil(DXCoilNum)%RatedHPWHCondWaterFlow/DXCoil(DXCoilNum)%RatedTotCap2),10))//'].')
    END IF
  END IF

  IF(SameString(Alphas(2),'Yes') .OR. SameString(Alphas(2),'No')) THEN
!  initialized to TRUE on allocate
   IF(SameString(Alphas(2),'No'))DXCoil(DXCoilNum)%FanPowerIncludedInCOP  = .FALSE.
  ELSE
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError(',,,invalid choice for '//TRIM(cAlphaFields(2))//'.  Entered choice = '//TRIM(Alphas(2)))
    CALL ShowContinueError('Valid choices are Yes or No.')
    ErrorsFound=.TRUE.
  END IF

  IF(SameString(Alphas(3),'Yes') .OR. SameString(Alphas(3),'No')) THEN
!  initialized to FALSE on allocate
    IF(SameString(Alphas(3),'Yes'))DXCoil(DXCoilNum)%CondPumpPowerInCOP     = .TRUE.
  ELSE
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError(',,,invalid choice for '//TRIM(cAlphaFields(3))//'.  Entered choice = '//TRIM(Alphas(3)))
    CALL ShowContinueError('Valid choices are Yes or No.')
    ErrorsFound=.TRUE.
  END IF

  IF(SameString(Alphas(4),'Yes') .OR. SameString(Alphas(4),'No')) THEN
!  initialized to FALSE on allocate
    IF(SameString(Alphas(4),'Yes'))DXCoil(DXCoilNum)%CondPumpHeatInCapacity = .TRUE.
  ELSE
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError(',,,invalid choice for '//TRIM(cAlphaFields(4))//'.  Entered choice = '//TRIM(Alphas(4)))
    CALL ShowContinueError('Valid choices are Yes or No.')
    ErrorsFound=.TRUE.
  END IF

  DXCoil(DXCoilNum)%HPWHCondPumpElecNomPower    = Numbers(9)
  IF(DXCoil(DXCoilNum)%HPWHCondPumpElecNomPower .LT. 0.0d0) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//TRIM(cNumericFields(9))//' must be >= 0.0 '//  &
                ' entered value=['//trim(TrimSigDigits(Numbers(9),3))//'].')
    ErrorsFound=.TRUE.
  END IF

  DXCoil(DXCoilNum)%HPWHCondPumpFracToWater     = Numbers(10)
  IF(DXCoil(DXCoilNum)%HPWHCondPumpFracToWater .LE. 0.0d0 .OR. DXCoil(DXCoilNum)%HPWHCondPumpFracToWater .GT. 1.0d0) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//TRIM(cNumericFields(10))//' must be >= 0 and <= 1. '//  &
                ' entered value=['//trim(TrimSigDigits(Numbers(10),3))//'].')
    ErrorsFound=.TRUE.
  END IF

  DXCoil(DXCoilNum)%AirInNode = &
               GetOnlySingleNode(Alphas(5),ErrorsFound,TRIM(CurrentModuleObject),Alphas(1), &
               NodeType_Air,NodeConnectionType_Inlet,1,ObjectIsNotParent)

  DXCoil(DXCoilNum)%AirOutNode = &
               GetOnlySingleNode(Alphas(6),ErrorsFound,TRIM(CurrentModuleObject),Alphas(1), &
               NodeType_Air,NodeConnectionType_Outlet,1,ObjectIsNotParent)

  CALL TestCompSet(TRIM(CurrentModuleObject),Alphas(1),Alphas(5),Alphas(6),'Air Nodes')

  DXCoil(DXCoilNum)%WaterInNode = &
               GetOnlySingleNode(Alphas(7),ErrorsFound,TRIM(CurrentModuleObject),Alphas(1), &
               NodeType_Water,NodeConnectionType_Inlet,2,ObjectIsNotParent)

  DXCoil(DXCoilNum)%WaterOutNode = &
               GetOnlySingleNode(Alphas(8),ErrorsFound,TRIM(CurrentModuleObject),Alphas(1), &
               NodeType_Water,NodeConnectionType_Outlet,2,ObjectIsNotParent)

  CALL TestCompSet(TRIM(CurrentModuleObject),Alphas(1),Alphas(7),Alphas(8),'Water Nodes')

  DXCoil(DXCoilNum)%CrankcaseHeaterCapacity     = Numbers(11)
  IF(DXCoil(DXCoilNum)%CrankcaseHeaterCapacity .LT. 0.0d0) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//TRIM(cNumericFields(11))//' must be >= 0.0 '//  &
                ' entered value=['//trim(TrimSigDigits(Numbers(11),1))//'].')
    ErrorsFound=.TRUE.
  END IF

  DXCoil(DXCoilNum)%MaxOATCrankcaseHeater       = Numbers(12)
  IF(DXCoil(DXCoilNum)%MaxOATCrankcaseHeater .LT. 0.0d0) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//TRIM(cNumericFields(12))//' must be >= 0 {C}. '//  &
                ' entered value=['//trim(TrimSigDigits(Numbers(12),1))//'].')
    ErrorsFound=.TRUE.
  END IF

  IF(SameString(Alphas(9),'DryBulbTemperature'))THEN
    DXCoil(DXCoilNum)%InletAirTemperatureType = DryBulbIndicator
  ELSEIF(SameString(Alphas(9),'WetBulbTemperature'))THEN
    DXCoil(DXCoilNum)%InletAirTemperatureType = WetBulbIndicator
  ELSE
!   wrong temperature type selection
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//TRIM(cAlphaFields(9))//' must be DryBulbTemperature or WetBulbTemperature.')
    CALL ShowContinueError('...entered value="'//trim(Alphas(9))//'".')
    ErrorsFound = .TRUE.
  END IF

! set rated inlet air temperature for curve object verification
  IF(DXCoil(DXCoilNum)%InletAirTemperatureType .EQ. WetBulbIndicator) THEN
    InletAirTemp = DXCoil(DXCoilNum)%RatedInletWBTemp
  ELSE
    InletAirTemp = DXCoil(DXCoilNum)%RatedInletDBTemp
  END IF
! set rated water temperature for curve object verification
  InletWaterTemp = DXCoil(DXCoilNum)%RatedInletWaterTemp

  IF (.NOT. lAlphaBlanks(10)) THEN
    DXCoil(DXCoilNum)%HCapFTemp = GetCurveIndex(Alphas(10))
    IF (DXCoil(DXCoilNum)%HCapFTemp .EQ. 0) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(10))//'="'//TRIM(Alphas(10))//'".')
      ErrorsFound = .TRUE.
    ELSE
      ! Verify Curve Object, only legal types are BiQuadratic or Cubic
      SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%HCapFTemp))

      CASE('BIQUADRATIC')
        DXCoil(DXCoilNum)%HCapFTempCurveType = Biquadratic
        HeatCapFTemp = CurveValue(DXCoil(DXCoilNum)%HCapFTemp,InletAirTemp,InletWaterTemp)

      CASE('CUBIC')
        DXCoil(DXCoilNum)%HCapFTempCurveType = Cubic
        HeatCapFTemp = CurveValue(DXCoil(DXCoilNum)%HCapFTemp,InletAirTemp)

      CASE DEFAULT
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(10))//' type for this object = '// &
                             TRIM(GetCurveType(DXCoil(DXCoilNum)%HCapFTemp)))
        CALL ShowContinueError('Curve type must be BiQuadratic or Cubic.')
        ErrorsFound=.TRUE.
        HeatCapFTemp = 1.0d0
      END SELECT

      IF(ABS(HeatCapFTemp - 1.0d0) .GT. 0.05d0)THEN
        CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'":')
        ! could remove name from the field for output
        CALL ShowContinueError('...The '//TRIM(cAlphaFields(10))//' should be normalized to 1.0 at the rating point.')
        CALL ShowContinueError('...Curve output at the rating point = '//TrimSigDigits(HeatCapFTemp,3))
        CALL ShowContinueError('...The simulation continues using the user-specified curve.')
      END IF

    END IF
  END IF

  IF (.NOT. lAlphaBlanks(11)) THEN
    DXCoil(DXCoilNum)%HCapFAirFlow = GetCurveIndex(Alphas(11))
    IF (DXCoil(DXCoilNum)%HCapFAirFlow .EQ. 0) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(11))//'="'//TRIM(Alphas(11))//'".')
      ErrorsFound = .TRUE.
    ELSE
      ! Verify Curve Object, only legal types are Cubic or Quadratic
      SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%HCapFAirFlow))

      CASE('CUBIC')
        HeatCapFAirFlow  = CurveValue(DXCoil(DXCoilNum)%HCapFAirFlow,1.0d0)

      CASE('QUADRATIC')
        HeatCapFAirFlow  = CurveValue(DXCoil(DXCoilNum)%HCapFAirFlow,1.0d0)

      CASE DEFAULT
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(11))//' type for this object = '// &
                             TRIM(GetCurveType(DXCoil(DXCoilNum)%HCapFAirFlow)))
        CALL ShowContinueError('Curve type must be Quadratic or Cubic.')
        ErrorsFound=.TRUE.
        HeatCapFAirFlow  = 1.0d0
      END SELECT

      IF(ABS(HeatCapFAirFlow - 1.0d0) .GT. 0.05d0)THEN
        CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'":')
        ! could remove name from the field for output
        CALL ShowContinueError('...The '//TRIM(cAlphaFields(11))//' should be normalized to 1.0 at the rating point.')
        CALL ShowContinueError('...Curve output at an air flow fraction of 1 = ' &
                               //TrimSigDigits(HeatCapFAirFlow,3))
        CALL ShowContinueError('...The simulation continues using the user-specified curve.')
      END IF

    END IF
  END IF

  IF (.NOT. lAlphaBlanks(12)) THEN
    DXCoil(DXCoilNum)%HCapFWaterFlow = GetCurveIndex(Alphas(12))
    IF (DXCoil(DXCoilNum)%HCapFWaterFlow .EQ. 0) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(12))//'="'//TRIM(Alphas(12))//'".')
      ErrorsFound = .TRUE.
    ELSE
      ! Verify Curve Object, only legal types are Cubic or Quadratic
      SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%HCapFWaterFlow))

      CASE('CUBIC')
        HeatCapFWaterFlow  = CurveValue(DXCoil(DXCoilNum)%HCapFWaterFlow,1.0d0)

      CASE('QUADRATIC')
        HeatCapFWaterFlow  = CurveValue(DXCoil(DXCoilNum)%HCapFWaterFlow,1.0d0)

      CASE DEFAULT
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(12))//' type for this object = '// &
                             TRIM(GetCurveType(DXCoil(DXCoilNum)%HCapFWaterFlow)))
        CALL ShowContinueError('Curve type must be Quadratic or Cubic.')
        ErrorsFound=.TRUE.
        HeatCapFWaterFlow  = 1.0d0
      END SELECT

      IF(ABS(HeatCapFWaterFlow - 1.0d0) .GT. 0.05d0)THEN
        CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'":')
        ! could remove name from the field for output
        CALL ShowContinueError('...The '//TRIM(cAlphaFields(11))//' should be normalized to 1.0 at the rating point.')
        CALL ShowContinueError('...Curve output at an air flow fraction of 1 = ' &
                             //TrimSigDigits(HeatCapFWaterFlow,3))
        CALL ShowContinueError('...The simulation continues using the user-specified curve.')
      END IF

    END IF
  END IF

  IF (.NOT. lAlphaBlanks(13)) THEN
    DXCoil(DXCoilNum)%HCOPFTemp = GetCurveIndex(Alphas(13))
    IF (DXCoil(DXCoilNum)%HCOPFTemp .EQ. 0) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(13))//'="'//TRIM(Alphas(13))//'".')
      ErrorsFound = .TRUE.
    ELSE
      ! Verify Curve Object, only legal types are BiQuadratic or Cubic
      SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%HCOPFTemp))

      CASE('BIQUADRATIC')
        DXCoil(DXCoilNum)%HCOPFTempCurveType = Biquadratic
        HeatCOPFTemp = CurveValue(DXCoil(DXCoilNum)%HCOPFTemp,InletAirTemp,InletWaterTemp)

      CASE('CUBIC')
        DXCoil(DXCoilNum)%HCOPFTempCurveType = Cubic
        HeatCOPFTemp = CurveValue(DXCoil(DXCoilNum)%HCOPFTemp,InletAirTemp)

      CASE DEFAULT
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(13))//' type for this object = '// &
                             TRIM(GetCurveType(DXCoil(DXCoilNum)%HCOPFTemp)))
        CALL ShowContinueError('Curve type must be BiQuadratic or Cubic.')
        ErrorsFound=.TRUE.
        HeatCOPFTemp = 1.0d0
      END SELECT

      IF(ABS(HeatCOPFTemp - 1.0d0) .GT. 0.05d0)THEN
        CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'":')
        ! could remove name from the field for output
        CALL ShowContinueError('...The '//TRIM(cAlphaFields(13))//' should be normalized to 1.0 at the rating point.')
        CALL ShowContinueError('...Curve output at an air flow fraction of 1 = ' &
                             //TrimSigDigits(HeatCOPFTemp,3))
        CALL ShowContinueError('...The simulation continues using the user-specified curve.')
      END IF

    END IF
  END IF

  IF (.NOT. lAlphaBlanks(14)) THEN
    DXCoil(DXCoilNum)%HCOPFAirFlow = GetCurveIndex(Alphas(14))
    IF (DXCoil(DXCoilNum)%HCOPFAirFlow .EQ. 0) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(14))//'="'//TRIM(Alphas(14))//'".')
      ErrorsFound = .TRUE.
    ELSE
      ! Verify Curve Object, only legal types are Cubic or Quadratic
      SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%HCOPFAirFlow))

      CASE('CUBIC')
        HeatCOPFAirFlow  = CurveValue(DXCoil(DXCoilNum)%HCOPFAirFlow,1.0d0)

      CASE('QUADRATIC')
        HeatCOPFAirFlow  = CurveValue(DXCoil(DXCoilNum)%HCOPFAirFlow,1.0d0)

      CASE DEFAULT
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(14))//' type for this object = '// &
                             TRIM(GetCurveType(DXCoil(DXCoilNum)%HCOPFAirFlow)))
        CALL ShowContinueError('Curve type must be Quadratic or Cubic.')
        ErrorsFound=.TRUE.
        HeatCOPFAirFlow  = 1.0d0
      END SELECT

      IF(ABS(HeatCOPFAirFlow - 1.0d0) .GT. 0.05d0)THEN
        CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'":')
        ! could remove name from the field for output
        CALL ShowContinueError('...The '//TRIM(cAlphaFields(14))//' should be normalized to 1.0 at the rating point.')
        CALL ShowContinueError('...Curve output at an air flow fraction of 1 = ' &
                                    //TrimSigDigits(HeatCOPFAirFlow,3))
        CALL ShowContinueError('...The simulation continues using the user-specified curve.')
      END IF

    END IF
  END IF

  IF (.NOT. lAlphaBlanks(15)) THEN
    DXCoil(DXCoilNum)%HCOPFWaterFlow = GetCurveIndex(Alphas(15))
    IF (DXCoil(DXCoilNum)%HCOPFWaterFlow .EQ. 0) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(15))//'="'//TRIM(Alphas(15))//'".')
      ErrorsFound = .TRUE.
    ELSE
      ! Verify Curve Object, only legal types are Cubic or Quadratic
      SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%HCOPFWaterFlow))

      CASE('CUBIC')
        HeatCOPFWaterFlow  = CurveValue(DXCoil(DXCoilNum)%HCOPFWaterFlow,1.0d0)

      CASE('QUADRATIC')
        HeatCOPFWaterFlow  = CurveValue(DXCoil(DXCoilNum)%HCOPFWaterFlow,1.0d0)

      CASE DEFAULT
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(15))//' type for this object = '// &
                             TRIM(GetCurveType(DXCoil(DXCoilNum)%HCOPFWaterFlow)))
        CALL ShowContinueError('Curve type must be Quadratic or Cubic.')
        ErrorsFound=.TRUE.
        HeatCOPFWaterFlow  = 1.0d0
      END SELECT

      IF(ABS(HeatCOPFWaterFlow - 1.0d0) .GT. 0.05d0)THEN
        CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'":')
        ! could remove name from the field for output
        CALL ShowContinueError('...The '//TRIM(cAlphaFields(15))//' should be normalized to 1.0 at the rating point.')
        CALL ShowContinueError('...Curve output at a water flow fraction of 1 = ' &
                               //TrimSigDigits(HeatCOPFWaterFlow,3))
        CALL ShowContinueError('...The simulation continues using the user-specified curve.')
      END IF

    END IF
  END IF

  IF (.NOT. lAlphaBlanks(16)) THEN
    DXCoil(DXCoilNum)%PLFFPLR(1) = GetCurveIndex(Alphas(16))
    IF (DXCoil(DXCoilNum)%PLFFPLR(1) .EQ. 0) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(16))//'="'//TRIM(Alphas(16))//'".')
      ErrorsFound = .TRUE.
    ELSE
      ! Verify Curve Object, only legal types are Cubic or Quadratic
      SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%PLFFPLR(1)))

      CASE('CUBIC')

      CASE('QUADRATIC')

      CASE DEFAULT
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(16))//' type for this object = '// &
                             TRIM(GetCurveType(DXCoil(DXCoilNum)%PLFFPLR(1))))
        CALL ShowContinueError('Curve type must be Quadratic or Cubic.')
        ErrorsFound=.true.
      END SELECT

      IF(.NOT. ErrorsFound)THEN
!       Test PLF curve minimum and maximum. Cap if less than 0.7 or greater than 1.0.
        MinCurveVal = 999.0d0
        MaxCurveVal = -999.0d0
        CurveInput = 0.0d0
        DO WHILE (CurveInput <= 1.0d0)
          CurveVal = CurveValue(DXCoil(DXCoilNum)%PLFFPLR(1),CurveInput)
          IF(CurveVal .LT. MinCurveVal)THEN
            MinCurveVal = CurveVal
            MinCurvePLR = CurveInput
          END IF
          IF(CurveVal .GT. MaxCurveVal)THEN
            MaxCurveVal = CurveVal
            MaxCurvePLR = CurveInput
          END IF
          CurveInput=CurveInput+0.01d0
        END DO
        IF(MinCurveVal .LT. 0.7d0)THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
          CALL ShowContinueError('...'//TRIM(cAlphaFields(16))//' = '//TRIM(Alphas(16))//' has out of range value.')
          CALL ShowContinueError('...Curve minimum must be >= 0.7, '// &
                        'curve min at PLR = '//TRIM(TrimSigDigits(MinCurvePLR,2))//' is '//TRIM(TrimSigDigits(MinCurveVal,3)))
          CALL ShowContinueError('...Setting curve minimum to 0.7 and simulation continues.')
          CALL SetCurveOutputMinMaxValues(DXCoil(DXCoilNum)%PLFFPLR(1),ErrorsFound,CurveMin=0.7d0)
        END IF

        IF(MaxCurveVal .GT. 1.0d0)THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
          CALL ShowContinueError('...'//TRIM(cAlphaFields(16))//' = '//TRIM(Alphas(16))//' has out of range value.')
          CALL ShowContinueError('...Curve maximum must be <= 1.0, '// &
                        'curve max at PLR = '//TRIM(TrimSigDigits(MaxCurvePLR,2))//' is '//TRIM(TrimSigDigits(MaxCurveVal,3)))
          CALL ShowContinueError('...Setting curve maximum to 1.0 and simulation continues.')
          CALL SetCurveOutputMinMaxValues(DXCoil(DXCoilNum)%PLFFPLR(1),ErrorsFound,CurveMax=1.0d0)
        END IF

      END IF

    END IF
  END IF

! assume compressor resides at the inlet to the DX Coil
  DXCoil(DXCoilNum)%CondenserInletNodeNum(1) = DXCoil(DXCoilNum)%AirInNode

! set condenser type as HPWH
  DXCoil(DXCoilNum)%CondenserType(1) = WaterHeater

END DO ! end of the DX water heater coil loop

IF (ErrorsFound) THEN
  CALL ShowFatalError(RoutineName//'Errors found in getting '//TRIM(CurrentModuleObject)//' input. '// &
                      'Preceding condition(s) causes termination.')
END IF

! DX Multispeed cooling coil
CurrentModuleObject='Coil:Cooling:DX:MultiSpeed'
DO DXCoilIndex = 1,NumDXMulSpeedCoolCoils

  DXCoilNum = DXCoilNum+1

  CALL GetObjectItem(TRIM(CurrentModuleObject),DXCoilIndex,Alphas,NumAlphas,Numbers,NumNumbers,IOStatus, &
                     NumBlank=lNumericBlanks,AlphaBlank=lAlphaBlanks, &
                     AlphaFieldNames=cAlphaFields,NumericFieldNames=cNumericFields)

  IsNotOK=.FALSE.
  IsBlank=.FALSE.
  CALL VerifyName(Alphas(1),DXCoil%Name,DXCoilNum-1,IsNotOK,IsBlank,TRIM(CurrentModuleObject)//' Name')
  IF (IsNotOK) THEN
    ErrorsFound=.true.
    IF (IsBlank) Alphas(1)='xxxxx'
  ENDIF
  DXCoil(DXCoilNum)%Name = Alphas(1)
! Initialize DataHeatBalance heat reclaim variable name for use by heat reclaim coils
  HeatReclaimDXCoil(DXCoilNum)%Name = DXCoil(DXCoilNum)%Name
  HeatReclaimDXCoil(DXCoilNum)%SourceType = TRIM(CurrentModuleObject)
  DXCoil(DXCoilNum)%DXCoilType = TRIM(CurrentModuleObject)
  DXCoil(DXCoilNum)%DXCoilType_Num = CoilDX_MultiSpeedCooling
  DXCoil(DXCoilNum)%Schedule = Alphas(2)
  DXCoil(DXCoilNum)%SchedPtr = GetScheduleIndex(Alphas(2))  ! convert schedule name to pointer
  IF (DXCoil(DXCoilNum)%SchedPtr .EQ. 0) THEN
    IF (lAlphaBlanks(2)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(2))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...'//TRIM(cAlphaFields(2))//'="'//TRIM(Alphas(2))//'".')
    END IF
    ErrorsFound=.TRUE.
  END IF

  DXCoil(DXCoilNum)%AirInNode = &
               GetOnlySingleNode(Alphas(3),ErrorsFound,TRIM(CurrentModuleObject),Alphas(1), &
               NodeType_Air,NodeConnectionType_Inlet,1,ObjectIsNotParent)

  DXCoil(DXCoilNum)%AirOutNode = &
               GetOnlySingleNode(Alphas(4),ErrorsFound,TRIM(CurrentModuleObject),Alphas(1), &
               NodeType_Air,NodeConnectionType_Outlet,1,ObjectIsNotParent)

  CALL TestCompSet(TRIM(CurrentModuleObject),Alphas(1),Alphas(3),Alphas(4),'Air Nodes')

! outdoor condenser node
  IF (lAlphaBlanks(5)) THEN
    DXCoil(DXCoilNum)%CondenserInletNodeNum(1) = 0
  ELSE
    DXCoil(DXCoilNum)%CondenserInletNodeNum(1) = &
       GetOnlySingleNode(Alphas(5),ErrorsFound,TRIM(CurrentModuleObject),DXCoil(DXCoilNum)%Name, &
                         NodeType_Air,NodeConnectionType_OutsideAirReference,1,ObjectIsNotParent)
    IF (.not. CheckOutAirNodeNumber(DXCoil(DXCoilNum)%CondenserInletNodeNum(1))) THEN
      CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", may be invalid')
      CALL ShowContinueError(TRIM(cAlphaFields(5))//'="'//TRIM(Alphas(5))//&
                               '", node does not appear in an OutdoorAir:NodeList or as an OutdoorAir:Node.')
      CALL ShowContinueError('This node needs to be included in an air system or the coil model will not be valid' &
                             //', and the simulation continues')
    END IF
  ENDIF

  IF ((SameString(Alphas(6),'AirCooled')) .OR. lAlphaBlanks(6))THEN
    DXCoil(DXCoilNum)%CondenserType(1) = AirCooled
  ELSEIF (SameString(Alphas(6),'EvaporativelyCooled')) THEN
    DXCoil(DXCoilNum)%CondenserType(1) = EvapCooled
    DXCoil(DXCoilNum)%ReportEvapCondVars = .TRUE.
  ELSE
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//TRIM(cAlphaFields(6))//'="'//TRIM(Alphas(6))//'":')
    CALL ShowContinueError('...must be AirCooled or EvaporativelyCooled.')
    ErrorsFound = .TRUE.
  END IF

  ! Get Water System tank connections
  !  A8, \field Name of Water Storage Tank for Supply
  DXCoil(DXCoilNum)%EvapWaterSupplyName = Alphas(7)
  IF (lAlphaBlanks(7)) THEN
    DXCoil(DXCoilNum)%EvapWaterSupplyMode = WaterSupplyFromMains
  ELSE
    DXCoil(DXCoilNum)%EvapWaterSupplyMode = WaterSupplyFromTank
    CALL SetupTankDemandComponent(DXCoil(DXCoilNum)%Name,TRIM(CurrentModuleObject), &
                 DXCoil(DXCoilNum)%EvapWaterSupplyName, ErrorsFound, DXCoil(DXCoilNum)%EvapWaterSupTankID, &
                 DXCoil(DXCoilNum)%EvapWaterTankDemandARRID )
  END IF

  !A9; \field Name of Water Storage Tank for Condensate Collection
  DXCoil(DXCoilNum)%CondensateCollectName = Alphas(8)
  IF (lAlphaBlanks(8)) THEN
    DXCoil(DXCoilNum)%CondensateCollectMode = CondensateDiscarded
  ELSE
    DXCoil(DxCoilNum)%CondensateCollectMode = CondensateToTank
    CALL SetupTankSupplyComponent(DXCoil(DXCoilNum)%Name, TRIM(CurrentModuleObject), &
                 DXCoil(DXCoilNum)%CondensateCollectName, ErrorsFound, DXCoil(DXCoilNum)%CondensateTankID, &
                 DXCoil(DXCoilNum)%CondensateTankSupplyARRID )
  END IF

  !Set crankcase heater capacity
  DXCoil(DXCoilNum)%CrankcaseHeaterCapacity = Numbers(1)
  IF (DXCoil(DXCoilNum)%CrankcaseHeaterCapacity .LT. 0.0d0) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//trim(cNumericFields(1))//' cannot be < 0.0.')
    CALL ShowContinueError('...entered value=['//trim(TrimSigDigits(Numbers(1),2))//'].')
    ErrorsFound = .TRUE.
  END IF

  !Set crankcase heater cutout temperature
  DXCoil(DXCoilNum)%MaxOATCrankcaseHeater = Numbers(2)

  IF (SameString(Alphas(9),'Yes')) THEN
    DXCoil(DXCoilNum)%PLRImpact = .TRUE.
  ELSE IF (SameString(Alphas(9),'No')) THEN
    DXCoil(DXCoilNum)%PLRImpact = .FALSE.
  ELSE
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError(',,,invalid choice for '//TRIM(cAlphaFields(9))//'.  Entered choice = '//TRIM(Alphas(9)))
    CALL ShowContinueError('The allowed choices are Yes or No.')
    ErrorsFound=.TRUE.
  END IF

  IF (SameString(Alphas(10),'Yes')) THEN
    DXCoil(DXCoilNum)%LatentImpact = .TRUE.
  ELSE IF (SameString(Alphas(10),'No')) THEN
    DXCoil(DXCoilNum)%LatentImpact = .FALSE.
  ELSE
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError(',,,invalid choice for '//TRIM(cAlphaFields(10))//'.  Entered choice = '//TRIM(Alphas(10)))
    CALL ShowContinueError('The allowed choices are Yes or No.')
    ErrorsFound=.TRUE.
  END IF

  !   Basin heater power as a function of temperature must be greater than or equal to 0
  DXCoil(DxCoilNum)%BasinHeaterPowerFTempDiff = Numbers(3)
  IF(Numbers(3) .LT. 0.0d0) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//TRIM(cNumericFields(3))//' must be >= 0.0, '//  &
                ' entered value=['//trim(TrimSigDigits(Numbers(3),3))//'].')
    ErrorsFound = .TRUE.
  END IF

  DXCoil(DxCoilNum)%BasinHeaterSetPointTemp = Numbers(4)
  IF(DXCoil(DxCoilNum)%BasinHeaterPowerFTempDiff .GT. 0.0d0) THEN
    IF(NumNumbers .LT. 4) THEN
      DXCoil(DxCoilNum)%BasinHeaterSetPointTemp = 2.0d0
    ENDIF
    IF(DXCoil(DxCoilNum)%BasinHeaterSetPointTemp < 2.0d0) THEN
      CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", freeze possible')
      CALL ShowContinueError('...'//trim(cNumericFields(4))//' is less than 2 {C}. Freezing could occur.')
      CALL ShowContinueError('...entered value=['//trim(TrimSigDigits(Numbers(4),2))//'].')
    END IF
  END IF

  IF(.NOT. lAlphaBlanks(11))THEN
    DXCoil(DxCoilNum)%BasinHeaterSchedulePtr   = GetScheduleIndex(Alphas(11))
    IF(DXCoil(DxCoilNum)%BasinHeaterSchedulePtr .EQ. 0)THEN
      CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(11))//'="'//TRIM(Alphas(11))//'".')
      CALL ShowContinueError('Basin heater will be available to operate throughout the simulation.')
    END IF
  END IF

  !A12; \field Fuel type
  IF (SameString(Alphas(12),'Electricity')) THEN
    DXCoil(DXCoilNum)%FuelType = FuelTypeElectricity
  ELSE IF (SameString(Alphas(12),'NaturalGas')) THEN
    DXCoil(DXCoilNum)%FuelType = FuelTypeNaturalGas
  ELSE IF (SameString(Alphas(12),'PropaneGas')) THEN
    DXCoil(DXCoilNum)%FuelType = FuelTypePropaneGas
  ELSE IF (SameString(Alphas(12),'Diesel')) THEN
    DXCoil(DXCoilNum)%FuelType = FuelTypeDiesel
  ELSE IF (SameString(Alphas(12),'Gasoline')) THEN
    DXCoil(DXCoilNum)%FuelType = FuelTypeGasoline
  ELSE IF (SameString(Alphas(12),'FuelOil#1')) THEN
    DXCoil(DXCoilNum)%FuelType = FuelTypeFuelOil1
  ELSE IF (SameString(Alphas(12),'FuelOil#2')) THEN
    DXCoil(DXCoilNum)%FuelType = FuelTypeFuelOil2
  ELSE
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError(',,,invalid choice for '//TRIM(cAlphaFields(12))//'.  Entered choice = '//TRIM(Alphas(12)))
    CALL ShowContinueError('Valid choices are Electricity, NaturalGas, PropaneGas, Diesel, Gasoline, FuelOil#1 or FuelOil#2')
    ErrorsFound=.TRUE.
  END IF

  DXCoil(DXCoilNum)%NumOfSpeeds = Numbers(5)       ! Number of speeds
  If (DXCoil(DXCoilNum)%NumOfSpeeds .LT. 2) Then
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//TRIM(cNumericFields(5))//' must be >= 2.'//  &
                                         ' entered number is '//TRIM(TrimSigDigits(Numbers(5),0)))
    ErrorsFound=.TRUE.
  End If

  ! Allocate arrays based on the number of speeds
  ALLOCATE(DXCoil(DXCoilNum)%MSErrIndex(DXCoil(DXCoilNum)%NumOfSpeeds))
  DXCoil(DXCoilNum)%MSErrIndex=0
  ALLOCATE(DXCoil(DXCoilNum)%MSRatedTotCap(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSRatedSHR(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSRatedCOP(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSRatedAirMassFlowRate(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSCCapFTemp(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSCCapFFlow(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSEIRFTemp(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSEIRFFlow(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSWasteHeat(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSEvapCondEffect(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSEvapCondAirFlow(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSEvapCondPumpElecNomPower(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSRatedCBF(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSWasteHeatFrac(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSPLFFPLR(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSTwet_Rated(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSGamma_Rated(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSMaxONOFFCyclesperHour(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSLatentCapacityTimeConstant(DXCoil(DXCoilNum)%NumOfSpeeds))

  Do I=1,DXCoil(DXCoilNum)%NumOfSpeeds
    DXCoil(DXCoilNum)%MSRatedTotCap(I) = Numbers(6+(I-1)*12)
    DXCoil(DXCoilNum)%MSRatedSHR(I)    = Numbers(7+(I-1)*12)
    DXCoil(DXCoilNum)%MSRatedCOP(I)    = Numbers(8+(I-1)*12)
    DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(I) = Numbers(9+(I-1)*12)

    DXCoil(DXCoilNum)%MSCCapFTemp(I) = GetCurveIndex(Alphas(13+(I-1)*6)) ! convert curve name to number
    IF (DXCoil(DXCoilNum)%MSCCapFTemp(I) .EQ. 0) THEN
      IF (lAlphaBlanks(13+(I-1)*6)) THEN
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
        CALL ShowContinueError('...required '//trim(cAlphaFields(13+(I-1)*6))//' is blank.')
      ELSE
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...not found '//TRIM(cAlphaFields(13+(I-1)*6))//'="'//TRIM(Alphas(13+(I-1)*6))//'".')
      END IF
      ErrorsFound = .TRUE.
    ELSE
      ! Verify Curve Object, only legal type is BiQuadratic
      SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%MSCCapFTemp(I)))

      CASE('BIQUADRATIC')
        CurveVal = CurveValue(DXCoil(DXCoilNum)%MSCCapFTemp(I),RatedInletWetbulbTemp,RatedOutdoorAirTemp)
        IF(CurveVal .GT. 1.10d0 .OR. CurveVal .LT. 0.90d0) THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", curve values')
          CALL ShowContinueError('...'//TRIM(cAlphaFields(13+(I-1)*6))//' output is not equal to 1.0 '//  &
                                             '(+ or - 10%) at rated conditions.')
          CALL ShowContinueError('...Curve output at rated conditions = '//TRIM(TrimSigDigits(CurveVal,3)))
        END IF

      CASE DEFAULT
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(13+(I-1)*6))//' type for this object = '// &
                             TRIM(GetCurveType(DXCoil(DXCoilNum)%MSCCapFTemp(I))))
        CALL ShowContinueError('Curve type must be BiQuadratic.')
        ErrorsFound=.TRUE.
      END SELECT
    END IF

    DXCoil(DXCoilNum)%MSCCapFFlow(I) = GetCurveIndex(Alphas(14+(I-1)*6)) ! convert curve name to number
    IF (DXCoil(DXCoilNum)%MSCCapFFlow(I) .EQ. 0) THEN
      IF (lAlphaBlanks(14+(I-1)*6)) THEN
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
        CALL ShowContinueError('...required '//trim(cAlphaFields(14+(I-1)*6))//' is blank.')
      ELSE
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...not found '//TRIM(cAlphaFields(14+(I-1)*6))//'="'//TRIM(Alphas(14+(I-1)*6))//'".')
      END IF
      ErrorsFound = .TRUE.
    ELSE
      ! Verify Curve Object, only legal type is Quadratic
      SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%MSCCapFFlow(I)))

      CASE('QUADRATIC')
        CurveVal = CurveValue(DXCoil(DXCoilNum)%MSCCapFFlow(I),1.0d0)
        IF(CurveVal .GT. 1.10d0 .OR. CurveVal .LT. 0.90d0) THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", curve values')
          CALL ShowContinueError('...'//TRIM(cAlphaFields(14+(I-1)*6))//' output is not equal to 1.0 '//  &
                                             '(+ or - 10%) at rated conditions.')
          CALL ShowContinueError('...Curve output at rated conditions = '//TRIM(TrimSigDigits(CurveVal,3)))
        END IF

      CASE('CUBIC')
        CurveVal = CurveValue(DXCoil(DXCoilNum)%MSCCapFFlow(I),1.0d0)
        IF(CurveVal .GT. 1.10d0 .OR. CurveVal .LT. 0.90d0) THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", curve values')
          CALL ShowContinueError('...'//TRIM(cAlphaFields(14+(I-1)*6))//' output is not equal to 1.0 '//  &
                                             '(+ or - 10%) at rated conditions.')
          CALL ShowContinueError('...Curve output at rated conditions = '//TRIM(TrimSigDigits(CurveVal,3)))
        END IF

      CASE DEFAULT
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(14+(I-1)*6))//' type for this object = '// &
                             TRIM(GetCurveType(DXCoil(DXCoilNum)%MSCCapFFlow(I))))
        CALL ShowContinueError('Curve type must be Quadratic or Cubic.')
        ErrorsFound=.TRUE.
      END SELECT
    END IF

    DXCoil(DXCoilNum)%MSEIRFTemp(I) = GetCurveIndex(Alphas(15+(I-1)*6)) ! convert curve name to number
    IF (DXCoil(DXCoilNum)%MSEIRFTemp(I) .EQ. 0) THEN
      IF (lAlphaBlanks(15+(I-1)*6)) THEN
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
        CALL ShowContinueError('...required '//trim(cAlphaFields(15+(I-1)*6))//' is blank.')
      ELSE
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...not found '//TRIM(cAlphaFields(15+(I-1)*6))//'="'//TRIM(Alphas(15+(I-1)*6))//'".')
      END IF
      ErrorsFound = .TRUE.
    ELSE
      ! Verify Curve Object, only legal type is Biquadratic
      SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%MSEIRFTemp(I)))

      CASE('BIQUADRATIC')
        CurveVal = CurveValue(DXCoil(DXCoilNum)%MSEIRFTemp(I),RatedInletWetbulbTemp,RatedOutdoorAirTemp)
        IF(CurveVal .GT. 1.10d0 .OR. CurveVal .LT. 0.90d0) THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", curve values')
          CALL ShowContinueError('...'//TRIM(cAlphaFields(15+(I-1)*6))//' output is not equal to 1.0 '//  &
                                             '(+ or - 10%) at rated conditions.')
          CALL ShowContinueError('...Curve output at rated conditions = '//TRIM(TrimSigDigits(CurveVal,3)))
        END IF

      CASE DEFAULT
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(15+(I-1)*6))//' type for this object = '// &
                             TRIM(GetCurveType(DXCoil(DXCoilNum)%MSEIRFTemp(1))))
        CALL ShowContinueError('Curve type must be BiQuadratic.')
        ErrorsFound=.TRUE.
      END SELECT
    END IF

    DXCoil(DXCoilNum)%MSEIRFFlow(I) = GetCurveIndex(Alphas(16+(I-1)*6)) ! convert curve name to number
    IF (DXCoil(DXCoilNum)%MSEIRFFlow(I) .EQ. 0) THEN
      IF (lAlphaBlanks(16+(I-1)*6)) THEN
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
        CALL ShowContinueError('...required '//trim(cAlphaFields(16+(I-1)*6))//' is blank.')
      ELSE
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...not found '//TRIM(cAlphaFields(16+(I-1)*6))//'="'//TRIM(Alphas(16+(I-1)*6))//'".')
      END IF
      ErrorsFound = .TRUE.
    ELSE
      ! Verify Curve Object, only legal type is Quadratic
      SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%MSEIRFFlow(I)))

      CASE('QUADRATIC')
        CurveVal = CurveValue(DXCoil(DXCoilNum)%MSEIRFFlow(I),1.0d0)
        IF(CurveVal .GT. 1.10d0 .OR. CurveVal .LT. 0.90d0) THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", curve values')
          CALL ShowContinueError('...'//TRIM(cAlphaFields(16+(I-1)*6))//' output is not equal to 1.0 '//  &
                                             '(+ or - 10%) at rated conditions.')
          CALL ShowContinueError('...Curve output at rated conditions = '//TRIM(TrimSigDigits(CurveVal,3)))
        END IF

      CASE('CUBIC')
        CurveVal = CurveValue(DXCoil(DXCoilNum)%MSEIRFFlow(I),1.0d0)
        IF(CurveVal .GT. 1.10d0 .OR. CurveVal .LT. 0.90d0) THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", curve values')
          CALL ShowContinueError('...'//TRIM(cAlphaFields(16+(I-1)*6))//' output is not equal to 1.0 '//  &
                                             '(+ or - 10%) at rated conditions.')
          CALL ShowContinueError('...Curve output at rated conditions = '//TRIM(TrimSigDigits(CurveVal,3)))
        END IF

      CASE DEFAULT
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(16+(I-1)*6))//' type for this object = '// &
              TRIM(GetCurveType(DXCoil(DXCoilNum)%MSEIRFFlow(I))))
        CALL ShowContinueError('Curve type must be Quadratic or Cubic.')
        ErrorsFound=.TRUE.
      END SELECT
    END IF

    DXCoil(DXCoilNum)%MSPLFFPLR(I) = GetCurveIndex(Alphas(17+(I-1)*6)) ! convert curve name to number
    IF (DXCoil(DXCoilNum)%MSPLFFPLR(I) .EQ. 0) THEN
      IF (lAlphaBlanks(17+(I-1)*6)) THEN
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
        CALL ShowContinueError('...required '//trim(cAlphaFields(17+(I-1)*6))//' is blank.')
      ELSE
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...not found '//TRIM(cAlphaFields(17+(I-1)*6))//'="'//TRIM(Alphas(16+(I-1)*6))//'".')
      END IF
      ErrorsFound = .TRUE.
    ELSE
    ! Verify Curve Object, only legal types are Quadratic or Cubic
      SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%MSPLFFPLR(I)))

      CASE('QUADRATIC')

      CASE('CUBIC')

      CASE DEFAULT
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(17+(I-1)*6))//' type for this object = '// &
                             TRIM(GetCurveType(DXCoil(DXCoilNum)%MSPLFFPLR(I))))
        CALL ShowContinueError('Curve type must be Quadratic or Cubic.')
        ErrorsFound=.TRUE.
      END SELECT

      IF(.NOT. ErrorsFound)THEN
!       Test PLF curve minimum and maximum. Cap if less than 0.7 or greater than 1.0.
        MinCurveVal = 999.0d0
        MaxCurveVal = -999.0d0
        CurveInput = 0.0d0
        DO WHILE (CurveInput <= 1.0d0)
          CurveVal = CurveValue(DXCoil(DXCoilNum)%MSPLFFPLR(I),CurveInput)
          IF(CurveVal .LT. MinCurveVal)THEN
            MinCurveVal = CurveVal
            MinCurvePLR = CurveInput
          END IF
          IF(CurveVal .GT. MaxCurveVal)THEN
            MaxCurveVal = CurveVal
            MaxCurvePLR = CurveInput
          END IF
          CurveInput=CurveInput+0.01d0
        END DO
        IF(MinCurveVal .LT. 0.7d0)THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
          CALL ShowContinueError('...'//TRIM(cAlphaFields2(17+(I-1)*6))//' = '//TRIM(Alphas2(17+(I-1)*6))//  &
             ' has out of range value.')
          CALL ShowContinueError('...Curve minimum must be >= 0.7, '// &
                   'curve min at PLR = '//TRIM(TrimSigDigits(MinCurvePLR,2))//' is '//TRIM(TrimSigDigits(MinCurveVal,3)))
          CALL ShowContinueError('...Setting curve minimum to 0.7 and simulation continues.')
          CALL SetCurveOutputMinMaxValues(DXCoil(DXCoilNum)%PLFFPLR(PerfModeNum),ErrorsFound,CurveMin=0.7d0)
        END IF

        IF(MaxCurveVal .GT. 1.0d0)THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
          CALL ShowContinueError('...'//TRIM(cAlphaFields2(17+(I-1)*6))//' = '//TRIM(Alphas2(17+(I-1)*6))//  &
             ' has out of range value.')
          CALL ShowContinueError('...Curve maximum must be <= 1.0, '// &
                        'curve max at PLR = '//TRIM(TrimSigDigits(MaxCurvePLR,2))//' is '//TRIM(TrimSigDigits(MaxCurveVal,3)))
          CALL ShowContinueError('...Setting curve maximum to 1.0 and simulation continues.')
          CALL SetCurveOutputMinMaxValues(DXCoil(DXCoilNum)%MSPLFFPLR(I),ErrorsFound,CurveMax=1.0d0)
        END IF

      END IF

    END IF

    ! read data for latent degradation
    DXCoil(DXCoilNum)%MSTwet_Rated(I) = Numbers(10+(I-1)*12)
    IF (DXCoil(DXCoilNum)%MSTwet_Rated(I) .LT. 0.0d0) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...'//TRIM(cNumericFields(10+(I-1)*12))//' cannot be < 0.0, '//  &
         'entered value=['//trim(TrimSigDigits(DXCoil(DXCoilNum)%MSTwet_Rated(I),4))//'].')
      ErrorsFound=.TRUE.
    END IF
    DXCoil(DXCoilNum)%MSGamma_Rated(I) = Numbers(11+(I-1)*12)
    IF (DXCoil(DXCoilNum)%MSGamma_Rated(I) .LT. 0.0d0) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...'//TRIM(cNumericFields(11+(I-1)*12))//' cannot be < 0.0, '//  &
         'entered value=['//trim(TrimSigDigits(DXCoil(DXCoilNum)%MSGamma_Rated(I),4))//'].')
      ErrorsFound=.TRUE.
    END IF
    DXCoil(DXCoilNum)%MSMaxONOFFCyclesperHour(I) = Numbers(12+(I-1)*12)
    IF (DXCoil(DXCoilNum)%Gamma_Rated(I) .LT. 0.0d0) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...'//TRIM(cNumericFields(12+(I-1)*12))//' cannot be < 0.0, '//  &
         'entered value=['//trim(TrimSigDigits(DXCoil(DXCoilNum)%MSMaxONOFFCyclesperHour(I),2))//'].')
      ErrorsFound=.TRUE.
    END IF
    DXCoil(DXCoilNum)%MSLatentCapacityTimeConstant(I) = Numbers(13+(I-1)*12)
    IF (DXCoil(DXCoilNum)%Gamma_Rated(I) .LT. 0.0d0) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...'//TRIM(cNumericFields(13+(I-1)*12))//' cannot be < 0.0, '//  &
         'entered value=['//trim(TrimSigDigits(DXCoil(DXCoilNum)%MSLatentCapacityTimeConstant(I),2))//'].')
      ErrorsFound=.TRUE.
    END IF

    DXCoil(DXCoilNum)%MSWasteHeatFrac(I) = Numbers(14+(I-1)*12)

    ! Read waste heat modifier curve name
    DXCoil(DXCoilNum)%MSWasteHeat(I) = GetCurveIndex(Alphas(18+(I-1)*6)) ! convert curve name to number
    IF (DXCoil(DXCoilNum)%MSWasteHeat(I) .EQ. 0) THEN
      IF (lAlphaBlanks(18+(I-1)*6)) THEN
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
        CALL ShowContinueError('...required '//trim(cAlphaFields(18+(I-1)*6))//' is blank.')
      ELSE
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...not found '//TRIM(cAlphaFields(18+(I-1)*6))//'="'//TRIM(Alphas(18+(I-1)*6))//'".')
      END IF
      ErrorsFound = .TRUE.
    ELSE
      ! Verify Curve Object, only legal types are BiQuadratic
      SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%MSWasteHeat(I)))

      CASE('BIQUADRATIC')
        CurveVal = CurveValue(DXCoil(DXCoilNum)%MSWasteHeat(I),RatedOutdoorAirTemp,RatedInletAirTemp)
        IF(CurveVal .GT. 1.10d0 .OR. CurveVal .LT. 0.90d0) THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", curve values')
          CALL ShowContinueError('...'//TRIM(cAlphaFields(18+(I-1)*6))//' output is not equal to 1.0 '//  &
                                             '(+ or - 10%) at rated conditions.')
          CALL ShowContinueError('...Curve output at rated conditions = '//TRIM(TrimSigDigits(CurveVal,3)))
        END IF

      CASE DEFAULT
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(18+(I-1)*6))//' type for this object = '// &
                             TRIM(GetCurveType(DXCoil(DXCoilNum)%MSWasteHeat(I))))
        CALL ShowContinueError('Curve type must be BiQuadratic.')
        ErrorsFound=.TRUE.
      END SELECT
    END IF

    DXCoil(DXCoilNum)%MSEvapCondEffect(I) = Numbers(15+(I-1)*12)
    IF (DXCoil(DXCoilNum)%MSEvapCondEffect(I) .LT. 0.0d0 .OR. DXCoil(DXCoilNum)%MSEvapCondEffect(I) .GT. 1.0d0) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...'//TRIM(cNumericFields(15+(I-1)*12))//' cannot be < 0.0 or > 1.0, '//  &
         'entered value=['//trim(TrimSigDigits(Numbers(15+(I-1)*12),3))//'].')
      ErrorsFound = .TRUE.
    END IF

    DXCoil(DXCoilNum)%MSEvapCondAirFlow(I) = Numbers(16+(I-1)*12)
    IF (DXCoil(DXCoilNum)%MSEvapCondAirFlow(I) .LT. 0.0 .AND. DXCoil(DXCoilNum)%MSEvapCondAirFlow(I) /= AutoSize) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...'//TRIM(cNumericFields(16+(I-1)*12))//' cannot be < 0.0, '//  &
         'entered value=['//trim(TrimSigDigits(Numbers(16+(I-1)*12),3))//'].')
      ErrorsFound = .TRUE.
    END IF

    DXCoil(DXCoilNum)%MSEvapCondPumpElecNomPower(I) = Numbers(17+(I-1)*12)
    IF (DXCoil(DXCoilNum)%MSEvapCondPumpElecNomPower(I) .LT. 0.0 .AND. &
        DXCoil(DXCoilNum)%MSEvapCondPumpElecNomPower(I) /= AutoSize) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...'//TRIM(cNumericFields(17+(I-1)*12))//' cannot be < 0.0, '//  &
         'entered value=['//trim(TrimSigDigits(Numbers(17+(I-1)*12),3))//'].')
      ErrorsFound = .TRUE.
    END IF

  END DO

END DO

IF (ErrorsFound) THEN
  CALL ShowFatalError(RoutineName//'Errors found in getting '//TRIM(CurrentModuleObject)//' input.  '//&
                      'Preceding condition(s) causes termination.')
END IF

! DX multispeed heating coil
CurrentModuleObject='Coil:Heating:DX:MultiSpeed'
DO DXCoilIndex = 1,NumDXMulSpeedHeatCoils

  DXCoilNum = DXCoilNum+1

  CALL GetObjectItem(TRIM(CurrentModuleObject),DXCoilIndex,Alphas,NumAlphas,Numbers,NumNumbers,IOStatus, &
                     NumBlank=lNumericBlanks,AlphaBlank=lAlphaBlanks, &
                     AlphaFieldNames=cAlphaFields,NumericFieldNames=cNumericFields)

  IsNotOK=.FALSE.
  IsBlank=.FALSE.
  CALL VerifyName(Alphas(1),DXCoil%Name,DXCoilNum-1,IsNotOK,IsBlank,TRIM(CurrentModuleObject)//' Name')
  IF (IsNotOK) THEN
    ErrorsFound=.true.
    IF (IsBlank) Alphas(1)='xxxxx'
  ENDIF
  DXCoil(DXCoilNum)%Name = Alphas(1)
! Initialize DataHeatBalance heat reclaim variable name for use by heat reclaim coils
  HeatReclaimDXCoil(DXCoilNum)%Name = DXCoil(DXCoilNum)%Name
  HeatReclaimDXCoil(DXCoilNum)%SourceType = TRIM(CurrentModuleObject)
  DXCoil(DXCoilNum)%DXCoilType = TRIM(CurrentModuleObject)
  DXCoil(DXCoilNum)%DXCoilType_Num = CoilDX_MultiSpeedHeating
  DXCoil(DXCoilNum)%Schedule = Alphas(2)
  DXCoil(DXCoilNum)%SchedPtr = GetScheduleIndex(Alphas(2))  ! convert schedule name to pointer
  IF (DXCoil(DXCoilNum)%SchedPtr .EQ. 0) THEN
    IF (lAlphaBlanks(2)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(2))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...'//TRIM(cAlphaFields(2))//'="'//TRIM(Alphas(2))//'".')
    END IF
    ErrorsFound=.TRUE.
  END IF

  DXCoil(DXCoilNum)%AirInNode = &
               GetOnlySingleNode(Alphas(3),ErrorsFound,TRIM(CurrentModuleObject),Alphas(1), &
               NodeType_Air,NodeConnectionType_Inlet,1,ObjectIsNotParent)

  DXCoil(DXCoilNum)%AirOutNode = &
               GetOnlySingleNode(Alphas(4),ErrorsFound,TRIM(CurrentModuleObject),Alphas(1), &
               NodeType_Air,NodeConnectionType_Outlet,1,ObjectIsNotParent)

  CALL TestCompSet(TRIM(CurrentModuleObject),Alphas(1),Alphas(3),Alphas(4),'Air Nodes')

  !Set minimum OAT for heat pump compressor operation
  DXCoil(DXCoilNum)%MinOATCompressor = Numbers(1)

  !Set crankcase heater capacity
  DXCoil(DXCoilNum)%CrankcaseHeaterCapacity = Numbers(2)
  IF (DXCoil(DXCoilNum)%CrankcaseHeaterCapacity .LT. 0.0d0) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//TRIM(cNumericFields(2))//' cannot be < 0.0,'//  &
        ' entered value=['//trim(TrimSigDigits(Numbers(2),2))//'].')
    ErrorsFound = .TRUE.
  END IF

  !Set crankcase heater cutout temperature
  DXCoil(DXCoilNum)%MaxOATCrankcaseHeater = Numbers(3)

! Only required for reverse cycle heat pumps
  DXCoil(DXCoilNum)%DefrostEIRFT = GetCurveIndex(Alphas(5)) ! convert curve name to number
  IF (SameString(Alphas(6),'ReverseCycle')) THEN
    IF (DXCoil(DXCoilNum)%DefrostEIRFT .EQ. 0) THEN
      IF (lAlphaBlanks(5)) THEN
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
        CALL ShowContinueError('...required '//trim(cAlphaFields(5))//' is blank.')
      ELSE
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...not found '//TRIM(cAlphaFields(5))//'="'//TRIM(Alphas(5))//'".')
      END IF
      ErrorsFound = .TRUE.
    ELSE
      ! Verify Curve Object, only legal type is BiQuadratic
      SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%DefrostEIRFT))
        CASE('BIQUADRATIC')
        CurveVal = CurveValue(DXCoil(DXCoilNum)%DefrostEIRFT,RatedInletWetbulbTempHeat,RatedOutdoorAirTempHeat)
        IF(CurveVal .GT. 1.10d0 .OR. CurveVal .LT. 0.90d0) THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", curve values')
          CALL ShowContinueError('...'//TRIM(cAlphaFields(5))//' output is not equal to 1.0 '//  &
                                             '(+ or - 10%) at rated conditions.')
          CALL ShowContinueError('...Curve output at rated conditions = '//TRIM(TrimSigDigits(CurveVal,3)))
        END IF

        CASE DEFAULT
          CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
          CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(5))//' type for this object = '// &
                             TRIM(GetCurveType(DXCoil(DXCoilNum)%DefrostEIRFT)))
        CALL ShowContinueError('Curve type must be BiQuadratic.')
        ErrorsFound=.TRUE.
      END SELECT
    END IF
  END IF

  IF (SameString(Alphas(6),'ReverseCycle'))  DXCoil(DXCoilNum)%DefrostStrategy = ReverseCycle
  IF (SameString(Alphas(6),'Resistive')) DXCoil(DXCoilNum)%DefrostStrategy = Resistive
  IF (DXCoil(DXCoilNum)%DefrostStrategy .EQ.0) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(6))//'="'//TRIM(Alphas(6))//'".')
    CALL ShowContinueError('...valid values for this field are ReverseCycle or Resistive.')
    ErrorsFound = .TRUE.
  END IF

  IF (SameString(Alphas(7),'Timed'))  DXCoil(DXCoilNum)%DefrostControl = Timed
  IF (SameString(Alphas(7),'OnDemand')) DXCoil(DXCoilNum)%DefrostControl = OnDemand
  IF (DXCoil(DXCoilNum)%DefrostControl .EQ.0) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(7))//'="'//TRIM(Alphas(7))//'".')
    CALL ShowContinueError('...valid values for this field are Timed or OnDemand.')
    ErrorsFound = .TRUE.
  END IF

  !Set maximum outdoor temp for defrost to occur
  DXCoil(DXCoilNum)%MaxOATDefrost = Numbers(4)

  !Set defrost time period
  DXCoil(DXCoilNum)%DefrostTime = Numbers(5)
  IF(DXCoil(DXCoilNum)%DefrostTime .EQ. 0.0d0 .AND. DXCoil(DXCoilNum)%DefrostControl .EQ. 1) THEN
    CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", ')
    CALL ShowContinueError('...'//TRIM(cNumericFields(5))//' = 0.0 for defrost control = TIMED.')
  END IF

  !Set defrost capacity (for resistive defrost)
  DXCoil(DXCoilNum)%DefrostCapacity = Numbers(6)
  IF(DXCoil(DXCoilNum)%DefrostCapacity .EQ. 0.0d0 .AND. DXCoil(DXCoilNum)%DefrostStrategy .EQ. 2) THEN
    CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", ')
    CALL ShowContinueError('...'//TRIM(cNumericFields(6))//' = 0.0 for defrost strategy = RESISTIVE.')
  END IF

  IF (SameString(Alphas(8),'Yes')) THEN
    DXCoil(DXCoilNum)%PLRImpact = .TRUE.
  ELSE IF (SameString(Alphas(8),'No')) THEN
    DXCoil(DXCoilNum)%PLRImpact = .FALSE.
  ELSE
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError(',,,invalid choice for '//TRIM(cAlphaFields(8))//'.  Entered choice = '//TRIM(Alphas(8)))
    CALL ShowContinueError('The allowed choices are Yes or No.')
    ErrorsFound=.TRUE.
  END IF

  !A10; \field Fuel type
  IF (SameString(Alphas(9),'Electricity')) THEN
    DXCoil(DXCoilNum)%FuelType = FuelTypeElectricity
  ELSE IF (SameString(Alphas(9),'NaturalGas')) THEN
    DXCoil(DXCoilNum)%FuelType = FuelTypeNaturalGas
  ELSE IF (SameString(Alphas(9),'PropaneGas')) THEN
    DXCoil(DXCoilNum)%FuelType = FuelTypePropaneGas
  ELSE IF (SameString(Alphas(9),'Diesel')) THEN
    DXCoil(DXCoilNum)%FuelType = FuelTypeDiesel
  ELSE IF (SameString(Alphas(9),'Gasoline')) THEN
    DXCoil(DXCoilNum)%FuelType = FuelTypeGasoline
  ELSE IF (SameString(Alphas(9),'FuelOil#1')) THEN
    DXCoil(DXCoilNum)%FuelType = FuelTypeFuelOil1
  ELSE IF (SameString(Alphas(9),'FuelOil#2')) THEN
    DXCoil(DXCoilNum)%FuelType = FuelTypeFuelOil2
  ELSE
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError(',,,invalid choice for '//TRIM(cAlphaFields(9))//'.  Entered choice = '//TRIM(Alphas(9)))
    CALL ShowContinueError('Valid choices are Electricity, NaturalGas, PropaneGas, Diesel, Gasoline, FuelOil#1 or FuelOil#2')
    ErrorsFound=.TRUE.
  END IF

  DXCoil(DXCoilNum)%NumOfSpeeds = Numbers(7)       ! Number of speeds
  IF (DXCoil(DXCoilNum)%NumOfSpeeds .LT. 2) THEN
    CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
    CALL ShowContinueError('...'//TRIM(cNumericFields(7))//' must be >= 2.'//  &
                                         ' entered number is '//TRIM(TrimSigDigits(Numbers(7),0)))
    ErrorsFound=.TRUE.
  END IF

  ! Allocate arrays based on the number of speeds
  ALLOCATE(DXCoil(DXCoilNum)%MSErrIndex(DXCoil(DXCoilNum)%NumOfSpeeds))
  DXCoil(DXCoilNum)%MSErrIndex=0
  ALLOCATE(DXCoil(DXCoilNum)%MSRatedTotCap(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSRatedCOP(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSRatedAirMassFlowRate(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSCCapFTemp(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSCCapFFlow(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSEIRFTemp(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSEIRFFlow(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSWasteHeat(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSPLFFPLR(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSRatedCBF(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSWasteHeatFrac(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSTotCapTempModFacCurveType(DXCoil(DXCoilNum)%NumOfSpeeds))
  ALLOCATE(DXCoil(DXCoilNum)%MSEIRTempModFacCurveType(DXCoil(DXCoilNum)%NumOfSpeeds))

  DXCoil(DXCoilNum)%RatedSHR(1) = 1.0d0

  DO I=1,DXCoil(DXCoilNum)%NumOfSpeeds

    DXCoil(DXCoilNum)%MSRatedTotCap(I) = Numbers(8+(I-1)*4)
    DXCoil(DXCoilNum)%MSRatedCOP(I) = Numbers(9+(I-1)*4)
    DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(I) = Numbers(10+(I-1)*4)
    DXCoil(DXCoilNum)%MSWasteHeatFrac(I) = Numbers(11+(I-1)*4)

    DXCoil(DXCoilNum)%MSCCapFTemp(I) = GetCurveIndex(Alphas(10+(I-1)*6)) ! convert curve name to number
    IF (DXCoil(DXCoilNum)%MSCCapFTemp(I) .EQ. 0) THEN
      CALL ShowSevereError(TRIM(CurrentModuleObject)//', "'//TRIM(DXCoil(DXCoilNum)%Name)//&
                           '" '//TRIM(cAlphaFields(10+(I-1)*6))//' not found:'//TRIM(Alphas(10+(I-1)*6)))
      ErrorsFound = .TRUE.
    ELSE
      ! only legal types are Quadratic, Biquadratic and Cubic
      SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%MSCCapFTemp(I)))

      CASE('QUADRATIC')
        DXCoil(DXCoilNum)%MSTotCapTempModFacCurveType(I) = Quadratic
        CurveVal = CurveValue(DXCoil(DXCoilNum)%MSCCapFTemp(I),RatedOutdoorAirTempHeat)
        IF(CurveVal .GT. 1.10d0 .OR. CurveVal .LT. 0.90d0) THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", curve values')
          CALL ShowContinueError('...'//TRIM(cAlphaFields(10+(I-1)*6))//' output is not equal to 1.0 '//  &
                                             '(+ or - 10%) at rated conditions.')
          CALL ShowContinueError('...Curve output at rated conditions = '//TRIM(TrimSigDigits(CurveVal,3)))
        END IF

      CASE('BIQUADRATIC')
        DXCoil(DXCoilNum)%MSTotCapTempModFacCurveType(I) = Biquadratic
        CurveVal = CurveValue(DXCoil(DXCoilNum)%MSCCapFTemp(I),RatedInletAirTempHeat,RatedOutdoorAirTempHeat)
        IF(CurveVal .GT. 1.10d0 .OR. CurveVal .LT. 0.90d0) THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", curve values')
          CALL ShowContinueError('...'//TRIM(cAlphaFields(10+(I-1)*6))//' output is not equal to 1.0 '//  &
                                             '(+ or - 10%) at rated conditions.')
          CALL ShowContinueError('...Curve output at rated conditions = '//TRIM(TrimSigDigits(CurveVal,3)))
        END IF

      CASE('CUBIC')
        DXCoil(DXCoilNum)%MSTotCapTempModFacCurveType(I) = Cubic
        CurveVal = CurveValue(DXCoil(DXCoilNum)%MSCCapFTemp(I),RatedOutdoorAirTempHeat)
        IF(CurveVal .GT. 1.10d0 .OR. CurveVal .LT. 0.90d0) THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", curve values')
          CALL ShowContinueError('...'//TRIM(cAlphaFields(10+(I-1)*6))//' output is not equal to 1.0 '//  &
                                             '(+ or - 10%) at rated conditions.')
          CALL ShowContinueError('...Curve output at rated conditions = '//TRIM(TrimSigDigits(CurveVal,3)))
        END IF

      CASE DEFAULT
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(10+(I-1)*6))//' type for this object = '// &
                             TRIM(GetCurveType(DXCoil(DXCoilNum)%MSCCapFTemp(I))))
        CALL ShowContinueError('Curve type must be BiQuadratic, Quadratic or Cubic.')
        ErrorsFound=.TRUE.
      END SELECT
    END IF

    DXCoil(DXCoilNum)%MSCCapFFlow(I) = GetCurveIndex(Alphas(11+(I-1)*6)) ! convert curve name to number
    IF (DXCoil(DXCoilNum)%MSCCapFFlow(I) .EQ. 0) THEN
      IF (lAlphaBlanks(11+(I-1)*6)) THEN
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
        CALL ShowContinueError('...required '//trim(cAlphaFields(11+(I-1)*6))//' is blank.')
      ELSE
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...not found '//TRIM(cAlphaFields(11+(I-1)*6))//'="'//TRIM(Alphas(11+(I-1)*6))//'".')
      END IF
      ErrorsFound = .TRUE.
    ELSE
      ! Verify Curve Object, only legal type is Quadratic
      SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%MSCCapFFlow(I)))

      CASE('QUADRATIC')
        CurveVal = CurveValue(DXCoil(DXCoilNum)%MSCCapFFlow(I),1.0d0)
        IF(CurveVal .GT. 1.10d0 .OR. CurveVal .LT. 0.90d0) THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", curve values')
          CALL ShowContinueError('...'//TRIM(cAlphaFields(11+(I-1)*6))//' output is not equal to 1.0 '//  &
                                             '(+ or - 10%) at rated conditions.')
          CALL ShowContinueError('...Curve output at rated conditions = '//TRIM(TrimSigDigits(CurveVal,3)))
        END IF

      CASE('CUBIC')
        CurveVal = CurveValue(DXCoil(DXCoilNum)%MSCCapFFlow(I),1.0d0)
        IF(CurveVal .GT. 1.10d0 .OR. CurveVal .LT. 0.90d0) THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", curve values')
          CALL ShowContinueError('...'//TRIM(cAlphaFields(11+(I-1)*6))//' output is not equal to 1.0 '//  &
                                             '(+ or - 10%) at rated conditions.')
          CALL ShowContinueError('...Curve output at rated conditions = '//TRIM(TrimSigDigits(CurveVal,3)))
        END IF

      CASE DEFAULT
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(11+(I-1)*6))//' type for this object = '// &
                             TRIM(GetCurveType(DXCoil(DXCoilNum)%MSCCapFFlow(I))))
        CALL ShowContinueError('Curve type must be Quadratic or Cubic.')
        ErrorsFound=.TRUE.
      END SELECT
    END IF

    DXCoil(DXCoilNum)%MSEIRFTemp(I) = GetCurveIndex(Alphas(12+(I-1)*6)) ! convert curve name to number
    IF (DXCoil(DXCoilNum)%MSEIRFTemp(I) .EQ. 0) THEN
      IF (lAlphaBlanks(12+(I-1)*6)) THEN
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
        CALL ShowContinueError('...required '//trim(cAlphaFields(12+(I-1)*6))//' is blank.')
      ELSE
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...not found '//TRIM(cAlphaFields(12+(I-1)*6))//'="'//TRIM(Alphas(15+(I-1)*6))//'".')
      END IF
      ErrorsFound = .TRUE.
    ELSE
      ! only legal types are Quadratic, Biquadratic and Cubic
      SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%MSEIRFTemp(I)))

      CASE('QUADRATIC')
        DXCoil(DXCoilNum)%MSEIRTempModFacCurveType(I) = Quadratic
        CurveVal = CurveValue(DXCoil(DXCoilNum)%MSEIRFTemp(I),RatedOutdoorAirTempHeat)
        IF(CurveVal .GT. 1.10d0 .OR. CurveVal .LT. 0.90d0) THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", curve values')
          CALL ShowContinueError('...'//TRIM(cAlphaFields(12+(I-1)*6))//' output is not equal to 1.0 '//  &
                                             '(+ or - 10%) at rated conditions.')
          CALL ShowContinueError('...Curve output at rated conditions = '//TRIM(TrimSigDigits(CurveVal,3)))
        END IF

      CASE('BIQUADRATIC')
        DXCoil(DXCoilNum)%MSEIRTempModFacCurveType(I) = BiQuadratic
        CurveVal = CurveValue(DXCoil(DXCoilNum)%MSEIRFTemp(I),RatedInletAirTempHeat,RatedOutdoorAirTempHeat)
        IF(CurveVal .GT. 1.10d0 .OR. CurveVal .LT. 0.90d0) THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", curve values')
          CALL ShowContinueError('...'//TRIM(cAlphaFields(12+(I-1)*6))//' output is not equal to 1.0 '//  &
                                             '(+ or - 10%) at rated conditions.')
          CALL ShowContinueError('...Curve output at rated conditions = '//TRIM(TrimSigDigits(CurveVal,3)))
        END IF

      CASE('CUBIC')
        DXCoil(DXCoilNum)%MSEIRTempModFacCurveType(I) = Cubic
        CurveVal = CurveValue(DXCoil(DXCoilNum)%MSEIRFTemp(I),RatedOutdoorAirTempHeat)
        IF(CurveVal .GT. 1.10d0 .OR. CurveVal .LT. 0.90d0) THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", curve values')
          CALL ShowContinueError('...'//TRIM(cAlphaFields(12+(I-1)*6))//' output is not equal to 1.0 '//  &
                                             '(+ or - 10%) at rated conditions.')
          CALL ShowContinueError('...Curve output at rated conditions = '//TRIM(TrimSigDigits(CurveVal,3)))
        END IF

      CASE DEFAULT
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(12+(I-1)*6))//' type for this object = '// &
                             TRIM(GetCurveType(DXCoil(DXCoilNum)%MSEIRFTemp(I))))
        CALL ShowContinueError('Curve type must be BiQuadratic, Quadratic or Cubic.')
        ErrorsFound=.TRUE.
      END SELECT
    END IF

    DXCoil(DXCoilNum)%MSEIRFFlow(I) = GetCurveIndex(Alphas(13+(I-1)*6)) ! convert curve name to number
    IF (DXCoil(DXCoilNum)%MSEIRFFlow(I) .EQ. 0) THEN
      IF (lAlphaBlanks(13+(I-1)*6)) THEN
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
        CALL ShowContinueError('...required '//trim(cAlphaFields(13+(I-1)*6))//' is blank.')
      ELSE
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...not found '//TRIM(cAlphaFields(13+(I-1)*6))//'="'//TRIM(Alphas(13+(I-1)*6))//'".')
      END IF
      ErrorsFound = .TRUE.
    ELSE
      ! Verify Curve Object, only legal type is Quadratic
      SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%MSEIRFFlow(I)))

      CASE('QUADRATIC')
        CurveVal = CurveValue(DXCoil(DXCoilNum)%MSEIRFFlow(I),1.0d0)
        IF(CurveVal .GT. 1.10d0 .OR. CurveVal .LT. 0.90d0) THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", curve values')
          CALL ShowContinueError('...'//TRIM(cAlphaFields(13+(I-1)*6))//' output is not equal to 1.0 '//  &
                                             '(+ or - 10%) at rated conditions.')
          CALL ShowContinueError('...Curve output at rated conditions = '//TRIM(TrimSigDigits(CurveVal,3)))
        END IF

      CASE('CUBIC')
        CurveVal = CurveValue(DXCoil(DXCoilNum)%MSEIRFFlow(I),1.0d0)
        IF(CurveVal .GT. 1.10d0 .OR. CurveVal .LT. 0.90d0) THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", curve values')
          CALL ShowContinueError('...'//TRIM(cAlphaFields(13+(I-1)*6))//' output is not equal to 1.0 '//  &
                                             '(+ or - 10%) at rated conditions.')
          CALL ShowContinueError('...Curve output at rated conditions = '//TRIM(TrimSigDigits(CurveVal,3)))
        END IF

      CASE DEFAULT
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(13+(I-1)*6))//' type for this object = '// &
                             TRIM(GetCurveType(DXCoil(DXCoilNum)%MSEIRFFlow(I))))
        CALL ShowContinueError('Curve type must be Quadratic or Cubic.')
        ErrorsFound=.TRUE.
      END SELECT
    END IF

    DXCoil(DXCoilNum)%MSPLFFPLR(I) = GetCurveIndex(Alphas(14+(I-1)*6)) ! convert curve name to number
    IF (DXCoil(DXCoilNum)%MSPLFFPLR(I) .EQ. 0) THEN
      IF (lAlphaBlanks(14+(I-1)*6)) THEN
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
        CALL ShowContinueError('...required '//trim(cAlphaFields(14+(I-1)*6))//' is blank.')
      ELSE
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...not found '//TRIM(cAlphaFields(14+(I-1)*6))//'="'//TRIM(Alphas(14+(I-1)*6))//'".')
      END IF
      ErrorsFound = .TRUE.
    ELSE
      ! Verify Curve Object, only legal types are Quadratic or Cubic
      SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%MSPLFFPLR(I)))

      CASE('QUADRATIC')

      CASE('CUBIC')

      CASE DEFAULT
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(14+(I-1)*6))//' type for this object = '// &
                             TRIM(GetCurveType(DXCoil(DXCoilNum)%MSPLFFPLR(I))))
        CALL ShowContinueError('Curve type must be Quadratic or Cubic.')
        ErrorsFound=.TRUE.
      END SELECT

      IF(.NOT. ErrorsFound)THEN
!       Test PLF curve minimum and maximum. Cap if less than 0.7 or greater than 1.0.
        MinCurveVal = 999.0d0
        MaxCurveVal = -999.0d0
        CurveInput = 0.0d0
        DO WHILE (CurveInput <= 1.0d0)
          CurveVal = CurveValue(DXCoil(DXCoilNum)%MSPLFFPLR(I),CurveInput)
          IF(CurveVal .LT. MinCurveVal)THEN
            MinCurveVal = CurveVal
            MinCurvePLR = CurveInput
          END IF
          IF(CurveVal .GT. MaxCurveVal)THEN
            MaxCurveVal = CurveVal
            MaxCurvePLR = CurveInput
          END IF
          CurveInput=CurveInput+0.01d0
        END DO
        IF(MinCurveVal .LT. 0.7d0)THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
          CALL ShowContinueError('...'//TRIM(cAlphaFields(14+(I-1)*6))//' = '//TRIM(Alphas(14+(I-1)*6))//' has out of range value.')
          CALL ShowContinueError('...Curve minimum must be >= 0.7, '// &
                        'curve min at PLR = '//TRIM(TrimSigDigits(MinCurvePLR,2))//' is '//TRIM(TrimSigDigits(MinCurveVal,3)))
          CALL ShowContinueError('...Setting curve minimum to 0.7 and simulation continues.')
          CALL SetCurveOutputMinMaxValues(DXCoil(DXCoilNum)%PLFFPLR(1),ErrorsFound,CurveMin=0.7d0)
          CALL SetCurveOutputMinMaxValues(DXCoil(DXCoilNum)%MSPLFFPLR(I),ErrorsFound,CurveMin=0.7d0)
        END IF

        IF(MaxCurveVal .GT. 1.0d0)THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
          CALL ShowContinueError('...'//TRIM(cAlphaFields(14+(I-1)*6))//' = '//TRIM(Alphas(14+(I-1)*6))//' has out of range value.')
          CALL ShowContinueError('...Curve maximum must be <= 1.0, '// &
                        'curve max at PLR = '//TRIM(TrimSigDigits(MaxCurvePLR,2))//' is '//TRIM(TrimSigDigits(MaxCurveVal,3)))
          CALL ShowContinueError('...Setting curve maximum to 1.0 and simulation continues.')
          CALL SetCurveOutputMinMaxValues(DXCoil(DXCoilNum)%MSPLFFPLR(I),ErrorsFound,CurveMax=1.0d0)
        END IF

      END IF

    END IF

    ! Read waste heat modifier curve name
    DXCoil(DXCoilNum)%MSWasteHeat(I) = GetCurveIndex(Alphas(15+(I-1)*6)) ! convert curve name to number
    IF (DXCoil(DXCoilNum)%MSWasteHeat(I) .EQ. 0) THEN
      IF (lAlphaBlanks(11)) THEN
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
        CALL ShowContinueError('...required '//trim(cAlphaFields(15+(I-1)*6))//' is blank.')
      ELSE
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...not found '//TRIM(cAlphaFields(15+(I-1)*6))//'="'//TRIM(Alphas(15+(I-1)*6))//'".')
      END IF
      ErrorsFound = .TRUE.
    ELSE
      ! Verify Curve Object, only legal types are BiQuadratic
      SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%MSWasteHeat(I)))

      CASE('BIQUADRATIC')
        CurveVal = CurveValue(DXCoil(DXCoilNum)%MSWasteHeat(I),RatedOutdoorAirTempHeat,RatedInletAirTempHeat)
        IF(CurveVal .GT. 1.10d0 .OR. CurveVal .LT. 0.90d0) THEN
          CALL ShowWarningError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", curve values')
          CALL ShowContinueError('...'//TRIM(cAlphaFields(15+(I-1)*6))//' output is not equal to 1.0 '//  &
                                             '(+ or - 10%) at rated conditions.')
          CALL ShowContinueError('...Curve output at rated conditions = '//TRIM(TrimSigDigits(CurveVal,3)))
        END IF

      CASE DEFAULT
        CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
        CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(15+(I-1)*6))//' type for this object = '// &
                             TRIM(GetCurveType(DXCoil(DXCoilNum)%MSWasteHeat(I))))
        CALL ShowContinueError('Curve type must be BiQuadratic.')
        ErrorsFound=.TRUE.
      END SELECT
    END IF

  END DO
END DO

! Loop over the VRF Cooling Coils and get & load the data
CurrentModuleObject=cAllCoilTypes(CoilVRF_Cooling)
DO DXCoilIndex = 1, NumVRFCoolingCoils

  CALL GetObjectItem(TRIM(CurrentModuleObject),DXCoilIndex,Alphas,NumAlphas,Numbers,NumNumbers,IOStatus, &
                     NumBlank=lNumericBlanks,AlphaBlank=lAlphaBlanks, &
                     AlphaFieldNames=cAlphaFields,NumericFieldNames=cNumericFields)

  DXCoilNum = DXCoilNum+1
  IsNotOK=.FALSE.
  IsBlank=.FALSE.
  CALL VerifyName(Alphas(1),DXCoil%Name,DXCoilNum-1,IsNotOK,IsBlank,TRIM(CurrentModuleObject)//' Name')
  IF (IsNotOK) THEN
    ErrorsFound=.true.
    IF (IsBlank) Alphas(1)='xxxxx'
  ENDIF
  DXCoil(DXCoilNum)%Name = Alphas(1)
  DXCoil(DXCoilNum)%DXCoilType = TRIM(CurrentModuleObject)
  DXCoil(DXCoilNum)%DXCoilType_Num = CoilVRF_Cooling
  DXCoil(DXCoilNum)%Schedule = Alphas(2)
  DXCoil(DXCoilNum)%SchedPtr = GetScheduleIndex(Alphas(2))  ! convert schedule name to pointer
  IF (DXCoil(DXCoilNum)%SchedPtr .EQ. 0) THEN
    IF (lAlphaBlanks(2)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(2))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...'//TRIM(cAlphaFields(2))//'="'//TRIM(Alphas(2))//'".')
    END IF
    ErrorsFound=.TRUE.
  END IF
  DXCoil(DXCoilNum)%RatedTotCap(1) = Numbers(1)
  DXCoil(DXCoilNum)%RatedSHR(1)    = Numbers(2)
  DXCoil(DXCoilNum)%RatedAirVolFlowRate(1) = Numbers(3)


  DXCoil(DXCoilNum)%CCapFTemp(1)     = GetCurveIndex(Alphas(3))
  ! Verify Curve Object, only legal type is Linear, Quadratic, Cubic, or BiQuadratic
    SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%CCapFTemp(1)))

    CASE('LINEAR')
      DXCoil(DXCoilNum)%TotCapTempModFacCurveType(1) = Linear

    CASE('QUADRATIC')
      DXCoil(DXCoilNum)%TotCapTempModFacCurveType(1) = Quadratic

    CASE('CUBIC')
      DXCoil(DXCoilNum)%TotCapTempModFacCurveType(1) = Cubic

    CASE('BIQUADRATIC')
      DXCoil(DXCoilNum)%TotCapTempModFacCurveType(1) = BiQuadratic

    CASE DEFAULT
          CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
          CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(3))//' type for this object = '// &
                          TRIM(GetCurveType(DXCoil(DXCoilNum)%CCapFTemp(1))))
      CALL ShowContinueError('... Curve type must be Linear, Quadratic, Cubic, or BiQuadratic.')
      ErrorsFound=.TRUE.
    END SELECT

  DXCoil(DXCoilNum)%CCapFFlow(1) = GetCurveIndex(Alphas(4)) ! convert curve name to number
  IF (DXCoil(DXCoilNum)%CCapFFlow(1) .EQ. 0) THEN
    IF (lAlphaBlanks(4)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(4))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(4))//'="'//TRIM(Alphas(4))//'".')
    END IF
    ErrorsFound = .TRUE.
  ELSE
    ! Verify Curve Object, only legal type is Linear, Quadratic or Cubic
    SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%CCapFFlow(1)))

    CASE('LINEAR', 'QUADRATIC', 'CUBIC')

    CASE DEFAULT
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(4))//' type for this object = '// &
                            TRIM(GetCurveType(DXCoil(DXCoilNum)%CCapFFlow(1))))
      CALL ShowContinueError('... Curve type must be Linear, Quadratic or Cubic.')
      ErrorsFound=.TRUE.
    END SELECT
  END IF

  DXCoil(DXCoilNum)%AirInNode = &
               GetOnlySingleNode(Alphas(5),ErrorsFound,TRIM(CurrentModuleObject),Alphas(1), &
               NodeType_Air,NodeConnectionType_Inlet,1,ObjectIsNotParent)

  DXCoil(DXCoilNum)%AirOutNode = &
               GetOnlySingleNode(Alphas(6),ErrorsFound,TRIM(CurrentModuleObject),Alphas(1), &
               NodeType_Air,NodeConnectionType_Outlet,1,ObjectIsNotParent)

  CALL TestCompSet(TRIM(CurrentModuleObject),Alphas(1),Alphas(5),Alphas(6),'Air Nodes')

  DXCoil(DXCoilNum)%CondensateCollectName = Alphas(7)
  IF (lAlphaBlanks(7)) THEN
    DXCoil(DXCoilNum)%CondensateCollectMode = CondensateDiscarded
  ELSE
    DXCoil(DxCoilNum)%CondensateCollectMode = CondensateToTank
    CALL SetupTankSupplyComponent(DXCoil(DXCoilNum)%Name, TRIM(CurrentModuleObject), &
                 DXCoil(DXCoilNum)%CondensateCollectName, ErrorsFound, DXCoil(DXCoilNum)%CondensateTankID, &
                 DXCoil(DXCoilNum)%CondensateTankSupplyARRID )
  END IF

END DO

! Loop over the VRF Heating Coils and get & load the data
CurrentModuleObject='COIL:Heating:DX:VariableRefrigerantFlow'
DO DXCoilIndex = 1, NumVRFHeatingCoils

  CALL GetObjectItem(TRIM(CurrentModuleObject),DXCoilIndex,Alphas,NumAlphas,Numbers,NumNumbers,IOStatus, &
                     NumBlank=lNumericBlanks,AlphaBlank=lAlphaBlanks, &
                     AlphaFieldNames=cAlphaFields,NumericFieldNames=cNumericFields)

  DXCoilNum = DXCoilNum+1
  IsNotOK=.FALSE.
  IsBlank=.FALSE.
  CALL VerifyName(Alphas(1),DXCoil%Name,DXCoilNum-1,IsNotOK,IsBlank,TRIM(CurrentModuleObject)//' Name')
  IF (IsNotOK) THEN
    ErrorsFound=.true.
    IF (IsBlank) Alphas(1)='xxxxx'
  ENDIF
  DXCoil(DXCoilNum)%Name = Alphas(1)
  DXCoil(DXCoilNum)%DXCoilType = TRIM(CurrentModuleObject)
  DXCoil(DXCoilNum)%DXCoilType_Num = CoilVRF_Heating
  DXCoil(DXCoilNum)%Schedule = Alphas(2)
  DXCoil(DXCoilNum)%SchedPtr = GetScheduleIndex(Alphas(2))  ! convert schedule name to pointer
  IF (DXCoil(DXCoilNum)%SchedPtr .EQ. 0) THEN
    IF (lAlphaBlanks(2)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(2))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...'//TRIM(cAlphaFields(2))//'="'//TRIM(Alphas(2))//'".')
    END IF
    ErrorsFound=.TRUE.
  END IF
  DXCoil(DXCoilNum)%RatedTotCap(1) = Numbers(1)
  DXCoil(DXCoilNum)%RatedAirVolFlowRate(1) = Numbers(2)

  DXCoil(DXCoilNum)%AirInNode = &
               GetOnlySingleNode(Alphas(3),ErrorsFound,TRIM(CurrentModuleObject),Alphas(1), &
               NodeType_Air,NodeConnectionType_Inlet,1,ObjectIsNotParent)

  DXCoil(DXCoilNum)%AirOutNode = &
               GetOnlySingleNode(Alphas(4),ErrorsFound,TRIM(CurrentModuleObject),Alphas(1), &
               NodeType_Air,NodeConnectionType_Outlet,1,ObjectIsNotParent)

  CALL TestCompSet(TRIM(CurrentModuleObject),Alphas(1),Alphas(3),Alphas(4),'Air Nodes')

  DXCoil(DXCoilNum)%CCapFTemp = GetCurveIndex(Alphas(5))
  IF (DXCoil(DXCoilNum)%CCapFTemp(1) .EQ. 0) THEN
    IF (lAlphaBlanks(5)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(5))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(5))//'="'//TRIM(Alphas(5))//'".')
    END IF
    ErrorsFound = .TRUE.
  ELSE
    ! Verify Curve Object, only legal type is Quadratic
    SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%CCapFTemp(1)))

    CASE('LINEAR')
      DXCoil(DXCoilNum)%TotCapTempModFacCurveType(1) = Linear
    CASE('QUADRATIC')
      DXCoil(DXCoilNum)%TotCapTempModFacCurveType(1) = Quadratic
    CASE('CUBIC')
      DXCoil(DXCoilNum)%TotCapTempModFacCurveType(1) = Cubic
    CASE('BIQUADRATIC')
      DXCoil(DXCoilNum)%TotCapTempModFacCurveType(1) = Biquadratic
    CASE DEFAULT
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(5))//' type for this object = '// &
                            TRIM(GetCurveType(DXCoil(DXCoilNum)%CCapFTemp(1))))
      CALL ShowContinueError('... Curve type must be Linear, Quadratic, Cubic or BiQuadratic.')
      ErrorsFound=.TRUE.
    END SELECT
  END IF

  DXCoil(DXCoilNum)%CCapFFlow(1) = GetCurveIndex(Alphas(6)) ! convert curve name to number
  IF (DXCoil(DXCoilNum)%CCapFFlow(1) .EQ. 0) THEN
    IF (lAlphaBlanks(6)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(6))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...not found '//TRIM(cAlphaFields(6))//'="'//TRIM(Alphas(6))//'".')
    END IF
    ErrorsFound = .TRUE.
  ELSE
    ! Verify Curve Object, only legal type is Quadratic
    SELECT CASE(GetCurveType(DXCoil(DXCoilNum)%CCapFFlow(1)))

    CASE('LINEAR', 'QUADRATIC', 'CUBIC')

    CASE DEFAULT
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...illegal '//TRIM(cAlphaFields(5))//' type for this object = '// &
                            TRIM(GetCurveType(DXCoil(DXCoilNum)%CCapFFlow(1))))
      CALL ShowContinueError('... Curve type must be linear, Quadratic or Cubic.')
      ErrorsFound=.TRUE.
    END SELECT
  END IF

END DO

!RS: Comment: Looping over HPSim coil 'data'
CurrentModuleObject='ZoneHVAC:HPSIM'
DO DXCoilIndex = 1, DXCoilHPSimNum2 !RS: Debugging: Moving to an air loop instead of packaged (10/31/14)

  CALL GetObjectItem(TRIM(CurrentModuleObject),DXCoilIndex,Alphas,NumAlphas,Numbers,NumNumbers,IOStatus, &
                     NumBlank=lNumericBlanks,AlphaBlank=lAlphaBlanks, &
                     AlphaFieldNames=cAlphaFields,NumericFieldNames=cNumericFields)

  DXCoilNum = DXCoilNum+1
  IsNotOK=.FALSE.
  IsBlank=.FALSE.
  CALL VerifyName(Alphas(1),DXCoil%Name,DXCoilNum-1,IsNotOK,IsBlank,TRIM(CurrentModuleObject)//' Name')
  IF (IsNotOK) THEN
    ErrorsFound=.true.
    IF (IsBlank) Alphas(1)='xxxxx'
  ENDIF
  DXCoil(DXCoilNum)%Name = Alphas(1)
! Initialize DataHeatBalance heat reclaim variable name for use by heat reclaim coils
  !HeatReclaimDXCoil(DXCoilNum)%Name = DXCoil(DXCoilNum)%Name
  !HeatReclaimDXCoil(DXCoilNum)%SourceType = TRIM(CurrentModuleObject)
  !DXCoil(DXCoilNum)%DXCoilType = TRIM(CurrentModuleObject)
  !DXCoil(DXCoilNum)%DXCoilType_Num = CoilDX_CoolingTwoStageWHumControl
  DXCoil(DXCoilNum)%DXCoilType = Alphas(9) !TRIM(CurrentModuleObject)
  DXCoil(DXCoilNum)%DXCoilType_Num = CoilDX_HPSim !CoolingSingleSpeed
  DXCoil(DXCoilNum)%Schedule = Alphas(2)
  DXCoil(DXCoilNum)%SchedPtr = GetScheduleIndex(Alphas(2))  ! convert schedule name to pointer
  DXCoil(DXCoilNum)%RatedTotCap(1) = Numbers(1)
  !DXCoil(DXCoilNum)%MinOATCompressor=-5 !RS: Debugging: Setting this here to test rather than adding to IDF (8/30/14)
  IF (DXCoil(DXCoilNum)%SchedPtr .EQ. 0) THEN
    IF (lAlphaBlanks(2)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(2))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...'//TRIM(cAlphaFields(2))//'="'//TRIM(Alphas(2))//'".')
    END IF
    ErrorsFound=.TRUE.
  END IF

  DXCoil(DXCoilNum)%AirInNode = &
               GetOnlySingleNode(Alphas(3),ErrorsFound,TRIM(CurrentModuleObject),Alphas(1), &
               NodeType_Air,NodeConnectionType_Inlet,1,ObjectIsParent) !ObjectIsNotParent)  !RS: Debugging: Trying this...

  DXCoil(DXCoilNum)%AirOutNode = &
               GetOnlySingleNode(Alphas(4),ErrorsFound,TRIM(CurrentModuleObject),Alphas(1), &
               NodeType_Air,NodeConnectionType_Outlet,1,ObjectIsParent) !ObjectIsNotParent)  !RS: Debugging: Trying this...
END DO

CurrentModuleObject='Coil:HPSim' !ZoneHVAC:HPSIM'
DO DXCoilIndex = 1, DXCoilHPSimNum  !RS: Debugging: Moving to an air loop instead of packaged (10/31/14)

  CALL GetObjectItem(TRIM(CurrentModuleObject),DXCoilIndex,Alphas,NumAlphas,Numbers,NumNumbers,IOStatus, &
                     NumBlank=lNumericBlanks,AlphaBlank=lAlphaBlanks, &
                     AlphaFieldNames=cAlphaFields,NumericFieldNames=cNumericFields)

  DXCoilNum = DXCoilNum+1
  IsNotOK=.FALSE.
  IsBlank=.FALSE.
  CALL VerifyName(Alphas(1),DXCoil%Name,DXCoilNum-1,IsNotOK,IsBlank,TRIM(CurrentModuleObject)//' Name')
  IF (IsNotOK) THEN
    ErrorsFound=.true.
    IF (IsBlank) Alphas(1)='xxxxx'
  ENDIF
  DXCoil(DXCoilNum)%Name = Alphas(6) !5) !RS: Debugging: We actually want the "coil" name (4/28/14) !(1)
! Initialize DataHeatBalance heat reclaim variable name for use by heat reclaim coils
  !HeatReclaimDXCoil(DXCoilNum)%Name = DXCoil(DXCoilNum)%Name
  !HeatReclaimDXCoil(DXCoilNum)%SourceType = TRIM(CurrentModuleObject)
  DXCoil(DXCoilNum)%DXCoilType = Alphas(6) !5) !TRIM(CurrentModuleObject)
  DXCoil(DXCoilNum)%DXCoilType_Num = CoilDX_HPSim !CoolingSingleSpeed
  DXCoil(DXCoilNum)%Schedule = Alphas(5) !4) !2)
  DXCoil(DXCoilNum)%SchedPtr = GetScheduleIndex(Alphas(5)) !4)) !2))  ! convert schedule name to pointer
  DXCoil(DXCoilNum)%RatedTotCap(1) = Numbers(1)
  !DXCoil(DXCoilNum)%MinOATCompressor=-5 !RS: Debugging: Setting this here to test rather than adding to IDF (8/30/14)
  IF (DXCoil(DXCoilNum)%SchedPtr .EQ. 0) THEN
    IF (lAlphaBlanks(2)) THEN
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", missing')
      CALL ShowContinueError('...required '//trim(cAlphaFields(2))//' is blank.')
    ELSE
      CALL ShowSevereError(RoutineName//trim(CurrentModuleObject)//'="'//trim(DXCoil(DXCoilNum)%Name)//'", invalid')
      CALL ShowContinueError('...'//TRIM(cAlphaFields(2))//'="'//TRIM(Alphas(2))//'".')
    END IF
    ErrorsFound=.TRUE.
  END IF

  DXCoil(DXCoilNum)%AirInNode = &
               GetOnlySingleNode(Alphas(2),ErrorsFound,TRIM(CurrentModuleObject),Alphas(1), &
               NodeType_Air,NodeConnectionType_Inlet,1,ObjectIsParent) !ObjectIsNotParent)  !RS: Debugging: Trying this...

  DXCoil(DXCoilNum)%AirOutNode = &
               GetOnlySingleNode(Alphas(3),ErrorsFound,TRIM(CurrentModuleObject),Alphas(1), &
               NodeType_Air,NodeConnectionType_Outlet,1,ObjectIsParent) !ObjectIsNotParent)  !RS: Debugging: Trying this...
END DO

IF (ErrorsFound) THEN
  CALL ShowFatalError(RoutineName//'Errors found in getting '//TRIM(CurrentModuleObject)//' input.  '//&
                      'Preceding condition(s) causes termination.')
END IF

DO DXCoilNum=1,NumDOE2DXCoils+NumDXMulModeCoils
  ! Setup Report Variables for Cooling Equipment
  ! CurrentModuleObject='Coil:Cooling:DX:SingleSpeed/Coil:Cooling:DX:TwoStageWithHumidityControlMode'
  CALL SetupOutputVariable('DX Coil Total Cooling Rate[W]',DXCoil(DXCoilNum)%TotalCoolingEnergyRate,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Total Cooling Energy[J]',DXCoil(DXCoilNum)%TotalCoolingEnergy,'System','Sum',&
                           DXCoil(DXCoilNum)%Name, &
                           ResourceTypeKey='ENERGYTRANSFER',EndUseKey='COOLINGCOILS',GroupKey='System')
  CALL SetupOutputVariable('DX Coil Sensible Cooling Rate[W]',DXCoil(DXCoilNum)%SensCoolingEnergyRate,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Sensible Cooling Energy[J]',DXCoil(DXCoilNum)%SensCoolingEnergy,'System','Sum',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Latent Cooling Rate[W]',DXCoil(DXCoilNum)%LatCoolingEnergyRate,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Latent Cooling Energy[J]',DXCoil(DXCoilNum)%LatCoolingEnergy,'System','Sum',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Cooling Coil Electric Power[W]',DXCoil(DXCoilNum)%ElecCoolingPower,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Cooling Coil Electric Consumption[J]',DXCoil(DXCoilNum)%ElecCoolingConsumption,'System','Sum',&
                           DXCoil(DXCoilNum)%Name, &
                           ResourceTypeKey='Electric',EndUseKey='COOLING',GroupKey='System')
  CALL SetupOutputVariable('DX Cooling Coil Runtime Fraction []',DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction,'System','Average',&
                           DXCoil(DXCoilNum)%Name)

! do we report these even if no storage tank?
  IF (DXCoil(DXCoilNum)%CondensateCollectMode == CondensateToTank) THEN
    CALL SetupOutputVariable('DX Cooling Coil Condensate Volumetric Flow Rate [m3/s]',DXCoil(DXCoilNum)%CondensateVdot,&
                           'System','Average', DXCoil(DXCoilNum)%Name)
    CALL SetupOutputVariable('DX Cooling Coil Condensate Volume [m3]',DXCoil(DXCoilNum)%CondensateVol,&
                           'System','Sum', DXCoil(DXCoilNum)%Name,  &
                           ResourceTypeKey='OnSiteWater', &
                           EndUseKey='Condensate', GroupKey='System')
  ENDIF

! Moved to Init
!  IF (DXCoil(DXCoilNum)%ReportCoolingCoilCrankcasePower) THEN
!    CALL SetupOutputVariable('DX Cooling Coil Crankcase Heater Power[W]',DXCoil(DXCoilNum)%CrankcaseHeaterPower,'System',&
!                             'Average',DXCoil(DXCoilNum)%Name)
!    CALL SetupOutputVariable('DX Cooling Coil Crankcase Heater Consumption[J]',DXCoil(DXCoilNum)%CrankcaseHeaterConsumption,&
!                             'System','Sum',DXCoil(DXCoilNum)%Name, &
!                              ResourceTypeKey='Electric',EndUseKey='COOLING',GroupKey='System')
!  END IF

  IF (DXCoil(DXCoilNum)%ReportEvapCondVars) THEN
    CALL SetupOutputVariable('DX Cooling Coil Condenser Inlet Temp [C]', &
                              DXCoil(DXCoilNum)%CondInletTemp,'System','Average', DXCoil(DXCoilNum)%Name)
    CALL SetupOutputVariable('DX Cooling Coil Evap Condenser Water Consumption [m3]',DXCoil(DXCoilNum)%EvapWaterConsump, &
                             'System','Sum',DXCoil(DXCoilNum)%Name, &
                              ResourceTypeKey='Water',EndUseKey='Cooling',GroupKey='System')
    CALL SetupOutputVariable('Mains Water Supply for DX Cooling Coil Evap Condenser [m3]',DXCoil(DXCoilNum)%EvapWaterConsump, &
                             'System','Sum',DXCoil(DXCoilNum)%Name, &
                              ResourceTypeKey='MainsWater',EndUseKey='Cooling',GroupKey='System')
    CALL SetupOutputVariable('DX Cooling Coil Evap Condenser Pump Electric Power[W]',DXCoil(DXCoilNum)%EvapCondPumpElecPower, &
                             'System','Average',DXCoil(DXCoilNum)%Name)
    CALL SetupOutputVariable('DX Cooling Coil Evap Condenser Pump Electric Consumption[J]', &
                              DXCoil(DXCoilNum)%EvapCondPumpElecConsumption,'System','Sum',DXCoil(DXCoilNum)%Name, &
                              ResourceTypeKey='Electric',EndUseKey='COOLING',GroupKey='System')
    IF(DXCoil(DXCoilNum)%BasinHeaterPowerFTempDiff .GT. 0.0d0)THEN
      CALL SetupOutputVariable('DX Cooling Coil Basin Heater Electric Power [W]', &
        DXCoil(DXCoilNum)%BasinHeaterPower,'System','Average',DXCoil(DXCoilNum)%Name)
      CALL SetupOutputVariable('DX Cooling Coil Basin Heater Electric Consumption [J]', &
        DXCoil(DXCoilNum)%BasinHeaterConsumption,'System','Sum',DXCoil(DXCoilNum)%Name, &
        ResourceTypeKey='Electric',EndUseKey='COOLING',GroupKey='System')
    END IF
  END IF

END DO

DO DXCoilNum=NumDOE2DXCoils+1,NumDOE2DXCoils+NumDXMulModeCoils
  ! Setup Report Variables for Cooling Equipment
  ! CurrentModuleObject='Cooling:DX:TwoStageWithHumidityControlMode'
  CALL SetupOutputVariable('DX Cooling Coil Stage 2 Runtime Fraction []',DXCoil(DXCoilNum)%CoolingCoilStg2RuntimeFrac,&
                           'System','Average',DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Cooling Coil Dehumidification Mode []',DXCoil(DXCoilNum)%DehumidificationMode,&
                           'System','Average',DXCoil(DXCoilNum)%Name)
END DO

DO DXCoilNum=NumDOE2DXCoils+NumDXMulModeCoils+1,NumDXHeatingCoils+NumDOE2DXCoils+NumDXMulModeCoils
  ! Setup Report Variables for Heating Equipment
  ! CurrentModuleObject='Coil:Heating:DX:SingleSpeed'
  CALL SetupOutputVariable('DX Coil Total Heating Rate[W]',DXCoil(DXCoilNum)%TotalHeatingEnergyRate,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Total Heating Energy[J]',DXCoil(DXCoilNum)%TotalHeatingEnergy,'System','Sum',&
                           DXCoil(DXCoilNum)%Name, &
                           ResourceTypeKey='ENERGYTRANSFER',EndUseKey='HEATINGCOILS',GroupKey='System')
  CALL SetupOutputVariable('DX Heating Coil Electric Power[W]',DXCoil(DXCoilNum)%ElecHeatingPower,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Heating Coil Electric Consumption[J]',DXCoil(DXCoilNum)%ElecHeatingConsumption,'System','Sum',&
                           DXCoil(DXCoilNum)%Name, &
                           ResourceTypeKey='Electric',EndUseKey='HEATING',GroupKey='System')
  CALL SetupOutputVariable('DX Heating Coil Electric Defrost Power[W]',DXCoil(DXCoilNum)%DefrostPower,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Heating Coil Electric Defrost Consumption[J]',DXCoil(DXCoilNum)%DefrostConsumption,'System',&
                           'Sum',DXCoil(DXCoilNum)%Name, &
                           ResourceTypeKey='Electric',EndUseKey='HEATING',GroupKey='System')
  CALL SetupOutputVariable('DX Heating Coil Crankcase Heater Power[W]',DXCoil(DXCoilNum)%CrankcaseHeaterPower,'System',&
                           'Average',DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Heating Coil Crankcase Heater Consumption[J]',DXCoil(DXCoilNum)%CrankcaseHeaterConsumption,&
                           'System','Sum',DXCoil(DXCoilNum)%Name, &
                           ResourceTypeKey='Electric',EndUseKey='HEATING',GroupKey='System')
  CALL SetupOutputVariable('DX Heating Coil Runtime Fraction []',DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
END DO

DO DXCoilNum=NumDOE2DXCoils+NumDXMulModeCoils+NumDXHeatingCoils+1, &
             NumDOE2DXCoils+NumDXMulModeCoils+NumDXHeatingCoils+NumDXMulSpeedCoils
  ! Setup Report Variables for Cooling Equipment
  ! CurrentModuleObject='Coil:Cooling:DX:TwoSpeed'
  CALL SetupOutputVariable('DX Coil Total Cooling Rate[W]',DXCoil(DXCoilNum)%TotalCoolingEnergyRate,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Total Cooling Energy[J]',DXCoil(DXCoilNum)%TotalCoolingEnergy,'System','Sum',&
                           DXCoil(DXCoilNum)%Name, &
                           ResourceTypeKey='ENERGYTRANSFER',EndUseKey='COOLINGCOILS',GroupKey='System')
  CALL SetupOutputVariable('DX Coil Sensible Cooling Rate[W]',DXCoil(DXCoilNum)%SensCoolingEnergyRate,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Sensible Cooling Energy[J]',DXCoil(DXCoilNum)%SensCoolingEnergy,'System','Sum',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Latent Cooling Rate[W]',DXCoil(DXCoilNum)%LatCoolingEnergyRate,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Latent Cooling Energy[J]',DXCoil(DXCoilNum)%LatCoolingEnergy,'System','Sum',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Cooling Coil Electric Power[W]',DXCoil(DXCoilNum)%ElecCoolingPower,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Cooling Coil Electric Consumption[J]',DXCoil(DXCoilNum)%ElecCoolingConsumption,'System','Sum',&
                           DXCoil(DXCoilNum)%Name, &
                           ResourceTypeKey='Electric',EndUseKey='COOLING',GroupKey='System')
  CALL SetupOutputVariable('DX Cooling Coil Runtime Fraction []',DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction,'System','Average',&
                           DXCoil(DXCoilNum)%Name)

  IF (DXCoil(DXCoilNum)%ReportEvapCondVars) THEN
    CALL SetupOutputVariable('DX Cooling Coil Condenser Inlet Temp [C]', &
                              DXCoil(DXCoilNum)%CondInletTemp,'System','Average', DXCoil(DXCoilNum)%Name)
    CALL SetupOutputVariable('DX Cooling Coil Evap Condenser Water Consumption [m3]',DXCoil(DXCoilNum)%EvapWaterConsump, &
                             'System','Sum',DXCoil(DXCoilNum)%Name, &
                              ResourceTypeKey='Water',EndUseKey='Cooling',GroupKey='System')
    CALL SetupOutputVariable('Mains water for DX Cooling Coil Evap Condenser [m3]',DXCoil(DXCoilNum)%EvapWaterConsump, &
                             'System','Sum',DXCoil(DXCoilNum)%Name, &
                              ResourceTypeKey='MainsWater',EndUseKey='Cooling',GroupKey='System')
    CALL SetupOutputVariable('DX Cooling Coil Evap Condenser Pump Electric Power[W]',DXCoil(DXCoilNum)%EvapCondPumpElecPower, &
                             'System','Average',DXCoil(DXCoilNum)%Name)
    CALL SetupOutputVariable('DX Cooling Coil Evap Condenser Pump Electric Consumption[J]', &
                              DXCoil(DXCoilNum)%EvapCondPumpElecConsumption,'System','Sum',DXCoil(DXCoilNum)%Name, &
                              ResourceTypeKey='Electric',EndUseKey='COOLING',GroupKey='System')
    IF(DXCoil(DXCoilNum)%BasinHeaterPowerFTempDiff .GT. 0.0d0)THEN
      CALL SetupOutputVariable('DX Cooling Coil Basin Heater Electric Power [W]', &
        DXCoil(DXCoilNum)%BasinHeaterPower,'System','Average',DXCoil(DXCoilNum)%Name)
      CALL SetupOutputVariable('DX Cooling Coil Basin Heater Electric Consumption [J]', &
        DXCoil(DXCoilNum)%BasinHeaterConsumption,'System','Sum',DXCoil(DXCoilNum)%Name, &
        ResourceTypeKey='Electric',EndUseKey='COOLING',GroupKey='System')
    END IF

  END IF

END DO

DO DXCoilNum=NumDOE2DXCoils+NumDXMulModeCoils+NumDXHeatingCoils+NumDXMulSpeedCoils+1, &
             NumDOE2DXCoils+NumDXMulModeCoils+NumDXHeatingCoils+NumDXMulSpeedCoils+NumDXHeatPumpWaterHeaterCoils
  ! Setup Report Variables for Cooling Equipment
  ! CurrentModuleObject='Coil:WaterHeating:AirToWaterHeatPump'
  CALL SetupOutputVariable('DX Coil Total Cooling Rate[W]',DXCoil(DXCoilNum)%TotalCoolingEnergyRate,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Total Cooling Energy[J]',DXCoil(DXCoilNum)%TotalCoolingEnergy,'System','Sum',&
                           DXCoil(DXCoilNum)%Name) !, &
!                           ResourceTypeKey='ENERGYTRANSFER',EndUseKey='COOLING',GroupKey='Plant')
  CALL SetupOutputVariable('DX Coil Sensible Cooling Rate[W]',DXCoil(DXCoilNum)%SensCoolingEnergyRate,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Sensible Cooling Energy[J]',DXCoil(DXCoilNum)%SensCoolingEnergy,'System','Sum',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Latent Cooling Rate[W]',DXCoil(DXCoilNum)%LatCoolingEnergyRate,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Latent Cooling Energy[J]',DXCoil(DXCoilNum)%LatCoolingEnergy,'System','Sum',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Runtime Fraction []',DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction,'System','Average',&
                           DXCoil(DXCoilNum)%Name)

  IF (DXCoil(DXCoilNum)%ReportCoolingCoilCrankcasePower) THEN
    CALL SetupOutputVariable('DX Coil Crankcase Heater Power[W]',DXCoil(DXCoilNum)%CrankcaseHeaterPower,'System',&
                             'Average',DXCoil(DXCoilNum)%Name)
    CALL SetupOutputVariable('DX Coil Crankcase Heater Consumption[J]',DXCoil(DXCoilNum)%CrankcaseHeaterConsumption,&
                             'System','Sum',DXCoil(DXCoilNum)%Name, &
                              ResourceTypeKey='Electric',EndUseKey='DHW',GroupKey='Plant')
  END IF

! new report variables for a HP water heater DX coil
  CALL SetupOutputVariable('DX Coil Total Water Heating Rate[W]',DXCoil(DXCoilNum)%TotalHeatingEnergyRate,'System', &
                           'Average',DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Total Water Heating Energy[J]',DXCoil(DXCoilNum)%TotalHeatingEnergy,'System', &
                           'Sum',DXCoil(DXCoilNum)%Name) !, &
!                           ResourceTypeKey='ENERGYTRANSFER',EndUseKey='HEATING',GroupKey='Plant')
  CALL SetupOutputVariable('DX Coil Water Heating Electric Power[W]',DXCoil(DXCoilNum)%ElecWaterHeatingPower,'System', &
                           'Average',DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Water Heating Electric Consumption[J]',DXCoil(DXCoilNum)%ElecWaterHeatingConsumption,'System', &
                           'Sum',DXCoil(DXCoilNum)%Name, &
                           ResourceTypeKey='Electric',EndUseKey='DHW',GroupKey='Plant')
END DO

DO DXCoilNum=NumDoe2DXCoils + NumDXHeatingCoils + NumDXMulSpeedCoils + NumDXMulModeCoils + NumDXHeatPumpWaterHeaterCoils+1, &
   NumDoe2DXCoils + NumDXHeatingCoils + NumDXMulSpeedCoils + NumDXMulModeCoils + NumDXHeatPumpWaterHeaterCoils &
             + NumDXMulSpeedCoolCoils
  ! Setup Report Variables for Cooling Equipment:
  ! CurrentModuleObject='Coil:Cooling:DX:MultiSpeed'
  CALL SetupOutputVariable('DX Coil Total Cooling Rate[W]',DXCoil(DXCoilNum)%TotalCoolingEnergyRate,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Total Cooling Energy[J]',DXCoil(DXCoilNum)%TotalCoolingEnergy,'System','Sum',&
                           DXCoil(DXCoilNum)%Name, &
                           ResourceTypeKey='ENERGYTRANSFER',EndUseKey='COOLINGCOILS',GroupKey='System')
  CALL SetupOutputVariable('DX Coil Sensible Cooling Rate[W]',DXCoil(DXCoilNum)%SensCoolingEnergyRate,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Sensible Cooling Energy[J]',DXCoil(DXCoilNum)%SensCoolingEnergy,'System','Sum',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Latent Cooling Rate[W]',DXCoil(DXCoilNum)%LatCoolingEnergyRate,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Latent Cooling Energy[J]',DXCoil(DXCoilNum)%LatCoolingEnergy,'System','Sum',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Cooling Coil Electric Power[W]',DXCoil(DXCoilNum)%ElecCoolingPower,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Cooling Coil Electric Consumption[J]',DXCoil(DXCoilNum)%ElecCoolingConsumption,'System','Sum',&
                           DXCoil(DXCoilNum)%Name, &
                           ResourceTypeKey='Electric',EndUseKey='COOLING',GroupKey='System')

  IF (DXCoil(DXCoilNum)%FuelType .NE. FuelTypeElectricity) THEN
    CALL SetupOutputVariable('DX Cooling Coil '//TRIM(cValidFuelTypes(DXCoil(DXCoilNum)%FuelType))//' Power[W]', &
       DXCoil(DXCoilNum)%FuelUsed,'System','Average',DXCoil(DXCoilNum)%Name)
    CALL SetupOutputVariable('DX Cooling Coil '//TRIM(cValidFuelTypes(DXCoil(DXCoilNum)%FuelType))//' Consumption[J]', &
       DXCoil(DXCoilNum)%FuelConsumed,'System','Sum',DXCoil(DXCoilNum)%Name, &
       ResourceTypeKey=TRIM(cValidFuelTypes(DXCoil(DXCoilNum)%FuelType)),EndUseKey='COOLING',GroupKey='System')
  END IF

  CALL SetupOutputVariable('DX Cooling Coil Runtime Fraction []',DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction,'System','Average',&
                           DXCoil(DXCoilNum)%Name)

  IF (DXCoil(DXCoilNum)%ReportEvapCondVars) THEN
    CALL SetupOutputVariable('DX Cooling Coil Condenser Inlet Temp [C]', &
                              DXCoil(DXCoilNum)%CondInletTemp,'System','Average', DXCoil(DXCoilNum)%Name)
    CALL SetupOutputVariable('DX Cooling Coil Evap Condenser Water Consumption [m3]',DXCoil(DXCoilNum)%EvapWaterConsump, &
                             'System','Sum',DXCoil(DXCoilNum)%Name, &
                              ResourceTypeKey='Water',EndUseKey='Cooling',GroupKey='System')
    CALL SetupOutputVariable('Mains water for DX Cooling Coil Evap Condenser [m3]',DXCoil(DXCoilNum)%EvapWaterConsump, &
                             'System','Sum',DXCoil(DXCoilNum)%Name, &
                              ResourceTypeKey='MainsWater',EndUseKey='Cooling',GroupKey='System')
    CALL SetupOutputVariable('DX Cooling Coil Evap Condenser Pump Electric Power[W]',DXCoil(DXCoilNum)%EvapCondPumpElecPower, &
                             'System','Average',DXCoil(DXCoilNum)%Name)
    CALL SetupOutputVariable('DX Cooling Coil Evap Condenser Pump Electric Consumption[J]', &
                              DXCoil(DXCoilNum)%EvapCondPumpElecConsumption,'System','Sum',DXCoil(DXCoilNum)%Name, &
                              ResourceTypeKey='Electric',EndUseKey='COOLING',GroupKey='System')
    IF(DXCoil(DXCoilNum)%BasinHeaterPowerFTempDiff .GT. 0.0d0)THEN
      CALL SetupOutputVariable('DX Cooling Coil Basin Heater Electric Power [W]', &
        DXCoil(DXCoilNum)%BasinHeaterPower,'System','Average',DXCoil(DXCoilNum)%Name)
      CALL SetupOutputVariable('DX Cooling Coil Basin Heater Electric Consumption [J]', &
        DXCoil(DXCoilNum)%BasinHeaterConsumption,'System','Sum',DXCoil(DXCoilNum)%Name, &
        ResourceTypeKey='Electric',EndUseKey='COOLING',GroupKey='System')
    END IF
  END IF

END DO

DO DXCoilNum=NumDoe2DXCoils + NumDXHeatingCoils + NumDXMulSpeedCoils + NumDXMulModeCoils + NumDXHeatPumpWaterHeaterCoils+ &
             NumDXMulSpeedCoolCoils+1, &
   NumDoe2DXCoils + NumDXHeatingCoils + NumDXMulSpeedCoils + NumDXMulModeCoils + NumDXHeatPumpWaterHeaterCoils &
             + NumDXMulSpeedCoolCoils + NumDXMulSpeedHeatCoils
  ! Setup Report Variables for Heating Equipment:
  ! CurrentModuleObject='Coil:Heating:DX:MultiSpeed'
  CALL SetupOutputVariable('DX Coil Total Heating Rate[W]',DXCoil(DXCoilNum)%TotalHeatingEnergyRate,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Total Heating Energy[J]',DXCoil(DXCoilNum)%TotalHeatingEnergy,'System','Sum',&
                           DXCoil(DXCoilNum)%Name, &
                           ResourceTypeKey='ENERGYTRANSFER',EndUseKey='HEATINGCOILS',GroupKey='System')
  CALL SetupOutputVariable('DX Heating Coil Electric Power[W]',DXCoil(DXCoilNum)%ElecHeatingPower,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Heating Coil Electric Consumption[J]',DXCoil(DXCoilNum)%ElecHeatingConsumption,'System','Sum',&
                           DXCoil(DXCoilNum)%Name, &
                           ResourceTypeKey='Electric',EndUseKey='HEATING',GroupKey='System')

  IF (DXCoil(DXCoilNum)%FuelType .NE. FuelTypeElectricity) THEN
    CALL SetupOutputVariable('DX Heating Coil '//TRIM(cValidFuelTypes(DXCoil(DXCoilNum)%FuelType))//' Power[W]', &
       DXCoil(DXCoilNum)%FuelUsed,'System','Average',DXCoil(DXCoilNum)%Name)
    CALL SetupOutputVariable('DX Heating Coil '//TRIM(cValidFuelTypes(DXCoil(DXCoilNum)%FuelType))//' Consumption[J]', &
       DXCoil(DXCoilNum)%FuelConsumed,'System','Sum',DXCoil(DXCoilNum)%Name, &
       ResourceTypeKey=TRIM(cValidFuelTypes(DXCoil(DXCoilNum)%FuelType)),EndUseKey='HEATING',GroupKey='System')
  END IF

  IF (DXCoil(DXCoilNum)%FuelType .NE. FuelTypeElectricity .AND. DXCoil(DXCoilNum)%DefrostStrategy .EQ. ReverseCycle) THEN
    CALL SetupOutputVariable('DX Heating Coil '//TRIM(cValidFuelTypes(DXCoil(DXCoilNum)%FuelType))//' Defrost Power[W]', &
         DXCoil(DXCoilNum)%DefrostPower,'System','Average',DXCoil(DXCoilNum)%Name)
    CALL SetupOutputVariable('DX Heating Coil '//TRIM(cValidFuelTypes(DXCoil(DXCoilNum)%FuelType))//' Defrost Consumption[J]', &
         DXCoil(DXCoilNum)%DefrostConsumption,'System','Sum',DXCoil(DXCoilNum)%Name, &
         ResourceTypeKey=TRIM(cValidFuelTypes(DXCoil(DXCoilNum)%FuelType)),EndUseKey='HEATING',GroupKey='System')
  ELSE
    CALL SetupOutputVariable('DX Heating Coil Electric Defrost Power[W]',DXCoil(DXCoilNum)%DefrostPower,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
    CALL SetupOutputVariable('DX Heating Coil Electric Defrost Consumption[J]',DXCoil(DXCoilNum)%DefrostConsumption,'System',&
                           'Sum',DXCoil(DXCoilNum)%Name, &
                           ResourceTypeKey='Electric',EndUseKey='HEATING',GroupKey='System')
  END IF

  CALL SetupOutputVariable('DX Heating Coil Crankcase Heater Power[W]',DXCoil(DXCoilNum)%CrankcaseHeaterPower,'System',&
                           'Average',DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Heating Coil Crankcase Heater Consumption[J]',DXCoil(DXCoilNum)%CrankcaseHeaterConsumption,&
                           'System','Sum',DXCoil(DXCoilNum)%Name, &
                           ResourceTypeKey='Electric',EndUseKey='HEATING',GroupKey='System')
  CALL SetupOutputVariable('DX Heating Coil Runtime Fraction []',DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction,'System','Average',&
                           DXCoil(DXCoilNum)%Name)

END DO

! VRF cooling coil report variables
DO DXCoilNum=NumDoe2DXCoils + NumDXHeatingCoils + NumDXMulSpeedCoils + NumDXMulModeCoils + NumDXHeatPumpWaterHeaterCoils+ &
             NumDXMulSpeedCoolCoils + NumDXMulSpeedHeatCoils + 1, &
   NumDoe2DXCoils + NumDXHeatingCoils + NumDXMulSpeedCoils + NumDXMulModeCoils + NumDXHeatPumpWaterHeaterCoils &
             + NumDXMulSpeedCoolCoils + NumDXMulSpeedHeatCoils + NumVRFCoolingCoils
  ! Setup Report Variables for Cooling Equipment:
  ! CurrentModuleObject='Coil:Cooling:DX:VariableRefrigerantFlow
  CALL SetupOutputVariable('DX Coil Total Cooling Rate[W]',DXCoil(DXCoilNum)%TotalCoolingEnergyRate,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Total Cooling Energy[J]',DXCoil(DXCoilNum)%TotalCoolingEnergy,'System','Sum',&
                           DXCoil(DXCoilNum)%Name, &
                           ResourceTypeKey='ENERGYTRANSFER',EndUseKey='COOLINGCOILS',GroupKey='System')
  CALL SetupOutputVariable('DX Coil Sensible Cooling Rate[W]',DXCoil(DXCoilNum)%SensCoolingEnergyRate,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Sensible Cooling Energy[J]',DXCoil(DXCoilNum)%SensCoolingEnergy,'System','Sum',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Latent Cooling Rate[W]',DXCoil(DXCoilNum)%LatCoolingEnergyRate,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Latent Cooling Energy[J]',DXCoil(DXCoilNum)%LatCoolingEnergy,'System','Sum',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Cooling Coil Runtime Fraction[]',DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  IF (DXCoil(DXCoilNum)%CondensateCollectMode == CondensateToTank) THEN
    CALL SetupOutputVariable('DX Cooling Coil Condensate Volumetric Flow Rate [m3/s]',DXCoil(DXCoilNum)%CondensateVdot,&
                           'System','Average', DXCoil(DXCoilNum)%Name)
    CALL SetupOutputVariable('DX Cooling Coil Condensate Volume [m3]',DXCoil(DXCoilNum)%CondensateVol,&
                           'System','Sum', DXCoil(DXCoilNum)%Name,  &
                           ResourceTypeKey='OnSiteWater', &
                           EndUseKey='Condensate', GroupKey='System')
  ENDIF
END DO

! VRF heating coil report variables
DO DXCoilNum=NumDoe2DXCoils + NumDXHeatingCoils + NumDXMulSpeedCoils + NumDXMulModeCoils + NumDXHeatPumpWaterHeaterCoils+ &
             NumDXMulSpeedCoolCoils + NumDXMulSpeedHeatCoils + NumVRFCoolingCoils + 1, &
   NumDoe2DXCoils + NumDXHeatingCoils + NumDXMulSpeedCoils + NumDXMulModeCoils + NumDXHeatPumpWaterHeaterCoils &
             + NumDXMulSpeedCoolCoils + NumDXMulSpeedHeatCoils + NumVRFCoolingCoils + NumVRFHeatingCoils
  ! Setup Report Variables for Heating Equipment:
  ! CurrentModuleObject='Coil:Heating:DX:VariableRefrigerantFlow
  CALL SetupOutputVariable('DX Coil Total Heating Rate[W]',DXCoil(DXCoilNum)%TotalHeatingEnergyRate,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Total Heating Energy[J]',DXCoil(DXCoilNum)%TotalHeatingEnergy,'System','Sum',&
                           DXCoil(DXCoilNum)%Name, &
                           ResourceTypeKey='ENERGYTRANSFER',EndUseKey='HEATINGCOILS',GroupKey='System')
  CALL SetupOutputVariable('DX Heating Coil Runtime Fraction',DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
END DO

!RS: Debugging: Adding in HPSim energy reporting (12/15/14)
DO DXCoilNum=NumDoe2DXCoils + NumDXHeatingCoils + NumDXMulSpeedCoils + NumDXMulModeCoils + NumDXHeatPumpWaterHeaterCoils+ &
             NumDXMulSpeedCoolCoils + NumDXMulSpeedHeatCoils + NumVRFCoolingCoils + NumVRFHeatingCoils + 1, &
   NumDoe2DXCoils + NumDXHeatingCoils + NumDXMulSpeedCoils + NumDXMulModeCoils + NumDXHeatPumpWaterHeaterCoils &
             + NumDXMulSpeedCoolCoils + NumDXMulSpeedHeatCoils + NumVRFCoolingCoils + NumVRFHeatingCoils + DXCoilHPSimNum !2
    
  CALL SetupOutputVariable('DX Coil Total Cooling Rate[W]',DXCoil(DXCoilNum)%TotalCoolingEnergyRate,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Total Cooling Energy[J]',DXCoil(DXCoilNum)%TotalCoolingEnergy,'System','Sum',&
                           DXCoil(DXCoilNum)%Name, &
                           ResourceTypeKey='ENERGYTRANSFER',EndUseKey='COOLINGCOILS',GroupKey='System')
  CALL SetupOutputVariable('DX Coil Sensible Cooling Rate[W]',DXCoil(DXCoilNum)%SensCoolingEnergyRate,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Sensible Cooling Energy[J]',DXCoil(DXCoilNum)%SensCoolingEnergy,'System','Sum',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Latent Cooling Rate[W]',DXCoil(DXCoilNum)%LatCoolingEnergyRate,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Coil Latent Cooling Energy[J]',DXCoil(DXCoilNum)%LatCoolingEnergy,'System','Sum',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Cooling Coil Electric Power[W]',DXCoil(DXCoilNum)%ElecCoolingPower,'System','Average',&
                           DXCoil(DXCoilNum)%Name)
  CALL SetupOutputVariable('DX Cooling Coil Electric Consumption[J]',DXCoil(DXCoilNum)%ElecCoolingConsumption,'System','Sum',&
                           DXCoil(DXCoilNum)%Name, &
                           ResourceTypeKey='Electric',EndUseKey='COOLING',GroupKey='System')

END DO
             
IF (AnyEnergyManagementSystemInModel) THEN
  ! setup EMS sizing actuators for single speed DX
  DO DXCoilNum=1,NumDOE2DXCoils
    CALL SetupEMSActuator('Coil:Cooling:DX:SingleSpeed', DXCoil(DXCoilNum)%Name, &
                          'Autosized Rated Air Flow Rate' , '[m3/s]', &
                           DXCoil(DXCoilNum)%RatedAirVolFlowRateEMSOverrideON(1),  &
                           DXCoil(DXCoilNum)%RatedAirVolFlowRateEMSOverrideValue(1) )

    CALL SetupEMSActuator('Coil:Cooling:DX:SingleSpeed', DXCoil(DXCoilNum)%Name, &
                          'Autosized Rated Sensible Heat Ratio' , '[W/W]', &
                           DXCoil(DXCoilNum)%RatedSHREMSOverrideOn(1),  &
                           DXCoil(DXCoilNum)%RatedSHREMSOverrideValue(1) )

    CALL SetupEMSActuator('Coil:Cooling:DX:SingleSpeed', DXCoil(DXCoilNum)%Name, &
                          'Autosized Rated Total Cooling Capacity' , '[W]', &
                           DXCoil(DXCoilNum)%RatedTotCapEMSOverrideOn(1),  &
                           DXCoil(DXCoilNum)%RatedTotCapEMSOverrideValue(1) )
  ENDDO

ENDIF
DEALLOCATE(Alphas)
DEALLOCATE(cAlphaFields)
DEALLOCATE(cNumericFields)
DEALLOCATE(Numbers)
DEALLOCATE(lAlphaBlanks)
DEALLOCATE(lNumericBlanks)

DEALLOCATE(Alphas2)
DEALLOCATE(cAlphaFields2)
DEALLOCATE(cNumericFields2)
DEALLOCATE(Numbers2)
DEALLOCATE(lAlphaBlanks2)
DEALLOCATE(lNumericBlanks2)
CALL ManageEMS(emsCallFromComponentGetInput)


RETURN

END SUBROUTINE GetDXCoils

SUBROUTINE InitDXCoil(DXCoilNum)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Fred Buhl
          !       DATE WRITTEN   May 2000
          !                      Feb 2005 M. J. Witte, GARD Analytics, Inc.
          !                        Add new coil type COIL:DX:MultiMode:CoolingEmpirical:
          !                      July 2005 R. Raustad, FSEC
          !                        Add new coil type COIL:DX:HEATPUMPWATERHEATER
          !                      June 2007 L. Gu, FSEC
          !                        Add new coil type COIL:DX:MULTISPEED:COOLING and HEATING
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! This subroutine is for initializations of DX Coil Components.

          ! METHODOLOGY EMPLOYED:
          ! Uses the status flags to trigger initializations.

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE DataHVACGlobals, ONLY: FanElecPower
  USE DataAirLoop,     ONLY: AirLoopInputsFilled
  USE General,         ONLY: TrimSigDigits
  USE ReportSizingManager, ONLY: ReportSizingOutput

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  INTEGER, INTENT (IN) :: DXCoilNum ! number of the current DX coil unit being simulated

          ! SUBROUTINE PARAMETER DEFINITIONS:
  REAL(r64) :: SmallDifferenceTest=0.00000001d0
  CHARACTER(len=*), PARAMETER :: RoutineName='InitDXCoil'

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  LOGICAL, ALLOCATABLE,SAVE, DIMENSION(:) :: MyEnvrnFlag  ! One time environment flag
  LOGICAL, ALLOCATABLE,SAVE, DIMENSION(:) :: MySizeFlag   ! One time sizing flag
  LOGICAL,SAVE :: MyOneTimeFlag = .TRUE.                  ! One time flag used to allocate MyEnvrnFlag and MySizeFlag
  LOGICAL,SAVE :: CrankcaseHeaterReportVarFlag = .TRUE.   ! One time flag used to report crankcase heater power for non-HP coils
  REAL(r64) :: RatedHeatPumpIndoorAirTemp ! Indoor dry-bulb temperature to heat pump evaporator at rated conditions [C]
  REAL(r64) :: RatedHeatPumpIndoorHumRat  ! Inlet humidity ratio to heat pump evaporator at rated conditions [kgWater/kgDryAir]
  REAL(r64) :: RatedVolFlowPerRatedTotCap ! Rated Air Volume Flow Rate divided by Rated Total Capacity [m3/s-W)
  REAL(r64) :: HPInletAirHumRat           ! Rated inlet air humidity ratio for heat pump water heater [kgWater/kgDryAir]
  LOGICAL :: ErrorsFound=.FALSE.        ! TRUE when errors found
  INTEGER :: CapacityStageNum           ! Loop index for 1,Number of capacity stages
  INTEGER :: DehumidModeNum             ! Loop index for 1,Number of enhanced dehumidification modes
  INTEGER :: Mode                       ! Performance mode for MultiMode DX coil; Always 1 for other coil types
  INTEGER :: DXCoilNumTemp              ! Counter for crankcase heater report variable DO loop
  INTEGER :: AirInletNode               ! Air inlet node number

  IF (MyOneTimeFlag) THEN
    ! initialize the environment and sizing flags
    ALLOCATE(MyEnvrnFlag(NumDXCoils))
    ALLOCATE(MySizeFlag(NumDXCoils))
    MyEnvrnFlag = .TRUE.
    MySizeFlag = .TRUE.
    MyOneTimeFlag = .FALSE.

  END IF

  IF(DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_HeatPumpWaterHeater .AND. MyEnvrnFlag(DXCoilNum)) THEN

    CALL SizeDXCoil(DXCoilNum)

    RatedVolFlowPerRatedTotCap = DXCoil(DXCoilNum)%RatedAirVolFlowRate(1)/DXCoil(DXCoilNum)%RatedTotCap2
    IF (((MinRatedVolFlowPerRatedTotCap - RatedVolFlowPerRatedTotCap) > SmallDifferenceTest).OR. &
        ((RatedVolFlowPerRatedTotCap - MaxHeatVolFlowPerRatedTotCap) > SmallDifferenceTest)) THEN
     CALL ShowSevereError (TRIM(DXCoil(DXCoilNum)%DXCoilType) // ' "'//TRIM(DXCoil(DXCoilNum)%Name)//  &
                           '": Rated air volume flow rate per watt of rated total water '// &
                           'heating capacity is out of range.')
     CALL ShowContinueError('Min Rated Vol Flow Per Watt=['//TRIM(TrimSigDigits(MinRatedVolFlowPerRatedTotCap,3))//'], '//  &
           'Rated Vol Flow Per Watt=['//TRIM(TrimSigDigits(RatedVolFlowPerRatedTotCap,3))//'], Max Rated Vol Flow Per Watt=['// &
           TRIM(TrimSigDigits(MaxHeatVolFlowPerRatedTotCap,3))//']. See Input-Output Reference Manual for valid range.')
    END IF
    HPInletAirHumRat = PsyWFnTdbTwbPb(DXCoil(DXCoilNum)%RatedInletDBTemp,DXCoil(DXCoilNum)%RatedInletWBTemp, &
                       StdBaroPress,RoutineName)
    DXCoil(DXCoilNum)%RatedAirMassFlowRate(1) = DXCoil(DXCoilNum)%RatedAirVolFlowRate(1)* &
                       PsyRhoAirFnPbTdbW(StdBaroPress,DXCoil(DXCoilNum)%RatedInletDBTemp,HPInletAirHumRat,RoutineName)
!   get rated coil bypass factor excluding fan heat
    FanElecPower = 0.0
!   call CalcHPWHDXCoil to determine DXCoil%RatedTotCap(1) for rated CBF calculation below
    CALL CalcHPWHDXCoil(DXCoilNum,1.0d0)
    IF(MySizeFlag(DXCoilNum))THEN
      CALL SizeDXCoil(DXCoilNum)
      MySizeFlag(DXCoilNum) = .FALSE.
    END IF

    DXCoil(DXCoilNum)%RatedCBF(1) = CalcCBF(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name,&
                                    DXCoil(DXCoilNum)%RatedInletDBTemp,HPInletAirHumRat,DXCoil(DXCoilNum)%RatedTotCap(1), &
                                    DXCoil(DXCoilNum)%RatedAirMassFlowRate(1),DXCoil(DXCoilNum)%RatedSHR(1))
    MyEnvrnFlag(DXCoilNum) = .FALSE.
  END IF

! Find the companion upstream coil (DX cooling coil) that is used with DX heating coils (HP AC units only)
  IF(DXCoil(DXCoilNum)%FindCompanionUpStreamCoil)THEN
    IF(DXCoil(DXCoilNum)%DXCoilType_Num .EQ. CoilDX_HeatingEmpirical .OR. &
       DXCoil(DXCoilNum)%DXCoilType_Num .EQ. CoilDX_MultiSpeedHeating)THEN
      DXCoil(DXCoilNum)%CompanionUpstreamDXCoil = &
                  GetHPCoolingCoilIndex(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, DXCoilNum)
      IF(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil .GT. 0)THEN
        DXCoil(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil)%ReportCoolingCoilCrankcasePower = .FALSE.
        DXCoil(DXCoilNum)%FindCompanionUpStreamCoil = .FALSE.
!       Copy condenser node number from DX cooling coil when used with a companion DX heating coil
        DO Mode = 1, MaxModes
          DXCoil(DXCoilNum)%CondenserInletNodeNum(Mode) = &
                          DXCoil(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil)%CondenserInletNodeNum(Mode)
        END DO
      END IF
    ELSE
      DXCoil(DXCoilNum)%FindCompanionUpStreamCoil = .FALSE.
    END IF
  END IF !IF(DXCoil(DXCoilNum)%FindCompanionUpStreamCoil)THEN

! CR7308 - Wait for zone and air loop equipment to be simulated, then print out report variables
  IF(CrankcaseHeaterReportVarFlag)THEN
    IF(AirLoopInputsFilled)THEN
!     Set report variables for DX cooling coils that will have a crankcase heater (all DX coils not used in a HP AC unit)
      DO DXCoilNumTemp=1,NumDOE2DXCoils+NumDXMulModeCoils
        IF (DXCoil(DXCoilNumTemp)%ReportCoolingCoilCrankcasePower) THEN
          CALL SetupOutputVariable('DX Cooling Coil Crankcase Heater Power[W]', &
                                    DXCoil(DXCoilNumTemp)%CrankcaseHeaterPower,'System', 'Average',DXCoil(DXCoilNumTemp)%Name)
          CALL SetupOutputVariable('DX Cooling Coil Crankcase Heater Consumption[J]', &
                                    DXCoil(DXCoilNumTemp)%CrankcaseHeaterConsumption, 'System','Sum',DXCoil(DXCoilNumTemp)%Name, &
                                    ResourceTypeKey='Electric',EndUseKey='COOLING',GroupKey='Plant')
          DXCoil(DXCoilNumTemp)%ReportCoolingCoilCrankcasePower = .FALSE.
        END IF
      END DO
      CrankcaseHeaterReportVarFlag = .FALSE.
    END IF !(AirLoopInputsFilled)THEN
  END IF !(CrankcaseHeaterReportVarFlag)THEN

! Find the companion upstream coil (DX cooling coil) that is used with DX heating coils (Multispeed HP units only)
  IF(DXCoil(DXCoilNum)%FindCompanionUpStreamCoil)THEN
    IF(DXCoil(DXCoilNum)%DXCoilType_Num .EQ. CoilDX_MultiSpeedHeating)THEN
      DXCoil(DXCoilNum)%CompanionUpstreamDXCoil = &
                  GetHPCoolingCoilIndex(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, DXCoilNum)
      IF(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil .GT. 0)THEN
        DXCoil(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil)%ReportCoolingCoilCrankcasePower = .FALSE.
        DXCoil(DXCoilNum)%FindCompanionUpStreamCoil = .FALSE.
      END IF
    ELSE
      DXCoil(DXCoilNum)%FindCompanionUpStreamCoil = .FALSE.
    END IF
  END IF !IF(DXCoil(DXCoilNum)%FindCompanionUpStreamCoil)THEN

  IF ( .NOT. SysSizingCalc .AND. MySizeFlag(DXCoilNum)) THEN
    ! for each coil, do the sizing once.
    CALL SizeDXCoil(DXCoilNum)
    MySizeFlag(DXCoilNum) = .FALSE.

    IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingSingleSpeed .OR. &
        DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoSpeed .OR. &
        DXCoil(DXCoilNum)%DXCoilType_Num == CoilVRF_Cooling) THEN

      Mode = 1
      ! Check for zero capacity or zero max flow rate
      IF (DXCoil(DXCoilNum)%RatedTotCap(Mode) <= 0.0) THEN
        CALL ShowSevereError('Sizing: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name)// &
                            ' has zero rated total capacity')
        ErrorsFound=.TRUE.
      END IF
      IF (DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) <= 0.0) THEN
        CALL ShowSevereError('Sizing: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name)// &
                            ' has zero rated air flow rate')
        ErrorsFound=.TRUE.
      END IF
      IF (ErrorsFound) THEN
        CALL ShowFatalError('Preceding condition causes termination.')
      ENDIF
      !
      ! Check for valid range of (Rated Air Volume Flow Rate / Rated Total Capacity)
      !
      RatedVolFlowPerRatedTotCap = DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode)/DXCoil(DXCoilNum)%RatedTotCap(Mode)
      IF (((MinRatedVolFlowPerRatedTotCap - RatedVolFlowPerRatedTotCap) > SmallDifferenceTest).OR. &
          ((RatedVolFlowPerRatedTotCap - MaxRatedVolFlowPerRatedTotCap) > SmallDifferenceTest)) THEN
       CALL ShowSevereError ('Sizing: '//TRIM(DXCoil(DXCoilNum)%DXCoilType) // ' "'//TRIM(DXCoil(DXCoilNum)%Name)//  &
                             '": Rated air volume flow rate per watt of rated total '// &
                             'cooling capacity is out of range.')
       CALL ShowContinueError('Min Rated Vol Flow Per Watt=['//TRIM(TrimSigDigits(MinRatedVolFlowPerRatedTotCap,3))//'], '//  &
           'Rated Vol Flow Per Watt=['//TRIM(TrimSigDigits(RatedVolFlowPerRatedTotCap,3))//'], Max Rated Vol Flow Per Watt=['// &
           TRIM(TrimSigDigits(MaxRatedVolFlowPerRatedTotCap,3))//']. See Input-Output Reference Manual for valid range.')
      END IF
      DXCoil(DXCoilNum)%RatedAirMassFlowRate(Mode) = DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode)* &
        PsyRhoAirFnPbTdbW(StdBaroPress,RatedInletAirTemp,RatedInletAirHumRat,RoutineName)
      ! get high speed rated coil bypass factor
      DXCoil(DXCoilNum)%RatedCBF(Mode) = CalcCBF(DXCoil(DXCoilNum)%DXCoilType,DXCoil(DXCoilNum)%Name,&
                                           RatedInletAirTemp,RatedInletAirHumRat,DXCoil(DXCoilNum)%RatedTotCap(Mode),&
                                           DXCoil(DXCoilNum)%RatedAirMassFlowRate(Mode),DXCoil(DXCoilNum)%RatedSHR(Mode))
    END IF

    IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoStageWHumControl) THEN
      DO DehumidModeNum = 0, DXCoil(DXCoilNum)%NumDehumidModes
        DO CapacityStageNum = 1, DXCoil(DXCoilNum)%NumCapacityStages
          Mode = DehumidModeNum*2 + CapacityStageNum
          ! Check for zero capacity or zero max flow rate
          IF (DXCoil(DXCoilNum)%RatedTotCap(Mode) <= 0.0) THEN
            CALL ShowSevereError('Sizing: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name)// &
                                ' has zero rated total capacity')
            CALL ShowContinueError('for CoilPerformance:DX:Cooling mode: '// &
                                   TRIM(DXCoil(DXCoilNum)%CoilPerformanceName(Mode)))
            ErrorsFound=.TRUE.
          END IF
          IF (DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) <= 0.0) THEN
            CALL ShowSevereError('Sizing: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name)// &
                                ' has zero rated air flow rate')
            CALL ShowContinueError('for CoilPerformance:DX:Cooling mode: '// &
                                   TRIM(DXCoil(DXCoilNum)%CoilPerformanceName(Mode)))
            ErrorsFound=.TRUE.
          END IF
          IF (ErrorsFound) THEN
            CALL ShowFatalError('Preceding condition causes termination.')
          ENDIF
          !
          ! Check for valid range of (Rated Air Volume Flow Rate / Rated Total Capacity)
          !
          RatedVolFlowPerRatedTotCap = DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode)/DXCoil(DXCoilNum)%RatedTotCap(Mode)
          IF (((MinRatedVolFlowPerRatedTotCap - RatedVolFlowPerRatedTotCap) > SmallDifferenceTest).OR. &
              ((RatedVolFlowPerRatedTotCap - MaxRatedVolFlowPerRatedTotCap) > SmallDifferenceTest)) THEN
            CALL ShowSevereError ('Sizing: '//TRIM(DXCoil(DXCoilNum)%DXCoilType) // ' "'//TRIM(DXCoil(DXCoilNum)%Name)//  &
                                 '": Rated air volume flow rate per watt of rated total '// &
                                 'cooling capacity is out of range.')
            CALL ShowContinueError('Min Rated Vol Flow Per Watt=['//TRIM(TrimSigDigits(MinRatedVolFlowPerRatedTotCap,3))//'], '// &
             'Rated Vol Flow Per Watt=['//TRIM(TrimSigDigits(RatedVolFlowPerRatedTotCap,3))//'], Max Rated Vol Flow Per Watt=['// &
             TRIM(TrimSigDigits(MaxRatedVolFlowPerRatedTotCap,3))//']. See Input-Output Reference Manual for valid range.')
            CALL ShowContinueError('for CoilPerformance:DX:Cooling mode: '// &
                                   TRIM(DXCoil(DXCoilNum)%CoilPerformanceName(Mode)))
          END IF
          DXCoil(DXCoilNum)%RatedAirMassFlowRate(Mode) = DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode)* &
            PsyRhoAirFnPbTdbW(StdBaroPress,RatedInletAirTemp,RatedInletAirHumRat,RoutineName)
          ! get rated coil bypass factor
          DXCoil(DXCoilNum)%RatedCBF(Mode) = CalcCBF(DXCoil(DXCoilNum)%CoilPerformanceType(Mode), &
                                               DXCoil(DXCoilNum)%CoilPerformanceName(Mode),&
                                               RatedInletAirTemp,RatedInletAirHumRat,DXCoil(DXCoilNum)%RatedTotCap(Mode),&
                                               DXCoil(DXCoilNum)%RatedAirMassFlowRate(Mode),DXCoil(DXCoilNum)%RatedSHR(Mode))
        END DO ! End capacity stages loop
      END DO ! End dehumidification modes loop

    END IF


    IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_HeatingEmpirical .OR. &
        DXCoil(DXCoilNum)%DXCoilType_Num == CoilVRF_Heating) THEN

      Mode = 1
      IF (DXCoil(DXCoilNum)%RatedTotCap(Mode) <= 0.0) THEN
        CALL ShowSevereError('Sizing: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name)// &
                            ' has zero rated total capacity')
        ErrorsFound=.TRUE.
      END IF
      IF (DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) <= 0.0) THEN
        CALL ShowSevereError('Sizing: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name)// &
                            ' has zero rated air flow rate')
        ErrorsFound=.TRUE.
      END IF
      IF (ErrorsFound) THEN
        CALL ShowFatalError('Preceding condition causes termination.')
      ENDIF
      RatedHeatPumpIndoorAirTemp = 21.11d0  ! 21.11C or 70F
      RatedHeatPumpIndoorHumRat = 0.00881d0 ! Humidity ratio corresponding to 70F dry bulb/60F wet bulb
      DXCoil(DXCoilNum)%RatedAirMassFlowRate(Mode) = DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode)* &
         PsyRhoAirFnPbTdbW(StdBaroPress,RatedHeatPumpIndoorAirTemp,RatedHeatPumpIndoorHumRat,RoutineName)
      ! Check for valid range of (Rated Air Volume Flow Rate / Rated Total Capacity)
      !
      RatedVolFlowPerRatedTotCap = DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode)/DXCoil(DXCoilNum)%RatedTotCap(Mode)
      IF (((MinRatedVolFlowPerRatedTotCap - RatedVolFlowPerRatedTotCap) > SmallDifferenceTest).OR. &
          ((RatedVolFlowperRatedTotCap - MaxRatedVolFlowPerRatedTotCap) > SmallDifferenceTest)) THEN
        CALL ShowSevereError ('Sizing: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name)//  &
                              ': Rated air volume flow rate per watt of rated total '// &
                              'heating capacity is out of range.')
        CALL ShowContinueError('Min Rated Vol Flow Per Watt=['//TRIM(TrimSigDigits(MinRatedVolFlowPerRatedTotCap,3))//'], '//  &
           'Rated Vol Flow Per Watt=['//TRIM(TrimSigDigits(RatedVolFlowPerRatedTotCap,3))//'], Max Rated Vol Flow Per Watt=['// &
           TRIM(TrimSigDigits(MaxRatedVolFlowPerRatedTotCap,3))//']. See Input-Output Reference Manual for valid range.')
      END IF

    END IF


    IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoSpeed) THEN
      !
      ! Check for valid range of (Rated Air Volume Flow Rate / Rated Total Capacity)
      RatedVolFlowPerRatedTotCap = DXCoil(DXCoilNum)%RatedAirVolFlowRate2/DXCoil(DXCoilNum)%RatedTotCap2
      IF (((MinRatedVolFlowPerRatedTotCap - RatedVolFlowPerRatedTotCap) > SmallDifferenceTest).OR. &
          ((RatedVolFlowPerRatedTotCap - MaxRatedVolFlowPerRatedTotCap) > SmallDifferenceTest)) THEN
       CALL ShowSevereError ('Coil:Cooling:DX:TwoSpeed "'//TRIM(DXCoil(DXCoilNum)%Name)//  &
                             '": At low speed rated air volume flow rate per watt of rated total '// &
                             'cooling capacity is out of range.')
       CALL ShowContinueError('Min Rated Vol Flow Per Watt=['//TRIM(TrimSigDigits(MinRatedVolFlowPerRatedTotCap,3))//'], '// &
           'Rated Vol Flow Per Watt=['//TRIM(TrimSigDigits(RatedVolFlowPerRatedTotCap,3))//'], Max Rated Vol Flow Per Watt=['// &
           TRIM(TrimSigDigits(MaxRatedVolFlowPerRatedTotCap,3))//']. See Input-Output Reference Manual for valid range.')
      END IF

      DXCoil(DXCoilNum)%RatedAirMassFlowRate2 = DXCoil(DXCoilNum)%RatedAirVolFlowRate2* &
        PsyRhoAirFnPbTdbW(StdBaroPress,RatedInletAirTemp,RatedInletAirHumRat,RoutineName)
      ! get low speed rated coil bypass factor
      DXCoil(DXCoilNum)%RatedCBF2 = CalcCBF(DXCoil(DXCoilNum)%DXCoilType,DXCoil(DXCoilNum)%Name,&
                                           RatedInletAirTemp,RatedInletAirHumRat,DXCoil(DXCoilNum)%RatedTotCap2,&
                                           DXCoil(DXCoilNum)%RatedAirMassFlowRate2,DXCoil(DXCoilNum)%RatedSHR2)
    END IF

!   Autosizing is completed in Size routine, however, the HPWH disrupts the flow of the eio and reporting
!   is done here while all other coils are sized and reported.
    IF(DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_HeatPumpWaterHeater .AND. DXCoil(DXCoilNum)%AirVolFlowAutoSized)THEN
      CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                             'Rated Air Volume Flow Rate [m3/s]', DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode))
      DXCoil(DXCoilNum)%AirVolFlowAutoSized = .FALSE.
    END IF
    IF(DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_HeatPumpWaterHeater .AND. DXCoil(DXCoilNum)%WaterVolFlowAutoSized)THEN
       CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                              'Rated Condenser Water Volume Flow Rate [m3/s]', DXCoil(DXCoilNum)%RatedHPWHCondWaterFlow)
       DXCoil(DXCoilNum)%WaterVolFlowAutoSized = .FALSE.
    END IF

    ! Multispeed Cooling
    IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_MultiSpeedCooling) THEN
      Do Mode = 1, DXCoil(DXCoilNum)%NumOfSpeeds
        ! Check for zero capacity or zero max flow rate
        IF (DXCoil(DXCoilNum)%MSRatedTotCap(Mode) <= 0.0) THEN
          CALL ShowSevereError('Sizing: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name)// &
                            ' has zero rated total capacity at speed '//Trim(TrimSigDigits(Mode)))
          ErrorsFound=.TRUE.
        END IF
        IF (DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode) <= 0.0) THEN
          CALL ShowSevereError('Sizing: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name)// &
                            ' has zero rated air flow rate at speed '//Trim(TrimSigDigits(Mode)))
          ErrorsFound=.TRUE.
        END IF
        IF (ErrorsFound) THEN
          CALL ShowFatalError('Preceding condition causes termination.')
        ENDIF
        !
        ! Check for valid range of (Rated Air Volume Flow Rate / Rated Total Capacity)
        !
        RatedVolFlowPerRatedTotCap = DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode)/  &
                                                 DXCoil(DXCoilNum)%MSRatedTotCap(Mode)
        IF (((MinRatedVolFlowPerRatedTotCap - RatedVolFlowPerRatedTotCap) > SmallDifferenceTest).OR. &
           ((RatedVolFlowPerRatedTotCap - MaxRatedVolFlowPerRatedTotCap) > SmallDifferenceTest)) THEN
          CALL ShowSevereError ('Sizing: '//TRIM(DXCoil(DXCoilNum)%DXCoilType) // ' "'//TRIM(DXCoil(DXCoilNum)%Name)//  &
                '": Rated air volume flow rate per watt of rated total '// &
                'cooling capacity is out of range at speed '//TRIM(TrimSigDigits(Mode)))
          CALL ShowContinueError('Min Rated Vol Flow Per Watt=['//TRIM(TrimSigDigits(MinRatedVolFlowPerRatedTotCap,3))//'], '// &
           'Rated Vol Flow Per Watt=['//TRIM(TrimSigDigits(RatedVolFlowPerRatedTotCap,3))//'], Max Rated Vol Flow Per Watt=['// &
           TRIM(TrimSigDigits(MaxRatedVolFlowPerRatedTotCap,3))//']. See Input-Output Reference Manual for valid range.')
        END IF
        DXCoil(DXCoilNum)%MSRatedAirMassFlowRate(Mode) = DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode)* &
          PsyRhoAirFnPbTdbW(StdBaroPress,RatedInletAirTemp,RatedInletAirHumRat,RoutineName)
        ! get high speed rated coil bypass factor
        DXCoil(DXCoilNum)%MSRatedCBF(Mode) = CalcCBF(DXCoil(DXCoilNum)%DXCoilType,DXCoil(DXCoilNum)%Name,&
                                           RatedInletAirTemp,RatedInletAirHumRat,DXCoil(DXCoilNum)%MSRatedTotCap(Mode),&
                                           DXCoil(DXCoilNum)%MSRatedAirMassFlowRate(Mode),DXCoil(DXCoilNum)%MSRatedSHR(Mode))
      END DO
    END IF

    ! Multispeed Heating
    IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_MultiSpeedHeating) THEN
      RatedHeatPumpIndoorAirTemp = 21.11d0  ! 21.11C or 70F
      RatedHeatPumpIndoorHumRat = 0.00881d0 ! Humidity ratio corresponding to 70F dry bulb/60F wet bulb
      Do Mode = 1, DXCoil(DXCoilNum)%NumOfSpeeds

        DXCoil(DXCoilNum)%MSRatedAirMassFlowRate(Mode) = DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode)* &
         PsyRhoAirFnPbTdbW(StdBaroPress,RatedHeatPumpIndoorAirTemp,RatedHeatPumpIndoorHumRat,RoutineName)
        ! Check for valid range of (Rated Air Volume Flow Rate / Rated Total Capacity)
        !
        RatedVolFlowPerRatedTotCap = DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode)/  &
                                                  DXCoil(DXCoilNum)%MSRatedTotCap(Mode)
        IF (((MinRatedVolFlowPerRatedTotCap - RatedVolFlowPerRatedTotCap) > SmallDifferenceTest).OR. &
            ((RatedVolFlowperRatedTotCap - MaxRatedVolFlowPerRatedTotCap) > SmallDifferenceTest)) THEN
          CALL ShowSevereError ('Coil:Heating:DX:MultiSpeed '//TRIM(DXCoil(DXCoilNum)%Name)//  &
                              ': Rated air volume flow rate per watt of rated total '// &
                'heating capacity is out of range at speed '//TRIM(TrimSigDigits(Mode)))
          CALL ShowContinueError('Min Rated Vol Flow Per Watt=['//TRIM(TrimSigDigits(MinRatedVolFlowPerRatedTotCap,3))//'], '// &
           'Rated Vol Flow Per Watt=['//TRIM(TrimSigDigits(RatedVolFlowPerRatedTotCap,3))//'], Max Rated Vol Flow Per Watt=['// &
           TRIM(TrimSigDigits(MaxRatedVolFlowPerRatedTotCap,3))//']. See Input-Output Reference Manual for valid range.')
        END IF
      End Do
    END IF

  END IF

  AirInletNode = DXCoil(DXCoilNum)%AirInNode

  ! Each iteration, load the coil data structure with the inlet conditions

  DXCoil(DXCoilNum)%InletAirMassFlowRate    = Node(AirInletNode)%MassFlowRate
  DXCoil(DXCoilNum)%InletAirMassFlowRateMax = MAX(Node(AirInletNode)%MassFlowRateMax,Node(AirInletNode)%MassFlowRate)
  DXCoil(DXCoilNum)%InletAirTemp            = Node(AirInletNode)%Temp
  DXCoil(DXCoilNum)%InletAirHumRat          = Node(AirInletNode)%HumRat
  DXCoil(DXCoilNum)%InletAirEnthalpy        = Node(AirInletNode)%Enthalpy
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!  DXCoil(DXCoilNum)%InletAirPressure        = Node(AirInletNode)%Press

  IF(DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_HeatPumpWaterHeater) THEN
    DXCoil(DXCoilNum)%TotalHeatingEnergyRate = 0.0
    DXCoil(DXCoilNum)%ElecWaterHeatingPower  = 0.0
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!  DXCoil(DXCoilNum)%InletAirPressure         = StdBaroPress

!   HPWH's that use an inlet air temperature schedule also need to have a valid barometric pressure
!   The DX Coil used in HPWH's does not know if it is using a scheduled inlet temperature so check the node pressure
    IF (DXCoil(DXCoilNum)%CondenserInletNodeNum(1) > 0) THEN
      IF(Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(1))%Press == 0.0)THEN
        Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(1))%Press = StdBaroPress
      END IF
    END IF

  END IF
  DXCoil(DXCoilNum)%BasinHeaterPower = 0.0d0

  RETURN

END SUBROUTINE InitDXCoil

SUBROUTINE SizeDXCoil(DXCoilNum)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Fred Buhl
          !       DATE WRITTEN   January 2002
          !                      Feb 2005 M. J. Witte, GARD Analytics, Inc.
          !                        Add new coil type COIL:DX:MultiMode:CoolingEmpirical:
          !                      July 2005 R. Raustad, FSEC
          !                        Add new coil type COIL:DX:HEATPUMPWATERHEATER
          !                      June 2007 L. Gu, FSEC
          !                        Add new coil type COIL:DX:MULTISPEED:COOLING and HEATING
          !                      January 2011, B. Griffithn NREL. add EMS overrides for autosized fields. 1`
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! This subroutine is for sizing DX Coil components for which nominal capacity and air flow rate
          ! have not been specified in the input.

          ! METHODOLOGY EMPLOYED:
          ! Obtains cooling capacities and air flow rates from the zone or system sizing arrays.

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE DataSizing
  USE CurveManager, ONLY: CurveValue
  USE General,      ONLY: RoundSigDigits, TrimSigDigits
  USE ReportSizingManager, ONLY: ReportSizingOutput
  USE OutputReportPredefined
  USE DataAirSystems, ONLY: PrimaryAirSystem
  USE StandardRatings,   ONLY: CalcDXCoilStandardRating

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  Integer, Intent(IN) :: DXCoilNum

          ! SUBROUTINE PARAMETER DEFINITIONS:
  CHARACTER(len=*), PARAMETER ::  RoutineName='SizeDXCoil'

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  REAL(r64) :: rhoair
  REAL(r64) :: CpAir
  REAL(r64) :: MixTemp
  REAL(r64) :: MixHumRat
  REAL(r64) :: MixEnth
  REAL(r64) :: MixWetBulb
  REAL(r64) :: SupTemp
  REAL(r64) :: SupHumRat
  REAL(r64) :: SupEnth
  REAL(r64) :: OutTemp
  REAL(r64) :: OutAirFrac
  REAL(r64) :: VolFlowRate
  REAL(r64) :: CoolCapAtPeak
  REAL(r64) :: TotCapTempModFac
  REAL(r64) :: RatedVolFlowPerRatedTotCap !Rated Air Volume Flow Rate divided by Rated Total Capacity[m3/s-W)
  INTEGER :: TimeStepNumAtMax
  INTEGER :: DDNum
  INTEGER :: CapacityStageNum     ! Loop index for 1,Number of capacity stages
  INTEGER :: DehumidModeNum       ! Loop index for 1,Number of enhanced dehumidification modes
  INTEGER :: Mode  ! Operating mode for MultiMode DX coil; Always 1 for other coil types
  INTEGER :: NumOfSpeedCompanion  ! Number of speed for a companion cooling coil (Multispeed HO heating coil only
  CHARACTER(len=MaxNameLength) :: equipName
!  EXTERNAL ReportSizingOutput

  ! NOTE: we are sizing COIL:DX:HeatingEmpirical on the COOLING load. Thus the cooling and
  ! and heating capacities of a DX heat pump system will be identical. In real life the ARI
  ! heating and cooling capacities are close but not identical.

  DO DehumidModeNum = 0, DXCoil(DXCoilNum)%NumDehumidModes
    DO CapacityStageNum = 1, DXCoil(DXCoilNum)%NumCapacityStages
      Mode = DehumidModeNum*2 + CapacityStageNum

      IF (DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) == AutoSize) THEN
        IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingSingleSpeed .OR. &
            DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_HeatingEmpirical .OR. &
            DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoSpeed .OR. &
            DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoStageWHumControl .OR. &
            DXCoil(DXCoilNum)%DXCoilType_Num == CoilVRF_Cooling .OR. &
            DXCoil(DXCoilNum)%DXCoilType_Num == CoilVRF_Heating) THEN

          IF (CurSysNum > 0) THEN

            CALL CheckSysSizing(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name)
            IF (CurOASysNum > 0) THEN
              DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) = FinalSysSizing(CurSysNum)%DesOutAirVolFlow
            ELSE
              DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) = FinalSysSizing(CurSysNum)%DesMainVolFlow
            END IF

          ELSE IF (CurZoneEqNum > 0) THEN

            CALL CheckZoneSizing(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name)
            DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) =   &
               MAX(FinalZoneSizing(CurZoneEqNum)%DesCoolVolFlow,FinalZoneSizing(CurZoneEqNum)%DesHeatVolFlow)

          END IF

          IF (DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) < SmallAirVolFlow) THEN
            DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) = 0.0
          END IF

          IF (DXCoil(DXCoilNum)%RatedAirVolFlowRateEMSOverrideON(Mode)) THEN
            DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) = DXCoil(DXCoilNum)%RatedAirVolFlowRateEMSOverrideValue(Mode)
        END IF

        END IF

        IF(DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_HeatPumpWaterHeater) THEN

          IF(DXCoil(DXCoilNum)%RatedAirVolFlowRate(1) == AutoCalculate)THEN
            DXCoil(DXCoilNum)%RatedAirVolFlowRate(1) = DXCoil(DXCoilNum)%RatedTotCap2 * 0.00005035d0
            DXCoil(DXCoilNum)%AirVolFlowAutoSized = .TRUE.
          END IF

          IF(DXCoil(DXCoilNum)%RatedHPWHCondWaterFlow == AutoCalculate)THEN
            DXCoil(DXCoilNum)%RatedHPWHCondWaterFlow = DXCoil(DXCoilNum)%RatedTotCap2 * 0.00000004487d0
            DXCoil(DXCoilNum)%WaterVolFlowAutoSized = .TRUE.
!            Reporting autosize info for DX coils used with HPWHs will list the info out of order in the eio, report it later
!            CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
!                                  'Rated Condenser Water Volume Flow Rate [m3/s]', DXCoil(DXCoilNum)%RatedHPWHCondWaterFlow)
          END IF

        ELSE

          IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoStageWHumControl) THEN
            CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                                  'Coil Performance Mode='//TRIM(DXCoil(DXCoilNum)%CoilPerformanceName(Mode))//': '// &
                                  'Rated Air Flow Rate [m3/s]', DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode))
               ! For Multimode Coil, Rated flow must be adjusted for bypass fraction
            DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode)= &
               DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) * (1-DXCoil(DXCoilNum)%BypassedFlowFrac(Mode))
            CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                                  'Coil Performance Mode='//TRIM(DXCoil(DXCoilNum)%CoilPerformanceName(Mode))//': '// &
                                  'Rated Air Flow Rate (non-bypassed) [m3/s]', DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode))
          ELSEIF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoSpeed) THEN
            CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                                  'Rated High Speed Air Flow Rate [m3/s]', DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode))
          ELSE
            CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                                  'Rated Air Flow Rate [m3/s]', DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode))
          END IF

        END IF

      END IF

      IF (DXCoil(DXCoilNum)%RatedTotCap(Mode) == AutoSize) THEN

        IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingSingleSpeed .OR. &
            DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoSpeed .OR. &
            DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoStageWHumControl .OR. &
            DXCoil(DXCoilNum)%DXCoilType_Num == CoilVRF_Cooling) THEN

          IF (CurSysNum > 0) THEN

            CALL CheckSysSizing(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name)
            VolFlowRate = DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode)
            IF (VolFlowRate >= SmallAirVolFlow) THEN
              IF (CurOASysNum > 0) THEN ! coil is in the OA stream
                MixTemp = FinalSysSizing(CurSysNum)%CoolOutTemp
                MixHumRat = FinalSysSizing(CurSysNum)%CoolOutHumRat
                SupTemp = FinalSysSizing(CurSysNum)%PrecoolTemp
                SupHumRat = FinalSysSizing(CurSysNum)%PrecoolHumRat
              ELSE ! coil is on the main air loop
           !     MixTemp = FinalSysSizing(CurSysNum)%CoolMixTemp
           !     MixHumRat = FinalSysSizing(CurSysNum)%CoolMixHumRat
                SupTemp = FinalSysSizing(CurSysNum)%CoolSupTemp
                SupHumRat = FinalSysSizing(CurSysNum)%CoolSupHumRat
                IF (PrimaryAirSystem(CurSysNum)%NumOACoolCoils == 0) THEN ! there is no precooling of the OA stream
                  MixTemp = FinalSysSizing(CurSysNum)%CoolMixTemp
                  MixHumRat = FinalSysSizing(CurSysNum)%CoolMixHumRat
                ELSE ! there is precooling of OA stream
                  IF (VolFlowRate > 0.0) THEN
                    OutAirFrac = FinalSysSizing(CurSysNum)%DesOutAirVolFlow / VolFlowRate
                  ELSE
                    OutAirFrac = 1.0d0
                  END IF
                  OutAirFrac = MIN(1.0d0,MAX(0.0d0,OutAirFrac))
                  MixTemp = OutAirFrac*FinalSysSizing(CurSysNum)%PrecoolTemp + &
                              (1.0d0-OutAirFrac)*FinalSysSizing(CurSysNum)%CoolRetTemp
                  MixHumRat = OutAirFrac*FinalSysSizing(CurSysNum)%PrecoolHumRat + &
                                (1.0d0-OutAirFrac)*FinalSysSizing(CurSysNum)%CoolRetHumRat
                END IF
              END IF
              OutTemp = FinalSysSizing(CurSysNum)%CoolOutTemp
              rhoair = PsyRhoAirFnPbTdbW(StdBaroPress,MixTemp,MixHumRat,RoutineName)
              MixEnth = PsyHFnTdbW(MixTemp,MixHumRat,RoutineName)
              MixWetBulb = PsyTwbFnTdbWPb(MixTemp,MixHumRat,StdBaroPress,RoutineName)
              SupEnth = PsyHFnTdbW(SupTemp,SupHumRat,RoutineName)
              TotCapTempModFac = CurveValue(DXCoil(DXCoilNum)%CCapFTemp(Mode),MixWetBulb,OutTemp)
              CoolCapAtPeak = rhoair * VolFlowRate * (MixEnth-SupEnth)
              DXCoil(DXCoilNum)%RatedTotCap(Mode) = CoolCapAtPeak / TotCapTempModFac
              RatedVolFlowPerRatedTotCap = DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) / DXCoil(DXCoilNum)%RatedTotCap(Mode)
              ! check capacity to make sure design volume flow per total capacity is within range
              IF (RatedVolFlowPerRatedTotCap .LT. MinRatedVolFlowPerRatedTotCap) THEN
                IF (.NOT. DXCoil(DXCoilNum)%RatedTotCapEMSOverrideOn(Mode) .AND. DisplayExtraWarnings)THEN
                  CALL ShowWarningError('SizeDXCoil: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name))
                  CALL ShowContinueError('...Rated Total Cooling Capacity will be limited by the minimum rated volume flow per'// &
                                       ' rated total capacity ratio.')
                  CALL ShowContinueError('...DX coil volume flow rate (m3/s) = '// &
                                        TrimSigDigits(DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode),6))
                  CALL ShowContinueError('...Requested capacity (W) = '//TrimSigDigits(DXCoil(DXCoilNum)%RatedTotCap(Mode),3))
                  CALL ShowContinueError('...Requested flow/capacity ratio (m3/s/W) = '// &
                                          TrimSigDigits(RatedVolFlowPerRatedTotCap,3))
                  CALL ShowContinueError('...Minimum flow/capacity ratio (m3/s/W) = '// &
                                          TrimSigDigits(MinRatedVolFlowPerRatedTotCap,3))
                END IF
                DXCoil(DXCoilNum)%RatedTotCap(Mode) = DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) / MinRatedVolFlowPerRatedTotCap
                IF (.NOT. DXCoil(DXCoilNum)%RatedTotCapEMSOverrideOn(Mode) .AND. DisplayExtraWarnings)THEN
                  CALL ShowContinueError('...Adjusted capacity (W) = '//TrimSigDigits(DXCoil(DXCoilNum)%RatedTotCap(Mode),3))
                END IF
              ELSEIF (RatedVolFlowPerRatedTotCap .GT. MaxRatedVolFlowPerRatedTotCap) THEN
                IF (.NOT. DXCoil(DXCoilNum)%RatedTotCapEMSOverrideOn(Mode) .AND. DisplayExtraWarnings)THEN
                  CALL ShowWarningError('SizeDXCoil: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name))
                  CALL ShowContinueError('...Rated Total Cooling Capacity will be limited by the maximum rated volume flow per'// &
                                       ' rated total capacity ratio.')
                  CALL ShowContinueError('...DX coil volume flow rate (m3/s) = '// &
                                        TrimSigDigits(DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode),6))
                  CALL ShowContinueError('...Requested capacity (W) = '//TrimSigDigits(DXCoil(DXCoilNum)%RatedTotCap(Mode),3))
                  CALL ShowContinueError('...Requested flow/capacity ratio (m3/s/W) = '// &
                                          TrimSigDigits(RatedVolFlowPerRatedTotCap,3))
                  CALL ShowContinueError('...Maximum flow/capacity ratio (m3/s/W) = '// &
                                          TrimSigDigits(MaxRatedVolFlowPerRatedTotCap,3))
                END IF
                DXCoil(DXCoilNum)%RatedTotCap(Mode) = DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) / MaxRatedVolFlowPerRatedTotCap
                IF (.NOT. DXCoil(DXCoilNum)%RatedTotCapEMSOverrideOn(Mode) .AND. DisplayExtraWarnings)THEN
                  CALL ShowContinueError('...Adjusted capacity (W) = '//TrimSigDigits(DXCoil(DXCoilNum)%RatedTotCap(Mode),3))
                END IF
              END IF
            ELSE
              DXCoil(DXCoilNum)%RatedTotCap(Mode) = 0.0
            END IF

          ELSE IF (CurZoneEqNum > 0) THEN

            CALL CheckZoneSizing(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name)
            VolFlowRate = DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode)
            IF (VolFlowRate >= SmallAirVolFlow) THEN
              IF(ZoneEqDXCoil)THEN
                IF (ZoneEqSizing(CurZoneEqNum)%OAVolFlow > 0.0) THEN
                  MixTemp = FinalZoneSizing(CurZoneEqNum)%DesCoolCoilInTemp
                  MixHumRat = FinalZoneSizing(CurZoneEqNum)%DesCoolCoilInHumRat
                ELSE
                  MixTemp = FinalZoneSizing(CurZoneEqNum)%ZoneRetTempAtCoolPeak
                  MixHumRat = FinalZoneSizing(CurZoneEqNum)%ZoneHumRatAtCoolPeak
                END IF
              ELSE
                MixTemp = FinalZoneSizing(CurZoneEqNum)%DesCoolCoilInTemp
                MixHumRat = FinalZoneSizing(CurZoneEqNum)%DesCoolCoilInHumRat
              END IF
              SupTemp = FinalZoneSizing(CurZoneEqNum)%CoolDesTemp
              SupHumRat = FinalZoneSizing(CurZoneEqNum)%CoolDesHumRat
              TimeStepNumAtMax = FinalZoneSizing(CurZoneEqNum)%TimeStepNumAtCoolMax
              DDNum = FinalZoneSizing(CurZoneEqNum)%CoolDDNum
              IF (DDNum > 0 .and. TimeStepNumAtMax > 0) THEN
                OutTemp = DesDayWeath(DDNum)%Temp(TimeStepNumAtMax)
              ELSE
                OutTemp = 0.0
              ENDIF
              rhoair = PsyRhoAirFnPbTdbW(StdBaroPress,MixTemp,MixHumRat,RoutineName)
              MixEnth = PsyHFnTdbW(MixTemp,MixHumRat,RoutineName)
              MixWetBulb = PsyTwbFnTdbWPb(MixTemp,MixHumRat,StdBaroPress,RoutineName)
              SupEnth = PsyHFnTdbW(SupTemp,SupHumRat,RoutineName)
              TotCapTempModFac = CurveValue(DXCoil(DXCoilNum)%CCapFTemp(Mode),MixWetBulb,OutTemp)
              CoolCapAtPeak = rhoair * VolFlowRate * (MixEnth-SupEnth)
              DXCoil(DXCoilNum)%RatedTotCap(Mode) = CoolCapAtPeak / TotCapTempModFac
              RatedVolFlowPerRatedTotCap = DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) / DXCoil(DXCoilNum)%RatedTotCap(Mode)
              ! check capacity to make sure design volume flow per total capacity is within range
              IF (RatedVolFlowPerRatedTotCap .LT. MinRatedVolFlowPerRatedTotCap) THEN
                IF (.NOT. DXCoil(DXCoilNum)%RatedTotCapEMSOverrideOn(Mode) .AND. DisplayExtraWarnings)THEN
                  CALL ShowWarningError('SizeDXCoil: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name))
                  CALL ShowContinueError('...Rated Total Cooling Capacity will be limited by the minimum rated volume flow per'// &
                                       ' rated total capacity ratio.')
                  CALL ShowContinueError('...DX coil volume flow rate (m3/s) = '// &
                                        TrimSigDigits(DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode),6))
                  CALL ShowContinueError('...Requested capacity (W) = '//TrimSigDigits(DXCoil(DXCoilNum)%RatedTotCap(Mode),3))
                  CALL ShowContinueError('...Requested flow/capacity ratio (m3/s/W) = '// &
                                          TrimSigDigits(RatedVolFlowPerRatedTotCap,3))
                  CALL ShowContinueError('...Minimum flow/capacity ratio (m3/s/W) = '// &
                                          TrimSigDigits(MinRatedVolFlowPerRatedTotCap,3))
                END IF
                 DXCoil(DXCoilNum)%RatedTotCap(Mode) = DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) / MinRatedVolFlowPerRatedTotCap
                IF (.NOT. DXCoil(DXCoilNum)%RatedTotCapEMSOverrideOn(Mode) .AND. DisplayExtraWarnings)THEN
                  CALL ShowContinueError('...Adjusted capacity (W) = '//TrimSigDigits(DXCoil(DXCoilNum)%RatedTotCap(Mode),3))
                END IF
              ELSEIF (RatedVolFlowPerRatedTotCap .GT. MaxRatedVolFlowPerRatedTotCap) THEN
                IF (.NOT. DXCoil(DXCoilNum)%RatedTotCapEMSOverrideOn(Mode) .AND. DisplayExtraWarnings)THEN
                  CALL ShowWarningError('SizeDXCoil: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name))
                  CALL ShowContinueError('...Rated Total Cooling Capacity will be limited by the maximum rated volume flow per'// &
                                       ' rated total capacity ratio.')
                  CALL ShowContinueError('...DX coil volume flow rate (m3/s) = '// &
                                        TrimSigDigits(DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode),6))
                  CALL ShowContinueError('...Requested capacity (W) = '//TrimSigDigits(DXCoil(DXCoilNum)%RatedTotCap(Mode),3))
                  CALL ShowContinueError('...Requested flow/capacity ratio (m3/s/W) = '// &
                                          TrimSigDigits(RatedVolFlowPerRatedTotCap,3))
                  CALL ShowContinueError('...Maximum flow/capacity ratio (m3/s/W) = '// &
                                          TrimSigDigits(MaxRatedVolFlowPerRatedTotCap,3))
                END IF
                 DXCoil(DXCoilNum)%RatedTotCap(Mode) = DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) / MaxRatedVolFlowPerRatedTotCap
                IF (.NOT. DXCoil(DXCoilNum)%RatedTotCapEMSOverrideOn(Mode) .AND. DisplayExtraWarnings)THEN
                  CALL ShowContinueError('...Adjusted capacity (W) = '//TrimSigDigits(DXCoil(DXCoilNum)%RatedTotCap(Mode),3))
                END IF
              END IF
            ELSE
              DXCoil(DXCoilNum)%RatedTotCap(Mode) = 0.0
            END IF

          END IF

          IF (DXCoil(DXCoilNum)%RatedTotCapEMSOverrideOn(Mode)) THEN
            DXCoil(DXCoilNum)%RatedTotCap(Mode) = DXCoil(DXCoilNum)%RatedTotCapEMSOverrideValue(Mode)
          ENDIF

          IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoStageWHumControl) THEN
            CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                                    'Coil Performance Mode='//TRIM(DXCoil(DXCoilNum)%CoilPerformanceName(Mode))//': '// &
                                    'Rated Total Cooling Capacity (gross) [W]', DXCoil(DXCoilNum)%RatedTotCap(Mode))
          ELSEIF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoSpeed) THEN
            CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                                    'Rated High Speed Total Cooling Capacity (gross) [W]', DXCoil(DXCoilNum)%RatedTotCap(Mode))
          ELSE
            CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                                    'Rated Total Cooling Capacity (gross) [W]', DXCoil(DXCoilNum)%RatedTotCap(Mode))
          END IF

        END IF

      END IF

      IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingSingleSpeed .OR. &
            DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoSpeed .OR. &
            DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoStageWHumControl .OR. &
            DXCoil(DXCoilNum)%DXCoilType_Num == CoilVRF_Cooling) THEN
        DXCoolCap = DXCoil(DXCoilNum)%RatedTotCap(Mode)
      END IF

      IF ((DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_HeatingEmpirical .OR. &
            DXCoil(DXCoilNum)%DXCoilType_Num == CoilVRF_Heating) .AND. &
          DXCoil(DXCoilNum)%RatedTotCap(Mode) == AutoSize) THEN
        IF(DXCoil(DXCoilNum)%CoolingCoilPresent)THEN
          DXCoil(DXCoilNum)%RatedTotCap(Mode) = DXCoolCap
        ELSE

          IF (CurSysNum > 0) THEN

            CALL CheckSysSizing(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name)
            VolFlowRate = DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode)
            IF (VolFlowRate >= SmallAirVolFlow) THEN
              IF (CurOASysNum > 0) THEN ! coil is in the OA stream
                MixTemp = FinalSysSizing(CurSysNum)%CoolOutTemp
                MixHumRat = FinalSysSizing(CurSysNum)%CoolOutHumRat
                SupTemp = FinalSysSizing(CurSysNum)%PrecoolTemp
                SupHumRat = FinalSysSizing(CurSysNum)%PrecoolHumRat
              ELSE ! coil is on the main air loop
           !     MixTemp = FinalSysSizing(CurSysNum)%CoolMixTemp
           !     MixHumRat = FinalSysSizing(CurSysNum)%CoolMixHumRat
                SupTemp = FinalSysSizing(CurSysNum)%CoolSupTemp
                SupHumRat = FinalSysSizing(CurSysNum)%CoolSupHumRat
                IF (PrimaryAirSystem(CurSysNum)%NumOACoolCoils == 0) THEN ! there is no precooling of the OA stream
                  MixTemp = FinalSysSizing(CurSysNum)%CoolMixTemp
                  MixHumRat = FinalSysSizing(CurSysNum)%CoolMixHumRat
                ELSE ! there is precooling of OA stream
                  IF (VolFlowRate > 0.0) THEN
                    OutAirFrac = FinalSysSizing(CurSysNum)%DesOutAirVolFlow / VolFlowRate
                  ELSE
                    OutAirFrac = 1.0d0
                  END IF
                  OutAirFrac = MIN(1.0d0,MAX(0.0d0,OutAirFrac))
                  MixTemp = OutAirFrac*FinalSysSizing(CurSysNum)%PrecoolTemp + &
                              (1.0d0-OutAirFrac)*FinalSysSizing(CurSysNum)%CoolRetTemp
                  MixHumRat = OutAirFrac*FinalSysSizing(CurSysNum)%PrecoolHumRat + &
                                (1.0d0-OutAirFrac)*FinalSysSizing(CurSysNum)%CoolRetHumRat
                END IF
              END IF
              OutTemp = FinalSysSizing(CurSysNum)%CoolOutTemp
              rhoair = PsyRhoAirFnPbTdbW(StdBaroPress,MixTemp,MixHumRat,RoutineName)
              MixEnth = PsyHFnTdbW(MixTemp,MixHumRat,RoutineName)
              MixWetBulb = PsyTwbFnTdbWPb(MixTemp,MixHumRat,StdBaroPress,RoutineName)
              SupEnth = PsyHFnTdbW(SupTemp,SupHumRat,RoutineName)
              TotCapTempModFac = CurveValue(DXCoil(DXCoilNum)%CCapFTemp(Mode),MixWetBulb,OutTemp)
              CoolCapAtPeak = rhoair * VolFlowRate * (MixEnth-SupEnth)
              DXCoil(DXCoilNum)%RatedTotCap(Mode) = CoolCapAtPeak / TotCapTempModFac
              RatedVolFlowPerRatedTotCap = DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) / DXCoil(DXCoilNum)%RatedTotCap(Mode)
              ! check capacity to make sure design volume flow per total capacity is within range
              IF (RatedVolFlowPerRatedTotCap .LT. MinRatedVolFlowPerRatedTotCap) THEN
                IF(DisplayExtraWarnings)THEN
                  CALL ShowWarningError('SizeDXCoil: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name))
                  CALL ShowContinueError('...Rated Total Cooling Capacity will be limited by the minimum rated volume flow per'// &
                                         ' rated total capacity ratio.')
                  CALL ShowContinueError('...DX coil volume flow rate (m3/s) = '// &
                                          TrimSigDigits(DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode),6))
                  CALL ShowContinueError('...Requested capacity (W) = '//TrimSigDigits(DXCoil(DXCoilNum)%RatedTotCap(Mode),3))
                  CALL ShowContinueError('...Requested flow/capacity ratio (m3/s/W) = '// &
                                          TrimSigDigits(RatedVolFlowPerRatedTotCap,3))
                  CALL ShowContinueError('...Minimum flow/capacity ratio (m3/s/W) = '// &
                                          TrimSigDigits(MinRatedVolFlowPerRatedTotCap,3))
                END IF
                DXCoil(DXCoilNum)%RatedTotCap(Mode) = DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) / MinRatedVolFlowPerRatedTotCap
                IF(DisplayExtraWarnings)THEN
                  CALL ShowContinueError('...Adjusted capacity (W) = '//TrimSigDigits(DXCoil(DXCoilNum)%RatedTotCap(Mode),3))
                END IF
              ELSEIF (RatedVolFlowPerRatedTotCap .GT. MaxRatedVolFlowPerRatedTotCap) THEN
                IF(DisplayExtraWarnings)THEN
                  CALL ShowWarningError('SizeDXCoil: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name))
                  CALL ShowContinueError('...Rated Total Cooling Capacity will be limited by the maximum rated volume flow per'// &
                                         ' rated total capacity ratio.')
                  CALL ShowContinueError('...DX coil volume flow rate (m3/s) = '// &
                                          TrimSigDigits(DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode),6))
                  CALL ShowContinueError('...Requested capacity (W) = '//TrimSigDigits(DXCoil(DXCoilNum)%RatedTotCap(Mode),3))
                  CALL ShowContinueError('...Requested flow/capacity ratio (m3/s/W) = '// &
                                          TrimSigDigits(RatedVolFlowPerRatedTotCap,3))
                  CALL ShowContinueError('...Maximum flow/capacity ratio (m3/s/W) = '// &
                                          TrimSigDigits(MaxRatedVolFlowPerRatedTotCap,3))
                END IF
                DXCoil(DXCoilNum)%RatedTotCap(Mode) = DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) / MaxRatedVolFlowPerRatedTotCap
                IF(DisplayExtraWarnings)THEN
                  CALL ShowContinueError('...Adjusted capacity (W) = '//TrimSigDigits(DXCoil(DXCoilNum)%RatedTotCap(Mode),3))
                END IF
              END IF
            ELSE
              DXCoil(DXCoilNum)%RatedTotCap(Mode) = 0.0
            END IF

          ELSE IF (CurZoneEqNum > 0) THEN

            CALL CheckZoneSizing(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name)
            VolFlowRate = DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode)
            IF (VolFlowRate >= SmallAirVolFlow) THEN
              IF(ZoneEqDXCoil)THEN
                IF (ZoneEqSizing(CurZoneEqNum)%OAVolFlow > 0.0) THEN
                  MixTemp = FinalZoneSizing(CurZoneEqNum)%DesCoolCoilInTemp
                  MixHumRat = FinalZoneSizing(CurZoneEqNum)%DesCoolCoilInHumRat
                ELSE
                  MixTemp = FinalZoneSizing(CurZoneEqNum)%ZoneRetTempAtCoolPeak
                  MixHumRat = FinalZoneSizing(CurZoneEqNum)%ZoneHumRatAtCoolPeak
                END IF
              ELSE
                MixTemp = FinalZoneSizing(CurZoneEqNum)%DesCoolCoilInTemp
                MixHumRat = FinalZoneSizing(CurZoneEqNum)%DesCoolCoilInHumRat
              END IF
              SupTemp = FinalZoneSizing(CurZoneEqNum)%CoolDesTemp
              SupHumRat = FinalZoneSizing(CurZoneEqNum)%CoolDesHumRat
              TimeStepNumAtMax = FinalZoneSizing(CurZoneEqNum)%TimeStepNumAtCoolMax
              DDNum = FinalZoneSizing(CurZoneEqNum)%CoolDDNum
              IF (DDNum > 0 .and. TimeStepNumAtMax > 0) THEN
                OutTemp = DesDayWeath(DDNum)%Temp(TimeStepNumAtMax)
              ELSE
                OutTemp = 0.0
              ENDIF
              rhoair = PsyRhoAirFnPbTdbW(StdBaroPress,MixTemp,MixHumRat,RoutineName)
              MixEnth = PsyHFnTdbW(MixTemp,MixHumRat,RoutineName)
              MixWetBulb = PsyTwbFnTdbWPb(MixTemp,MixHumRat,StdBaroPress,RoutineName)
              SupEnth = PsyHFnTdbW(SupTemp,SupHumRat,RoutineName)
              TotCapTempModFac = CurveValue(DXCoil(DXCoilNum)%CCapFTemp(Mode),MixWetBulb,OutTemp)
              CoolCapAtPeak = rhoair * VolFlowRate * (MixEnth-SupEnth)
              DXCoil(DXCoilNum)%RatedTotCap(Mode) = CoolCapAtPeak / TotCapTempModFac
              RatedVolFlowPerRatedTotCap = DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) / DXCoil(DXCoilNum)%RatedTotCap(Mode)
              ! check capacity to make sure design volume flow per total capacity is within range
              IF (RatedVolFlowPerRatedTotCap .LT. MinRatedVolFlowPerRatedTotCap) THEN
                IF(DisplayExtraWarnings)THEN
                  CALL ShowWarningError('SizeDXCoil: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name))
                  CALL ShowContinueError('...Rated Total Cooling Capacity will be limited by the minimum rated volume flow per'// &
                                         ' rated total capacity ratio.')
                  CALL ShowContinueError('...DX coil volume flow rate (m3/s) = '// &
                                          TrimSigDigits(DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode),6))
                  CALL ShowContinueError('...Requested capacity (W) = '//TrimSigDigits(DXCoil(DXCoilNum)%RatedTotCap(Mode),3))
                  CALL ShowContinueError('...Requested flow/capacity ratio (m3/s/W) = '// &
                                          TrimSigDigits(RatedVolFlowPerRatedTotCap,3))
                  CALL ShowContinueError('...Minimum flow/capacity ratio (m3/s/W) = '// &
                                          TrimSigDigits(MinRatedVolFlowPerRatedTotCap,3))
                END IF
                DXCoil(DXCoilNum)%RatedTotCap(Mode) = DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) / MinRatedVolFlowPerRatedTotCap
                IF(DisplayExtraWarnings)THEN
                  CALL ShowContinueError('...Adjusted capacity (W) = '//TrimSigDigits(DXCoil(DXCoilNum)%RatedTotCap(Mode),3))
                END IF
              ELSEIF (RatedVolFlowPerRatedTotCap .GT. MaxRatedVolFlowPerRatedTotCap) THEN
                IF(DisplayExtraWarnings)THEN
                  CALL ShowWarningError('SizeDXCoil: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name))
                  CALL ShowContinueError('...Rated Total Cooling Capacity will be limited by the maximum rated volume flow per'// &
                                         ' rated total capacity ratio.')
                  CALL ShowContinueError('...DX coil volume flow rate (m3/s) = '// &
                                          TrimSigDigits(DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode),6))
                  CALL ShowContinueError('...Requested capacity (W) = '//TrimSigDigits(DXCoil(DXCoilNum)%RatedTotCap(Mode),3))
                  CALL ShowContinueError('...Requested flow/capacity ratio (m3/s/W) = '// &
                                          TrimSigDigits(RatedVolFlowPerRatedTotCap,3))
                  CALL ShowContinueError('...Maximum flow/capacity ratio (m3/s/W) = '// &
                                          TrimSigDigits(MaxRatedVolFlowPerRatedTotCap,3))
                END IF
                DXCoil(DXCoilNum)%RatedTotCap(Mode) = DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) / MaxRatedVolFlowPerRatedTotCap
                IF(DisplayExtraWarnings)THEN
                  CALL ShowContinueError('...Adjusted capacity (W) = '//TrimSigDigits(DXCoil(DXCoilNum)%RatedTotCap(Mode),3))
                END IF
              END IF
            ELSE
              DXCoil(DXCoilNum)%RatedTotCap(Mode) = 0.0
            END IF

          END IF

        END IF
        CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
               'Rated Total Heating Capacity [W]', DXCoil(DXCoilNum)%RatedTotCap(Mode))
      END IF

      IF (DXCoil(DXCoilNum)%RatedSHR(Mode) == AutoSize) THEN

        IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingSingleSpeed .OR. &
            DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoSpeed .OR. &
            DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoStageWHumControl .OR. &
            DXCoil(DXCoilNum)%DXCoilType_Num == CoilVRF_Cooling) THEN

          CpAir = PsyCpAirFnWTdb(RatedInletAirHumRat,RatedInletAirTemp,RoutineName)

          IF (CurSysNum > 0) THEN

            CALL CheckSysSizing(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name)

          ELSE IF (CurZoneEqNum > 0) THEN

            CALL CheckZoneSizing(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name)

          END IF

          IF (DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) >= SmallAirVolFLow .AND. DXCoil(DXCoilNum)%RatedTotCap(Mode) > 0.0) THEN
    ! For autosizing the rated SHR, we set a minimum SHR of 0.676 and a maximum of 0.798. The min SHR occurs occurs at the
    ! minimum flow / capacity ratio = MinRatedVolFlowPerRatedTotCap = 0.00004027 [m3/s / W] = 300 [cfm/ton].
    ! The max SHR occurs at maximum flow / capacity ratio = MaxRatedVolFlowPerRatedTotCap = 0.00006041 [m3/s / W] = 450 [cfm/ton].
    ! For flow / capacity ratios between the min and max we linearly interpolate between min and max SHR. Thus rated SHR is a
    ! linear function of the rated flow / capacity ratio. This linear function (see below) is the result of a regression
    ! of flow/capacity ratio vs SHR for several actual coils.
            RatedVolFlowPerRatedTotCap = DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) / DXCoil(DXCoilNum)%RatedTotCap(Mode)
            IF (RatedVolFlowPerRatedTotCap > MaxRatedVolFlowPerRatedTotCap) THEN
              DXCoil(DXCoilNum)%RatedSHR(Mode) = 0.431d0 + 6086.d0*MaxRatedVolFlowPerRatedTotCap
            ELSE IF (RatedVolFlowPerRatedTotCap < MinRatedVolFlowPerRatedTotCap) THEN
              DXCoil(DXCoilNum)%RatedSHR(Mode) = 0.431d0 + 6086.d0*MinRatedVolFlowPerRatedTotCap
            ELSE
              DXCoil(DXCoilNum)%RatedSHR(Mode) = 0.431d0 + 6086.d0*RatedVolFlowPerRatedTotCap
            END IF
          ELSE
            DXCoil(DXCoilNum)%RatedSHR(Mode) = 1.0d0
          END IF

          IF (DXCoil(DXCoilNum)%RatedSHREMSOverrideOn(Mode)) THEN
            DXCoil(DXCoilNum)%RatedSHR(Mode) = DXCoil(DXCoilNum)%RatedSHREMSOverrideValue(mode)
          ENDIF

          IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoStageWHumControl) THEN
            CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                                    'Coil Performance Mode='//TRIM(DXCoil(DXCoilNum)%CoilPerformanceName(Mode))//': '// &
                                    'Rated Sensible Heat Ratio', DXCoil(DXCoilNum)%RatedSHR(Mode))
          ELSEIF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoSpeed) THEN
            CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                                    'Rated High Speed Sensible Heat Ratio', DXCoil(DXCoilNum)%RatedSHR(Mode))
          ELSE
            CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                                    'Rated Sensible Heat Ratio', DXCoil(DXCoilNum)%RatedSHR(Mode))
          END IF

        END IF

      END IF

      IF (DXCoil(DXCoilNum)%EvapCondAirFlow(Mode) == AutoSize) THEN

        IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingSingleSpeed .OR. &
            DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoSpeed .OR. &
            DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoStageWHumControl) THEN

    !     Auto size condenser air flow to Total Capacity * 0.000114 m3/s/w (850 cfm/ton)
          DXCoil(DXCoilNum)%EvapCondAirFlow(Mode) = DXCoil(DXCoilNum)%RatedTotCap(Mode)*0.000114d0

          IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoStageWHumControl) THEN
            CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                                    'Coil Performance Mode='//TRIM(DXCoil(DXCoilNum)%CoilPerformanceName(Mode))//': '// &
                                    'Evaporative Condenser Air Flow Rate [m3/s]', DXCoil(DXCoilNum)%EvapCondAirFlow(Mode))
          ELSEIF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoSpeed) THEN
            CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                            'High Speed Evaporative Condenser Air Flow Rate [m3/s]', DXCoil(DXCoilNum)%EvapCondAirFlow(Mode))
          ELSE
            CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                            'Evaporative Condenser Air Flow Rate [m3/s]', DXCoil(DXCoilNum)%EvapCondAirFlow(Mode))
          END IF

        END IF

      END IF

      IF (DXCoil(DXCoilNum)%EvapCondAirFlow2 == AutoSize) THEN

        IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoSpeed) THEN
    !     Auto size low speed condenser air flow to 1/3 Total Capacity * 0.000114 m3/s/w (850 cfm/ton)
          DXCoil(DXCoilNum)%EvapCondAirFlow2 = 0.3333d0*DXCoil(DXCoilNum)%RatedTotCap(Mode)*0.000114d0

          CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                            'Low Speed Evaporative Condenser Air Flow Rate [m3/s]', DXCoil(DXCoilNum)%EvapCondAirFlow2)
        END IF

      END IF


      IF (DXCoil(DXCoilNum)%EvapCondPumpElecNomPower(Mode) == AutoSize) THEN

        IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingSingleSpeed .OR. &
            DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoSpeed .OR. &
            DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoStageWHumControl) THEN

    !     Auto size high speed evap condenser pump power to Total Capacity * 0.004266 w/w (15 w/ton)
          DXCoil(DXCoilNum)%EvapCondPumpElecNomPower(Mode) = DXCoil(DXCoilNum)%RatedTotCap(Mode)*0.004266d0

          IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoStageWHumControl) THEN
            CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                                'Coil Performance Mode='//TRIM(DXCoil(DXCoilNum)%CoilPerformanceName(Mode))//': '// &
                                'Evaporative Condenser Pump Rated Power Consumption [W]',   &
                                 DXCoil(DXCoilNum)%EvapCondPumpElecNomPower(Mode))
          ELSEIF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoSpeed) THEN
            CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                            'High Speed Evaporative Condenser Pump Rated Power Consumption [W]',   &
                               DXCoil(DXCoilNum)%EvapCondAirFlow(Mode))
          ELSE
            CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                                'Evaporative Condenser Pump Rated Power Consumption [W]',   &
                                 DXCoil(DXCoilNum)%EvapCondPumpElecNomPower(Mode))
          END IF

        END IF

      END IF

      IF (DXCoil(DXCoilNum)%EvapCondPumpElecNomPower2 == AutoSize) THEN

         IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoSpeed) THEN

    !     Auto size low speed evap condenser pump power to 1/3 Total Capacity * 0.004266 w/w (15 w/ton)
          DXCoil(DXCoilNum)%EvapCondPumpElecNomPower2 = 0.3333d0*DXCoil(DXCoilNum)%RatedTotCap(Mode)*0.004266d0

          CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                              'Low Speed Evaporative Condenser Pump Rated Power Consumption [W]',   &
                               DXCoil(DXCoilNum)%EvapCondPumpElecNomPower2)
        END IF

      END IF


      IF (DXCoil(DXCoilNum)%RatedAirVolFlowRate2 == AutoSize) THEN

        IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoSpeed) THEN
          DXCoil(DXCoilNum)%RatedAirVolFlowRate2 = 0.3333d0*DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode)

          CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                                  'Rated Low Speed Air Flow Rate [m3/s]', DXCoil(DXCoilNum)%RatedAirVolFlowRate2)
        END IF

      END IF
      IF (DXCoil(DXCoilNum)%RatedTotCap2 == AutoSize) THEN

        IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoSpeed) THEN
          DXCoil(DXCoilNum)%RatedTotCap2 = 0.3333d0*DXCoil(DXCoilNum)%RatedTotCap(Mode)

          CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                                  'Rated Low Speed Total Cooling Capacity (gross) [W]', DXCoil(DXCoilNum)%RatedTotCap2)
        END IF

      END IF


      IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoSpeed) THEN

        IF(DXCoil(DXCoilNum)%EvapCondAirFlow2 .GT. DXCoil(DXCoilNum)%EvapCondAirFlow(Mode)) THEN
          CALL ShowSevereError('SizeDXCoil: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name)//', '// &
                  'Evaporative Condenser low speed air flow must be less than or equal to high speed air flow.')
          CALL ShowContinueError('Instead, '//TRIM(RoundSigDigits(DXCoil(DXCoilNum)%EvapCondAirFlow2,2))//' > '//  &
                  TRIM(RoundSigDigits(DXCoil(DXCoilNum)%EvapCondAirFlow(Mode),2)))
          CALL ShowFatalError('Preceding conditions cause termination.')
        END IF

        IF(DXCoil(DXCoilNum)%EvapCondPumpElecNomPower2 .GT. DXCoil(DXCoilNum)%EvapCondPumpElecNomPower(Mode))THEN
          CALL ShowSevereError('SizeDXCoil: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name)//', '// &
                  'Evaporative Condenser low speed pump power must be less than or equal to high speed pump power.')
          CALL ShowContinueError('Instead, '//TRIM(RoundSigDigits(DXCoil(DXCoilNum)%EvapCondPumpElecNomPower2,2))//' > '//  &
                  TRIM(RoundSigDigits(DXCoil(DXCoilNum)%EvapCondPumpElecNomPower(Mode),2)))
          CALL ShowFatalError('Preceding conditions cause termination.')
        END IF

        IF(DXCoil(DXCoilNum)%RatedTotCap2 .GT. DXCoil(DXCoilNum)%RatedTotCap(Mode))THEN
          CALL ShowSevereError('SizeDXCoil: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name)//', '// &
                  'Rated Total Cooling Capacity, Low Speed must be less than or equal to '// &
                  'Rated Total Cooling Capacity, High Speed.')
          CALL ShowContinueError('Instead, '//TRIM(RoundSigDigits(DXCoil(DXCoilNum)%RatedTotCap2,2))//' > '//  &
                  TRIM(RoundSigDigits(DXCoil(DXCoilNum)%RatedTotCap(Mode),2)))
          CALL ShowFatalError('Preceding conditions cause termination.')
        END IF

        IF(DXCoil(DXCoilNum)%RatedAirVolFlowRate2 .GT. DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode))THEN
          CALL ShowFatalError('SizeDXCoil: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name)//', '// &
                  'Rated Air Volume Flow Rate, low speed must be less than or equal to '//&
                  'Rated Air Volume Flow Rate, high speed.')
          CALL ShowContinueError('Instead, '//TRIM(RoundSigDigits(DXCoil(DXCoilNum)%RatedAirVolFlowRate2,2))//' > '//  &
                  TRIM(RoundSigDigits(DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode),2)))
          CALL ShowFatalError('Preceding conditions cause termination.')
        END IF

      END IF

      IF (DXCoil(DXCoilNum)%RatedSHR2 == AutoSize) THEN

        IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoSpeed) THEN
          DXCoil(DXCoilNum)%RatedSHR2 = DXCoil(DXCoilNum)%RatedSHR(Mode)

          CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                                  'Rated Low Speed Sensible Heat Ratio', DXCoil(DXCoilNum)%RatedSHR2)
        END IF

      END IF

      IF (DXCoil(DXCoilNum)%DefrostCapacity == AutoSize) THEN

        IF (DXCoil(DXCoilNum)%DefrostStrategy == Resistive) THEN

          DXCoil(DXCoilNum)%DefrostCapacity = DXCoolCap

        ELSE

          DXCoil(DXCoilNum)%DefrostCapacity = 0.0

        END IF

        CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                                'Resistive Defrost Heater Capacity', DXCoil(DXCoilNum)%DefrostCapacity)

      END IF
    END DO ! End capacity stages loop
  END DO ! End dehumidification modes loop

  ! Autosizing for multispeed cooling coil
  IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_MultiSpeedCooling) THEN
    ! flow rate auto size
    Do Mode = DXCoil(DXCoilNum)%NumOfSpeeds,1,-1
      If (DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode) == AutoSize) Then
        If (Mode == DXCoil(DXCoilNum)%NumOfSpeeds) Then
          IF (CurSysNum > 0) THEN
            CALL CheckSysSizing(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name)
            IF (CurOASysNum > 0) THEN
              DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode) = FinalSysSizing(CurSysNum)%DesOutAirVolFlow
            ELSE
              DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode) = FinalSysSizing(CurSysNum)%DesMainVolFlow
            END IF
          ELSE IF (CurZoneEqNum > 0) THEN
            CALL CheckZoneSizing(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name)
            DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode) = MAX(FinalZoneSizing(CurZoneEqNum)%DesCoolVolFlow, &
                                                      FinalZoneSizing(CurZoneEqNum)%DesHeatVolFlow)
          END IF
          IF (DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) < SmallAirVolFlow) THEN
            DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) = 0.0
          END IF
        Else
          DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode) = DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(DXCoil(DXCoilNum)%NumOfSpeeds)* &
                                                          Mode/DXCoil(DXCoilNum)%NumOfSpeeds
        End If
        CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
             'Speed '//Trim(TrimSigDigits(Mode))//' Rated Air Flow Rate [m3/s]',  &
             DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode))
      End If
    End Do

   ! Ensure flow rate at lower speed must be lower or equal to the flow rate at higher speed. Otherwise, a severe error is isssued.
    Do Mode = 1,DXCoil(DXCoilNum)%NumOfSpeeds-1
      If (DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode) .GT. DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode+1)) Then
        CALL ShowWarningError('SizeDXCoil: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name)//', '// &
          'Speed '//Trim(TrimSigDigits(Mode))//' Rated Air Flow Rate must be less than or equal to '//&
          'Speed '//Trim(TrimSigDigits(Mode+1))//' Rated Air Flow Rate.')
        CALL ShowContinueError('Instead, '//TRIM(RoundSigDigits(DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode),2))//' > '//  &
                  TRIM(RoundSigDigits(DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode+1),2)))
        CALL ShowFatalError('Preceding conditions cause termination.')
      End If
    End Do

    Do Mode = DXCoil(DXCoilNum)%NumOfSpeeds,1,-1
      If (DXCoil(DXCoilNum)%MSRatedTotCap(Mode) == AutoSize) Then
        IF (Mode .eq. DXCoil(DXCoilNum)%NumOfSpeeds) Then
          IF (CurSysNum > 0) THEN
            CALL CheckSysSizing(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name)
            VolFlowRate = DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode)
            IF (VolFlowRate >= SmallAirVolFlow) THEN
             IF (CurOASysNum > 0) THEN ! coil is in the OA stream
                MixTemp = FinalSysSizing(CurSysNum)%CoolOutTemp
                MixHumRat = FinalSysSizing(CurSysNum)%CoolOutHumRat
                SupTemp = FinalSysSizing(CurSysNum)%PrecoolTemp
                SupHumRat = FinalSysSizing(CurSysNum)%PrecoolHumRat
              ELSE ! coil is on the main air loop
                SupTemp = FinalSysSizing(CurSysNum)%CoolSupTemp
                SupHumRat = FinalSysSizing(CurSysNum)%CoolSupHumRat
                IF (PrimaryAirSystem(CurSysNum)%NumOACoolCoils == 0) THEN ! there is no precooling of the OA stream
                  MixTemp = FinalSysSizing(CurSysNum)%CoolMixTemp
                  MixHumRat = FinalSysSizing(CurSysNum)%CoolMixHumRat
                ELSE ! there is precooling of OA stream
                  IF (VolFlowRate > 0.0) THEN
                    OutAirFrac = FinalSysSizing(CurSysNum)%DesOutAirVolFlow / VolFlowRate
                  ELSE
                    OutAirFrac = 1.0d0
                  END IF
                  OutAirFrac = MIN(1.0d0,MAX(0.0d0,OutAirFrac))
                  MixTemp = OutAirFrac*FinalSysSizing(CurSysNum)%PrecoolTemp + &
                              (1.0d0-OutAirFrac)*FinalSysSizing(CurSysNum)%CoolRetTemp
                  MixHumRat = OutAirFrac*FinalSysSizing(CurSysNum)%PrecoolHumRat + &
                                (1.0d0-OutAirFrac)*FinalSysSizing(CurSysNum)%CoolRetHumRat
                END IF
              END IF
              OutTemp = FinalSysSizing(CurSysNum)%CoolOutTemp
              rhoair = PsyRhoAirFnPbTdbW(StdBaroPress,MixTemp,MixHumRat,RoutineName)
              MixEnth = PsyHFnTdbW(MixTemp,MixHumRat,RoutineName)
              MixWetBulb = PsyTwbFnTdbWPb(MixTemp,MixHumRat,StdBaroPress,RoutineName)
              SupEnth = PsyHFnTdbW(SupTemp,SupHumRat,RoutineName)
              TotCapTempModFac = CurveValue(DXCoil(DXCoilNum)%MSCCapFTemp(Mode),MixWetBulb,OutTemp)
              CoolCapAtPeak = rhoair * VolFlowRate * (MixEnth-SupEnth)
              DXCoil(DXCoilNum)%MSRatedTotCap(Mode) = CoolCapAtPeak / TotCapTempModFac
              RatedVolFlowPerRatedTotCap = DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode) / DXCoil(DXCoilNum)%MSRatedTotCap(Mode)
              ! check capacity to make sure design volume flow per total capacity is within range
              IF (RatedVolFlowPerRatedTotCap .LT. MinRatedVolFlowPerRatedTotCap) THEN
                IF(DisplayExtraWarnings)THEN
                  CALL ShowWarningError('SizeDXCoil: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name))
                  CALL ShowContinueError('...Rated Total Cooling Capacity will be limited by the minimum rated volume flow per'// &
                                         ' rated total capacity ratio.')
                  CALL ShowContinueError('...DX coil speed = '//TrimSigDigits(Mode,0))
                  CALL ShowContinueError('...DX coil volume flow rate (m3/s) = '// &
                                          TrimSigDigits(DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode),6))
                  CALL ShowContinueError('...Requested capacity (W) = '//TrimSigDigits(DXCoil(DXCoilNum)%MSRatedTotCap(Mode),3))
                  CALL ShowContinueError('...Requested flow/capacity ratio (m3/s/W) = '// &
                                          TrimSigDigits(RatedVolFlowPerRatedTotCap,3))
                  CALL ShowContinueError('...Minimum flow/capacity ratio (m3/s/W) = '// &
                                          TrimSigDigits(MinRatedVolFlowPerRatedTotCap,3))
                END IF
                DXCoil(DXCoilNum)%MSRatedTotCap(Mode) =   &
                   DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode) / MinRatedVolFlowPerRatedTotCap
                IF(DisplayExtraWarnings)THEN
                  CALL ShowContinueError('...Adjusted capacity (W) = '//TrimSigDigits(DXCoil(DXCoilNum)%MSRatedTotCap(Mode),3))
                END IF
              ELSEIF (RatedVolFlowPerRatedTotCap .GT. MaxRatedVolFlowPerRatedTotCap) THEN
                IF(DisplayExtraWarnings)THEN
                  CALL ShowWarningError('SizeDXCoil: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name))
                  CALL ShowContinueError('...Rated Total Cooling Capacity will be limited by the maximum rated volume flow per'// &
                                         ' rated total capacity ratio.')
                  CALL ShowContinueError('...DX coil speed = '//TrimSigDigits(Mode,0))
                  CALL ShowContinueError('...DX coil volume flow rate (m3/s) = '// &
                                          TrimSigDigits(DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode),6))
                  CALL ShowContinueError('...Requested capacity (W) = '//TrimSigDigits(DXCoil(DXCoilNum)%MSRatedTotCap(Mode),3))
                  CALL ShowContinueError('...Requested flow/capacity ratio (m3/s/W) = '// &
                                          TrimSigDigits(RatedVolFlowPerRatedTotCap,3))
                  CALL ShowContinueError('...Maximum flow/capacity ratio (m3/s/W) = '// &
                                          TrimSigDigits(MaxRatedVolFlowPerRatedTotCap,3))
                END IF
                DXCoil(DXCoilNum)%MSRatedTotCap(Mode) =   &
                   DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode) / MaxRatedVolFlowPerRatedTotCap
                IF(DisplayExtraWarnings)THEN
                  CALL ShowContinueError('...Adjusted capacity (W) = '//TrimSigDigits(DXCoil(DXCoilNum)%MSRatedTotCap(Mode),3))
                END IF
              END IF
            ELSE
              DXCoil(DXCoilNum)%MSRatedTotCap(Mode) = 0.0
            END IF
          ELSE IF (CurZoneEqNum > 0) THEN
            CALL CheckZoneSizing(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name)
            VolFlowRate = DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode)
            IF (VolFlowRate >= SmallAirVolFlow) THEN
              MixTemp = FinalZoneSizing(CurZoneEqNum)%DesCoolCoilInTemp
              MixHumRat = FinalZoneSizing(CurZoneEqNum)%DesCoolCoilInHumRat
              SupTemp = FinalZoneSizing(CurZoneEqNum)%CoolDesTemp
              SupHumRat = FinalZoneSizing(CurZoneEqNum)%CoolDesHumRat
              TimeStepNumAtMax = FinalZoneSizing(CurZoneEqNum)%TimeStepNumAtCoolMax
              DDNum = FinalZoneSizing(CurZoneEqNum)%CoolDDNum
              IF (DDNum > 0 .and. TimeStepNumAtMax > 0) THEN
                OutTemp = DesDayWeath(DDNum)%Temp(TimeStepNumAtMax)
              ELSE
                OutTemp = 0.0
              ENDIF
              rhoair = PsyRhoAirFnPbTdbW(StdBaroPress,MixTemp,MixHumRat,RoutineName)
              MixEnth = PsyHFnTdbW(MixTemp,MixHumRat,RoutineName)
              MixWetBulb = PsyTwbFnTdbWPb(MixTemp,MixHumRat,StdBaroPress,RoutineName)
              SupEnth = PsyHFnTdbW(SupTemp,SupHumRat,RoutineName)
              TotCapTempModFac = CurveValue(DXCoil(DXCoilNum)%MSCCapFTemp(Mode),MixWetBulb,OutTemp)
              CoolCapAtPeak = rhoair * VolFlowRate * (MixEnth-SupEnth)
              DXCoil(DXCoilNum)%MSRatedTotCap(Mode) = CoolCapAtPeak / TotCapTempModFac
              RatedVolFlowPerRatedTotCap = DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode) / DXCoil(DXCoilNum)%MSRatedTotCap(Mode)
              ! check capacity to make sure design volume flow per total capacity is within range
              IF (RatedVolFlowPerRatedTotCap .LT. MinRatedVolFlowPerRatedTotCap) THEN
                IF(DisplayExtraWarnings)THEN
                  CALL ShowWarningError('SizeDXCoil: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name))
                  CALL ShowContinueError('...Rated Total Cooling Capacity will be limited by the minimum rated volume flow per'// &
                                         ' rated total capacity ratio.')
                  CALL ShowContinueError('...DX coil speed = '//TrimSigDigits(Mode,0))
                  CALL ShowContinueError('...DX coil volume flow rate (m3/s) = '// &
                                          TrimSigDigits(DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode),6))
                  CALL ShowContinueError('...Requested capacity (W) = '//TrimSigDigits(DXCoil(DXCoilNum)%MSRatedTotCap(Mode),3))
                  CALL ShowContinueError('...Requested flow/capacity ratio (m3/s/W) = '// &
                                          TrimSigDigits(RatedVolFlowPerRatedTotCap,3))
                  CALL ShowContinueError('...Minimum flow/capacity ratio (m3/s/W) = '// &
                                          TrimSigDigits(MinRatedVolFlowPerRatedTotCap,3))
                END IF
                DXCoil(DXCoilNum)%MSRatedTotCap(Mode) =DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode)/MinRatedVolFlowPerRatedTotCap
                IF(DisplayExtraWarnings)THEN
                  CALL ShowContinueError('...Adjusted capacity (W) = '//TrimSigDigits(DXCoil(DXCoilNum)%MSRatedTotCap(Mode),3))
                END IF
              ELSEIF (RatedVolFlowPerRatedTotCap .GT. MaxRatedVolFlowPerRatedTotCap) THEN
                IF(DisplayExtraWarnings)THEN
                  CALL ShowWarningError('SizeDXCoil: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name))
                  CALL ShowContinueError('...Rated Total Cooling Capacity will be limited by the maximum rated volume flow per'// &
                                         ' rated total capacity ratio.')
                  CALL ShowContinueError('...DX coil speed = '//TrimSigDigits(Mode,0))
                  CALL ShowContinueError('...DX coil volume flow rate (m3/s) = '// &
                                          TrimSigDigits(DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode),6))
                  CALL ShowContinueError('...Requested capacity (W) = '//TrimSigDigits(DXCoil(DXCoilNum)%MSRatedTotCap(Mode),3))
                  CALL ShowContinueError('...Requested flow/capacity ratio (m3/s/W) = '// &
                                          TrimSigDigits(RatedVolFlowPerRatedTotCap,3))
                  CALL ShowContinueError('...Maximum flow/capacity ratio (m3/s/W) = '// &
                                          TrimSigDigits(MaxRatedVolFlowPerRatedTotCap,3))
                END IF
                DXCoil(DXCoilNum)%MSRatedTotCap(Mode) =DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode)/MaxRatedVolFlowPerRatedTotCap
                IF(DisplayExtraWarnings)THEN
                  CALL ShowContinueError('...Adjusted capacity (W) = '//TrimSigDigits(DXCoil(DXCoilNum)%MSRatedTotCap(Mode),3))
                END IF
              END IF
            ELSE
              DXCoil(DXCoilNum)%MSRatedTotCap(Mode) = 0.0
            END IF
          END IF
        Else
          DXCoil(DXCoilNum)%MSRatedTotCap(Mode) = DXCoil(DXCoilNum)%MSRatedTotCap(DXCoil(DXCoilNum)%NumOfSpeeds)* &
                                                         Mode/DXCoil(DXCoilNum)%NumOfSpeeds
        End If
        CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                 'Speed '//Trim(TrimSigDigits(Mode))//' Rated Total Cooling Capacity [W]', DXCoil(DXCoilNum)%MSRatedTotCap(Mode))
      End If
    End Do

   ! Ensure capacity at lower speed must be lower or equal to the capacity at higher speed.
    Do Mode = 1,DXCoil(DXCoilNum)%NumOfSpeeds-1
      If (DXCoil(DXCoilNum)%MSRatedTotCap(Mode) .GT. DXCoil(DXCoilNum)%MSRatedTotCap(Mode+1)) Then
        CALL ShowWarningError('SizeDXCoil: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name)//', '// &
          'Speed '//Trim(TrimSigDigits(Mode))//' Rated Total Cooling Capacity must be less than or equal to '//&
          'Speed '//Trim(TrimSigDigits(Mode+1))//' Rated Total Cooling Capacity.')
        CALL ShowContinueError('Instead, '//TRIM(RoundSigDigits(DXCoil(DXCoilNum)%MSRatedTotCap(Mode),2))//' > '//  &
                  TRIM(RoundSigDigits(DXCoil(DXCoilNum)%MSRatedTotCap(Mode+1),2)))
        CALL ShowFatalError('Preceding conditions cause termination.')
      End If
    End Do

    ! Rated SHR
    Do Mode = DXCoil(DXCoilNum)%NumOfSpeeds,1,-1
      If (DXCoil(DXCoilNum)%MSRatedSHR(Mode) == AutoSize) Then
        IF (Mode .eq. DXCoil(DXCoilNum)%NumOfSpeeds) Then
          IF (CurSysNum > 0) THEN
            CALL CheckSysSizing(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name)
          ELSE IF (CurZoneEqNum > 0) THEN
            CALL CheckZoneSizing(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name)
          END IF
          IF (DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode) >= SmallAirVolFLow .AND.   &
              DXCoil(DXCoilNum)%MSRatedTotCap(Mode) > 0.0) THEN
    ! For autosizing the rated SHR, we set a minimum SHR of 0.676 and a maximum of 0.798. The min SHR occurs occurs at the
    ! minimum flow / capacity ratio = MinRatedVolFlowPerRatedTotCap = 0.00004027 [m3/s / W] = 300 [cfm/ton].
    ! The max SHR occurs at maximum flow / capacity ratio = MaxRatedVolFlowPerRatedTotCap = 0.00006041 [m3/s / W] = 450 [cfm/ton].
    ! For flow / capacity ratios between the min and max we linearly interpolate between min and max SHR. Thus rated SHR is a
    ! linear function of the rated flow / capacity ratio. This linear function (see below) is the result of a regression
    ! of flow/capacity ratio vs SHR for several actual coils.
            RatedVolFlowPerRatedTotCap = DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode) / DXCoil(DXCoilNum)%MSRatedTotCap(Mode)
            IF (RatedVolFlowPerRatedTotCap > MaxRatedVolFlowPerRatedTotCap) THEN
              DXCoil(DXCoilNum)%RatedSHR(Mode) = 0.431d0 + 6086.d0*MaxRatedVolFlowPerRatedTotCap
            ELSE IF (RatedVolFlowPerRatedTotCap < MinRatedVolFlowPerRatedTotCap) THEN
              DXCoil(DXCoilNum)%RatedSHR(Mode) = 0.431d0 + 6086.d0*MinRatedVolFlowPerRatedTotCap
            ELSE
              DXCoil(DXCoilNum)%MSRatedSHR(Mode) = 0.431d0 + 6086.d0*RatedVolFlowPerRatedTotCap
            END IF
          ELSE
            DXCoil(DXCoilNum)%MSRatedSHR(Mode) = 1.0d0
          END IF
        Else
          DXCoil(DXCoilNum)%MSRatedSHR(Mode) = DXCoil(DXCoilNum)%MSRatedSHR(Mode+1)
        End If
        CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
              'Speed '//Trim(TrimSigDigits(Mode))//' Rated Sensible Heat Ratio', DXCoil(DXCoilNum)%MSRatedSHR(Mode))
      End If
    End Do

    ! Rated Evapovative condenser airflow rates
    Do Mode = 1,DXCoil(DXCoilNum)%NumOfSpeeds
      If (DXCoil(DXCoilNum)%MSEvapCondAirFlow(Mode) == AutoSize) Then
        ! Auto size condenser air flow to Total Capacity * 0.000114 m3/s/w (850 cfm/ton)
        DXCoil(DXCoilNum)%MSEvapCondAirFlow(Mode) = DXCoil(DXCoilNum)%MSRatedTotCap(Mode)*0.000114d0
        CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
          'Speed '//TRIM(TrimSigDigits(Mode))//' Evaporative Condenser Air Flow Rate [m3/s]', &
          DXCoil(DXCoilNum)%MSEvapCondAirFlow(mode))
      End If
    End Do
   ! Ensure evaporative condesner airflow rate at lower speed must be lower or equal to one at higher speed.
    Do Mode = 1,DXCoil(DXCoilNum)%NumOfSpeeds-1
      If (DXCoil(DXCoilNum)%MSEvapCondAirFlow(Mode) .GT. DXCoil(DXCoilNum)%MSEvapCondAirFlow(Mode+1)) Then
        CALL ShowWarningError('SizeDXCoil: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name)//', '// &
          'Speed '//TRIM(TrimSigDigits(Mode))//' Evaporative Condenser Air Flow Rate must be less than or equal to '//&
          'Speed '//TRIM(TrimSigDigits(Mode+1))//' Evaporative Condenser Air Flow Rate.')
        CALL ShowContinueError('Instead, '//TRIM(RoundSigDigits(DXCoil(DXCoilNum)%MSEvapCondAirFlow(Mode),2))//' > '//  &
                  TRIM(RoundSigDigits(DXCoil(DXCoilNum)%MSEvapCondAirFlow(Mode+1),2)))
        CALL ShowFatalError('Preceding conditions cause termination.')
      End If
    End Do

    ! Rated Evapovative condenser pump power
    Do Mode = 1,DXCoil(DXCoilNum)%NumOfSpeeds
      IF (DXCoil(DXCoilNum)%MSEvapCondPumpElecNomPower(mode) == AutoSize) THEN
        ! Auto size low speed evap condenser pump power to 1/3 Total Capacity * 0.004266 w/w (15 w/ton)
        DXCoil(DXCoilNum)%MSEvapCondPumpElecNomPower(Mode) = DXCoil(DXCoilNum)%MSRatedTotCap(Mode)*0.004266d0
          CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
               'Speed '//TRIM(TrimSigDigits(Mode))//' Rated Evaporative Condenser Pump Power Consumption [W]', &
               DXCoil(DXCoilNum)%MSEvapCondPumpElecNomPower(mode))
      END IF
    End Do
   ! Ensure evaporative condesner pump power at lower speed must be lower or equal to one at higher speed.
    Do Mode = 1,DXCoil(DXCoilNum)%NumOfSpeeds-1
      If (DXCoil(DXCoilNum)%MSEvapCondPumpElecNomPower(Mode) .GT. DXCoil(DXCoilNum)%MSEvapCondPumpElecNomPower(Mode+1)) Then
        CALL ShowWarningError('SizeDXCoil: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name)//', '// &
          'Speed '//TRIM(TrimSigDigits(Mode))//  &
             ' Rated Evaporative Condenser Pump Power Consumption must be less than or equal to '//&
          'Speed '//TRIM(TrimSigDigits(Mode+1))//  &
             ' Rated Evaporative Condenser Pump Power Consumption.')
        CALL ShowContinueError('Instead, '//TRIM(RoundSigDigits(DXCoil(DXCoilNum)%MSEvapCondPumpElecNomPower(Mode),2))//' > '//  &
                  TRIM(RoundSigDigits(DXCoil(DXCoilNum)%MSEvapCondPumpElecNomPower(Mode+1),2)))
        CALL ShowFatalError('Preceding conditions cause termination.')
      End If
    End Do

  End If

  ! Autosizing for multispeed heating coil
  IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_MultiSpeedHeating) THEN
    ! flow rate auto size
    Do Mode = DXCoil(DXCoilNum)%NumOfSpeeds,1,-1
      If (DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode) == AutoSize) Then
        If (Mode == DXCoil(DXCoilNum)%NumOfSpeeds) Then
          IF (CurSysNum > 0) THEN
            CALL CheckSysSizing(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name)
            IF (CurOASysNum > 0) THEN
              DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode) = FinalSysSizing(CurSysNum)%DesOutAirVolFlow
            ELSE
              DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode) = FinalSysSizing(CurSysNum)%DesMainVolFlow
            END IF
          ELSE IF (CurZoneEqNum > 0) THEN
            CALL CheckZoneSizing(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name)
            DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode) = MAX(FinalZoneSizing(CurZoneEqNum)%DesCoolVolFlow, &
                                                      FinalZoneSizing(CurZoneEqNum)%DesHeatVolFlow)
          END IF
          IF (DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) < SmallAirVolFlow) THEN
            DXCoil(DXCoilNum)%RatedAirVolFlowRate(Mode) = 0.0
          END IF
        Else
          DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode) = DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(DXCoil(DXCoilNum)%NumOfSpeeds)* &
                                                          Mode/DXCoil(DXCoilNum)%NumOfSpeeds
        End If
        CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
             'Speed '//Trim(TrimSigDigits(Mode))//' Rated Air Flow Rate [m3/s]', &
              DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode))
      End If
    End Do

   ! Ensure flow rate at lower speed must be lower or equal to the flow rate at higher speed. Otherwise, a severe error is isssued.
    Do Mode = 1,DXCoil(DXCoilNum)%NumOfSpeeds-1
      If (DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode) .GT. DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode+1)) Then
        CALL ShowWarningError('SizeDXCoil: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name)//', '// &
          'Speed '//Trim(TrimSigDigits(Mode))//' Rated Air Flow Rate must be less than or equal to '//&
          'Speed '//Trim(TrimSigDigits(Mode+1))//' Rated Air Flow Rate.')
        CALL ShowContinueError('Instead, '//TRIM(RoundSigDigits(DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode),2))//' > '//  &
                  TRIM(RoundSigDigits(DXCoil(DXCoilNum)%MSRatedAirVolFlowRate(Mode+1),2)))
        CALL ShowFatalError('Preceding conditions cause termination.')
      End If
    End Do

    Do Mode = DXCoil(DXCoilNum)%NumOfSpeeds,1,-1
      If (DXCoil(DXCoilNum)%MSRatedTotCap(Mode) == AutoSize) Then
        IF (Mode .eq. DXCoil(DXCoilNum)%NumOfSpeeds) Then
          ! Heating capacity is assumed to be equal to the cooling capacity
          NumOfSpeedCompanion = DXCoil(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil)%NumOfSpeeds
          DXCoil(DXCoilNum)%MSRatedTotCap(Mode) = &
            DXCoil(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil)%MSRatedTotCap(NumOfSpeedCompanion)
        Else
          DXCoil(DXCoilNum)%MSRatedTotCap(Mode) = DXCoil(DXCoilNum)%MSRatedTotCap(DXCoil(DXCoilNum)%NumOfSpeeds)* &
                                                         Mode/DXCoil(DXCoilNum)%NumOfSpeeds
        End If
        CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
           'Speed '//Trim(TrimSigDigits(Mode))//' Rated Total Heating Capacity [W]', DXCoil(DXCoilNum)%MSRatedTotCap(Mode))
      End If
    End Do

   ! Ensure capacity at lower speed must be lower or equal to the capacity at higher speed.
    Do Mode = 1,DXCoil(DXCoilNum)%NumOfSpeeds-1
      If (DXCoil(DXCoilNum)%MSRatedTotCap(Mode) .GT. DXCoil(DXCoilNum)%MSRatedTotCap(Mode+1)) Then
        CALL ShowWarningError('SizeDXCoil: '//TRIM(DXCoil(DXCoilNum)%DXCoilType)//' '//TRIM(DXCoil(DXCoilNum)%Name)//', '// &
          'Speed '//Trim(TrimSigDigits(Mode))//' Rated Total Heating Capacity must be less than or equal to '//&
          'Speed '//Trim(TrimSigDigits(Mode+1))//' Rated Total Heating Capacity.')
        CALL ShowContinueError('Instead, '//TRIM(RoundSigDigits(DXCoil(DXCoilNum)%MSRatedTotCap(Mode),2))//' > '//  &
                  TRIM(RoundSigDigits(DXCoil(DXCoilNum)%MSRatedTotCap(Mode+1),2)))
        CALL ShowFatalError('Preceding conditions cause termination.')
      End If
    End Do

    ! Resistive Defrost Heater Capacity = capacity at the first stage
    IF (DXCoil(DXCoilNum)%DefrostCapacity == AutoSize) THEN
      IF (DXCoil(DXCoilNum)%DefrostStrategy == Resistive) THEN
        DXCoil(DXCoilNum)%DefrostCapacity = DXCoil(DXCoilNum)%MSRatedTotCap(1)
      ELSE
        DXCoil(DXCoilNum)%DefrostCapacity = 0.0
      END IF
        CALL ReportSizingOutput(DXCoil(DXCoilNum)%DXCoilType, DXCoil(DXCoilNum)%Name, &
                                'Resistive Defrost Heater Capacity', DXCoil(DXCoilNum)%DefrostCapacity)
    END IF
  End If

  ! Call routine that computes AHRI certified rating for single-speed DX Cooling Coils
  IF ((DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingSingleSpeed .AND. DXCoil(DXCoilNum)%CondenserType(1) == AirCooled) &
       .OR. DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_HeatingEmpirical) THEN
      CALL CalcDXCoilStandardRating(DXCoil(DXCoilNum)%Name,                          &
                                    DXCoil(DXCoilNum)%DXCoilType,                    &
                                    DXCoil(DXCoilNum)%DXCoilType_Num,                &
                                    DXCoil(DXCoilNum)%RatedTotCap(1),                &
                                    DXCoil(DXCoilNum)%RatedCOP(1),                   &
                                    DXCoil(DXCoilNum)%CCapFFlow(1),                  &
                                    DXCoil(DXCoilNum)%CCapFTemp(1),                  &
                                    DXCoil(DXCoilNum)%EIRFFlow(1),                   &
                                    DXCoil(DXCoilNum)%EIRFTemp(1),                   &
                                    DXCoil(DXCoilNum)%PLFFPLR(1),                    &
                                    DXCoil(DXCoilNum)%RatedAirVolFlowRate(1),        &
                                    DXCoil(DXCoilNum)%FanPowerPerEvapAirFlowRate(1), &
                                    DXCoil(DXCoilNum)%RegionNum,                     &
                                    DXCoil(DXCoilNum)%MinOATCompressor,              &
                                    DXCoil(DXCoilNum)%OATempCompressorOn,            &
                                    DXCoil(DXCoilNum)%OATempCompressorOnOffBlank,    &
                                    DXCoil(DXCoilNum)%DefrostControl)
  END IF

  !create predefined report entries
  equipName = DXCoil(DXCoilNum)%Name
  ! put tables for cooling and heating separate
  SELECT CASE (DXCoil(DXCoilNum)%DXCoilType_Num)
    CASE (CoilDX_CoolingSingleSpeed, &
          CoilDX_CoolingTwoSpeed, &
          CoilDX_CoolingTwoStageWHumControl, &
          CoilDX_MultiSpeedCooling)
      CALL PreDefTableEntry(pdchCoolCoilType,equipName,DXCoil(DXCoilNum)%DXCoilType)
      IF (DXCoil(DXCoilNum)%NumOfSpeeds .EQ. 0) Then
        IF (DXCoil(DXCoilNum)%NumCapacityStages .EQ. 1) THEN
          CALL PreDefTableEntry(pdchCoolCoilTotCap,equipName,DXCoil(DXCoilNum)%RatedTotCap(1))
          CALL PreDefTableEntry(pdchCoolCoilSensCap,equipName,DXCoil(DXCoilNum)%RatedTotCap(1) &
                    * DXCoil(DXCoilNum)%RatedSHR(1))
          CALL PreDefTableEntry(pdchCoolCoilLatCap,equipName,DXCoil(DXCoilNum)%RatedTotCap(1) &
                    - DXCoil(DXCoilNum)%RatedTotCap(1) * DXCoil(DXCoilNum)%RatedSHR(1))
          CALL PreDefTableEntry(pdchCoolCoilSHR,equipName,DXCoil(DXCoilNum)%RatedSHR(1))
          CALL PreDefTableEntry(pdchCoolCoilNomEff,equipName,DXCoil(DXCoilNum)%RatedCOP(1))
        ELSE
          CALL PreDefTableEntry(pdchCoolCoilTotCap,equipName,DXCoil(DXCoilNum)%RatedTotCap(2))
          CALL PreDefTableEntry(pdchCoolCoilSensCap,equipName,DXCoil(DXCoilNum)%RatedTotCap(2) &
                    * DXCoil(DXCoilNum)%RatedSHR(2))
          CALL PreDefTableEntry(pdchCoolCoilLatCap,equipName,DXCoil(DXCoilNum)%RatedTotCap(2) &
                    - DXCoil(DXCoilNum)%RatedTotCap(2) * DXCoil(DXCoilNum)%RatedSHR(2))
          CALL PreDefTableEntry(pdchCoolCoilSHR,equipName,DXCoil(DXCoilNum)%RatedSHR(2))
          CALL PreDefTableEntry(pdchCoolCoilNomEff,equipName,DXCoil(DXCoilNum)%RatedCOP(2))
        END IF
      ELSE
        DO Mode=1,DXCoil(DXCoilNum)%NumOfSpeeds
          CALL PreDefTableEntry(pdchCoolCoilTotCap,equipName,DXCoil(DXCoilNum)%MSRatedTotCap(Mode))
          CALL PreDefTableEntry(pdchCoolCoilSensCap,equipName,DXCoil(DXCoilNum)%MSRatedTotCap(Mode) &
                    * DXCoil(DXCoilNum)%MSRatedSHR(Mode))
          CALL PreDefTableEntry(pdchCoolCoilLatCap,equipName,DXCoil(DXCoilNum)%MSRatedTotCap(Mode) &
                    - DXCoil(DXCoilNum)%MSRatedTotCap(Mode) * DXCoil(DXCoilNum)%MSRatedSHR(Mode))
          CALL PreDefTableEntry(pdchCoolCoilSHR,equipName,DXCoil(DXCoilNum)%MSRatedSHR(Mode))
          CALL PreDefTableEntry(pdchCoolCoilNomEff,equipName,DXCoil(DXCoilNum)%MSRatedCOP(Mode))
        END DO
      END IF
    CASE (CoilDX_HeatingEmpirical, &
          CoilDX_MultiSpeedHeating, &
          CoilDX_HeatPumpWaterHeater)
      CALL PreDefTableEntry(pdchHeatCoilType,equipName,DXCoil(DXCoilNum)%DXCoilType)
      IF (DXCoil(DXCoilNum)%NumOfSpeeds .EQ. 0) Then
        IF (DXCoil(DXCoilNum)%NumCapacityStages .EQ. 1) THEN
          CALL PreDefTableEntry(pdchHeatCoilNomCap,equipName,DXCoil(DXCoilNum)%RatedTotCap(1))
          CALL PreDefTableEntry(pdchHeatCoilNomEff,equipName,DXCoil(DXCoilNum)%RatedCOP(1))
        ELSE
          CALL PreDefTableEntry(pdchHeatCoilNomCap,equipName,DXCoil(DXCoilNum)%RatedTotCap(2))
          CALL PreDefTableEntry(pdchHeatCoilNomEff,equipName,DXCoil(DXCoilNum)%RatedCOP(2))
        END IF
      ELSE
        DO Mode=1,DXCoil(DXCoilNum)%NumOfSpeeds
          CALL PreDefTableEntry(pdchHeatCoilNomCap,equipName,DXCoil(DXCoilNum)%MSRatedTotCap(Mode))
          CALL PreDefTableEntry(pdchHeatCoilNomEff,equipName,DXCoil(DXCoilNum)%MSRatedCOP(Mode))
        END DO
      END IF
  END SELECT
  RETURN

END SUBROUTINE SizeDXCoil

SUBROUTINE CalcHPWHDXCoil(DXCoilNum,PartLoadRatio)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Richard Raustad
          !       DATE WRITTEN   May 2005
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! Calculates the gross cooling capacity of a heat pump water heater evaporator and
          ! heating capacity of the condenser coil given the rated heating capacity and COP.

          ! METHODOLOGY EMPLOYED:
          ! The routine requires the user to enter the total heating capacity and COP for the
          ! heat pump water heater along with logicals defining if fan and condenser pump are included.
          ! Since manufacturer's can rate their HPWH equipment with or without including condenser
          ! pump heat, this information is required to accurately determine the condenser's leaving
          ! water temperature. In addition, knowledge of the fan heat is required to back into
          ! a compressor COP.

          ! USE STATEMENTS:
  USE CurveManager,    ONLY: CurveValue
  USE General,         ONLY: TrimSigDigits
  USE DataHVACGlobals, ONLY: FanElecPower, HPWHInletDBTemp, HPWHInletWBTemp, DXCoilTotalCapacity

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  INTEGER, INTENT(IN) :: DXCoilNum        ! the number of the DX coil to be simulated
  REAL(r64),    INTENT(IN) :: PartLoadRatio    ! sensible water heating load / full load sensible water heating capacity

          ! SUBROUTINE PARAMETER DEFINITIONS:
  CHARACTER(len=*), PARAMETER ::  RoutineName='CalcHPWHDXCoil'

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
REAL(r64) :: RatedHeatingCapacity       ! Water heating rated capacity with or without condenser water pump heat (W)
REAL(r64) :: RatedHeatingCOP            ! Water heating rated COP with or without evap fan and cond water pump heat (W/W)
REAL(r64) :: OperatingHeatingCapacity   ! Water heating operating capacity including the impact of capacity and COP curves (W)
REAL(r64) :: OperatingHeatingCOP        ! Water heating operating COP including the impact of capacity and COP curves (W/W)
REAL(r64) :: OperatingHeatingPower      ! Water heating operating Power (W)
REAL(r64) :: CompressorPower            ! Power consumed by compressor only (W)

REAL(r64) :: TotalTankHeatingCapacity   ! Water heating capacity corrected for condenser water pump heat (W)
REAL(r64) :: TankHeatingCOP             ! Water heating COP corrected for fan and condenser water pump power (W/W)
                                      ! (these previous 2 variables also include the impact of capacity and COP curves)
REAL(r64) :: EvapCoolingCapacity        ! Air cooling capacity corrected for evap fan and cond water pump heat (W)
REAL(r64) :: InletWaterTemp             ! Condenser water inlet temperature (C)
REAL(r64) :: OutletWaterTemp            ! Condenser water outlet temperature (C)
REAL(r64) :: EvapInletMassFlowRate      ! Evaporator air inlet mass flow rate (m3/s)
REAL(r64) :: CondInletMassFlowRate      ! Condenser water inlet mass flow rate (m3/s)
REAL(r64) :: CpWater                    ! Specific heat of condenser inlet water (J/Kg/k)
REAL(r64) :: InletAirTemp               ! HPWH inlet air temperature (dry-bulb or wet-bulb) (C)
REAL(r64) :: HeatCapFTemp               ! Output of HPWH Heating Capacity as a Function of Temperature curve
REAL(r64) :: HeatCapFAirFlow            ! Output of HPWH Heating Capacity as a Function of Air Flow Rate Ratio curve
REAL(r64) :: HeatCapFWaterFlow          ! Output of HPWH Heating Capacity as a Function of Water Flow Rate Ratio curve
REAL(r64) :: HeatCOPFTemp               ! Output of HPWH COP as a Function of Temperature curve
REAL(r64) :: HeatCOPFAirFlow            ! Output of HPWH COP as a Function of Air Flow Rate Ratio curve
REAL(r64) :: HeatCOPFWAterFlow          ! Output of HPWH COP as a Function of Water Flow Rate Ratio curve
REAL(r64) :: AirFlowRateRatio           ! Ratio of evaporator inlet air mass flow rate to rated mass flow rate
REAL(r64) :: WaterFlowRateRatio         ! Ratio of evaporator inlet water mass flow rate to rated mass flow rate
REAL(r64) :: PartLoadFraction             ! Output of Part Load Fraction as a Function of Part Load Ratio curve
REAL(r64) :: PumpHeatToWater            ! Amount of pump heat attributed to heating water
REAL(r64) :: HPRTF                      ! Heat pump run time fraction
INTEGER :: EvapInletNode              ! Evaporator air inlet node number
INTEGER :: EvapOutletNode             ! Evaporator air outlet node number
INTEGER :: CondInletNode              ! Condenser water inlet node number
INTEGER :: CondOutletNode             ! Condenser water outlet node number

  CondInletNode  = DXCoil(DXCoilNum)%WaterInNode
  CondOutletNode = DXCoil(DXCoilNum)%WaterOutNode

! If heat pump water heater is OFF, set outlet to inlet and RETURN
  IF(PartLoadRatio .EQ. 0.)THEN
    Node(CondOutletNode) = Node(CondInletNode)
    RETURN
  ELSE
    RatedHeatingCapacity   = DXCoil(DXCoilNum)%RatedTotCap2
    RatedHeatingCOP        = DXCoil(DXCoilNum)%RatedCOP(1)
    EvapInletNode          = DXCoil(DXCoilNum)%AirInNode
    EvapOutletNode         = DXCoil(DXCoilNum)%AirOutNode
    InletWaterTemp         = Node(CondInletNode)%Temp
    CondInletMassFlowRate  = Node(CondInletNode)%MassFlowRate / PartLoadRatio
    EvapInletMassFlowRate  = Node(EvapInletNode)%MassFlowRate / PartLoadRatio
    CpWater                = CPHW(InletWaterTemp)
    CompressorPower        = 0.0
    OperatingHeatingPower  = 0.0
    TankHeatingCOP         = 0.0
  END IF

! determine inlet air temperature type for curve objects
  IF(DXCoil(DXCoilNum)%InletAirTemperatureType .EQ. WetBulbIndicator) THEN
    InletAirTemp = HPWHInletWBTemp
  ELSE
    InletAirTemp = HPWHInletDBTemp
  END IF

! get output of Heating Capacity and Heating COP curves (curves default to 1 if user has not specified curve name)
  IF(DXCoil(DXCoilNum)%HCapFTemp .GT. 0)THEN
    IF(DXCoil(DXCoilNum)%HCapFTempCurveType .EQ. Cubic) THEN
      HeatCapFTemp = CurveValue(DXCoil(DXCoilNum)%HCapFTemp,InletAirTemp)
    ELSE
      HeatCapFTemp = CurveValue(DXCoil(DXCoilNum)%HCapFTemp,InletAirTemp,InletWaterTemp)
    END IF
!   Warn user if curve output goes negative
    IF(HeatCapFTemp .LT. 0.0)THEN
      IF(DXCoil(DXCoilNum)%HCapFTempErrorIndex == 0)THEN
        CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//'":')
        CALL ShowContinueError(' HPWH Heating Capacity Modifier curve (function of temperature) output is negative (' &
                           //TRIM(TrimSigDigits(HeatCapFTemp,3))//').')
        IF(DXCoil(DXCoilNum)%HCapFTempCurveType .EQ. Biquadratic) THEN
          CALL ShowContinueError(' Negative value occurs using an inlet air temperature of ' &
                            //TRIM(TrimSigDigits(InletAirTemp,1))// &
                               ' and an inlet water temperature of '//TRIM(TrimSigDigits(InletWaterTemp,1))//'.')
        ELSE
          CALL ShowContinueError(' Negative value occurs using an inlet air temperature of ' &
                            //TRIM(TrimSigDigits(InletAirTemp,1))//'.')
        END IF
        CALL ShowContinueErrorTimeStamp(' Resetting curve output to zero and continuing simulation.')
      END IF
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//'":'//&
          ' HPWH Heating Capacity Modifier curve (function of temperature) output is negative warning continues...' &
          , DXCoil(DXCoilNum)%HCapFTempErrorIndex, HeatCapFTemp, HeatCapFTemp,  &
            ReportMinUnits='[C]',ReportMaxUnits='[C]')
      HeatCapFTemp = 0.0
    END IF
  ELSE
    HeatCapFTemp = 1.0d0
  END IF

  IF(DXCoil(DXCoilNum)%HCOPFTemp .GT. 0)THEN
    IF(DXCoil(DXCoilNum)%HCOPFTempCurveType .EQ. Cubic) THEN
      HeatCOPFTemp = CurveValue(DXCoil(DXCoilNum)%HCOPFTemp,InletAirTemp)
    ELSE
      HeatCOPFTemp = CurveValue(DXCoil(DXCoilNum)%HCOPFTemp,InletAirTemp,InletWaterTemp)
    END IF
!   Warn user if curve output goes negative
    IF(HeatCOPFTemp .LT. 0.0)THEN
      IF(DXCoil(DXCoilNum)%HCOPFTempErrorIndex == 0)THEN
        CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//'":')
        CALL ShowContinueError(' HPWH Heating COP Modifier curve (function of temperature) output is negative (' &
                           //TRIM(TrimSigDigits(HeatCOPFTemp,3))//').')
        IF(DXCoil(DXCoilNum)%HCOPFTempCurveType .EQ. Biquadratic) THEN
          CALL ShowContinueError(' Negative value occurs using an inlet air temperature of ' &
                            //TRIM(TrimSigDigits(InletAirTemp,1))// &
                               ' and an inlet water temperature of '//TRIM(TrimSigDigits(InletWaterTemp,1))//'.')
        ELSE
          CALL ShowContinueError(' Negative value occurs using an inlet air temperature of ' &
                            //TRIM(TrimSigDigits(InletAirTemp,1))//'.')
        END IF
        CALL ShowContinueErrorTimeStamp(' Resetting curve output to zero and continuing simulation.')
      END IF
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//'":'//&
          ' HPWH Heating COP Modifier curve (function of temperature) output is negative warning continues...' &
          , DXCoil(DXCoilNum)%HCOPFTempErrorIndex, HeatCOPFTemp, HeatCOPFTemp,  &
            ReportMinUnits='[C]',ReportMaxUnits='[C]')
      HeatCOPFTemp = 0.0
    END IF
  ELSE
    HeatCOPFTemp = 1.0d0
  END IF

  IF(DXCoil(DXCoilNum)%HCapFAirFlow .GT. 0) THEN
    AirFlowRateRatio = EvapInletMassFlowRate / (DXCoil(DXCoilNum)%RatedAirMassFlowRate(1))
    HeatCapFAirFlow  = CurveValue(DXCoil(DXCoilNum)%HCapFAirFlow,AirFlowRateRatio)
!   Warn user if curve output goes negative
    IF(HeatCapFAirFlow .LT. 0.0)THEN
      IF(DXCoil(DXCoilNum)%HCapFAirFlowErrorIndex == 0)THEN
        CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//'":')
        CALL ShowContinueError(' HPWH Heating Capacity Modifier curve (function of air flow fraction) output is negative (' &
                           //TRIM(TrimSigDigits(HeatCapFAirFlow,3))//').')
        CALL ShowContinueError(' Negative value occurs using an air flow fraction of ' &
                            //TRIM(TrimSigDigits(AirFlowRateRatio,3))//'.')
        CALL ShowContinueErrorTimeStamp(' Resetting curve output to zero and continuing simulation.')
      END IF
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//'":'//&
          ' HPWH Heating Capacity Modifier curve (function of air flow fraction) output is negative warning continues...' &
          , DXCoil(DXCoilNum)%HCapFAirFlowErrorIndex, HeatCapFAirFlow, HeatCapFAirFlow)
      HeatCapFAirFlow = 0.0
    END IF
  ELSE
    HeatCapFAirFlow  = 1.0d0
  END IF

  IF(DXCoil(DXCoilNum)%HCOPFAirFlow .GT. 0) THEN
    AirFlowRateRatio = EvapInletMassFlowRate / (DXCoil(DXCoilNum)%RatedAirMassFlowRate(1))
    HeatCOPFAirFlow  = CurveValue(DXCoil(DXCoilNum)%HCOPFAirFlow,AirFlowRateRatio)
!   Warn user if curve output goes negative
    IF(HeatCOPFAirFlow .LT. 0.0)THEN
      IF(DXCoil(DXCoilNum)%HCOPFAirFlowErrorIndex == 0)THEN
        CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//'":')
        CALL ShowContinueError(' HPWH Heating COP Modifier curve (function of air flow fraction) output is negative (' &
                           //TRIM(TrimSigDigits(HeatCOPFAirFlow,3))//').')
        CALL ShowContinueError(' Negative value occurs using an air flow fraction of ' &
                            //TRIM(TrimSigDigits(AirFlowRateRatio,3))//'.')
        CALL ShowContinueErrorTimeStamp(' Resetting curve output to zero and continuing simulation.')
      END IF
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//'":'//&
          ' HPWH Heating COP Modifier curve (function of air flow fraction) output is negative warning continues...' &
          , DXCoil(DXCoilNum)%HCOPFAirFlowErrorIndex, HeatCOPFAirFlow, HeatCOPFAirFlow)
      HeatCOPFAirFlow = 0.0
    END IF
  ELSE
    HeatCOPFAirFlow  = 1.0d0
  END IF

  IF(DXCoil(DXCoilNum)%HCapFWaterFlow .GT. 0) THEN
    WaterFlowRateRatio = CondInletMassFlowRate / (DXCoil(DXCoilNum)%RatedHPWHCondWaterFlow*RhoH2O(InletWaterTemp))
    HeatCapFWaterFlow  = CurveValue(DXCoil(DXCoilNum)%HCapFWaterFlow,WaterFlowRateRatio)
!   Warn user if curve output goes negative
    IF(HeatCapFWaterFlow .LT. 0.0)THEN
      IF(DXCoil(DXCoilNum)%HCapFWaterFlowErrorIndex == 0)THEN
        CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//'":')
        CALL ShowContinueError(' HPWH Heating Capacity Modifier curve (function of water flow fraction) output is negative (' &
                           //TRIM(TrimSigDigits(HeatCapFWaterFlow,3))//').')
        CALL ShowContinueError(' Negative value occurs using a water flow fraction of ' &
                            //TRIM(TrimSigDigits(WaterFlowRateRatio,3))//'.')
        CALL ShowContinueErrorTimeStamp(' Resetting curve output to zero and continuing simulation.')
      END IF
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//'":'//&
          ' HPWH Heating Capacity Modifier curve (function of water flow fraction) output is negative warning continues...' &
          , DXCoil(DXCoilNum)%HCapFWaterFlowErrorIndex, HeatCapFWaterFlow, HeatCapFWaterFlow)
      HeatCapFWaterFlow = 0.0
    END IF
  ELSE
    HeatCapFWaterFlow  = 1.0d0
  END IF

  IF(DXCoil(DXCoilNum)%HCOPFWaterFlow .GT. 0) THEN
    WaterFlowRateRatio = CondInletMassFlowRate / (DXCoil(DXCoilNum)%RatedHPWHCondWaterFlow*RhoH2O(InletWaterTemp))
    HeatCOPFWaterFlow  = CurveValue(DXCoil(DXCoilNum)%HCOPFWaterFlow,WaterFlowRateRatio)
!   Warn user if curve output goes negative
    IF(HeatCOPFWaterFlow .LT. 0.0)THEN
      IF(DXCoil(DXCoilNum)%HCOPFWaterFlowErrorIndex == 0)THEN
        CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//'":')
        CALL ShowContinueError(' HPWH Heating COP Modifier curve (function of water flow fraction) output is negative (' &
                           //TRIM(TrimSigDigits(HeatCOPFWaterFlow,3))//').')
        CALL ShowContinueError(' Negative value occurs using a water flow fraction of ' &
                            //TRIM(TrimSigDigits(WaterFlowRateRatio,3))//'.')
        CALL ShowContinueErrorTimeStamp(' Resetting curve output to zero and continuing simulation.')
      END IF
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//'":'//&
          ' HPWH Heating COP Modifier curve (function of water flow fraction) output is negative warning continues...' &
          , DXCoil(DXCoilNum)%HCOPFWaterFlowErrorIndex, HeatCOPFWaterFlow, HeatCOPFWaterFlow)
      HeatCOPFWaterFlow = 0.0
    END IF
  ELSE
    HeatCOPFWaterFlow  = 1.0d0
  END IF

! adjust Heating Capacity and COP for off-design conditions
  OperatingHeatingCapacity = RatedHeatingCapacity * HeatCapFTemp * HeatCapFAirFlow * HeatCapFWaterFlow
  OperatingHeatingCOP      = RatedHeatingCOP      * HeatCOPFTemp * HeatCOPFAirFlow * HeatCOPFWaterFlow

  IF(OperatingHeatingCOP .GT. 0.0) OperatingHeatingPower = OperatingHeatingCapacity / OperatingHeatingCOP

  PumpHeatToWater          = DXCoil(DXCoilNum)%HPWHCondPumpElecNomPower * DXCoil(DXCoilNum)%HPWHCondPumpFracToWater
  TankHeatingCOP           = OperatingHeatingCOP

! account for pump heat if not included in total water heating capacity
  IF(DXCoil(DXCoilNum)%CondPumpHeatInCapacity)THEN
    TotalTankHeatingCapacity = OperatingHeatingCapacity
  ELSE
    TotalTankHeatingCapacity = OperatingHeatingCapacity + PumpHeatToWater
  END IF

! find part load fraction to calculate RTF
  IF(DXCoil(DXCoilNum)%PLFFPLR(1) .GT. 0) THEN
    PartLoadFraction = MAX(0.7d0,CurveValue(DXCoil(DXCoilNum)%PLFFPLR(1),PartLoadRatio))
  ELSE
    PartLoadFraction = 1.0d0
  END IF

  HPRTF = MIN(1.0d0,(PartLoadRatio / PartLoadFraction))

! calculate evaporator total cooling capacity
  IF(HPRTF .GT. 0.0)THEN
    IF(DXCoil(DXCoilNum)%FanPowerIncludedInCOP)THEN
      IF(DXCoil(DXCoilNum)%CondPumpPowerInCOP)THEN
!       make sure fan power is full load fan power
        CompressorPower = OperatingHeatingPower - FanElecPower/HPRTF - DXCoil(DXCoilNum)%HPWHCondPumpElecNomPower
        IF(OperatingHeatingPower .GT. 0.0)TankHeatingCOP = TotalTankHeatingCapacity / OperatingHeatingPower
      ELSE
        CompressorPower = OperatingHeatingPower - FanElecPower/HPRTF
        IF((OperatingHeatingPower + DXCoil(DXCoilNum)%HPWHCondPumpElecNomPower) .GT. 0.0) &
          TankHeatingCOP = TotalTankHeatingCapacity / (OperatingHeatingPower + DXCoil(DXCoilNum)%HPWHCondPumpElecNomPower)
      END IF
    ELSE
      IF(DXCoil(DXCoilNum)%CondPumpPowerInCOP)THEN
!       make sure fan power is full load fan power
        CompressorPower = OperatingHeatingPower - DXCoil(DXCoilNum)%HPWHCondPumpElecNomPower
        IF((OperatingHeatingPower + FanElecPower/HPRTF) .GT. 0.0) &
          TankHeatingCOP = TotalTankHeatingCapacity / (OperatingHeatingPower + FanElecPower/HPRTF)
      ELSE
        CompressorPower = OperatingHeatingPower
        IF((OperatingHeatingPower + FanElecPower/HPRTF + DXCoil(DXCoilNum)%HPWHCondPumpElecNomPower) .GT. 0.0) &
          TankHeatingCOP = TotalTankHeatingCapacity / &
                          (OperatingHeatingPower + FanElecPower/HPRTF + DXCoil(DXCoilNum)%HPWHCondPumpElecNomPower)
      END IF
    END IF
  END IF

  IF(DXCoil(DXCoilNum)%CondPumpHeatInCapacity)THEN
    EvapCoolingCapacity = TotalTankHeatingCapacity - PumpHeatToWater - CompressorPower
  ELSE
    EvapCoolingCapacity = TotalTankHeatingCapacity - CompressorPower
  END IF

! set evaporator total cooling capacity prior to CalcDOE2DXCoil subroutine
  DXCoil(DXCoilNum)%RatedTotCap(1) = EvapCoolingCapacity

! determine condenser water inlet/outlet condition at full capacity
  IF(CondInletMassFlowRate .EQ. 0.0)THEN
    OutletWaterTemp         = InletWaterTemp
  ELSE
    OutletWaterTemp         = InletWaterTemp + TotalTankHeatingCapacity/(CpWater * CondInletMassFlowRate)
  END IF

  Node(CondOutletNode)%Temp = OutletWaterTemp

  Node(CondOutletNode)%MassFlowRate = Node(CondInletNode)%MassFlowRate

! send heating capacity and COP to water heater module for standards rating calculation
! total heating capacity including condenser pump
  HPWHHeatingCapacity = TotalTankHeatingCapacity
! total heating COP including compressor, fan, and condenser pump
  HPWHHeatingCOP      = TankHeatingCOP

! send DX coil total cooling capacity to HPWH for reporting
  DXCoilTotalCapacity = EvapCoolingCapacity

  DXCoil(DXCoilNum)%TotalHeatingEnergyRate = TotalTankHeatingCapacity * PartLoadRatio

! calculate total compressor plus condenser pump power, fan power reported in fan module
  DXCoil(DXCoilNum)%ElecWaterHeatingPower = (CompressorPower + DXCoil(DXCoilNum)%HPWHCondPumpElecNomPower) * HPRTF

RETURN
END SUBROUTINE CalcHPWHDXCoil

SUBROUTINE CalcDoe2DXCoil(DXCoilNum,CompOp,FirstHVACIteration,PartLoadRatio,FanOpMode,PerfMode,OnOffAirFlowRatio, &
                          CoolingHeatingPLR)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Fred Buhl
          !       DATE WRITTEN   May 2000
          !       MODIFIED       Shirey, Feb/October 2001, Feb/Mar 2004
          !                      Feb 2005 M. J. Witte, GARD Analytics, Inc.
          !                        Add new coil type COIL:DX:MultiMode:CoolingEmpirical:
          !                      April 2010 Chandan Sharma, FSEC, Added basin heater
          !       RE-ENGINEERED  Don Shirey, Aug/Sept 2000

          ! PURPOSE OF THIS SUBROUTINE:
          ! Calculates the air-side performance and electrical energy use of a direct-
          ! expansion, air-cooled cooling unit.

          ! METHODOLOGY EMPLOYED:
          ! This routine simulates the performance of air-cooled DX cooling equipment.
          ! The routine requires the user to enter the total cooling capacity, sensible heat ratio,
          ! and COP for the unit at ARI 210/240 rating conditions (26.67C [80F] dry-bulb, 19.44C [67F]
          ! wet-bulb air entering the cooling coil, 35C [95F] dry-bulb air entering the outdoor
          ! condenser. Since different manufacturer's rate their equipment at different air flow rates,
          ! the supply air flow rate corresponding to the rated capacities and rated COP must also be
          ! entered (should be between 300 cfm/ton and 450 cfm/ton). The rated information entered by
          ! the user should NOT include the thermal or electrical impacts of the supply air fan, as
          ! this is addressed by another module.

          ! With the rated performance data entered by the user, the model employs some of the
          ! DOE-2.1E curve fits to adjust the capacity and efficiency of the unit as a function
          ! of entering air temperatures and supply air flow rate (actual vs rated flow). The model
          ! does NOT employ the exact same methodology to calculate performance as DOE-2, although
          ! some of the DOE-2 curve fits are employed by this model.

          ! The model checks for coil dryout conditions, and adjusts the calculated performance
          ! appropriately.

          ! REFERENCES:
          ! ASHRAE HVAC 2 Toolkit page 4-81.
          !
          ! Henderson, H.I. Jr., K. Rengarajan and D.B. Shirey, III. 1992.The impact of comfort
          ! control on air conditioner energy use in humid climates. ASHRAE Transactions 98(2):
          ! 104-113.
          !
          ! Henderson, H.I. Jr., Danny Parker and Y.J. Huang. 2000.Improving DOE-2's RESYS routine:
          ! User Defined Functions to Provide More Accurate Part Load Energy Use and Humidity
          ! Predictions. Proceedings of ACEEE Conference.


          ! USE STATEMENTS:
  USE CurveManager,    ONLY: CurveValue
  USE DataGlobals,     ONLY: CurrentTime
  USE DataHVACGlobals, ONLY: HPWHCrankcaseDBTemp, TimeStepSys, SysTimeElapsed
  USE General,         ONLY: TrimSigDigits, RoundSigDigits, CreateSysTimeIntervalString
  USE DataWater,       ONLY: WaterStorage

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  INTEGER,   INTENT(IN)           :: DXCoilNum           ! the number of the DX coil to be simulated
  INTEGER,   INTENT(IN)           :: CompOp              ! compressor operation; 1=on, 0=off
  LOGICAL,   INTENT(IN)           :: FirstHVACIteration  ! true if this is the first iteration of HVAC
  REAL(r64), INTENT(IN)           :: PartLoadRatio       ! sensible cooling load / full load sensible cooling capacity
  INTEGER,   INTENT(IN)           :: FanOpMode           ! Allows parent object to control fan operation
  INTEGER,   INTENT(IN), OPTIONAL :: PerfMode            ! Performance mode for MultiMode DX coil; Always 1 for other coil types
  REAL(r64), INTENT(IN), OPTIONAL :: OnOffAirFlowRatio   ! ratio of compressor on airflow to compressor off airflow
  REAL(r64), INTENT(IN), OPTIONAL :: CoolingHeatingPLR   ! used for cycling fan RH control

          ! SUBROUTINE PARAMETER DEFINITIONS:
CHARACTER(len=*), PARAMETER :: RoutineName='CalcDoe2DXCoil'

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
REAL(r64) :: AirMassFlow           ! dry air mass flow rate through coil [kg/s] (adjusted for bypass if any)
REAL(r64) :: AirMassFlowRatio      ! Ratio of actual air mass flow to rated air mass flow (adjusted for bypass if any)
REAL(r64) :: AirVolumeFlowRate     ! Air volume flow rate across the cooling coil [m3/s] (adjusted for bypass if any)
                                   ! (average flow if cycling fan, full flow if constant fan)
REAL(r64) :: VolFlowperRatedTotCap ! Air volume flow rate divided by rated total cooling capacity [m3/s-W] (adjusted for bypass)
REAL(r64) :: BypassFlowFraction    ! Fraction of total flow which is bypassed around the cooling coil
REAL(r64) :: TotCap                ! gross total cooling capacity at off-rated conditions [W]
REAL(r64) :: TotCapTempModFac      ! Total capacity modifier (function of entering wetbulb, outside drybulb)
REAL(r64) :: TotCapFlowModFac      ! Total capacity modifier (function of actual supply air flow vs rated flow)
REAL(r64) :: InletAirWetBulbC      ! wetbulb temperature of inlet air [C]
REAL(r64) :: InletAirDryBulbTemp   ! inlet air dry bulb temperature [C]
REAL(r64) :: InletAirEnthalpy      ! inlet air enthalpy [J/kg]
REAL(r64) :: InletAirHumRat        ! inlet air humidity ratio [kg/kg]
REAL(r64) :: InletAirHumRatTemp    ! inlet air humidity ratio used in ADP/BF loop [kg/kg]
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!REAL(r64) :: InletAirPressure      ! inlet air pressure [Pa]
REAL(r64) :: RatedCBF              ! coil bypass factor at rated conditions
REAL(r64) :: SHR                   ! Sensible Heat Ratio (sensible/total) of the cooling coil
REAL(r64) :: CBF                   ! coil bypass factor at off rated conditions
REAL(r64) :: A0                    ! NTU * air mass flow rate, used in CBF calculation
REAL(r64) :: hDelta                ! Change in air enthalpy across the cooling coil [J/kg]
REAL(r64) :: hADP                  ! Apparatus dew point enthalpy [J/kg]
REAL(r64) :: hTinwADP              ! Enthalpy at inlet dry-bulb and wADP [J/kg]
REAL(r64) :: hTinwout              ! Enthalpy at inlet dry-bulb and outlet humidity ratio [J/kg]
REAL(r64) :: tADP                  ! Apparatus dew point temperature [C]
REAL(r64) :: wADP                  ! Apparatus dew point humidity ratio [kg/kg]
REAL(r64) :: FullLoadOutAirEnth    ! outlet full load enthalpy [J/kg]
REAL(r64) :: FullLoadOutAirHumRat  ! outlet humidity ratio at full load
REAL(r64) :: FullLoadOutAirTemp    ! outlet air temperature at full load [C]
REAL(r64) :: EIRTempModFac         ! EIR modifier (function of entering wetbulb, outside drybulb)
REAL(r64) :: EIRFlowModFac         ! EIR modifier (function of actual supply air flow vs rated flow)
REAL(r64) :: EIR                   ! EIR at part load and off rated conditions
REAL(r64) :: PLF                   ! Part load factor, accounts for thermal lag at compressor startup, used in power calculation
REAL(r64) :: QLatActual            ! operating latent capacity of DX coil
REAL(r64) :: QLatRated             ! Rated latent capacity of DX coil
REAL(r64) :: SHRUnadjusted         ! SHR prior to latent degradation effective SHR calculation
INTEGER :: Counter                 ! Counter for dry evaporator iterations
INTEGER :: MaxIter                 ! Maximum number of iterations for dry evaporator calculations
REAL(r64) :: RF                    ! Relaxation factor for dry evaporator iterations
REAL(r64) :: Tolerance             ! Error tolerance for dry evaporator iterations
REAL(r64) :: werror                ! Deviation of humidity ratio in dry evaporator iteration loop
REAL(r64) :: CondInletTemp         ! Condenser inlet temperature (C). Outdoor dry-bulb temp for air-cooled condenser.
                                 ! Outdoor Wetbulb +(1 - effectiveness)*(outdoor drybulb - outdoor wetbulb) for evap condenser.
REAL(r64) :: CondInletHumrat       ! Condenser inlet humidity ratio (kg/kg). Zero for air-cooled condenser.
                                 ! For evap condenser, its the humidity ratio of the air leaving the evap cooling pads.
REAL(r64) :: CondAirMassFlow       ! Condenser air mass flow rate [kg/s]
REAL(r64) :: RhoAir                ! Density of air [kg/m3]
REAL(r64) :: RhoWater              ! Density of water [kg/m3]
REAL(r64) :: CrankcaseHeatingPower ! power due to crankcase heater
REAL(r64) :: CompAmbTemp = 0.0     ! Ambient temperature at compressor
REAL(r64) :: AirFlowRatio          ! ratio of compressor on airflow to average timestep airflow
                                 ! used when constant fan mode yields different air flow rates when compressor is ON and OFF
                                 ! (e.g. Packaged Terminal Heat Pump)
REAL(r64) :: OutdoorDryBulb        ! Outdoor dry-bulb temperature at condenser (C)
REAL(r64) :: OutdoorWetBulb        ! Outdoor wet-bulb temperature at condenser (C)
REAL(r64) :: OutdoorHumRat         ! Outdoor humidity ratio at condenser (kg/kg)
REAL(r64) :: OutdoorPressure       ! Outdoor barometric pressure at condenser (Pa)

REAL(r64) :: CurrentEndTime = 0.0  ! end time of time step for current simulation time step
REAL(r64) :: MinAirHumRat = 0.0    ! minimum of the inlet air humidity ratio and the outlet air humidity ratio
!CHARACTER(len=6) :: OutputChar = ' '     ! character string for warning messages
!INTEGER,SAVE     :: ErrCount3=0    ! Counter used to minimize the occurrence of output warnings
!INTEGER,SAVE     :: ErrCount4=0    ! Counter used to minimize the occurrence of output warnings
!CHARACTER(len=6) :: CharPLR        ! used in warning messages
!CHARACTER(len=6) :: CharPLF        ! used in warning messages
INTEGER          :: Mode           ! Performance mode for Multimode DX coil; Always 1 for other coil types
!CHARACTER(len=MaxNameLength) :: MinVol      ! character string used for error messages
!CHARACTER(len=MaxNameLength) :: MaxVol      ! character string used for error messages
!CHARACTER(len=MaxNameLength) :: VolFlowChar ! character string used for error messages
REAL(r64) :: OutletAirTemp           ! Supply air temperature (average value if constant fan, full output if cycling fan)
REAL(r64) :: OutletAirHumRat         ! Supply air humidity ratio (average value if constant fan, full output if cycling fan)
REAL(r64) :: OutletAirEnthalpy       ! Supply air enthalpy (average value if constant fan, full output if cycling fan)
REAL(r64) :: Adiff                   ! Used for exponential
REAL(r64) :: DXcoolToHeatPLRRatio    ! ratio of cooling PLR to heating PLR, used for cycling fan RH control
REAL(r64) :: HeatRTF                 ! heating coil part-load ratio, used for cycling fan RH control
REAL(r64) :: HeatingCoilPLF          ! heating coil PLF (function of PLR), used for cycling fan RH control

! If Performance mode not present, then set to 1.  Used only by Multimode/Multispeed DX coil (otherwise mode = 1)
IF (PRESENT(PerfMode)) THEN
  Mode = PerfMode
ELSE
  Mode = 1
END IF

! If AirFlowRatio not present, then set to 1. Used only by DX coils with different air flow
! during cooling and when no cooling is required (constant fan, fan speed changes)
IF (PRESENT(OnOffAirFlowRatio)) THEN
  AirFlowRatio = OnOffAirFlowRatio
ELSE
  AirFlowRatio = 1.0d0
END IF

! If CoolingHeatingPLR not present, then set to 1. Used for cycling fan systems where
! heating PLR is greater than cooling PLR, otherwise CoolingHeatingPLR = 1.
IF(PRESENT(CoolingHeatingPLR))THEN
  DXcoolToHeatPLRRatio = CoolingHeatingPLR
ELSE
  DXcoolToHeatPLRRatio = 1.0d0
END IF

MaxIter         = 30
RF              = 0.4d0
Counter         = 0
Tolerance       = 0.01d0
CondInletTemp   = 0.0d0
CondInletHumrat = 0.0d0
BypassFlowFraction  = DXCoil(DXCoilNum)%BypassedFlowFrac(Mode)
AirMassFlow         = DXCoil(DXCoilNum)%InletAirMassFlowRate * (1.d0-BypassFlowFraction)
InletAirDryBulbTemp = DXCoil(DXCoilNum)%InletAirTemp
InletAirEnthalpy    = DXCoil(DXCoilNum)%InletAirEnthalpy
InletAirHumRat      = DXCoil(DXCoilNum)%InletAirHumRat
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!InletAirPressure    = DXCoil(DXCoilNum)%InletAirPressure
HeatReclaimDXCoil(DXCoilNum)%AvailCapacity   = 0.0d0
DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction = 0.0d0
DXCoil(DXCoilNum)%PartLoadRatio              = 0.0d0
DXCoil(DXCoilNum)%BasinHeaterPower           = 0.0d0

IF (DXCoil(DXCoilNum)%CondenserType(Mode) /= WaterHeater) THEN
  IF (DXCoil(DXCoilNum)%CondenserInletNodeNum(Mode) /= 0) THEN
    OutdoorDryBulb  = Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(Mode))%Temp
    OutdoorHumRat   = Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(Mode))%HumRat
    OutdoorPressure = Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(Mode))%Press
    OutdoorWetBulb  = Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(Mode))%OutAirWetBulb
  ELSE
    OutdoorDryBulb  = OutDryBulbTemp
    OutdoorHumRat   = OutHumRat
    OutdoorPressure = OutBaroPress
    OutdoorWetBulb  = OutWetBulbTemp
  ENDIF
ELSE
  IF (DXCoil(DXCoilNum)%CondenserInletNodeNum(Mode) /= 0) THEN
    OutdoorPressure = Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(Mode))%Press
  ELSE
    OutdoorPressure = OutBaroPress
  ENDIF
END IF

IF (DXCoil(DXCoilNum)%CondenserType(Mode) == AirCooled) THEN
  CondInletTemp   = OutdoorDryBulb ! Outdoor dry-bulb temp
  CompAmbTemp     = OutdoorDryBulb
ELSEIF (DXCoil(DXCoilNum)%CondenserType(Mode) == EvapCooled) THEN
  RhoAir          = PsyRhoAirFnPbTdbW(OutdoorPressure,OutdoorDryBulb,OutdoorHumRat)
  CondAirMassFlow =  RhoAir * DXCoil(DXCoilNum)%EvapCondAirFlow(Mode)
 ! (Outdoor wet-bulb temp from DataEnvironment) + (1.0-EvapCondEffectiveness) * (drybulb - wetbulb)
  CondInletTemp   = OutdoorWetBulb + (OutdoorDryBulb-OutdoorWetBulb)*(1.0d0 - DXCoil(DXCoilNum)%EvapCondEffect(Mode))
  CondInletHumrat = PsyWFnTdbTwbPb(CondInletTemp,OutdoorWetBulb,OutdoorPressure)
  CompAmbTemp     = CondInletTemp
ELSEIF (DXCoil(DXCoilNum)%CondenserType(Mode) == WaterHeater) THEN
  CompAmbTemp     = HPWHCrankcaseDBTemp ! Temperature at HP water heater compressor
END IF

! Initialize crankcase heater, operates below OAT defined in input deck for HP DX cooling coil
! If used in a heat pump, the value of MaxOAT in the heating coil overrides that in the cooling coil (in GetInput)
IF (CompAmbTemp .LT. DXCoil(DXCoilNum)%MaxOATCrankcaseHeater)THEN
  CrankcaseHeatingPower = DXCoil(DXCoilNum)%CrankcaseHeaterCapacity
ELSE
  CrankcaseHeatingPower = 0.0
END IF

! calculate end time of current time step to determine if error messages should be printed
CurrentEndTime = CurrentTime + SysTimeElapsed

!   Print warning messages only when valid and only for the first ocurrance. Let summary provide statistics.
!   Wait for next time step to print warnings. If simulation iterates, print out
!   the warning for the last iteration only. Must wait for next time step to accomplish this.
!   If a warning occurs and the simulation down shifts, the warning is not valid.
IF(DXCoil(DXCoilNum)%PrintLowAmbMessage)THEN ! .AND. &
  IF(CurrentEndTime .GT. DXCoil(DXCoilNum)%CurrentEndTimeLast .AND. &
     TimeStepSys .GE. DXCoil(DXCoilNum)%TimeStepSysLast)THEN
    IF (DXCoil(DXCoilNum)%LowAmbErrIndex == 0) THEN
      CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%LowAmbBuffer1))
      CALL ShowContinueError(TRIM(DXCoil(DXCoilNum)%LowAmbBuffer2))
      CALL ShowContinueError('... Operation at low ambient temperatures may require special performance curves.')
    ENDIF
    IF (DXCoil(DXCoilNum)%CondenserType(Mode) .EQ. AirCooled) THEN
        CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'&
          //TRIM(DXCoil(DXCoilNum)%Name)//'" - Low condenser dry-bulb temperature error continues...' &
          ,DXCoil(DXCoilNum)%LowAmbErrIndex,DXCoil(DXCoilNum)%LowTempLast,DXCoil(DXCoilNum)%LowTempLast,  &
            ReportMinUnits='[C]',ReportMaxUnits='[C]')
    ELSE
        CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'&
          //TRIM(DXCoil(DXCoilNum)%Name)//'" - Low condenser wet-bulb temperature error continues...' &
          ,DXCoil(DXCoilNum)%LowAmbErrIndex,DXCoil(DXCoilNum)%LowTempLast,DXCoil(DXCoilNum)%LowTempLast,  &
            ReportMinUnits='[C]',ReportMaxUnits='[C]')
    END IF
  END IF
END IF

IF(DXCoil(DXCoilNum)%PrintLowOutTempMessage)THEN
  IF(CurrentEndTime .GT. DXCoil(DXCoilNum)%CurrentEndTimeLast .AND. &
       TimeStepSys .GE. DXCoil(DXCoilNum)%TimeStepSysLast)THEN
    IF(DXCoil(DXCoilNum)%LowOutletTempIndex == 0)THEN
      CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%LowOutTempBuffer1))
      CALL ShowContinueError(TRIM(DXCoil(DXCoilNum)%LowOutTempBuffer2))
      CALL ShowContinueError('... Possible reasons for low outlet air dry-bulb temperatures are: This DX coil')
      CALL ShowContinueError('   1) may have a low inlet air dry-bulb temperature. Inlet air temperature = '// &
                                    TRIM(TrimSigDigits(DXCoil(DXCoilNum)%FullLoadInletAirTempLast,3))//' C.')
      CALL ShowContinueError('   2) may have a low air flow rate per watt of cooling capacity. Check inputs.')
      CALL ShowContinueError('   3) is used as part of a HX assisted cooling coil which uses a high sensible'// &
                                   ' effectiveness. Check inputs.')
    END IF
    CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'&
               //TRIM(DXCoil(DXCoilNum)%Name)//'" - Full load outlet temperature'// &
          ' indicates a possibility of frost/freeze error continues. Outlet air temperature statistics follow:', &
          DXCoil(DXCoilNum)%LowOutletTempIndex, DXCoil(DXCoilNum)%FullLoadOutAirTempLast, &
             DXCoil(DXCoilNum)%FullLoadOutAirTempLast)
  END IF
END IF

! save last system time step and last end time of current time step (used to determine if warning is valid)
DXCoil(DXCoilNum)%TimeStepSysLast    = TimeStepSys
DXCoil(DXCoilNum)%CurrentEndTimeLast = CurrentEndTime
DXCoil(DXCoilNum)%PrintLowAmbMessage = .FALSE.
DXCoil(DXCoilNum)%PrintLowOutTempMessage = .FALSE.

IF((AirMassFlow .GT. 0.0) .AND. &
   (GetCurrentScheduleValue(DXCoil(DXCoilNum)%SchedPtr) .GT. 0.0 .OR. &
    DXCoil(DXCoilNum)%DXCoilType_Num .EQ. CoilDX_HeatPumpWaterHeater) .AND. &
   (PartLoadRatio .GT. 0.0) .AND. (CompOp == On)) THEN      ! for cycling fan, reset mass flow to full on rate
  IF (FanOpMode .EQ. CycFanCycCoil) THEN
    AirMassFlow = AirMassFlow / (PartLoadRatio/DXcoolToHeatPLRRatio)
  ELSE IF (FanOpMode .EQ. ContFanCycCoil .AND. &
           DXCoil(DXCoilNum)%DXCoilType_Num .NE. CoilDX_CoolingTwoSpeed) THEN
    AirMassFlow = AirMassFlow * AirFlowRatio
  ELSE
    AirMassFlow = DXCoil(DXCoilNum)%RatedAirMassFlowRate(Mode)
  END IF

! Check for valid air volume flow per rated total cooling capacity (200 - 500 cfm/ton)

! for some reason there are diff's when using coil inlet air pressure
! these lines (more to follow) are commented out for the time being

  InletAirWetbulbC = PsyTwbFnTdbWPb(InletAirDryBulbTemp,InletAirHumRat,OutdoorPressure)
  AirVolumeFlowRate = AirMassFlow/ PsyRhoAirFnPbTdbW(OutdoorPressure,InletAirDryBulbTemp, InletAirHumRat)
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!  InletAirWetbulbC = PsyTwbFnTdbWPb(InletAirDryBulbTemp,InletAirHumRat,InletAirPressure)
!  AirVolumeFlowRate = AirMassFlow/ PsyRhoAirFnPbTdbW(InletAirPressure,InletAirDryBulbTemp, InletAirHumRat)
  IF (DXCoil(DXCoilNum)%RatedTotCap(Mode) .LE. 0.0) THEN
      CALL ShowFatalError(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//&
        '" - Rated total cooling capacity is zero or less.')
  END IF
  IF(DXCoil(DXCoilNum)%DXCoilType_Num .EQ. CoilDX_HeatPumpWaterHeater)THEN
    VolFlowperRatedTotCap = AirVolumeFlowRate/DXCoil(DXCoilNum)%RatedTotCap2
  ELSE
    VolFlowperRatedTotCap = AirVolumeFlowRate/DXCoil(DXCoilNum)%RatedTotCap(Mode)
  END IF
  IF (.NOT. FirstHVACIteration .AND. &
      .NOT. WarmupFlag .AND. DXCoil(DXCoilNum)%DXCoilType_Num .NE. CoilDX_HeatPumpWaterHeater .AND. &
      ((VolFlowperRatedTotCap .LT. MinOperVolFlowPerRatedTotCap) .OR. &
       (VolFlowperRatedTotCap .GT. MaxCoolVolFlowPerRatedTotCap))) THEN
    IF (DXCoil(DXCoilNum)%ErrIndex1 == 0) THEN
      CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//&
        '" - Air volume flow rate per watt of rated total cooling capacity is out of range at '//  &
          TRIM(RoundSigDigits(VolFlowperRatedTotCap,3))//' m3/s/W.')
      CALL ShowContinueErrorTimeStamp(' ')
      CALL ShowContinueError('Expected range for VolumeFlowPerRatedTotalCapacity=['//  &
          TRIM(RoundSigDigits(MinOperVolFlowPerRatedTotCap,3))//'--'//  &
          TRIM(RoundSigDigits(MaxCoolVolFlowPerRatedTotCap,3))//']')
      CALL ShowContinueError('Possible causes include inconsistent air flow rates in system components,')
      CALL ShowContinueError('or variable air volume [VAV] system using incorrect coil type.')
    ENDIF
    CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//&
          '" - Air volume flow rate per watt of rated total cooling capacity is out ' //&
          'of range error continues...',DXCoil(DXCoilNum)%ErrIndex1,VolFlowperRatedTotCap,VolFlowperRatedTotCap)
  ELSEIF (.NOT. WarmupFlag .AND. DXCoil(DXCoilNum)%DXCoilType_Num .EQ. CoilDX_HeatPumpWaterHeater .AND. &
      ((VolFlowperRatedTotCap .LT. MinOperVolFlowPerRatedTotCap) .OR. &
       (VolFlowperRatedTotCap .GT. MaxHeatVolFlowPerRatedTotCap))) THEN
    IF (DXCoil(DXCoilNum)%ErrIndex1 == 0) THEN
      CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//&
             '" - Air volume flow rate per watt of rated total water heating capacity is out of range at '// &
          TRIM(RoundSigDigits(VolFlowperRatedTotCap,2))//' m3/s/W.')
      CALL ShowContinueErrorTimeStamp(' ')
      CALL ShowContinueError('Expected range for VolumeFlowPerRatedTotalCapacity=['//  &
          TRIM(RoundSigDigits(MinOperVolFlowPerRatedTotCap,3))//'--'//  &
          TRIM(RoundSigDigits(MaxHeatVolFlowPerRatedTotCap,3))//']')
      CALL ShowContinueError('Possible causes may be that the parent object is calling for an actual supply air flow'//&
                               ' rate that is much higher or lower than the DX coil rated supply air flow rate.')
    ENDIF
    CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//&
          '" - Air volume flow rate per watt of rated total water heating capacity is out ' //&
          'of range error continues...',DXCoil(DXCoilNum)%ErrIndex1,VolFlowperRatedTotCap,VolFlowperRatedTotCap)
  END IF
!
!    Adjust coil bypass factor for actual air flow rate. Use relation CBF = exp(-NTU) where
!    NTU = A0/(m*cp). Relationship models the cooling coil as a heat exchanger with Cmin/Cmax = 0.

  RatedCBF = DXCoil(DXCoilNum)%RatedCBF(Mode)
  IF (RatedCBF .gt. 0.0) THEN
     A0 = -log(RatedCBF)*DXCoil(DXCoilNum)%RatedAirMassFlowRate(Mode)
  ELSE
     A0 = 0.
  END IF
  ADiff=-A0/AirMassFlow
  IF (ADiff >= EXP_LowerLimit) THEN
     CBF = exp(ADiff)
  ELSE
     CBF = 0.0
  END IF

!   check boundary for low ambient temperature and post warnings to individual DX coil buffers to print at end of time step
    IF (DXCoil(DXCoilNum)%CondenserType(Mode) .EQ. AirCooled) THEN
      IF(OutdoorDryBulb .LT. 0.0 .AND. .NOT. WarmupFlag) THEN !Same threshold as for air-cooled electric chiller
        DXCoil(DXCoilNum)%PrintLowAmbMessage = .TRUE.
        DXCoil(DXCoilNum)%LowTempLast = OutdoorDryBulb
        IF(DXCoil(DXCoilNum)%LowAmbErrIndex == 0)THEN
          DXCoil(DXCoilNum)%LowAmbBuffer1 = TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)// &
           '" - Air-cooled condenser inlet dry-bulb temperature below 0 C. Outdoor dry-bulb temperature = '//  &
              TRIM(RoundSigDigits(OutdoorDryBulb,2))
          DXCoil(DXCoilNum)%LowAmbBuffer2 = ' '//'... Occurrence info = '//TRIM(EnvironmentName)//', '//Trim(CurMnDy)//' '&
                     //TRIM(CreateSysTimeIntervalString())
        END IF
      END IF
    ELSEIF (DXCoil(DXCoilNum)%CondenserType(Mode) == EvapCooled) THEN
      IF(OutdoorWetBulb .LT. 10.0d0 .AND. .NOT. WarmUpFlag) THEN !Same threshold as for evap-cooled electric chiller
        DXCoil(DXCoilNum)%PrintLowAmbMessage = .TRUE.
        DXCoil(DXCoilNum)%LowTempLast = OutdoorWetBulb
        IF(DXCoil(DXCoilNum)%LowAmbErrIndex == 0)THEN
          DXCoil(DXCoilNum)%LowAmbBuffer1 = TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)// &
           '" - Evap-cooled condenser inlet wet-bulb temperature below 10 C. Outdoor wet-bulb temperature = '//  &
            TRIM(RoundSigDigits(OutdoorWetBulb,2))
          DXCoil(DXCoilNum)%LowAmbBuffer2 = ' '//'... Occurrence info = '//TRIM(EnvironmentName)//', '//Trim(CurMnDy)//' '&
                     //TRIM(CreateSysTimeIntervalString())
        END IF
      END IF
    END IF

!  Get total capacity modifying factor (function of temperature) for off-rated conditions
!  InletAirHumRat may be modified in this ADP/BF loop, use temporary varible for calculations
   InletAirHumRatTemp = InletAirHumRat
50 IF(DXCoil(DXCoilNum)%DXCoilType_Num .EQ. CoilDX_HeatPumpWaterHeater) THEN
!    Coil:DX:HeatPumpWaterHeater does not have total cooling capacity as a function of temp or flow curve
     TotCapTempModFac = 1.0d0
     TotCapFlowModFac = 1.0d0
   ELSE
     IF(DXCoil(DXCoilNum)%TotCapTempModFacCurveType(Mode) .EQ. Biquadratic) THEN
       TotCapTempModFac = CurveValue(DXCoil(DXCoilNum)%CCapFTemp(Mode),InletAirWetbulbC,CondInletTemp)
     ELSE
       TotCapTempModFac = CurveValue(DXCoil(DXCoilNum)%CCapFTemp(Mode),CondInletTemp)
     END IF

!    Warn user if curve output goes negative
     IF(TotCapTempModFac .LT. 0.0)THEN
       IF(DXCoil(DXCoilNum)%CCapFTempErrorIndex == 0)THEN
         CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//'":')
         CALL ShowContinueError(' Total Cooling Capacity Modifier curve (function of temperature) output is negative (' &
                           //TRIM(TrimSigDigits(TotCapTempModFac,3))//').')
         IF(DXCoil(DXCoilNum)%TotCapTempModFacCurveType(Mode) .EQ. Biquadratic) THEN
           CALL ShowContinueError(' Negative value occurs using a condenser inlet air temperature of ' &
                            //TRIM(TrimSigDigits(CondInletTemp,1))// &
                               ' and an inlet air wet-bulb temperature of '//TRIM(TrimSigDigits(InletAirWetbulbC,1))//'.')
         ELSE
           CALL ShowContinueError(' Negative value occurs using a condenser inlet air temperature of ' &
                            //TRIM(TrimSigDigits(CondInletTemp,1))//'.')
         END IF
         IF(Mode .GT. 1)THEN
           Call ShowContinueError(' Negative output results from stage '//TRIM(TrimSigDigits(Mode))// &
                                  ' compressor operation.')
         END IF
         CALL ShowContinueErrorTimeStamp(' Resetting curve output to zero and continuing simulation.')
       END IF
       CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//'":'//&
           ' Total Cooling Capacity Modifier curve (function of temperature) output is negative warning continues...' &
           , DXCoil(DXCoilNum)%CCapFTempErrorIndex, TotCapTempModFac, TotCapTempModFac)
       TotCapTempModFac = 0.0
     END IF

!    Get total capacity modifying factor (function of mass flow) for off-rated conditions
     AirMassFlowRatio = AirMassFlow/DXCoil(DXCoilNum)%RatedAirMassFlowRate(Mode)
     TotCapFlowModFac = CurveValue(DXCoil(DXCoilNum)%CCapFFlow(Mode),AirMassFlowRatio)

!    Warn user if curve output goes negative
     IF(TotCapFlowModFac .LT. 0.0)THEN
       IF(DXCoil(DXCoilNum)%CCapFFlowErrorIndex == 0)THEN
         CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//'":')
         CALL ShowContinueError(' Total Cooling Capacity Modifier curve (function of flow fraction) output is negative (' &
                           //TRIM(TrimSigDigits(TotCapFlowModFac,3))//').')
         CALL ShowContinueError(' Negative value occurs using an air flow fraction of ' &
                            //TRIM(TrimSigDigits(AirMassFlowRatio,3))//'.')
         CALL ShowContinueErrorTimeStamp(' Resetting curve output to zero and continuing simulation.')
         IF(Mode .GT. 1)THEN
           Call ShowContinueError(' Negative output results from stage '//TRIM(TrimSigDigits(Mode))// &
                                  ' compressor operation.')
         END IF
       END IF
       CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//'":'//&
          ' Total Cooling Capacity Modifier curve (function of flow fraction) output is negative warning continues...' &
          , DXCoil(DXCoilNum)%CCapFFlowErrorIndex, TotCapFlowModFac, TotCapFlowModFac)
       TotCapFlowModFac = 0.0
     END IF
    END IF

  TotCap = DXCoil(DXCoilNum)%RatedTotCap(Mode) * TotCapFlowModFac * TotCapTempModFac

! Calculate apparatus dew point conditions using TotCap and CBF
  hDelta = TotCap/AirMassFlow
  hADP = InletAirEnthalpy - hDelta/(1.d0-CBF)
  tADP = PsyTsatFnHPb(hADP,OutdoorPressure,'CalcDoe2DXCoil')
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!  tADP = PsyTsatFnHPb(hADP,InletAirPressure)
  wADP = PsyWFnTdbH(tADP,hADP,'CalcDoe2DXCoil')
  hTinwADP = PsyHFnTdbW(InletAirDryBulbTemp,wADP,'CalcDoe2DXCoil')
  IF((InletAirEnthalpy-hADP) .NE. 0)THEN
    SHR = MIN((hTinwADP-hADP)/(InletAirEnthalpy-hADP),1.d0)
  ELSE
    SHR = 1.0d0
  END IF
!
! Check for dry evaporator conditions (win < wadp)
!
  IF (wADP .gt. InletAirHumRatTemp .or. (Counter .ge. 1 .and. Counter .lt. MaxIter)) THEN
     If(InletAirHumRatTemp == 0.0)InletAirHumRatTemp=0.00001d0
     werror = (InletAirHumRatTemp - wADP)/InletAirHumRatTemp
!
! Increase InletAirHumRatTemp at constant InletAirTemp to find coil dry-out point. Then use the
! capacity at the dry-out point to determine exiting conditions from coil. This is required
! since the TotCapTempModFac doesn't work properly with dry-coil conditions.
!
      InletAirHumRatTemp = RF*wADP + (1.d0-RF)*InletAirHumRatTemp
      InletAirWetbulbC = PsyTwbFnTdbWPb(InletAirDryBulbTemp,InletAirHumRatTemp,OutdoorPressure)
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!      InletAirWetbulbC = PsyTwbFnTdbWPb(InletAirDryBulbTemp,InletAirHumRatTemp,InletAirPressure)
      Counter = Counter + 1
      IF (ABS(werror) .gt. Tolerance) go to 50   ! Recalculate with modified inlet conditions

  END IF

  IF(DXCoil(DXCoilNum)%PLFFPLR(mode) .GT. 0)THEN
    PLF = CurveValue(DXCoil(DXCoilNum)%PLFFPLR(mode),PartLoadRatio) ! Calculate part-load factor
  ELSE
    PLF = 1.0d0
  END IF

  IF (PLF < 0.7d0) THEN
    IF (DXCoil(DXCoilNum)%ErrIndex2 == 0) THEN
      IF(DXCoil(DXCoilNum)%DXCoilType_Num .EQ. CoilDX_HeatPumpWaterHeater)THEN
        CALL ShowWarningMessage('The PLF curve value for the Heat Pump Water Heater DX coil '//TRIM(DXCoil(DXCoilNum)%Name)//&
                         ' ='//TRIM(RoundSigDigits(PLF,3))//  &
                         ' for part-load ratio ='//TRIM(RoundSigDigits(PartLoadRatio,3)))
        CALL ShowContinueErrorTimeStamp('PLF curve values must be >= 0.7. PLF has been reset to 0.7 and simulation is continuing.')
        CALL ShowContinueError('Check the IO reference manual for PLF curve guidance [Coil:WaterHeating:AirToWaterHeatPump].')
      ELSE
        CALL ShowWarningMessage('The PLF curve value for the DX cooling coil '//TRIM(DXCoil(DXCoilNum)%Name)//&
                         ' ='//TRIM(RoundSigDigits(PLF,3))//  &
                         ' for part-load ratio ='//TRIM(RoundSigDigits(PartLoadRatio,3)))
        CALL ShowContinueErrorTimeStamp('PLF curve values must be >= 0.7. PLF has been reset to 0.7 and simulation is continuing.')
        CALL ShowContinueError('Check the IO reference manual for PLF curve guidance [Coil:Cooling:DX:SingleSpeed].')
      END IF
    END IF
    IF(DXCoil(DXCoilNum)%DXCoilType_Num .EQ. CoilDX_HeatPumpWaterHeater)THEN
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%Name)// &
                                            ', Heat Pump Water Heater DX coil PLF curve < 0.7 warning continues...', &
        DXCoil(DXCoilNum)%ErrIndex2,PLF,PLF)
    ELSE
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%Name)// &
                                            ', DX cooling coil PLF curve < 0.7 warning continues...', &
        DXCoil(DXCoilNum)%ErrIndex2,PLF,PLF)
    END IF
    PLF = 0.7d0
  END IF


  DXCoil(DXCoilNum)%PartLoadRatio = PartLoadRatio
  DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction = PartLoadRatio / PLF
  IF (DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction > 1.0d0 .and.   &
      ABS(DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction-1.0d0) > .001d0 ) THEN
    IF (DXCoil(DXCoilNum)%ErrIndex3 == 0) THEN
      IF(DXCoil(DXCoilNum)%DXCoilType_Num .EQ. CoilDX_HeatPumpWaterHeater)THEN
        CALL ShowWarningMessage('The runtime fraction for Heat Pump Water Heater DX coil '//TRIM(DXCoil(DXCoilNum)%Name)//&
                           ' exceeded 1.0. ['//TRIM(RoundSigDigits(DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction,4))//'].')
        CALL ShowContinueError('Runtime fraction reset to 1 and the simulation will continue.')
        CALL ShowContinueError('Check the IO reference manual for PLF curve guidance [Coil:WaterHeating:AirToWaterHeatPump].')
        CALL ShowContinueErrorTimeStamp(' ')
      ELSE
        CALL ShowWarningMessage('The runtime fraction for DX cooling coil '//TRIM(DXCoil(DXCoilNum)%Name)//&
                           ' exceeded 1.0. ['//TRIM(RoundSigDigits(DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction,4))//'].')
        CALL ShowContinueError('Runtime fraction reset to 1 and the simulation will continue.')
        CALL ShowContinueError('Check the IO reference manual for PLF curve guidance [Coil:Cooling:DX:SingleSpeed].')
        CALL ShowContinueErrorTimeStamp(' ')
      END IF
    END IF
    IF(DXCoil(DXCoilNum)%DXCoilType_Num .EQ. CoilDX_HeatPumpWaterHeater)THEN
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%Name)//              &
                   ', Heat Pump Water Heater DX coil runtime fraction > 1.0 warning continues...', &
        DXCoil(DXCoilNum)%ErrIndex3,DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction,DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction)
    ELSE
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%Name)//              &
                     ', DX cooling coil runtime fraction > 1.0 warning continues...', &
        DXCoil(DXCoilNum)%ErrIndex3,DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction,DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction)
    END IF
    DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction = 1.0d0 ! Reset coil runtime fraction to 1.0
  ELSEIF (DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction > 1.0d0) THEN
    DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction = 1.0d0 ! Reset coil runtime fraction to 1.0
  END IF

  ! If cycling fan, send coil part-load fraction to on/off fan via HVACDataGlobals
  IF (FanOpMode .EQ. CycFanCycCoil) OnOffFanPartLoadFraction = PLF

  !  Calculate full load output conditions
  IF (SHR .gt. 1. .OR. Counter .gt. 0) SHR = 1.d0
  FullLoadOutAirEnth = InletAirEnthalpy - TotCap/AirMassFlow
  hTinwout = InletAirEnthalpy - (1.0d0-SHR)*hDelta
  IF (SHR < 1.0d0) THEN
    FullLoadOutAirHumRat = PsyWFnTdbH(InletAirDryBulbTemp,hTinwout)
  ELSE
    FullLoadOutAirHumRat = InletAirHumRat
  END IF
  FullLoadOutAirTemp = PsyTdbFnHW(FullLoadOutAirEnth,FullLoadOutAirHumRat)

! Check for saturation error and modify temperature at constant enthalpy
   IF(FullLoadOutAirTemp .LT. PsyTsatFnHPb(FullLoadOutAirEnth,OutdoorPressure)) THEN
    FullLoadOutAirTemp = PsyTsatFnHPb(FullLoadOutAirEnth,OutdoorPressure)
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!   IF(FullLoadOutAirTemp .LT. PsyTsatFnHPb(FullLoadOutAirEnth,InletAirPressure)) THEN
!    FullLoadOutAirTemp = PsyTsatFnHPb(FullLoadOutAirEnth,InletAirPressure)
    FullLoadOutAirHumRat  = PsyWFnTdbH(FullLoadOutAirTemp,FullLoadOutAirEnth)
   END IF

  ! Store actual outlet conditions when DX coil is ON for use in heat recovery module
  DXCoilFullLoadOutAirTemp(DXCoilNum) = FullLoadOutAirTemp
  DXCoilFullLoadOutAirHumRat(DXCoilNum) = FullLoadOutAirHumRat

! Add warning message for cold cooling coil (FullLoadOutAirTemp < 2 C)
  IF(FullLoadOutAirTemp .LT. 2.0d0 .AND. .NOT. FirstHVACIteration .AND. .NOT. WarmupFlag)THEN
    DXCoil(DXCoilNum)%PrintLowOutTempMessage = .TRUE.
    DXCoil(DXCoilNum)%FullLoadOutAirTempLast = FullLoadOutAirTemp
    IF(DXCoil(DXCoilNum)%LowOutletTempIndex == 0)THEN
      DXCoil(DXCoilNum)%FullLoadInletAirTempLast = InletAirDryBulbTemp
      DXCoil(DXCoilNum)%LowOutTempBuffer1= TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)// &
       '" - Full load outlet air dry-bulb temperature < 2C. This indicates the possibility of coil frost/freeze.'// &
       ' Outlet temperature = '//TRIM(RoundSigDigits(FullLoadOutAirTemp,2))//' C.'
      DXCoil(DXCoilNum)%LowOutTempBuffer2 = ' '//'...Occurrence info = '//TRIM(EnvironmentName)//', '//Trim(CurMnDy)//' '&
                   //TRIM(CreateSysTimeIntervalString())
    END IF
  END IF

  !  If constant fan with cycling compressor, call function to determine "effective SHR"
  !  which includes the part-load degradation on latent capacity
  IF (FanOpMode .EQ. ContFanCycCoil) THEN
     QLatRated = DXCoil(DXCoilNum)%RatedTotCap(Mode) * (1.d0 - DXCoil(DXCoilNum)%RatedSHR(Mode))
     QLatActual = TotCap * (1.d0 - SHR)
     SHRUnadjusted = SHR
     SHR = CalcEffectiveSHR(DXCoilNum, SHR, DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction, &
                          QLatRated, QLatActual, InletAirDryBulbTemp, InletAirWetbulbC, Mode)
    ! For multimode coil, if stage-2 operation (modes 2 or 4), adjust Stage1&2 SHR to account for
    ! Stage 1 operating at full load, so there is no degradation for that portion
    ! Use the stage 1 bypass fraction to allocate
    IF (Mode .EQ. 2 .OR. Mode .EQ. 4) THEN
      SHR = SHRUnadjusted*(1.d0-DXCoil(DXCoilNum)%BypassedFlowFrac(Mode-1))+ &
            SHR*DXCoil(DXCoilNum)%BypassedFlowFrac(Mode-1)
    END IF

  !  Calculate full load output conditions
    IF (SHR .gt. 1.d0 .OR. Counter .gt. 0) SHR = 1.d0
    FullLoadOutAirEnth = InletAirEnthalpy - TotCap/AirMassFlow
    hTinwout = InletAirEnthalpy - (1.0d0-SHR)*hDelta
    IF (SHR < 1.0d0) THEN
      FullLoadOutAirHumRat = PsyWFnTdbH(InletAirDryBulbTemp,hTinwout)
    ELSE
      FullLoadOutAirHumRat = InletAirHumRat
    END IF
    FullLoadOutAirTemp = PsyTdbFnHW(FullLoadOutAirEnth,FullLoadOutAirHumRat)

! apply latent degradation model to cycling fan when RH control is desired and heating coil operates
! longer than the cooling coil. DXcoolToHeatPLRRatio = Cooling coil PLR / Heating coil PLR.
  ELSE IF(FanOpMode .EQ. CycFanCycCoil) THEN
    IF(DXcoolToHeatPLRRatio .LT. 1.0d0)THEN
      QLatRated = DXCoil(DXCoilNum)%RatedTotCap(Mode) * (1.d0 - DXCoil(DXCoilNum)%RatedSHR(Mode))
      QLatActual = TotCap * (1.d0 - SHR)
      HeatRTF = PartLoadRatio/DXcoolToHeatPLRRatio
      IF(DXCoil(DXCoilNum)%HeatingCoilPLFCurvePTR .GT. 0)THEN
        HeatingCoilPLF = CurveValue(DXCoil(DXCoilNum)%HeatingCoilPLFCurvePTR,HeatRTF)
        IF(HeatingCoilPLF .GT. 0) HeatRTF = HeatRTF/HeatingCoilPLF
      END IF
      SHRUnadjusted = SHR
      SHR = CalcEffectiveSHR(DXCoilNum, SHR, DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction, &
                     QLatRated, QLatActual, InletAirDryBulbTemp, InletAirWetbulbC, Mode, HeatRTF)
  !   Calculate full load output conditions
      IF (SHR .gt. 1.d0 .OR. Counter .gt. 0) SHR = 1.d0
      FullLoadOutAirEnth = InletAirEnthalpy - TotCap/AirMassFlow
      hTinwout = InletAirEnthalpy - (1.0d0-SHR)*hDelta
      IF (SHR < 1.0d0) THEN
        FullLoadOutAirHumRat = PsyWFnTdbH(InletAirDryBulbTemp,hTinwout)
      ELSE
        FullLoadOutAirHumRat = InletAirHumRat
      END IF
      FullLoadOutAirTemp = PsyTdbFnHW(FullLoadOutAirEnth,FullLoadOutAirHumRat)
    END IF

  END IF

!  Calculate actual outlet conditions for the input part load ratio
!  Actual outlet conditions are "average" for time step

! For multimode coil, if stage-2 operation (modes 2 or 4), return "full load" outlet conditions
  IF ((FanOpMode .EQ. ContFanCycCoil) .AND. &
      (Mode .EQ. 1) .OR. (Mode .EQ. 3)) THEN
    ! Continuous fan, cycling compressor
    OutletAirEnthalpy = ((PartLoadRatio * AirFlowRatio)*FullLoadOutAirEnth + &
                                            (1.d0-(PartLoadRatio * AirFlowRatio))*InletAirEnthalpy)
    OutletAirHumRat = ((PartLoadRatio * AirFlowRatio)*FullLoadOutAirHumRat + &
                                            (1.d0-(PartLoadRatio * AirFlowRatio))*InletAirHumRat)
    OutletAirTemp = PsyTdbFnHW(OutletAirEnthalpy,OutletAirHumRat)
  ELSE
    ! Default to cycling fan, cycling compressor
    ! Also return this result for stage 2 operation of multimode coil
    ! Cycling fan typically provides full outlet conditions. When RH control is used, account for additional
    ! heating run time by using cooing/heating ratio the same as constant fan (otherwise PLRRatio = 1).
    OutletAirEnthalpy = FullLoadOutAirEnth * DXcoolToHeatPLRRatio + InletAirEnthalpy * (1.0d0 - DXcoolToHeatPLRRatio)
    OutletAirHumRat = FullLoadOutAirHumRat * DXcoolToHeatPLRRatio + InletAirHumRat * (1.0d0 - DXcoolToHeatPLRRatio)
    OutletAirTemp = FullLoadOutAirTemp * DXcoolToHeatPLRRatio + InletAirDryBulbTemp * (1.0d0 - DXcoolToHeatPLRRatio)
  END IF

! Check for saturation error and modify temperature at constant enthalpy
   IF(OutletAirTemp .LT. PsyTsatFnHPb(OutletAirEnthalpy,OutdoorPressure,'CalcDOE2DXCoil')) THEN
    OutletAirTemp = PsyTsatFnHPb(OutletAirEnthalpy,OutdoorPressure)
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!   IF(OutletAirTemp .LT. PsyTsatFnHPb(OutletAirEnthalpy,InletAirPressure)) THEN
!    OutletAirTemp = PsyTsatFnHPb(OutletAirEnthalpy,InletAirPressure)
    OutletAirHumRat  = PsyWFnTdbH(OutletAirTemp,OutletAirEnthalpy)
   END IF

! Mix with air that was bypassed around coil, if any
   IF(BypassFlowFraction .GT. 0.0) THEN
     OutletAirEnthalpy = (1.d0-BypassFlowFraction)*OutletAirEnthalpy + BypassFlowFraction*InletAirEnthalpy
     OutletAirHumRat = (1.d0-BypassFlowFraction)*OutletAirHumRat + BypassFlowFraction*InletAirHumRat
     OutletAirTemp = PsyTdbFnHW(OutletAirEnthalpy,OutletAirHumRat)
     ! Check for saturation error and modify temperature at constant enthalpy
     IF(OutletAirTemp .LT. PsyTsatFnHPb(OutletAirEnthalpy,OutdoorPressure)) THEN
       OutletAirTemp = PsyTsatFnHPb(OutletAirEnthalpy,OutdoorPressure)
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!     IF(OutletAirTemp .LT. PsyTsatFnHPb(OutletAirEnthalpy,InletAirPressure)) THEN
!       OutletAirTemp = PsyTsatFnHPb(OutletAirEnthalpy,InletAirPressure)
       OutletAirHumRat  = PsyWFnTdbH(OutletAirTemp,OutletAirEnthalpy)
     END IF
    END IF

! Calculate electricity consumed. First, get EIR modifying factors for off-rated conditions
  IF(DXCoil(DXCoilNum)%DXCoilType_Num .EQ. CoilDX_HeatPumpWaterHeater) THEN
!   Coil:DX:HeatPumpWaterHeater does not have EIR temp or flow curves
    EIRTempModFac = 1.0d0
    EIRFlowModFac = 1.0d0
  ELSE
    EIRTempModFac = CurveValue(DXCoil(DXCoilNum)%EIRFTemp(Mode),InletAirWetbulbC,CondInletTemp)

!   Warn user if curve output goes negative
    IF(EIRTempModFac .LT. 0.0)THEN
      IF(DXCoil(DXCoilNum)%EIRFTempErrorIndex == 0)THEN
        CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//'":')
        CALL ShowContinueError(' Energy Input Ratio Modifier curve (function of temperature) output is negative (' &
                          //TRIM(TrimSigDigits(EIRTempModFac,3))//').')
        IF(DXCoil(DXCoilNum)%EIRTempModFacCurveType(Mode) .EQ. Biquadratic)THEN
          CALL ShowContinueError(' Negative value occurs using a condenser inlet air temperature of ' &
                         //TRIM(TrimSigDigits(CondInletTemp,1))// &
                               ' and an inlet air wet-bulb temperature of '//TRIM(TrimSigDigits(InletAirWetbulbC,1))//'.')
        ELSE
          CALL ShowContinueError(' Negative value occurs using a condenser inlet air temperature of ' &
                         //TRIM(TrimSigDigits(CondInletTemp,1))//'.')
        END IF
        IF(Mode .GT. 1)THEN
          Call ShowContinueError(' Negative output results from stage '//TRIM(TrimSigDigits(Mode))// &
                                 ' compressor operation.')
        END IF
        CALL ShowContinueErrorTimeStamp(' Resetting curve output to zero and continuing simulation.')
      END IF
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//'":'//&
          ' Energy Input Ratio Modifier curve (function of temperature) output is negative warning continues...' &
          , DXCoil(DXCoilNum)%EIRFTempErrorIndex, EIRTempModFac, EIRTempModFac)
      EIRTempModFac = 0.0
    END IF

    EIRFlowModFac = CurveValue(DXCoil(DXCoilNum)%EIRFFlow(Mode),AirMassFlowRatio)

!   Warn user if curve output goes negative
    IF(EIRFlowModFac .LT. 0.0)THEN
      IF(DXCoil(DXCoilNum)%EIRFFlowErrorIndex == 0)THEN
        CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//'":')
        CALL ShowContinueError(' Energy Input Ratio Modifier curve (function of flow fraction) output is negative (' &
                          //TRIM(TrimSigDigits(EIRFlowModFac,3))//').')
        CALL ShowContinueError(' Negative value occurs using an air flow fraction of ' &
                           //TRIM(TrimSigDigits(AirMassFlowRatio,3))//'.')
        CALL ShowContinueErrorTimeStamp(' Resetting curve output to zero and continuing simulation.')
        IF(Mode .GT. 1)THEN
          Call ShowContinueError(' Negative output results from stage '//TRIM(TrimSigDigits(Mode))// &
                                 ' compressor operation.')
        END IF
      END IF
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//'":'//&
         ' Energy Input Ratio Modifier curve (function of flow fraction) output is negative warning continues...' &
         , DXCoil(DXCoilNum)%EIRFFlowErrorIndex, EIRFlowModFac, EIRFlowModFac)
      EIRFlowModFac = 0.0
    END IF
  END IF

  EIR = DXCoil(DXCoilNum)%RatedEIR(Mode) * EIRFlowModFac * EIRTempModFac

! For multimode coil, if stage-2 operation (Modes 2 or 4), return "full load" power adjusted for PLF
  IF (Mode .EQ. 1 .OR. Mode .EQ. 3) THEN
    DXCoil(DXCoilNum)%ElecCoolingPower = TotCap * EIR * DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction
  ELSE
    DXCoil(DXCoilNum)%ElecCoolingPower = TotCap * EIR * DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction / PartLoadRatio
  END IF

! Reset AirMassFlow to inlet node air mass flow for final total, sensible and latent calculations
! since AirMassFlow might have been modified above (in this subroutine):
!     IF (FanOpMode .EQ. CycFanCycCoil) AirMassFlow = AirMassFlow / PartLoadRatio
!
! For multimode coil, this should be full flow including bypassed fraction
  AirMassFlow = DXCoil(DXCoilNum)%InletAirMassFlowRate
  DXCoil(DXCoilNum)%TotalCoolingEnergyRate = AirMassFlow * (InletAirEnthalpy - OutletAirEnthalpy)

! Set DataHeatGlobal heat reclaim variable for use by heat reclaim coil (part load ratio is accounted for)
! Calculation for heat reclaim needs to be corrected to use compressor power (not including condenser fan power)
  HeatReclaimDXCoil(DXCoilNum)%AvailCapacity = DXCoil(DXCoilNum)%TotalCoolingEnergyRate + DXCoil(DXCoilNum)%ElecCoolingPower

  MinAirHumRat = MIN(InletAirHumRat,OutletAirHumRat)
  DXCoil(DXCoilNum)%SensCoolingEnergyRate = AirMassFlow * &
                                         (PsyHFnTdbW(InletAirDryBulbTemp,MinAirHumRat) - &
                                          PsyHFnTdbW(OutletAirTemp,MinAirHumRat))
!  Don't let sensible capacity be greater than total capacity
  IF (DXCoil(DXCoilNum)%SensCoolingEnergyRate .GT. DXCoil(DXCoilNum)%TotalCoolingEnergyRate) THEN
     DXCoil(DXCoilNum)%SensCoolingEnergyRate = DXCoil(DXCoilNum)%TotalCoolingEnergyRate
  END IF
!
  DXCoil(DXCoilNum)%LatCoolingEnergyRate = DXCoil(DXCoilNum)%TotalCoolingEnergyRate - DXCoil(DXCoilNum)%SensCoolingEnergyRate

! Calculate crankcase heater power using the runtime fraction for this DX cooling coil only if there is no companion DX coil.
! Else use the largest runtime fraction of this DX cooling coil and the companion DX heating coil.
  IF(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil .EQ. 0) THEN
    DXCoil(DXCoilNum)%CrankcaseHeaterPower = CrankcaseHeatingPower * &
                                           (1.0d0 - DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction)
  ELSE
    DXCoil(DXCoilNum)%CrankcaseHeaterPower = CrankcaseHeatingPower * &
                                           (1.0d0 - MAX(DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction, &
                                           DXCoil(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil)%HeatingCoilRuntimeFraction))
  END IF

  IF (DXCoil(DXCoilNum)%CondenserType(Mode) == EvapCooled) THEN
  !******************
  !             WATER CONSUMPTION IN m3 OF WATER FOR DIRECT
  !             H2O [m3/sec] = Delta W[KgH2O/Kg air]*Mass Flow Air[Kg air]
  !                                /RhoWater [kg H2O/m3 H2O]
  !******************
     RhoWater = RhoH2O(OutdoorDryBulb)
     DXCoil(DXCoilNum)%EvapWaterConsumpRate =  &
             (CondInletHumrat - OutdoorHumRat) *  &
              CondAirMassFlow/RhoWater * DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction
     DXCoil(DXCoilNum)%EvapCondPumpElecPower = DXCoil(DXCoilNum)%EvapCondPumpElecNomPower(Mode) * &
                                               DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction
  ! Calculate basin heater power
    CALL CalcBasinHeaterPower(DXCoil(DXCoilNum)%BasinHeaterPowerFTempDiff,&
                              DXCoil(DXCoilNum)%BasinHeaterSchedulePtr,&
                              DXCoil(DXCoilNum)%BasinHeaterSetPointTemp,DXCoil(DXCoilNum)%BasinHeaterPower)
    IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingSingleSpeed) THEN
      DXCoil(DXCoilNum)%BasinHeaterPower = DXCoil(DXCoilNum)%BasinHeaterPower * &
                                         (1.d0 - DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction)
    ENDIF
  END IF

  DXCoil(DXCoilNum)%OutletAirTemp     = OutletAirTemp
  DXCoil(DXCoilNum)%OutletAirHumRat   = OutletAirHumRat
  DXCoil(DXCoilNum)%OutletAirEnthalpy = OutletAirEnthalpy

ELSE

  ! DX coil is off; just pass through conditions
  DXCoil(DXCoilNum)%OutletAirEnthalpy = DXCoil(DXCoilNum)%InletAirEnthalpy
  DXCoil(DXCoilNum)%OutletAirHumRat   = DXCoil(DXCoilNum)%InletAirHumRat
  DXCoil(DXCoilNum)%OutletAirTemp     = DXCoil(DXCoilNum)%InletAirTemp

  DXCoil(DXCoilNum)%ElecCoolingPower = 0.0
  DXCoil(DXCoilNum)%TotalCoolingEnergyRate = 0.0
  DXCoil(DXCoilNum)%SensCoolingEnergyRate = 0.0
  DXCoil(DXCoilNum)%LatCoolingEnergyRate = 0.0
  DXCoil(DXCoilNum)%EvapCondPumpElecPower = 0.0
  DXCoil(DXCoilNum)%EvapWaterConsumpRate = 0.0

! Reset globals when DX coil is OFF for use in heat recovery module
  DXCoilFullLoadOutAirTemp(DXCoilNum) = 0.0
  DXCoilFullLoadOutAirHumRat(DXCoilNum) = 0.0

! Calculate crankcase heater power using the runtime fraction for this DX cooling coil (here DXCoolingCoilRTF=0) if
! there is no companion DX coil, or the runtime fraction of the companion DX heating coil (here DXHeatingCoilRTF>=0).
  IF(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil .EQ. 0) THEN
    DXCoil(DXCoilNum)%CrankcaseHeaterPower = CrankcaseHeatingPower
  ELSE
    DXCoil(DXCoilNum)%CrankcaseHeaterPower = CrankcaseHeatingPower * &
                                            (1.d0-DXCoil(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil)%HeatingCoilRuntimeFraction)
  END IF

! Calculate basin heater power
  IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingTwoStageWHumControl) THEN
    IF (ANY(DXCoil(DXCoilNum)%CondenserType == EvapCooled)) THEN
      CALL CalcBasinHeaterPower(DXCoil(DXCoilNum)%BasinHeaterPowerFTempDiff,&
                               DXCoil(DXCoilNum)%BasinHeaterSchedulePtr,&
                               DXCoil(DXCoilNum)%BasinHeaterSetPointTemp,DXCoil(DXCoilNum)%BasinHeaterPower)
    ENDIF
  ELSE
    IF (DXCoil(DXCoilNum)%CondenserType(Mode) == EvapCooled) THEN
      CALL CalcBasinHeaterPower(DXCoil(DXCoilNum)%BasinHeaterPowerFTempDiff,&
                              DXCoil(DXCoilNum)%BasinHeaterSchedulePtr,&
                              DXCoil(DXCoilNum)%BasinHeaterSetPointTemp,DXCoil(DXCoilNum)%BasinHeaterPower)
    ENDIF
  ENDIF

END IF ! end of on/off if - else

!set water system demand request (if needed)
IF ( DXCoil(DxCoilNum)%EvapWaterSupplyMode == WaterSupplyFromTank) THEN
   WaterStorage(DXCoil(DXCoilNum)%EvapWaterSupTankID)%VdotRequestDemand(DXCoil(DXCoilNum)%EvapWaterTankDemandARRID) &
       = DXCoil(DXCoilNum)%EvapWaterConsumpRate
ENDIF

DXCoilOutletTemp(DXCoilNum)         = DXCoil(DXCoilNum)%OutletAirTemp
DXCoilOutletHumRat(DXCoilNum)       = DXCoil(DXCoilNum)%OutletAirHumRat
DXCoilPartLoadRatio(DXCoilNum)      = DXCoil(DXCoilNum)%PartLoadRatio
DXCoilFanOpMode(DXCoilNum)          = FanOpMode
DXCoil(DXCoilNum)%CondInletTemp     = CondInletTemp

RETURN
END SUBROUTINE CalcDoe2DXCoil

SUBROUTINE CalcHPSimDXCoil(CompOp,FirstHVACIteration,PartLoadRatio,FanOpMode,PerfMode,OnOffAirFlowRatio, &
                          CoolingHeatingPLR)    !RS: Implementation: Trying to create a new subroutine to call HPSim

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Fred Buhl
          !       DATE WRITTEN   May 2000
          !       MODIFIED       Shirey, Feb/October 2001, Feb/Mar 2004
          !                      Feb 2005 M. J. Witte, GARD Analytics, Inc.
          !                        Add new coil type COIL:DX:MultiMode:CoolingEmpirical:
          !                      April 2010 Chandan Sharma, FSEC, Added basin heater
          !       RE-ENGINEERED  Don Shirey, Aug/Sept 2000

          ! PURPOSE OF THIS SUBROUTINE:
          ! Calculates the air-side performance and electrical energy use of a direct-
          ! expansion, air-cooled cooling unit.

          ! METHODOLOGY EMPLOYED:
          ! This routine simulates the performance of air-cooled DX cooling equipment.
          ! The routine requires the user to enter the total cooling capacity, sensible heat ratio,
          ! and COP for the unit at ARI 210/240 rating conditions (26.67C [80F] dry-bulb, 19.44C [67F]
          ! wet-bulb air entering the cooling coil, 35C [95F] dry-bulb air entering the outdoor
          ! condenser. Since different manufacturer's rate their equipment at different air flow rates,
          ! the supply air flow rate corresponding to the rated capacities and rated COP must also be
          ! entered (should be between 300 cfm/ton and 450 cfm/ton). The rated information entered by
          ! the user should NOT include the thermal or electrical impacts of the supply air fan, as
          ! this is addressed by another module.

          ! With the rated performance data entered by the user, the model employs some of the
          ! DOE-2.1E curve fits to adjust the capacity and efficiency of the unit as a function
          ! of entering air temperatures and supply air flow rate (actual vs rated flow). The model
          ! does NOT employ the exact same methodology to calculate performance as DOE-2, although
          ! some of the DOE-2 curve fits are employed by this model.

          ! The model checks for coil dryout conditions, and adjusts the calculated performance
          ! appropriately.

          ! REFERENCES:
          ! ASHRAE HVAC 2 Toolkit page 4-81.
          !
          ! Henderson, H.I. Jr., K. Rengarajan and D.B. Shirey, III. 1992.The impact of comfort
          ! control on air conditioner energy use in humid climates. ASHRAE Transactions 98(2):
          ! 104-113.
          !
          ! Henderson, H.I. Jr., Danny Parker and Y.J. Huang. 2000.Improving DOE-2's RESYS routine:
          ! User Defined Functions to Provide More Accurate Part Load Energy Use and Humidity
          ! Predictions. Proceedings of ACEEE Conference.


          ! USE STATEMENTS:
  USE CurveManager,    ONLY: CurveValue
  USE DataGlobals,     ONLY: CurrentTime
  USE DataHVACGlobals, ONLY: HPWHCrankcaseDBTemp, TimeStepSys, SysTimeElapsed
  USE General,         ONLY: TrimSigDigits, RoundSigDigits, CreateSysTimeIntervalString
  USE DataWater,       ONLY: WaterStorage
  USE HVACInterfaceManager  !RS: Debugging: Testing to see if UpdateHVACInterface needs to be called (8/8/14)
  USE PackagedTerminalHeatPump !, ONLY: PTUnitHPSimNum    !RS: Debugging: Trying to bring this over so Abhijit's code works (10/15/14)
  USE WeatherManager    !RS: Debugging: Trying to bring this over so Abhijit's code works (10/15/14)
  USE HVACDXSystem  !RS: Debugging: Trying to bring over the coil inputs and outputs (11/1/14)
  USE DataZoneEnergyDemands !RS: Debugging: Trying to get it only call HPSim and optimization when there's a load (11/6/14)
  USE DataHeatBalFanSys, ONLY: ZoneThermostatSetPointHi, ZoneThermostatSetPointLo   !RS: Debugging: Trying to keep a constant setpoint (12/5/15)
  USE ZoneTempPredictorCorrector, ONLY: LoadToCoolingSetPoint_HPSim !RS: Debugging: Trying to bring this in (11/3/18)

  USE dfwin !RS: Debugging: Trying to get the code to wait until the exe has finished running before moving on (11/24/14)
  USE dflib !RS: Debugging: Trying to get the code to wait until the exe has finished running before moving on (11/24/14)
  
  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  !INTEGER,   INTENT(IN)           :: DXCoilNum           ! the number of the DX coil to be simulated
  INTEGER,   INTENT(IN)           :: CompOp              ! compressor operation; 1=on, 0=off
  LOGICAL,   INTENT(IN)           :: FirstHVACIteration  ! true if this is the first iteration of HVAC
  REAL(r64), INTENT(IN)           :: PartLoadRatio       ! sensible cooling load / full load sensible cooling capacity
  INTEGER,   INTENT(IN)           :: FanOpMode           ! Allows parent object to control fan operation
  INTEGER,   INTENT(IN), OPTIONAL :: PerfMode            ! Performance mode for MultiMode DX coil; Always 1 for other coil types
  REAL(r64), INTENT(IN), OPTIONAL :: OnOffAirFlowRatio   ! ratio of compressor on airflow to compressor off airflow
  REAL(r64), INTENT(IN), OPTIONAL :: CoolingHeatingPLR   ! used for cycling fan RH control

          ! SUBROUTINE PARAMETER DEFINITIONS:
CHARACTER(len=*), PARAMETER :: RoutineName='CalcHPSimDXCoil'

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
REAL(r64) :: AirMassFlow           ! dry air mass flow rate through coil [kg/s] (adjusted for bypass if any)
REAL(r64) :: AirMassFlowRatio      ! Ratio of actual air mass flow to rated air mass flow (adjusted for bypass if any)
REAL(r64) :: AirVolumeFlowRate     ! Air volume flow rate across the cooling coil [m3/s] (adjusted for bypass if any)
                                   ! (average flow if cycling fan, full flow if constant fan)
REAL(r64) :: VolFlowperRatedTotCap ! Air volume flow rate divided by rated total cooling capacity [m3/s-W] (adjusted for bypass)
REAL(r64) :: BypassFlowFraction    ! Fraction of total flow which is bypassed around the cooling coil
REAL(r64) :: TotCap                ! gross total cooling capacity at off-rated conditions [W]
REAL(r64) :: InletAirWetBulbC      ! wetbulb temperature of inlet air [C]
REAL(r64) :: InletAirDryBulbTemp   ! inlet air dry bulb temperature [C]
REAL(r64) :: InletAirEnthalpy      ! inlet air enthalpy [J/kg]
REAL(r64) :: InletAirHumRat        ! inlet air humidity ratio [kg/kg]
REAL(r64) :: InletAirHumRatTemp    ! inlet air humidity ratio used in ADP/BF loop [kg/kg]
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
REAL(r64) :: SHR                   ! Sensible Heat Ratio (sensible/total) of the cooling coil
REAL(r64) :: CBF                   ! coil bypass factor at off rated conditions
REAL(r64) :: hDelta                ! Change in air enthalpy across the cooling coil [J/kg]
REAL(r64) :: hADP                  ! Apparatus dew point enthalpy [J/kg]
REAL(r64) :: hTinwADP              ! Enthalpy at inlet dry-bulb and wADP [J/kg]
REAL(r64) :: hTinwout              ! Enthalpy at inlet dry-bulb and outlet humidity ratio [J/kg]
REAL(r64) :: tADP                  ! Apparatus dew point temperature [C]
REAL(r64) :: wADP                  ! Apparatus dew point humidity ratio [kg/kg]
REAL(r64) :: FullLoadOutAirEnth    ! outlet full load enthalpy [J/kg]
REAL(r64) :: FullLoadOutAirHumRat  ! outlet humidity ratio at full load
REAL(r64) :: FullLoadOutAirTemp    ! outlet air temperature at full load [C]
REAL(r64) :: EIRTempModFac         ! EIR modifier (function of entering wetbulb, outside drybulb)
REAL(r64) :: EIRFlowModFac         ! EIR modifier (function of actual supply air flow vs rated flow)
REAL(r64) :: EIR                   ! EIR at part load and off rated conditions
REAL(r64) :: PLF                   ! Part load factor, accounts for thermal lag at compressor startup, used in power calculation
REAL(r64) :: QSens                 ! actual sensible capacity
!REAL(r64) :: SysLat                ! latent add/removal from system (kg/s) !RS: Never used currently (1/16/14)
REAL(r64) :: QLatActual            ! operating latent capacity of DX coil
REAL(r64) :: QLatRated             ! Rated latent capacity of DX coil
REAL(r64) :: SHRUnadjusted         ! SHR prior to latent degradation effective SHR calculation
INTEGER :: Counter                 ! Counter for dry evaporator iterations
INTEGER :: MaxIter                 ! Maximum number of iterations for dry evaporator calculations
REAL(r64) :: RF                    ! Relaxation factor for dry evaporator iterations
REAL(r64) :: Tolerance             ! Error tolerance for dry evaporator iterations
REAL(r64) :: werror                ! Deviation of humidity ratio in dry evaporator iteration loop
REAL(r64) :: CondInletTemp         ! Condenser inlet temperature (C). Outdoor dry-bulb temp for air-cooled condenser.
                                 ! Outdoor Wetbulb +(1 - effectiveness)*(outdoor drybulb - outdoor wetbulb) for evap condenser.
REAL(r64) :: CondInletHumrat       ! Condenser inlet humidity ratio (kg/kg). Zero for air-cooled condenser.
                                 ! For evap condenser, its the humidity ratio of the air leaving the evap cooling pads.
!REAL(r64) :: CondAirMassFlow       ! Condenser air mass flow rate [kg/s]
!REAL(r64) :: RhoAir                ! Density of air [kg/m3]
!REAL(r64) :: RhoWater              ! Density of water [kg/m3]
REAL(r64) :: CrankcaseHeatingPower ! power due to crankcase heater
!REAL(r64) :: CompAmbTemp = 0.0     ! Ambient temperature at compressor
REAL(r64) :: AirFlowRatio          ! ratio of compressor on airflow to average timestep airflow
                                 ! used when constant fan mode yields different air flow rates when compressor is ON and OFF
                                 ! (e.g. Packaged Terminal Heat Pump)
REAL(r64) :: OutdoorDryBulb        ! Outdoor dry-bulb temperature at condenser (C)
REAL(r64) :: OutdoorWetBulb        ! Outdoor wet-bulb temperature at condenser (C)
REAL(r64) :: OutdoorHumRat         ! Outdoor humidity ratio at condenser (kg/kg)
REAL(r64) :: OutdoorPressure       ! Outdoor barometric pressure at condenser (Pa)

REAL(r64) :: CurrentEndTime = 0.0  ! end time of time step for current simulation time step
!REAL(r64) :: MinAirHumRat = 0.0    ! minimum of the inlet air humidity ratio and the outlet air humidity ratio
INTEGER          :: Mode           ! Performance mode for Multimode DX coil; Always 1 for other coil types
REAL(r64) :: OutletAirTemp           ! Supply air temperature (average value if constant fan, full output if cycling fan)
REAL(r64) :: OutletAirHumRat         ! Supply air humidity ratio (average value if constant fan, full output if cycling fan)
REAL(r64) :: OutletAirEnthalpy       ! Supply air enthalpy (average value if constant fan, full output if cycling fan)
!REAL(r64) :: Adiff                   ! Used for exponential
REAL(r64) :: DXcoolToHeatPLRRatio    ! ratio of cooling PLR to heating PLR, used for cycling fan RH control
REAL(r64) :: HeatRTF                 ! heating coil part-load ratio, used for cycling fan RH control
REAL(r64) :: HeatingCoilPLF          ! heating coil PLF (function of PLR), used for cycling fan RH control
INTEGER :: DebugFile       =150 !RS: Debugging file denotion, hopefully this works.
INTEGER :: InletNode, OutletNode    !RS: Debugging: UpdateHVACInterface (8/8/14)

LOGICAL :: file_exists  !RS: Debugging: Setting this in for AT's modification below (10/6/14)
REAL :: RhoOutAir   !RS: Debugging: Evaporator outlet air density

!RS: Comment: The following definitions were put in by Abhijit (10/15/14)
INTEGER :: MixedNode !=0
    INTEGER :: OutsideNode !=0    
    REAL :: RHiCAct !=0.0      
    REAL :: TaiC  ! =0.0!RS: Drybulb Temperature, Outdoor Entering (C)
    REAL :: TaiE  ! =0.0!RS: Drybulb Temperature, Indoor Entering (C)
    REAL :: OutDryBulbTemp !=0.0
    REAL :: OutWetBulBTemp !=0.0 !RS: Wetbulb temperature, Outdoor (C)
    REAL :: DummyHROutSideNode !=0.0 !RS: Debugging
    REAL :: DummyHR !=0.0!RS: Debugging
    REAL(r64) :: OutBaroPressAT
    REAL :: XMaE
    integer :: ok
    real :: x
    character(len=90) :: name
    REAL :: RHiE
    REAL TWiC
    REAL TWiE
    CHARACTER(LEN=*), PARAMETER :: & XMaEL = ",! XMaE"
    CHARACTER(LEN=*), PARAMETER :: & OutDryBulbTempL = ",! OutDryBulbTemp"
    CHARACTER(LEN=*), PARAMETER :: & OutWetBulBTempL = ",! OutWetBulBTemp"
    CHARACTER(LEN=*), PARAMETER :: & OutBaroPressATL = ",! OutBaroPress"
    CHARACTER(LEN=*), PARAMETER :: & RHiCActL = ",! RHiCAct"
    CHARACTER(LEN=*), PARAMETER :: & RHiEL = ",! RHiE"
    CHARACTER(LEN=*), PARAMETER :: & TaiCL = ",! TaiC"
    CHARACTER(LEN=*), PARAMETER :: & TaiEL = ";! TaiE"
    CHARACTER(LEN=*), PARAMETER :: & TWiCL = ",! TWiC"
    CHARACTER(LEN=*), PARAMETER :: & TWiEL = ",! TWiE"
    CHARACTER(LEN=*), PARAMETER :: & DummyHROutSideNodeL = ",! DummyHROutSideNode"
    CHARACTER(LEN=*), PARAMETER :: & DummyHRMixedL = ";! DummyHRMixed"
    CHARACTER(LEN=*), PARAMETER :: & VariablesToPass = "VariablesToPass,"
    CHARACTER(LEN=500) :: & FolderPath
    integer :: i
    CHARACTER(Len=200), dimension(30) :: arr
    
    integer,dimension(8) :: t ! arguments for date_and_time
    integer :: s1,s2,ms1,ms2  ! start and end times [ms]
    real :: dt                ! desired sleep interval [ms]
    
    REAL ::    totalPowerVal
    REAL ::    compressorPowerVal
    REAL ::    volumetricFlowrateVal
    REAL ::    dryBulbTemperatureVal
    REAL ::    enthalpyVal
    REAL ::    humidityRatioVal
    REAL ::    QTotalVal
    REAL ::    QSensVal
    integer :: ios = 0
    integer :: pos
    character(len=200)           :: newDirPath
    character(len=256)           :: mkdirCmd
    logical                      :: dirExists
    INTEGER :: temp !RS: Debugging: Trying to get the code to wait until the exe has finished running before moving on (11/24/14)
    LOGICAL :: result2  !RS: Debugging: Trying to get the code to wait until the exe has finished running before moving on (11/24/14)
    REAL :: Superheat_Comp  !RS: Debugging: Keeping track of the compressor superheat for HPSim (5/25/15)
    REAL :: Subcool_ExpDev  !RS: Debugging: Keeping track of the expansion device subcool for HPSim (5/25/15)
    REAL :: Subcool_Cond    !RS: Debugging: Keeping track of the condenser subcool for HPSim (5/25/15)
    REAL :: Superheat_Evap  !RS: Debugging: Keeping track of the evaporator superheat for HPSim (5/25/15)
    REAL :: TsiCmp  !RS: Debugging: Keeping track of the compressor suction saturation temperature for HPSim (5/25/15)
    REAL :: TsoCmp  !RS: Debugging: Keeping track of the compressor discharge saturation temperature for HPSim (5/25/15)
    REAL :: Qevap   !RS: Debugging: The cooling needed to cool the zone (total) (10/10/15)
    REAL :: Qzone   !RS: Debugging: The cooling needed to bring the zone air to deck temperature (10/10/15)
    REAL :: QtoSP   !RS: Debugging: The cooling needed to bring the zone air to set point (10/10/15)
    REAL :: cpair
    CHARACTER(LEN=14),PARAMETER :: FMT_101 = "(7(F10.3,','))"
    CHARACTER(LEN=14),PARAMETER :: FMT_102 = "(2(F10.5,','))"
    REAL :: CompRatio   !RS: Debugging: Keeping track of the compressor ratio for HPSim (5/25/15)
    REAL :: OutTemp_Current !RS: Debugging: Current timestep's outdoor temperature (DB) (3/10/18)
    REAL :: OutTemp_Previous  !RS: Debugging: Previous timestep's outdoor temperature (DB) (3/10/18)
    REAL :: InDTemp_Current !RS: Debugging: Current timestep's indoor temperature (DB) (3/10/18)
    REAL :: InDTemp_Previous    !RS: Debugging: Previous timestep's indoor temperature (DB) (3/10/18)
    REAL :: InWTemp_Current !RS: Debugging: Current timestep's indoor temperature (WB) (3/10/18)
    REAL :: InWTemp_Previous    !RS: Debugging: Previous timestep's indoor temperature (WB) (3/10/18)
    INTEGER, SAVE :: HPSimCounter   !RS: Debugging: Counter for how often HPSim is called in a run (4/28/18)
       
!, only: getcharqq, nargs
    !implicit none
    CHARACTER*(MAX_PATH) path   !RS: Debugging: Trying to get the code to wait until the exe has finished running before moving on (11/24/14)
    TYPE (T_STARTUPINFO) :: StartupInfo !RS: Debugging: Trying to get the code to wait until the exe has finished running before moving on (11/24/14)
    TYPE (T_PROCESS_INFORMATION) :: ProcessInfo !RS: Debugging: Trying to get the code to wait until the exe has finished running before moving on (11/24/14)
    INTEGER ret, path_len   !RS: Debugging: Trying to get the code to wait until the exe has finished running before moving on (11/24/14)
    CHARACTER*1 cret    !RS: Debugging: Trying to get the code to wait until the exe has finished running before moving on (11/24/14)
        
    CHARACTER(100) ::  path_1, path_2, StrScalar, fullname !path,
REAL a,b2,c,path_2_len,  path_1_len !path_len,i,
    
    OPEN(unit=107,file='C:\Users\lab303user\Desktop\HPSimProject\HPSimBuild\GenOpt\Qevap_E+.txt') !C:\Users\lab303user\Documents\betsrg_dual\GenOpt\Qevap_E+.txt')  !RS: Debugging: Printing out the Qevap to a file that can be read by TestProgram (10/10/15)
    OPEN(unit=DebugFile,file='Debug.txt')    !RS: Debugging
    OPEN(unit=151,file='HPSim_temps.csv')   !RS: Debugging: Keeping track of the temperatures for HPSim (5/25/15)
    OPEN(unit=152,file='HPSim_comprat.csv')   !RS: Debugging: Keeping track of the compressor ratio for HPSim (9/12/15)


! If Performance mode not present, then set to 1.  Used only by Multimode/Multispeed DX coil (otherwise mode = 1)
IF (PRESENT(PerfMode)) THEN
  Mode = PerfMode
ELSE
  Mode = 1
END IF

! If AirFlowRatio not present, then set to 1. Used only by DX coils with different air flow
! during cooling and when no cooling is required (constant fan, fan speed changes)
IF (PRESENT(OnOffAirFlowRatio)) THEN
  AirFlowRatio = OnOffAirFlowRatio
ELSE
  AirFlowRatio = 1.0d0
END IF

! If CoolingHeatingPLR not present, then set to 1. Used for cycling fan systems where
! heating PLR is greater than cooling PLR, otherwise CoolingHeatingPLR = 1.
IF(PRESENT(CoolingHeatingPLR))THEN
  DXcoolToHeatPLRRatio = CoolingHeatingPLR
ELSE
  DXcoolToHeatPLRRatio = 1.0d0
END IF

MaxIter         = 30
RF              = 0.4d0
Counter         = 0
Tolerance       = 0.01d0
CondInletTemp   = 0.0d0
CondInletHumrat = 0.0d0
BypassFlowFraction  = DXCoil(DXCoilHPSimNum)%BypassedFlowFrac(Mode)
AirMassFlow         = DXCoil(DXCoilHPSimNum)%InletAirMassFlowRate * (1.d0-BypassFlowFraction)
InletAirDryBulbTemp = DXCoil(DXCoilHPSimNum)%InletAirTemp   !RS: Debugging: Will be set in HPSim
InletAirEnthalpy    = DXCoil(DXCoilHPSimNum)%InletAirEnthalpy
InletAirHumRat      = DXCoil(DXCoilHPSimNum)%InletAirHumRat
HeatReclaimDXCoil(DXCoilHPSimNum)%AvailCapacity   = 0.0d0
DXCoil(DXCoilHPSimNum)%CoolingCoilRuntimeFraction = 0.0d0
DXCoil(DXCoilHPSimNum)%PartLoadRatio              = 0.0d0
DXCoil(DXCoilHPSimNum)%BasinHeaterPower           = 0.0d0

OutdoorDryBulb  = OutDryBulbTemp
OutdoorHumRat   = OutHumRat
OutdoorPressure = OutBaroPress
OutdoorWetBulb  = OutWetBulbTemp

!RS: Debugging: Trying to have HPSim only called when it needs to be (3/17/18)
OutTemp_Current = WeathOutDryBulb
InDTemp_Current = DXCoil(DXCoilHPSimNum)%InletAirTemp

! calculate end time of current time step to determine if error messages should be printed
CurrentEndTime = CurrentTime + SysTimeElapsed

!   Print warning messages only when valid and only for the first ocurrance. Let summary provide statistics.
!   Wait for next time step to print warnings. If simulation iterates, print out
!   the warning for the last iteration only. Must wait for next time step to accomplish this.
!   If a warning occurs and the simulation down shifts, the warning is not valid.
IF(DXCoil(DXCoilHPSimNum)%PrintLowAmbMessage)THEN ! .AND. &
  IF(CurrentEndTime .GT. DXCoil(DXCoilHPSimNum)%CurrentEndTimeLast .AND. &
     TimeStepSys .GE. DXCoil(DXCoilHPSimNum)%TimeStepSysLast)THEN
    IF (DXCoil(DXCoilHPSimNum)%LowAmbErrIndex == 0) THEN
      CALL ShowWarningMessage(TRIM(DXCoil(DXCoilHPSimNum)%LowAmbBuffer1))
      CALL ShowContinueError(TRIM(DXCoil(DXCoilHPSimNum)%LowAmbBuffer2))
      CALL ShowContinueError('... Operation at low ambient temperatures may require special performance curves.')
    ENDIF
    IF (DXCoil(DXCoilHPSimNum)%CondenserType(Mode) .EQ. AirCooled) THEN
        CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilHPSimNum)%DXCoilType)//' "'&
          //TRIM(DXCoil(DXCoilHPSimNum)%Name)//'" - Low condenser dry-bulb temperature error continues...' &
          ,DXCoil(DXCoilHPSimNum)%LowAmbErrIndex,DXCoil(DXCoilHPSimNum)%LowTempLast,DXCoil(DXCoilHPSimNum)%LowTempLast,  &
            ReportMinUnits='[C]',ReportMaxUnits='[C]')
    ELSE
        CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilHPSimNum)%DXCoilType)//' "'&
          //TRIM(DXCoil(DXCoilHPSimNum)%Name)//'" - Low condenser wet-bulb temperature error continues...' &
          ,DXCoil(DXCoilHPSimNum)%LowAmbErrIndex,DXCoil(DXCoilHPSimNum)%LowTempLast,DXCoil(DXCoilHPSimNum)%LowTempLast,  &
            ReportMinUnits='[C]',ReportMaxUnits='[C]')
    END IF
  END IF
END IF

IF(DXCoil(DXCoilHPSimNum)%PrintLowOutTempMessage)THEN
  IF(CurrentEndTime .GT. DXCoil(DXCoilHPSimNum)%CurrentEndTimeLast .AND. &
       TimeStepSys .GE. DXCoil(DXCoilHPSimNum)%TimeStepSysLast)THEN
    IF(DXCoil(DXCoilHPSimNum)%LowOutletTempIndex == 0)THEN
      CALL ShowWarningMessage(TRIM(DXCoil(DXCoilHPSimNum)%LowOutTempBuffer1))
      CALL ShowContinueError(TRIM(DXCoil(DXCoilHPSimNum)%LowOutTempBuffer2))
      CALL ShowContinueError('... Possible reasons for low outlet air dry-bulb temperatures are: This DX coil')
      CALL ShowContinueError('   1) may have a low inlet air dry-bulb temperature. Inlet air temperature = '// &
                                    TRIM(TrimSigDigits(DXCoil(DXCoilHPSimNum)%FullLoadInletAirTempLast,3))//' C.')
      CALL ShowContinueError('   2) may have a low air flow rate per watt of cooling capacity. Check inputs.')
      CALL ShowContinueError('   3) is used as part of a HX assisted cooling coil which uses a high sensible'// &
                                   ' effectiveness. Check inputs.')
    END IF
    CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilHPSimNum)%DXCoilType)//' "'&
               //TRIM(DXCoil(DXCoilHPSimNum)%Name)//'" - Full load outlet temperature'// &
          ' indicates a possibility of frost/freeze error continues. Outlet air temperature statistics follow:', &
          DXCoil(DXCoilHPSimNum)%LowOutletTempIndex, DXCoil(DXCoilHPSimNum)%FullLoadOutAirTempLast, &
             DXCoil(DXCoilHPSimNum)%FullLoadOutAirTempLast)
  END IF
END IF

! save last system time step and last end time of current time step (used to determine if warning is valid)
DXCoil(DXCoilHPSimNum)%TimeStepSysLast    = TimeStepSys
DXCoil(DXCoilHPSimNum)%CurrentEndTimeLast = CurrentEndTime
DXCoil(DXCoilHPSimNum)%PrintLowAmbMessage = .FALSE.
DXCoil(DXCoilHPSimNum)%PrintLowOutTempMessage = .FALSE.

!RS: Debugging: Resetting AirMassFlow because EMS case has it resetting to 0 when line 8465 is executed (7/16/16)
AirMassFlow         = DXCoil(DXCoilHPSimNum)%InletAirMassFlowRate * (1.d0-BypassFlowFraction)

!RS: Debugging: Trying to have HPSim only called when it needs to be (3/17/18)
IF(OutTemp_Current .EQ. OutTemp_Previous .AND. InDTemp_Current .EQ. InDTemp_Previous) THEN
    InDTemp_Previous = InDTemp_Current  !RS: Debugging: This doesn't do anything but this IF statement keeps the Python program (3/17/18)
    !RS: Debugging: con. from being called if we don't need it. (3/17/18)

!IF((AirMassFlow .GT. 0.0) .AND. &
!   (GetCurrentScheduleValue(DXCoil(DXCoilHPSimNum)%SchedPtr) .GT. 0.0 .OR. &
!    DXCoil(DXCoilHPSimNum)%DXCoilType_Num .EQ. CoilDX_HeatPumpWaterHeater) .AND. &
!   (PartLoadRatio .GT. 0.0) .AND. (CompOp == On)) THEN      ! for cycling fan, reset mass flow to full on rate
!RS: Debugging: Trying to reflect only what HPSim actually needs to be considered on (11/21/15)
!RS: Debugging: con. 
!IF((AirMassFlow .GT. 0.0) .AND. &  !RS: Debugging: Trying to have HPSim only called when it needs to be - turning IF to ELSEIF (3/17/18)
ELSEIF((AirMassFlow .GT. 0.0) .AND. &
   (GetCurrentScheduleValue(DXCoil(DXCoilHPSimNum)%SchedPtr) .GT. 0.0) .AND. &
   !((Node(9)%TempSetPoint .LT. Node(12)%Temp))) THEN !InletAirDryBulbTemp))) THEN      ! for cycling fan, reset mass flow to full on rate
   ((ZoneThermostatSetPointHi(1) .LT. Node(12)%Temp))) THEN  !RS: Debugging: Trying to keep a constant setpoint (12/5/15)
  IF (FanOpMode .EQ. CycFanCycCoil) THEN
    AirMassFlow = AirMassFlow / (PartLoadRatio/DXcoolToHeatPLRRatio)
  ELSE IF (FanOpMode .EQ. ContFanCycCoil .AND. &
           DXCoil(DXCoilHPSimNum)%DXCoilType_Num .NE. CoilDX_CoolingTwoSpeed) THEN
    AirMassFlow = AirMassFlow * AirFlowRatio
  ELSE
    AirMassFlow = DXCoil(DXCoilHPSimNum)%RatedAirMassFlowRate(Mode)
  END IF

! Check for valid air volume flow per rated total cooling capacity (200 - 500 cfm/ton)

! for some reason there are diff's when using coil inlet air pressure
! these lines (more to follow) are commented out for the time being

  InletAirWetbulbC = PsyTwbFnTdbWPb(InletAirDryBulbTemp,InletAirHumRat,OutdoorPressure)
  AirVolumeFlowRate = AirMassFlow/ PsyRhoAirFnPbTdbW(OutdoorPressure,InletAirDryBulbTemp, InletAirHumRat)
  !IF (DXCoil(DXCoilHPSimNum)%RatedTotCap(Mode) .LE. 0.0) THEN
  !    CALL ShowFatalError(TRIM(DXCoil(DXCoilHPSimNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilHPSimNum)%Name)//&
  !      '" - Rated total cooling capacity is zero or less.')
  !END IF
  IF(DXCoil(DXCoilHPSimNum)%DXCoilType_Num .EQ. CoilDX_HeatPumpWaterHeater)THEN
    VolFlowperRatedTotCap = AirVolumeFlowRate/DXCoil(DXCoilHPSimNum)%RatedTotCap2
  ELSE
    VolFlowperRatedTotCap = AirVolumeFlowRate/DXCoil(DXCoilHPSimNum)%RatedTotCap(Mode)
  END IF
 
!
!    Adjust coil bypass factor for actual air flow rate. Use relation CBF = exp(-NTU) where
!    NTU = A0/(m*cp). Relationship models the cooling coil as a heat exchanger with Cmin/Cmax = 0.
     CBF = 0.0

!   check boundary for low ambient temperature and post warnings to individual DX coil buffers to print at end of time step
    !IF (DXCoil(DXCoilNum)%CondenserType(Mode) .EQ. AirCooled) THEN
      IF(OutdoorDryBulb .LT. 0.0 .AND. .NOT. WarmupFlag) THEN !Same threshold as for air-cooled electric chiller
        DXCoil(DXCoilHPSimNum)%PrintLowAmbMessage = .TRUE.
        DXCoil(DXCoilHPSimNum)%LowTempLast = OutdoorDryBulb
        IF(DXCoil(DXCoilHPSimNum)%LowAmbErrIndex == 0)THEN
          DXCoil(DXCoilHPSimNum)%LowAmbBuffer1 = TRIM(DXCoil(DXCoilHPSimNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilHPSimNum)%Name)// &
           '" - Air-cooled condenser inlet dry-bulb temperature below 0 C. Outdoor dry-bulb temperature = '//  &
              TRIM(RoundSigDigits(OutdoorDryBulb,2))
          DXCoil(DXCoilHPSimNum)%LowAmbBuffer2 = ' '//'... Occurrence info = '//TRIM(EnvironmentName)//', '//Trim(CurMnDy)//' '&
                     //TRIM(CreateSysTimeIntervalString())
        END IF
      END IF

!  Get total capacity modifying factor (function of temperature) for off-rated conditions
!  InletAirHumRat may be modified in this ADP/BF loop, use temporary varible for calculations
   InletAirHumRatTemp = InletAirHumRat


    XMaE=DXCoil(DXCoilHPSimNum)%InletAirMassFlowRate
    PTUnitHPSimNum=1    !RS: Debugging: Setting this here (11/3/14)
    CALL HPSimNodes(PTUnitHPSimNum,MixedNode) !,OutsideNode) !RS: Debugging: Bringing in the node numbers
    !MixedNode=DXCoolingSystem(NumDXSystem)%DXCoolingCoilInletNodeNum
    !CALL GetTempsOut(OutDryBulbTemp, OutWetBulbTemp, OutBaroPressAT, RHiCAct)
    
    TaiC=WeathOutDryBulb !Node(OutsideNode)%Temp
    DummyHR=WeathOutHumRat !Node(OutsideNode)%HumRat
    CALL PsyTwbFnTdbWPb2(TaiC,DummyHR,WeathOutBaroPress,TWiC)
    !CALL PsyRhFnTdbWPb2(TaiC,DummyHR,OutBaroPressAT,RHiC)  !RS: Debugging: Converting from humidity ratio to relative humidity
    !RHiC=RHiC*100 
    
    TaiE = DXCoil(DXCoilHPSimNum)%InletAirTemp
    DummyHR=DXCoil(DXCoilHPSimNum)%InletAirHumRat 
    CALL PsyTwbFnTdbWPb2(TaiE,DummyHR,WeathOutBaroPress,TWiE)
    !CALL PsyRhFnTdbWPb2(TaiE,DummyHR,OutBaroPressAT,RHiE)  !RS: Debugging: Converting from humidity ratio to relative humidity
    !RHiE=RHiE*100 
    OutBaroPressAT=WeathOutBaroPress
    
    !RS: Debugging: Calculating the total Q needed by the zones (10/10/15)
    cpair=PsyCpAirFnWTdb(DXCoil(DXCoilHPSimNum)%InletAirHumRat,DXCoil(DXCoilHPSimNum)%InletAirTemp)
    Qzone=DXCoil(DXCoilHPSimNum)%InletAirMassFlowRate *cpair *(Node(ControlPNode)%Temp-DXCoil(DXCoilHPSimNum)%InletAirTemp) !F
    QtoSP=cpair*DXCoil(DXCoilHPSimNum)%InletAirMassFlowRate*(Node(9)%TempSetPoint - DXCoil(DXCoilHPSimNum)%InletAirTemp) !23C
    Qevap=Qzone+QtoSP !(WATTS)
    
    Qevap=-Qevap/0.2927  !RS: Debugging: Converting to Btu/h for Buffer Program comparison and turning negative (1/16/16)
    
    Qevap=-ZoneSysEnergyDemand(1)%TotalOutputRequired/0.2927    !RS: Debugging: This is a much greater value than the previous Qevap (9/8/18)
    
    Qevap=-LoadToCoolingSetPoint_HPSim/0.2927   !RS: Debugging: Trying to bring this in; converting to Btu/h (11/3/18) Adding a "-" to try to correct the case (11/10/18)
    
    WRITE(107,*)Qevap,Node(9)%TempSetPoint,DXCoil(DXCoilHPSimNum)%InletAirTemp   !RS: Debugging: Printing out the Qevap to a file that can be read by TestProgram (10/10/15)
    CLOSE(107)
    
    ! INQUIRE(FILE="FilePathBufferProgram.txt", EXIST=file_exists)
    !if (file_exists) THEN
    !    OPEN (UNIT=580, FILE="FilePathBufferProgram.txt", STATUS="OLD")   ! Current directory
    !    read(580,'(A)', iostat=ios) FolderPath
    !    read(580,'(A)', iostat=ios) FolderPath
    !    CLOSE (UNIT=580)
    !    FolderPath=TRIM(FolderPath)
    !    pos = scan(FolderPath, '\')
    !    FolderPath = TRIM(ADJUSTL(FolderPath(1:pos-1)))
    !end if
        
        path_len = GetModuleFileName (NULL, path, &
  len(path))

path_1=TRIM(path)
StrScalar='C:\Users\lab303user\Desktop\HPSimProject\HPSimBuild\GenOpt\tmp-genopt-run-6' !'C:\Users\lab303user\Documents\betsrg_dual\GenOpt\tmp-genopt-run-6'


path_1_len=LEN_TRIM(path_1)
a= INDEX(path, 'Debug',BACK=.FALSE.)

b2=path_len-a

c=path_len-b2

path_2_len=i-7
path_2=path_1(1:c)
path_2_len=len(StrScalar)
a=INDEX(StrScalar, 'tmp-genopt-run')
path_2=StrScalar(1:(a-2))

FolderPath=path_2
    
    !name=TRIM(ADJUSTL(FolderPath(1:pos-1))) // "\\HPSim_Variables_Template.idf"
    name=TRIM(FolderPath) // "\HPSim_Variables_Template.idf"
    open(130,file=name,status='old',iostat=ok,action='write',form='formatted')
    !print '("file ",A8," opened ok=",I3)', name,ok
    
    open(120,file=name,status='old')
    do i=1,7
      read(120,'(A)')arr(i)
    end do
    close(120)
    
    do i=1,7
      write(130,'(A)')arr(i)
    end do
    
    !RS: Debugging: Hardcoding the following values ONLY to test the E+/HPSim data flow (10/31/15)
    !XMaE=0.4532593
    !WeathOutBaroPress=97900.00
    !TWiC=22.67
    !TaiC=34.72
    !TWiE=19.44
    !TaiE=26.67
    !RS: Debugging: The test showed the inputs and outputs passed correctly through (10/31/15)
    
    !RS: Debugging: Hardcoding the following values ONLY to test that a constant cooling energy/rate is carried through (11/13/15)
    !XMaE=0.4532593
    !WeathOutBaroPress=97900.00
    !TWiC=17.68352
    !TaiC=21.03333
    !TWiE=13.1747
    !TaiE=24.15519
    !RS: Debugging: The test showed the constant cooling/temp was correctly brought in/carried through the code (11/13/15 \ 11/14/15)
    
    write(130,*) VariablesToPass
    write(130,*) XMaE , XMaEL
    !write(1,*) OutWetBulbTemp , OutWetBulbTempL
    write(130,*) WeathOutBaroPress , OutBaroPressATL
    write(130,*) TWiC, TWiCL
    write(130,*) TaiC, TaiCL
    !write(1,*) DummyHROutSideNode, DummyHROutSideNodeL
    write(130,*) TWiE, TWiEL
    write(130,*) TaiE, TaiEL
    !write(1,*) DummyHRMixed, DummyHRMixedL
    close(130)
    
    !newDirPath=TRIM(ADJUSTL(FolderPath(1:pos-1))) // "\\XMLData"
    ! inquire( directory=newDirPath, exist=dirExists )         ! Works with ifort, but not gfortran
    !
    !
    !if (dirExists) then
    !else
    !    mkdirCmd = 'mkdir "'//trim(newDirPath)//'"'
    !    write(*,'(a)') "Creating new directory: '"//trim(mkdirCmd)//"'"
    !    call system( mkdirCmd )
    !endif
    !
    !
    !
    !INQUIRE(FILE=TRIM(ADJUSTL(FolderPath(1:pos-1))) // "\\XMLData\\inputValue.txt", EXIST=file_exists)
    !if (file_exists) THEN
    !    OPEN (UNIT=555, FILE=TRIM(ADJUSTL(FolderPath(1:pos-1))) // "\\XMLData\\inputValue.txt", STATUS="OLD")   ! Current directory
    !    CLOSE (UNIT=555, STATUS="DELETE")
    !end if
    !
    !open(120,file=TRIM(ADJUSTL(FolderPath(1:pos-1))) // "\\XMLData\\inputValue.txt",status='NEW')   
    !write(120,*) TaiC
    !write(120,*) TWiC
    !!write(1,*) DummyHROutSideNode, DummyHROutSideNodeL
    !write(120,*) TaiE
    !write(120,*) TWiE
    !close(120)
    
!IF (ZoneSysEnergyDemand(1)%TotalOutputRequired .NE. 0) THEN !RS: Debugging: Trying to only call optimization and HPSim when there's a load (11/6/14)

!AT : Added in code to run optimization    
    !INQUIRE(FILE=TRIM(ADJUSTL(FolderPath(1:pos-1))) // "\\GenOpt.log", EXIST=file_exists)
    !if (file_exists) THEN
    !    OPEN (UNIT=545, FILE=TRIM(ADJUSTL(FolderPath(1:pos-1))) // "\\GenOpt.log", STATUS="OLD")   ! Current directory
    !    CLOSE (UNIT=545, STATUS="DELETE")
    !end if

    !INQUIRE(FILE=TRIM(ADJUSTL(FolderPath(1:pos-1))) // "\\XMLData\\outputVal.txt", EXIST=file_exists)
        INQUIRE(FILE=TRIM(FolderPath) // "\\GenOpt\\data_Temp_python\\output.txt", EXIST=file_exists)
    if (file_exists) THEN
        !OPEN (UNIT=556, FILE=TRIM(ADJUSTL(FolderPath(1:pos-1))) // "\\XMLData\\outputVal.txt", STATUS="OLD")   ! Current directory
        OPEN (UNIT=556, FILE=TRIM(FolderPath) // "\\GenOpt\\data_Temp_python\\output.txt", STATUS="OLD")   ! Current directory
        CLOSE (UNIT=556, STATUS="DELETE")
    end if
    
    ! INQUIRE(FILE=TRIM(ADJUSTL(FolderPath(1:pos-1))) // "\\XMLData\\exitLoop.txt", EXIST=file_exists)
    !if (file_exists) THEN
    !    OPEN (UNIT=556, FILE=TRIM(ADJUSTL(FolderPath(1:pos-1))) // "\\XMLData\\exitLoop.txt", STATUS="OLD")   ! Current directory
    !    CLOSE (UNIT=556, STATUS="DELETE")
    !end if
    
    !RS: Debugging: Trying to get the code to wait until the exe has finished running before moving on (11/24/14)
    
!    StartupInfo = T_STARTUPINFO(SIZEOF(StartupInfo),&
!  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
!
!! Create a new process to run this executable, passing the 
!! complete executable path and a "-child" switch
!! as the command line (anything that has two or more tokens), 
!! and specifying that it should inherit the handles
!! (including console) of this process.
!!
!!ret = CreateProcess (NULL_CHARACTER,  & ! Application Name
!!	 '"'//path(1:path_len)//'" -child'C, &  ! Command line
!!	 NULL_SECURITY_ATTRIBUTES, &
!!	 NULL_SECURITY_ATTRIBUTES, &
!!	 TRUE, &	! InheritHandles
!!	 0, &	    ! CreationFlags
!!	 NULL, & ! Environment variables
!!	 NULL_CHARACTER, & ! Current directory
!!	 StartupInfo, &
!!	 ProcessInfo)
!     
!ret = CreateProcess (NULL_CHARACTER,  & ! Application Name
!	 !'c:\\windows\\system32\\notepad.exe', &  ! Command line
!     'GenOptExe.exe',& !'C:\Users\Rajesh\Desktop\Coupled HP\\genoptEXE.exe',&
!     !'C:\\Users\\lab303user\\Documents\\betsrg_dual\\Buffer Program\\HeatPumpSimulator_debug.exe', &
!	 NULL_SECURITY_ATTRIBUTES, &
!	 NULL_SECURITY_ATTRIBUTES, &
!	 TRUE, &	! InheritHandles
!	 0, &	    ! CreationFlags
!	 NULL, & ! Environment variables
!	 NULL_CHARACTER, & ! Current directory
!	 StartupInfo, &
!	 ProcessInfo)      
!
!if (ret == 0) then
!  ret = GetLastError ()
!  write (*,'(A,I0)') "Create process failed with error ",ret
!else
!  ! CreateProcess succeeded.  Wait for the process to finish
!  !
!  ret = WaitForSingleObject (ProcessInfo%hProcess, INFINITE)
!
!  ! Close handles, otherwise resources are lost
!  !
!  ret = CloseHandle (ProcessInfo%hThread)
!  ret = CloseHandle (ProcessInfo%hProcess)
!  end if
        
    !CALL System('GenOpt_utility_updated_shared_old.exe')
    CALL System('Genopt_utility_updated_shared_V3.py') !('Genopt_utility_updated_shared.exe')
    HPSimCounter = HPSimCounter + 1 !RS: Debugging: Putting in a counter to determine how often HPSim is called (4/28/18)
    
    !  do
    !     INQUIRE(FILE=TRIM(ADJUSTL(FolderPath(1:pos-1))) // "\\GenOpt.log", EXIST=file_exists)
    !     if (file_exists) exit
    !     INQUIRE(FILE=TRIM(ADJUSTL(FolderPath(1:pos-1))) // "\\XMLData\\outputVal.txt", EXIST=file_exists)
    !     if (file_exists) exit
    !  end do
    !  
    !call date_and_time(values=t)
    !  ms1=(t(5)*3600+t(6)*60+t(7))*1000+t(8)
    !
    !  dt = 1000
    !
    !  do
    !    call date_and_time(values=t)
    !    ms2=(t(5)*3600+t(6)*60+t(7))*1000+t(8)
    !    if(abs(ms2-ms1)>=dt)exit
    !    if (ms2==0) exit
    !  enddo
      
    !  CALL SLEEP@(10)

        !INQUIRE(FILE=TRIM(ADJUSTL(FolderPath(1:pos-1))) // "\\XMLData\\outputVal.txt", EXIST=file_exists)
        INQUIRE(FILE=TRIM(FolderPath)//'\data_Temp_python\output.txt', EXIST=file_exists)
         if (file_exists) THEN
              !open(122,file=TRIM(ADJUSTL(FolderPath(1:pos-1))) // "\\XMLData\\outputVal.txt",status='old')       
              open(122,file=TRIM(FolderPath)//"\data_Temp_python\output.txt",status='old')          
              read(122,*)totalPowerVal
              read(122,*)compressorPowerVal
              read(122,*)volumetricFlowrateVal
              read(122,*)dryBulbTemperatureVal
              read(122,*)enthalpyVal
              read(122,*)humidityRatioVal
              read(122,*)QTotalVal
              read(122,*)QSensVal
              !read(122,*)Subcool_Cond    !RS: Debugging: Keeping track of the condenser subcool for HPSim (5/25/15)
              !read(122,*)Subcool_ExpDev  !RS: Debugging: Keeping track of the expansion device subcool for HPSim (5/25/15)
              !read(122,*)Superheat_Evap  !RS: Debugging: Keeping track of the evaporator superheat for HPSim (5/25/15)
              !read(122,*)Superheat_Comp  !RS: Debugging: Keeping track of the compressor superheat for HPSim (5/25/15
              !read(122,*)TsiCmp  !RS: Debugging: Keeping track of the compressor suction saturation temperature for HPSim (5/25/15)
              !read(122,*)TsoCmp  !RS: Debugging: Keeping track of the compressor discharge saturation temperature for HPSim (5/25/15)
              read(122,*)CompRatio  !RS: Debugging: Keeping track of the compressor ratio for HPSim (9/12/15)
            close(122)
            DXCoil(DXCoilHPSimNum)%OutletAirTemp=(dryBulbTemperatureVal-32)/1.8
            DXCoil(DXCoilHPSimNum)%OutletAirEnthalpy=enthalpyVal*1000
            DXCoil(DXCoilHPSimNum)%OutletAirHumRat=humidityRatioVal
            RhoOutAir=PsyRhoAirFnPbTdbW(OutBaroPressAT,DXCoil(DXCoilHPSimNum)%OutletAirTemp,DXCoil(DXCoilHPSimNum)%OutletAirHumRat)
            DXCoil(DXCoilHPSimNum)%InletAirMassFlowRate=volumetricFlowrateVal*RhoOutAir*0.472E-3   !RS: Debugging: Actually mass flow rate
            DXCoil(DXCoilHPSimNum)%ElecCoolingPower=totalPowerVal
            TotCap=0.2927*QTotalVal
            QSens=0.2927*QSensVal
            !WRITE(151,FMT_101)CurrentEndTime,Subcool_Cond,Subcool_ExpDev,Superheat_Evap,Superheat_Comp,TsiCmp,TsoCmp    !RS: Debugging: Keeping track of the temperatures for HPSim (5/25/15)
            WRITE(152,FMT_101)CurrentEndTime,CompRatio
         end if
      
!     !AT : Added in code to run optimization ---- END
!ELSE    !RS: Debugging: Setting it so that it doesn't crash if there's no required load (11/15/14)
!    TotCap=0
!    Qsens=0
!END IF  !RS: Debugging: Trying to only call optimization and HPSim when there's a load (11/6/14)

  !CALL SimulationCycle(QSens,TotCap,CondInletTemp) !SysLat,TotCap,CondInletTemp)    !RS: SysLat isn't used (1/16/14)
  
  AirMassFlow = DXCoil(DXCoilHPSimNum)%InletAirMassFlowRate !RS: Debugging: Setting this equal immediately after
  !RS: Debugging: I think that the apparatus properties are calculated with entering conditions
  !RS: Debugging: con. The above can only be set here if the flow is constant throughout the HP
  Node(4)%MassFlowRate=AirMassFlow  !RS: Debugging: Attempting to set system air flow to the HPSim airflow (8/27/16)

! Calculate apparatus dew point conditions using TotCap and CBF
  hDelta = TotCap/AirMassFlow
  hADP = InletAirEnthalpy - hDelta/(1.d0-CBF)
  IF (AirMassFlow .GT. 0) THEN  !RS: Debugging: Does the air flow ever get above 0 for this node?
    tADP = PsyTsatFnHPb(hADP,OutdoorPressure,'CalcHPSimDXCoil')
    wADP = PsyWFnTdbH(tADP,hADP,'CalcHPSimDXCoil')
    hTinwADP = PsyHFnTdbW(InletAirDryBulbTemp,wADP,'CalcHPSimDXCoil')
  END IF
  !IF((InletAirEnthalpy-hADP) .NE. 0)THEN
  !  SHR = MIN((hTinwADP-hADP)/(InletAirEnthalpy-hADP),1.d0)
  !ELSE
  !  SHR = 1.0d0
  !END IF
  SHR=QSens/TotCap  !RS: Implementation: This should be the definition of SHR?
!
! Check for dry evaporator conditions (win < wadp)
!
  IF (AirMassFlow .GT. 0) THEN  !RS: Debugging: Does the air flow ever get above 0 for this node?
  IF (wADP .gt. InletAirHumRatTemp .or. (Counter .ge. 1 .and. Counter .lt. MaxIter)) THEN
     If(InletAirHumRatTemp == 0.0)InletAirHumRatTemp=0.00001d0
     werror = (InletAirHumRatTemp - wADP)/InletAirHumRatTemp
!!
!! Increase InletAirHumRatTemp at constant InletAirTemp to find coil dry-out point. Then use the
!! capacity at the dry-out point to determine exiting conditions from coil. This is required
!! since the TotCapTempModFac doesn't work properly with dry-coil conditions.
!!
      InletAirHumRatTemp = RF*wADP + (1.d0-RF)*InletAirHumRatTemp
      InletAirWetbulbC = PsyTwbFnTdbWPb(InletAirDryBulbTemp,InletAirHumRatTemp,OutdoorPressure)
      Counter = Counter + 1
!      IF (ABS(werror) .gt. Tolerance) go to 50   ! Recalculate with modified inlet conditions

  END IF
  END IF

  IF(DXCoil(DXCoilHPSimNum)%PLFFPLR(mode) .GT. 0)THEN
    PLF = CurveValue(DXCoil(DXCoilHPSimNum)%PLFFPLR(mode),PartLoadRatio) ! Calculate part-load factor
  ELSE
    PLF = 1.0d0
  END IF

  IF (PLF < 0.7d0) THEN
    IF (DXCoil(DXCoilHPSimNum)%ErrIndex2 == 0) THEN
      IF(DXCoil(DXCoilHPSimNum)%DXCoilType_Num .EQ. CoilDX_HeatPumpWaterHeater)THEN
        CALL ShowWarningMessage('The PLF curve value for the Heat Pump Water Heater DX coil '//TRIM(DXCoil(DXCoilHPSimNum)%Name)//&
                         ' ='//TRIM(RoundSigDigits(PLF,3))//  &
                         ' for part-load ratio ='//TRIM(RoundSigDigits(PartLoadRatio,3)))
        CALL ShowContinueErrorTimeStamp('PLF curve values must be >= 0.7. PLF has been reset to 0.7 and simulation is continuing.')
        CALL ShowContinueError('Check the IO reference manual for PLF curve guidance [Coil:WaterHeating:AirToWaterHeatPump].')
      ELSE
        CALL ShowWarningMessage('The PLF curve value for the DX cooling coil '//TRIM(DXCoil(DXCoilHPSimNum)%Name)//&
                         ' ='//TRIM(RoundSigDigits(PLF,3))//  &
                         ' for part-load ratio ='//TRIM(RoundSigDigits(PartLoadRatio,3)))
        CALL ShowContinueErrorTimeStamp('PLF curve values must be >= 0.7. PLF has been reset to 0.7 and simulation is continuing.')
        CALL ShowContinueError('Check the IO reference manual for PLF curve guidance [Coil:Cooling:DX:SingleSpeed].')
      END IF
    END IF
    IF(DXCoil(DXCoilHPSimNum)%DXCoilType_Num .EQ. CoilDX_HeatPumpWaterHeater)THEN
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilHPSimNum)%Name)// &
                                            ', Heat Pump Water Heater DX coil PLF curve < 0.7 warning continues...', &
        DXCoil(DXCoilHPSimNum)%ErrIndex2,PLF,PLF)
    ELSE
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilHPSimNum)%Name)// &
                                            ', DX cooling coil PLF curve < 0.7 warning continues...', &
        DXCoil(DXCoilHPSimNum)%ErrIndex2,PLF,PLF)
    END IF
    PLF = 0.7d0
  END IF

  DXCoil(DXCoilHPSimNum)%PartLoadRatio = PartLoadRatio
  DXCoil(DXCoilHPSimNum)%CoolingCoilRuntimeFraction = PartLoadRatio / PLF
  IF (DXCoil(DXCoilHPSimNum)%CoolingCoilRuntimeFraction > 1.0d0 .and.   &
      ABS(DXCoil(DXCoilHPSimNum)%CoolingCoilRuntimeFraction-1.0d0) > .001d0 ) THEN
    IF (DXCoil(DXCoilHPSimNum)%ErrIndex3 == 0) THEN
      IF(DXCoil(DXCoilHPSimNum)%DXCoilType_Num .EQ. CoilDX_HeatPumpWaterHeater)THEN
        CALL ShowWarningMessage('The runtime fraction for Heat Pump Water Heater DX coil '//TRIM(DXCoil(DXCoilHPSimNum)%Name)//&
                           ' exceeded 1.0. ['//TRIM(RoundSigDigits(DXCoil(DXCoilHPSimNum)%CoolingCoilRuntimeFraction,4))//'].')
        CALL ShowContinueError('Runtime fraction reset to 1 and the simulation will continue.')
        CALL ShowContinueError('Check the IO reference manual for PLF curve guidance [Coil:WaterHeating:AirToWaterHeatPump].')
        CALL ShowContinueErrorTimeStamp(' ')
      ELSE
        !CALL ShowWarningMessage('The runtime fraction for DX cooling coil '//TRIM(DXCoil(DXCoilHPSimNum)%Name)//&
        !                   ' exceeded 1.0. ['//TRIM(RoundSigDigits(DXCoil(DXCoilHPSimNum)%CoolingCoilRuntimeFraction,4))//'].')
        !CALL ShowContinueError('Runtime fraction reset to 1 and the simulation will continue.')
        !CALL ShowContinueError('Check the IO reference manual for PLF curve guidance [Coil:Cooling:DX:SingleSpeed].')
        !CALL ShowContinueErrorTimeStamp(' ')   !RS: Secret Search String
        WRITE(DebugFile,*) 'The runtime fraction for DX cooling coil ',TRIM(DXCoil(DXCoilHPSimNum)%Name),&
                           ' exceeded 1.0. [',TRIM(RoundSigDigits(DXCoil(DXCoilHPSimNum)%CoolingCoilRuntimeFraction,4)),'].'
        WRITE(DebugFile,*) 'Runtime fraction reset to 1 and the simulation will continue.'
        WRITE(DebugFile,*) 'Check the IO reference manual for PLF curve guidance [Coil:Cooling:DX:SingleSpeed]'
      END IF
    END IF
    IF(DXCoil(DXCoilHPSimNum)%DXCoilType_Num .EQ. CoilDX_HeatPumpWaterHeater)THEN
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilHPSimNum)%Name)//              &
                   ', Heat Pump Water Heater DX coil runtime fraction > 1.0 warning continues...', &
        DXCoil(DXCoilHPSimNum)%ErrIndex3,DXCoil(DXCoilHPSimNum)%CoolingCoilRuntimeFraction,DXCoil(DXCoilHPSimNum)%CoolingCoilRuntimeFraction)
    ELSE
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilHPSimNum)%Name)//              &
                     ', DX cooling coil runtime fraction > 1.0 warning continues...', &
        DXCoil(DXCoilHPSimNum)%ErrIndex3,DXCoil(DXCoilHPSimNum)%CoolingCoilRuntimeFraction,DXCoil(DXCoilHPSimNum)%CoolingCoilRuntimeFraction)
    END IF
    DXCoil(DXCoilHPSimNum)%CoolingCoilRuntimeFraction = 1.0d0 ! Reset coil runtime fraction to 1.0
  ELSEIF (DXCoil(DXCoilHPSimNum)%CoolingCoilRuntimeFraction > 1.0d0) THEN
    DXCoil(DXCoilHPSimNum)%CoolingCoilRuntimeFraction = 1.0d0 ! Reset coil runtime fraction to 1.0
  END IF

  ! If cycling fan, send coil part-load fraction to on/off fan via HVACDataGlobals
  IF (FanOpMode .EQ. CycFanCycCoil) OnOffFanPartLoadFraction = PLF

  !  Calculate full load output conditions
  IF (SHR .gt. 1. .OR. Counter .gt. 0) SHR = 1.d0
IF (AirMassFlow .GT. 0) THEN    !RS: Debugging: Keeping it from crashing if HPSim isn't called and there's no incoming flow
  FullLoadOutAirEnth = InletAirEnthalpy - TotCap/AirMassFlow
  hTinwout = InletAirEnthalpy - (1.0d0-SHR)*hDelta
  IF (SHR < 1.0d0) THEN
    FullLoadOutAirHumRat = PsyWFnTdbH(InletAirDryBulbTemp,hTinwout)
  ELSE
    FullLoadOutAirHumRat = InletAirHumRat
  END IF
  FullLoadOutAirTemp = PsyTdbFnHW(FullLoadOutAirEnth,FullLoadOutAirHumRat)

! Check for saturation error and modify temperature at constant enthalpy
   IF(FullLoadOutAirTemp .LT. PsyTsatFnHPb(FullLoadOutAirEnth,OutdoorPressure)) THEN
    FullLoadOutAirTemp = PsyTsatFnHPb(FullLoadOutAirEnth,OutdoorPressure)
    FullLoadOutAirHumRat  = PsyWFnTdbH(FullLoadOutAirTemp,FullLoadOutAirEnth)
   END IF

  ! Store actual outlet conditions when DX coil is ON for use in heat recovery module
  DXCoilFullLoadOutAirTemp(DXCoilHPSimNum) = FullLoadOutAirTemp
  DXCoilFullLoadOutAirHumRat(DXCoilHPSimNum) = FullLoadOutAirHumRat
ELSE    !RS: Debugging: Giving another case so it is set
    FullLoadOutAirTemp=DXCoil(DXCoilHPSimNum)%OutletAirTemp
    !OutletAirEnthalpy=DXCoil(DXCoilHPSimNum)%OutletAirEnthalpy
END IF

! Add warning message for cold cooling coil (FullLoadOutAirTemp < 2 C)
  IF(FullLoadOutAirTemp .LT. 2.0d0 .AND. .NOT. FirstHVACIteration .AND. .NOT. WarmupFlag)THEN
    DXCoil(DXCoilHPSimNum)%PrintLowOutTempMessage = .TRUE.
    DXCoil(DXCoilHPSimNum)%FullLoadOutAirTempLast = FullLoadOutAirTemp
    IF(DXCoil(DXCoilHPSimNum)%LowOutletTempIndex == 0)THEN
      DXCoil(DXCoilHPSimNum)%FullLoadInletAirTempLast = InletAirDryBulbTemp
      DXCoil(DXCoilHPSimNum)%LowOutTempBuffer1= TRIM(DXCoil(DXCoilHPSimNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilHPSimNum)%Name)// &
       '" - Full load outlet air dry-bulb temperature < 2C. This indicates the possibility of coil frost/freeze.'// &
       ' Outlet temperature = '//TRIM(RoundSigDigits(FullLoadOutAirTemp,2))//' C.'
      DXCoil(DXCoilHPSimNum)%LowOutTempBuffer2 = ' '//'...Occurrence info = '//TRIM(EnvironmentName)//', '//Trim(CurMnDy)//' '&
                   //TRIM(CreateSysTimeIntervalString())
    END IF
  END IF

  !  If constant fan with cycling compressor, call function to determine "effective SHR"
  !  which includes the part-load degradation on latent capacity
  IF (FanOpMode .EQ. ContFanCycCoil) THEN
     QLatRated = DXCoil(DXCoilHPSimNum)%RatedTotCap(Mode) * (1.d0 - DXCoil(DXCoilHPSimNum)%RatedSHR(Mode))
     QLatActual = TotCap * (1.d0 - SHR)
     SHRUnadjusted = SHR
     SHR = CalcEffectiveSHR(DXCoilHPSimNum, SHR, DXCoil(DXCoilHPSimNum)%CoolingCoilRuntimeFraction, &
                          QLatRated, QLatActual, InletAirDryBulbTemp, InletAirWetbulbC, Mode)
    ! For multimode coil, if stage-2 operation (modes 2 or 4), adjust Stage1&2 SHR to account for
    ! Stage 1 operating at full load, so there is no degradation for that portion
    ! Use the stage 1 bypass fraction to allocate
    IF (Mode .EQ. 2 .OR. Mode .EQ. 4) THEN
      SHR = SHRUnadjusted*(1.d0-DXCoil(DXCoilHPSimNum)%BypassedFlowFrac(Mode-1))+ &
            SHR*DXCoil(DXCoilHPSimNum)%BypassedFlowFrac(Mode-1)
    END IF

  !  Calculate full load output conditions
    IF (SHR .gt. 1.d0 .OR. Counter .gt. 0) SHR = 1.d0
    FullLoadOutAirEnth = InletAirEnthalpy - TotCap/AirMassFlow
    hTinwout = InletAirEnthalpy - (1.0d0-SHR)*hDelta
    IF (SHR < 1.0d0) THEN
      FullLoadOutAirHumRat = PsyWFnTdbH(InletAirDryBulbTemp,hTinwout)
    ELSE
      FullLoadOutAirHumRat = InletAirHumRat
    END IF
    FullLoadOutAirTemp = PsyTdbFnHW(FullLoadOutAirEnth,FullLoadOutAirHumRat)

! apply latent degradation model to cycling fan when RH control is desired and heating coil operates
! longer than the cooling coil. DXcoolToHeatPLRRatio = Cooling coil PLR / Heating coil PLR.
  ELSE IF(FanOpMode .EQ. CycFanCycCoil) THEN
    IF(DXcoolToHeatPLRRatio .LT. 1.0d0)THEN
      QLatRated = DXCoil(DXCoilHPSimNum)%RatedTotCap(Mode) * (1.d0 - DXCoil(DXCoilHPSimNum)%RatedSHR(Mode))
      QLatActual = TotCap * (1.d0 - SHR)
      HeatRTF = PartLoadRatio/DXcoolToHeatPLRRatio
      IF(DXCoil(DXCoilHPSimNum)%HeatingCoilPLFCurvePTR .GT. 0)THEN
        HeatingCoilPLF = CurveValue(DXCoil(DXCoilHPSimNum)%HeatingCoilPLFCurvePTR,HeatRTF)
        IF(HeatingCoilPLF .GT. 0) HeatRTF = HeatRTF/HeatingCoilPLF
      END IF
      SHRUnadjusted = SHR
      SHR = CalcEffectiveSHR(DXCoilHPSimNum, SHR, DXCoil(DXCoilHPSimNum)%CoolingCoilRuntimeFraction, &
                     QLatRated, QLatActual, InletAirDryBulbTemp, InletAirWetbulbC, Mode, HeatRTF)
  !   Calculate full load output conditions
      IF (SHR .gt. 1.d0 .OR. Counter .gt. 0) SHR = 1.d0
      FullLoadOutAirEnth = InletAirEnthalpy - TotCap/AirMassFlow
      hTinwout = InletAirEnthalpy - (1.0d0-SHR)*hDelta
      IF (SHR < 1.0d0) THEN
        FullLoadOutAirHumRat = PsyWFnTdbH(InletAirDryBulbTemp,hTinwout)
      ELSE
        FullLoadOutAirHumRat = InletAirHumRat
      END IF
      FullLoadOutAirTemp = PsyTdbFnHW(FullLoadOutAirEnth,FullLoadOutAirHumRat)
    END IF

  END IF

!  Calculate actual outlet conditions for the input part load ratio
!  Actual outlet conditions are "average" for time step

! For multimode coil, if stage-2 operation (modes 2 or 4), return "full load" outlet conditions
  IF ((FanOpMode .EQ. ContFanCycCoil) .AND. &
      (Mode .EQ. 1) .OR. (Mode .EQ. 3)) THEN
    ! Continuous fan, cycling compressor
    OutletAirEnthalpy = ((PartLoadRatio * AirFlowRatio)*DXCoil(DXCoilHPSimNum)%OutletAirEnthalpy + & !FullLoadOutAirEnth + &
                                            (1.d0-(PartLoadRatio * AirFlowRatio))*InletAirEnthalpy)
    OutletAirHumRat = ((PartLoadRatio * AirFlowRatio)*DXCoil(DXCoilHPSimNum)%OutletAirHumRat + & !FullLoadOutAirHumRat + &
                                            (1.d0-(PartLoadRatio * AirFlowRatio))*InletAirHumRat)
    OutletAirTemp = PsyTdbFnHW(OutletAirEnthalpy,OutletAirHumRat)
    Qsens=Qsens !*PartLoadRatio
  ELSE
    ! Default to cycling fan, cycling compressor
    ! Also return this result for stage 2 operation of multimode coil
    ! Cycling fan typically provides full outlet conditions. When RH control is used, account for additional
    ! heating run time by using cooling/heating ratio the same as constant fan (otherwise PLRRatio = 1).
    OutletAirEnthalpy = FullLoadOutAirEnth * DXcoolToHeatPLRRatio + InletAirEnthalpy * (1.0d0 - DXcoolToHeatPLRRatio)
    OutletAirHumRat = FullLoadOutAirHumRat * DXcoolToHeatPLRRatio + InletAirHumRat * (1.0d0 - DXcoolToHeatPLRRatio)
    OutletAirTemp = FullLoadOutAirTemp * DXcoolToHeatPLRRatio + InletAirDryBulbTemp * (1.0d0 - DXcoolToHeatPLRRatio)
  END IF  
  
  !RS: Implementation: Setting these here because otherwise they're not getting updated
  !OutletAirTemp=DXCoil(DXCoilHPSimNum)%OutletAirTemp
  !OutletAirEnthalpy=DXCoil(DXCoilHPSimNum)%OutletAirEnthalpy

! Check for saturation error and modify temperature at constant enthalpy
   IF(OutletAirTemp .LT. PsyTsatFnHPb(OutletAirEnthalpy,OutdoorPressure,'CalcHPSimDXCoil')) THEN
    OutletAirTemp = PsyTsatFnHPb(OutletAirEnthalpy,OutdoorPressure)
    OutletAirHumRat  = PsyWFnTdbH(OutletAirTemp,OutletAirEnthalpy)
   END IF

! Mix with air that was bypassed around coil, if any
   IF(BypassFlowFraction .GT. 0.0) THEN
     OutletAirEnthalpy = (1.d0-BypassFlowFraction)*OutletAirEnthalpy + BypassFlowFraction*InletAirEnthalpy
     OutletAirHumRat = (1.d0-BypassFlowFraction)*OutletAirHumRat + BypassFlowFraction*InletAirHumRat
     OutletAirTemp = PsyTdbFnHW(OutletAirEnthalpy,OutletAirHumRat)
     ! Check for saturation error and modify temperature at constant enthalpy
     IF(OutletAirTemp .LT. PsyTsatFnHPb(OutletAirEnthalpy,OutdoorPressure)) THEN
       OutletAirTemp = PsyTsatFnHPb(OutletAirEnthalpy,OutdoorPressure)
       OutletAirHumRat  = PsyWFnTdbH(OutletAirTemp,OutletAirEnthalpy)
     END IF
    END IF

! Calculate electricity consumed. First, get EIR modifying factors for off-rated conditions
  IF(DXCoil(DXCoilHPSimNum)%DXCoilType_Num .EQ. CoilDX_HeatPumpWaterHeater) THEN   !RS: Debugging: Pulling out the electricity for now
!!   Coil:DX:HeatPumpWaterHeater does not have EIR temp or flow curves
    EIRTempModFac = 1.0d0
    EIRFlowModFac = 1.0d0
  ELSE
!    EIRTempModFac = CurveValue(DXCoil(DXCoilHPSimNum)%EIRFTemp(Mode),InletAirWetbulbC,CondInletTemp)
!
!!   Warn user if curve output goes negative
!    IF(EIRTempModFac .LT. 0.0)THEN
!      IF(DXCoil(DXCoilHPSimNum)%EIRFTempErrorIndex == 0)THEN
!        CALL ShowWarningMessage(TRIM(DXCoil(DXCoilHPSimNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilHPSimNum)%Name)//'":')
!        CALL ShowContinueError(' Energy Input Ratio Modifier curve (function of temperature) output is negative (' &
!                          //TRIM(TrimSigDigits(EIRTempModFac,3))//').')
!        IF(DXCoil(DXCoilHPSimNum)%EIRTempModFacCurveType(Mode) .EQ. Biquadratic)THEN
!          CALL ShowContinueError(' Negative value occurs using a condenser inlet air temperature of ' &
!                         //TRIM(TrimSigDigits(CondInletTemp,1))// &
!                               ' and an inlet air wet-bulb temperature of '//TRIM(TrimSigDigits(InletAirWetbulbC,1))//'.')
!        ELSE
!          CALL ShowContinueError(' Negative value occurs using a condenser inlet air temperature of ' &
!                         //TRIM(TrimSigDigits(CondInletTemp,1))//'.')
!        END IF
!        IF(Mode .GT. 1)THEN
!          Call ShowContinueError(' Negative output results from stage '//TRIM(TrimSigDigits(Mode))// &
!                                 ' compressor operation.')
!        END IF
!        CALL ShowContinueErrorTimeStamp(' Resetting curve output to zero and continuing simulation.')
!      END IF
!      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilHPSimNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilHPSimNum)%Name)//'":'//&
!          ' Energy Input Ratio Modifier curve (function of temperature) output is negative warning continues...' &
!          , DXCoil(DXCoilHPSimNum)%EIRFTempErrorIndex, EIRTempModFac, EIRTempModFac)
!      EIRTempModFac = 0.0
!    END IF
!
!    EIRFlowModFac = CurveValue(DXCoil(DXCoilHPSimNum)%EIRFFlow(Mode),AirMassFlowRatio)
!
!!   Warn user if curve output goes negative
!    IF(EIRFlowModFac .LT. 0.0)THEN
!      IF(DXCoil(DXCoilHPSimNum)%EIRFFlowErrorIndex == 0)THEN
!        CALL ShowWarningMessage(TRIM(DXCoil(DXCoilHPSimNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilHPSimNum)%Name)//'":')
!        CALL ShowContinueError(' Energy Input Ratio Modifier curve (function of flow fraction) output is negative (' &
!                          //TRIM(TrimSigDigits(EIRFlowModFac,3))//').')
!        CALL ShowContinueError(' Negative value occurs using an air flow fraction of ' &
!                           //TRIM(TrimSigDigits(AirMassFlowRatio,3))//'.')
!        CALL ShowContinueErrorTimeStamp(' Resetting curve output to zero and continuing simulation.')
!        IF(Mode .GT. 1)THEN
!          Call ShowContinueError(' Negative output results from stage '//TRIM(TrimSigDigits(Mode))// &
!                                 ' compressor operation.')
!        END IF
!      END IF
!      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilHPSimNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilHPSimNum)%Name)//'":'//&
!         ' Energy Input Ratio Modifier curve (function of flow fraction) output is negative warning continues...' &
!         , DXCoil(DXCoilHPSimNum)%EIRFFlowErrorIndex, EIRFlowModFac, EIRFlowModFac)
!      EIRFlowModFac = 0.0
!    END IF
  END IF
!
  EIR = DXCoil(DXCoilHPSimNum)%RatedEIR(Mode) * EIRFlowModFac * EIRTempModFac

! For multimode coil, if stage-2 operation (Modes 2 or 4), return "full load" power adjusted for PLF
  !IF (Mode .EQ. 1 .OR. Mode .EQ. 3) THEN
  !  DXCoil(DXCoilHPSimNum)%ElecCoolingPower = TotCap * EIR * DXCoil(DXCoilHPSimNum)%CoolingCoilRuntimeFraction
  !ELSE
  !  DXCoil(DXCoilHPSimNum)%ElecCoolingPower = TotCap * EIR * DXCoil(DXCoilHPSimNum)%CoolingCoilRuntimeFraction / PartLoadRatio
  !END IF

! Reset AirMassFlow to inlet node air mass flow for final total, sensible and latent calculations
! since AirMassFlow might have been modified above (in this subroutine):
!     IF (FanOpMode .EQ. CycFanCycCoil) AirMassFlow = AirMassFlow / PartLoadRatio
!
! For multimode coil, this should be full flow including bypassed fraction
  AirMassFlow = DXCoil(DXCoilHPSimNum)%InletAirMassFlowRate
  DXCoil(DXCoilHPSimNum)%TotalCoolingEnergyRate = TotCap*(-LoadToCoolingSetPoint_HPSim/TotCap)  !RS:Debugging: Using this for part load fraction (3/16/19) !RS: Debugging: Which?

! Set DataHeatGlobal heat reclaim variable for use by heat reclaim coil (part load ratio is accounted for)
! Calculation for heat reclaim needs to be corrected to use compressor power (not including condenser fan power)
  HeatReclaimDXCoil(DXCoilHPSimNum)%AvailCapacity = DXCoil(DXCoilHPSimNum)%TotalCoolingEnergyRate + DXCoil(DXCoilHPSimNum)%ElecCoolingPower

  !MinAirHumRat = MIN(InletAirHumRat,OutletAirHumRat)   !RS: Debugging: Adding a "-" below so that it has the same convention (8/15/14)
  DXCoil(DXCoilHPSimNum)%SensCoolingEnergyRate = QSens !AirMassFlow * &
  !                                       (PsyHFnTdbW(InletAirDryBulbTemp,MinAirHumRat) - &
  !                                        PsyHFnTdbW(OutletAirTemp,MinAirHumRat))
!  Don't let sensible capacity be greater than total capacity
  !IF (DXCoil(DXCoilHPSimNum)%SensCoolingEnergyRate .GT. DXCoil(DXCoilHPSimNum)%TotalCoolingEnergyRate) THEN
  !   DXCoil(DXCoilHPSimNum)%SensCoolingEnergyRate = DXCoil(DXCoilHPSimNum)%TotalCoolingEnergyRate
  !END IF
  !RS: Debugging: Trying to handle the sign conventions for when cooling is "-" (8/18/14)
  IF (ABS(DXCoil(DXCoilHPSimNum)%SensCoolingEnergyRate) .GT. ABS(DXCoil(DXCoilHPSimNum)%TotalCoolingEnergyRate)) THEN
     DXCoil(DXCoilHPSimNum)%SensCoolingEnergyRate = DXCoil(DXCoilHPSimNum)%TotalCoolingEnergyRate
  END IF
!
  DXCoil(DXCoilHPSimNum)%LatCoolingEnergyRate = DXCoil(DXCoilHPSimNum)%TotalCoolingEnergyRate - DXCoil(DXCoilHPSimNum)%SensCoolingEnergyRate

! Calculate crankcase heater power using the runtime fraction for this DX cooling coil only if there is no companion DX coil.
! Else use the largest runtime fraction of this DX cooling coil and the companion DX heating coil.
  !IF(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil .EQ. 0) THEN
    DXCoil(DXCoilHPSimNum)%CrankcaseHeaterPower = CrankcaseHeatingPower * &
                                           (1.0d0 - DXCoil(DXCoilHPSimNum)%CoolingCoilRuntimeFraction)
  !ELSE
  !  DXCoil(DXCoilNum)%CrankcaseHeaterPower = CrankcaseHeatingPower * &
  !                                         (1.0d0 - MAX(DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction, &
  !                                         DXCoil(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil)%HeatingCoilRuntimeFraction))
  !END IF

  !IF (DXCoil(DXCoilNum)%CondenserType(Mode) == EvapCooled) THEN
  !!******************
  !!             WATER CONSUMPTION IN m3 OF WATER FOR DIRECT
  !!             H2O [m3/sec] = Delta W[KgH2O/Kg air]*Mass Flow Air[Kg air]
  !!                                /RhoWater [kg H2O/m3 H2O]
  !!******************
  !   RhoWater = RhoH2O(OutdoorDryBulb)
  !   DXCoil(DXCoilNum)%EvapWaterConsumpRate =  &
  !           (CondInletHumrat - OutdoorHumRat) *  &
  !            CondAirMassFlow/RhoWater * DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction
  !   DXCoil(DXCoilNum)%EvapCondPumpElecPower = DXCoil(DXCoilNum)%EvapCondPumpElecNomPower(Mode) * &
  !                                             DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction
  !! Calculate basin heater power
  !  CALL CalcBasinHeaterPower(DXCoil(DXCoilNum)%BasinHeaterPowerFTempDiff,&
  !                            DXCoil(DXCoilNum)%BasinHeaterSchedulePtr,&
  !                            DXCoil(DXCoilNum)%BasinHeaterSetPointTemp,DXCoil(DXCoilNum)%BasinHeaterPower)
  !  IF (DXCoil(DXCoilNum)%DXCoilType_Num == CoilDX_CoolingSingleSpeed) THEN
  !    DXCoil(DXCoilNum)%BasinHeaterPower = DXCoil(DXCoilNum)%BasinHeaterPower * &
  !                                       (1.d0 - DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction)
  !  ENDIF
  !END IF
    !RS: Debugging: Removing because they're already set and the recalculation is slightly off
  !DXCoil(DXCoilHPSimNum)%OutletAirTemp     = OutletAirTemp
  !DXCoil(DXCoilHPSimNum)%OutletAirHumRat   = OutletAirHumRat   !RS: Debugging: Unnecessary and if OutletAirHumRat not updated, a bad idea
  !DXCoil(DXCoilHPSimNum)%OutletAirEnthalpy = OutletAirEnthalpy

DXCoolingSystem(DXCoilHPSimNum)%PartLoadFrac    = (LoadToCoolingSetPoint_HPSim/TotCap)  !RS:Debugging: Seeing if this calculates the part load fraction (3/16/19)

ELSE

  ! DX coil is off; just pass through conditions
  DXCoil(DXCoilHPSimNum)%OutletAirEnthalpy = DXCoil(DXCoilHPSimNum)%InletAirEnthalpy
  DXCoil(DXCoilHPSimNum)%OutletAirHumRat   = DXCoil(DXCoilHPSimNum)%InletAirHumRat
  DXCoil(DXCoilHPSimNum)%OutletAirTemp     = DXCoil(DXCoilHPSimNum)%InletAirTemp

  DXCoil(DXCoilHPSimNum)%ElecCoolingPower = 0.0
  DXCoil(DXCoilHPSimNum)%TotalCoolingEnergyRate = 0.0
  DXCoil(DXCoilHPSimNum)%SensCoolingEnergyRate = 0.0
  DXCoil(DXCoilHPSimNum)%LatCoolingEnergyRate = 0.0
  DXCoil(DXCoilHPSimNum)%EvapCondPumpElecPower = 0.0
  DXCoil(DXCoilHPSimNum)%EvapWaterConsumpRate = 0.0

! Reset globals when DX coil is OFF for use in heat recovery module
  DXCoilFullLoadOutAirTemp(DXCoilHPSimNum) = 0.0
  DXCoilFullLoadOutAirHumRat(DXCoilHPSimNum) = 0.0

END IF ! end of on/off if - else

!set water system demand request (if needed)
IF ( DXCoil(DXCoilHPSimNum)%EvapWaterSupplyMode == WaterSupplyFromTank) THEN
   WaterStorage(DXCoil(DXCoilHPSimNum)%EvapWaterSupTankID)%VdotRequestDemand(DXCoil(DXCoilHPSimNum)%EvapWaterTankDemandARRID) &
       = DXCoil(DXCoilHPSimNum)%EvapWaterConsumpRate
ENDIF

DXCoilOutletTemp(DXCoilHPSimNum)         = DXCoil(DXCoilHPSimNum)%OutletAirTemp
DXCoilOutletHumRat(DXCoilHPSimNum)       = DXCoil(DXCoilHPSimNum)%OutletAirHumRat
DXCoilPartLoadRatio(DXCoilHPSimNum)      = DXCoil(DXCoilHPSimNum)%PartLoadRatio
DXCoilFanOpMode(DXCoilHPSimNum)          = FanOpMode
DXCoil(DXCoilHPSimNum)%CondInletTemp     = CondInletTemp

!RS: Debugging: Trying to have HPSim only called when it needs to be (3/17/18)
OutTemp_Previous = OutTemp_Current
InDTemp_Previous = InDTemp_Current
InWTemp_Previous = InWTemp_Current

!CALL HPSimPowerUpdate(DXCoilHPSimNum)   !RS: Debugging: Trying to pass through power from HPSim (9/3/14)

!ERRORFLAG=.FALSE.
!OutletNode=2
!InletNode=5 !4
!       CALL UpdateHVACInterface(OutletNode,InletNode,ERRORFLAG) !RS: Testing: ??? (8/8/14)
RETURN
END SUBROUTINE CalcHPSimDXCoil

SUBROUTINE CalcVRFCoolingCoil(DXCoilNum,CompOp,FirstHVACIteration,PartLoadRatio,FanOpMode,CompCycRatio, &
                              PerfMode,OnOffAirFlowRatio, MaxCoolCap)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Richard Raustad
          !       DATE WRITTEN   August 2010
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! Calculates the air-side performance of a direct-expansion, air-cooled
          ! VRF terminal unit cooling coil.
          ! A new subroutine was created in case this DX coil model is significantly
          ! different from the existing CalcDoe2DXCoil subroutine. The VRF heating coil
          ! uses the existing DX heating coil subroutine (CalcDXHeatingCoil).

          ! METHODOLOGY EMPLOYED:
          ! This routine simulates the performance of a variable refrigerant flow cooling coil.
          ! The routine requires the user to enter the total cooling capacity and sensible heat ratio.
          ! Since different manufacturer's rate their equipment at different air flow rates,
          ! the supply air flow rate corresponding to the rated capacities must also be
          ! entered (should be between 300 cfm/ton and 450 cfm/ton). The rated information entered by
          ! the user should NOT include the thermal or electrical impacts of the supply air fan, as
          ! this is addressed by another module.

          ! With the rated performance data entered by the user, the model employs some of the
          ! DOE-2.1E curve fits to adjust the capacity and efficiency of the unit as a function
          ! of entering air temperatures and supply air flow rate (actual vs rated flow). The model
          ! does NOT employ the exact same methodology to calculate performance as DOE-2.
          !
          ! This VRF cooling coil model adjusts the rated total cooling capacity by the CAPFT
          ! and CAP funciton of flow curve/model currently used by the existing DX coil model.
          ! The part-load ratio is then applied to the total operating capacity to find the capacity
          ! required to meet the load. This VRF model then uses the ADP/bypass method to find the
          ! SHR and resulting outlet conditions given that total capacity (or delta H).

          ! The model checks for coil dryout conditions, and adjusts the calculated performance
          ! appropriately.

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE CurveManager,    ONLY: CurveValue
  USE DataGlobals,     ONLY: CurrentTime
  USE DataHVACGlobals, ONLY: HPWHCrankcaseDBTemp, TimeStepSys, SysTimeElapsed
  USE General,         ONLY: TrimSigDigits, RoundSigDigits, CreateSysTimeIntervalString
  USE DataWater,       ONLY: WaterStorage

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  INTEGER,   INTENT(IN)           :: DXCoilNum           ! the number of the DX coil to be simulated
  INTEGER,   INTENT(IN)           :: CompOp              ! compressor operation; 1=on, 0=off
  LOGICAL,   INTENT(IN)           :: FirstHVACIteration  ! true if this is the first iteration of HVAC
  REAL(r64), INTENT(IN)           :: PartLoadRatio       ! sensible cooling load / full load sensible cooling capacity
  INTEGER,   INTENT(IN)           :: FanOpMode           ! Allows parent object to control fan operation
  REAL(r64), INTENT(IN)           :: CompCycRatio        ! cycling ratio of VRF condenser
  INTEGER,   INTENT(IN), OPTIONAL :: PerfMode            ! Performance mode for MultiMode DX coil; Always 1 for other coil types
  REAL(r64), INTENT(IN), OPTIONAL :: OnOffAirFlowRatio   ! ratio of compressor on airflow to compressor off airflow
!  REAL(r64), INTENT(IN), OPTIONAL :: CoolingHeatingPLR   ! used for cycling fan RH control
 REAL(r64), INTENT(IN), OPTIONAL :: MaxCoolCap           ! maximum capacity of DX coil

          ! SUBROUTINE PARAMETER DEFINITIONS:
CHARACTER(len=*), PARAMETER :: RoutineName='CalcVRFCoolingCoil'

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
REAL(r64) :: AirMassFlow           ! dry air mass flow rate through coil [kg/s] (adjusted for bypass if any)
REAL(r64) :: AirMassFlowRatio      ! Ratio of actual air mass flow to rated air mass flow (adjusted for bypass if any)
REAL(r64) :: AirVolumeFlowRate     ! Air volume flow rate across the cooling coil [m3/s] (adjusted for bypass if any)
                                   ! (average flow if cycling fan, full flow if constant fan)
REAL(r64) :: VolFlowperRatedTotCap ! Air volume flow rate divided by rated total cooling capacity [m3/s-W] (adjusted for bypass)
REAL(r64) :: TotCap                ! gross total cooling capacity at off-rated conditions [W]
REAL(r64) :: TotCapTempModFac      ! Total capacity modifier (function of entering wetbulb, outside drybulb)
REAL(r64) :: TotCapFlowModFac      ! Total capacity modifier (function of actual supply air flow vs rated flow)
REAL(r64) :: InletAirWetBulbC      ! wetbulb temperature of inlet air [C]
REAL(r64) :: InletAirDryBulbTemp   ! inlet air dry bulb temperature [C]
REAL(r64) :: InletAirEnthalpy      ! inlet air enthalpy [J/kg]
REAL(r64) :: InletAirHumRat        ! inlet air humidity ratio [kg/kg]
REAL(r64) :: InletAirHumRatTemp    ! inlet air humidity ratio used in ADP/BF loop [kg/kg]
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!REAL(r64) :: InletAirPressure      ! inlet air pressure [Pa]
REAL(r64) :: RatedCBF              ! coil bypass factor at rated conditions
REAL(r64) :: SHR                   ! Sensible Heat Ratio (sensible/total) of the cooling coil
REAL(r64) :: CBF                   ! coil bypass factor at off rated conditions
REAL(r64) :: A0                    ! NTU * air mass flow rate, used in CBF calculation
REAL(r64) :: hDelta                ! Change in air enthalpy across the cooling coil [J/kg]
REAL(r64) :: hADP                  ! Apparatus dew point enthalpy [J/kg]
REAL(r64) :: hTinwADP              ! Enthalpy at inlet dry-bulb and wADP [J/kg]
REAL(r64) :: hTinwout              ! Enthalpy at inlet dry-bulb and outlet humidity ratio [J/kg]
REAL(r64) :: tADP                  ! Apparatus dew point temperature [C]
REAL(r64) :: wADP                  ! Apparatus dew point humidity ratio [kg/kg]
REAL(r64) :: FullLoadOutAirEnth    ! outlet full load enthalpy [J/kg]
REAL(r64) :: FullLoadOutAirHumRat  ! outlet humidity ratio at full load
REAL(r64) :: FullLoadOutAirTemp    ! outlet air temperature at full load [C]
REAL(r64) :: PLF                   ! Part load factor, accounts for thermal lag at compressor startup, used in power calculation
REAL(r64) :: QLatActual            ! operating latent capacity of DX coil
REAL(r64) :: QLatRated             ! Rated latent capacity of DX coil
REAL(r64) :: SHRUnadjusted         ! SHR prior to latent degradation effective SHR calculation
INTEGER :: Counter                 ! Counter for dry evaporator iterations
INTEGER :: MaxIter                 ! Maximum number of iterations for dry evaporator calculations
REAL(r64) :: RF                    ! Relaxation factor for dry evaporator iterations
REAL(r64) :: Tolerance             ! Error tolerance for dry evaporator iterations
REAL(r64) :: werror                ! Deviation of humidity ratio in dry evaporator iteration loop
REAL(r64) :: CondInletTemp         ! Condenser inlet temperature (C). Outdoor dry-bulb temp for air-cooled condenser.
                                 ! Outdoor Wetbulb +(1 - effectiveness)*(outdoor drybulb - outdoor wetbulb) for evap condenser.
REAL(r64) :: CondInletHumrat       ! Condenser inlet humidity ratio (kg/kg). Zero for air-cooled condenser.
                                 ! For evap condenser, its the humidity ratio of the air leaving the evap cooling pads.
REAL(r64) :: CondAirMassFlow       ! Condenser air mass flow rate [kg/s]
REAL(r64) :: RhoAir                ! Density of air [kg/m3]
REAL(r64) :: CrankcaseHeatingPower ! power due to crankcase heater
REAL(r64) :: CompAmbTemp = 0.0     ! Ambient temperature at compressor
REAL(r64) :: AirFlowRatio          ! ratio of compressor on airflow to average timestep airflow
                                 ! used when constant fan mode yields different air flow rates when compressor is ON and OFF
                                 ! (e.g. Packaged Terminal Heat Pump)
REAL(r64) :: OutdoorDryBulb        ! Outdoor dry-bulb temperature at condenser (C)
REAL(r64) :: OutdoorWetBulb        ! Outdoor wet-bulb temperature at condenser (C)
REAL(r64) :: OutdoorHumRat         ! Outdoor humidity ratio at condenser (kg/kg)
REAL(r64) :: OutdoorPressure       ! Outdoor barometric pressure at condenser (Pa)

REAL(r64) :: CurrentEndTime = 0.0  ! end time of time step for current simulation time step
REAL(r64) :: MinAirHumRat = 0.0    ! minimum of the inlet air humidity ratio and the outlet air humidity ratio
INTEGER          :: Mode           ! Performance mode for Multimode DX coil; Always 1 for other coil types
REAL(r64) :: OutletAirTemp           ! Supply air temperature (average value if constant fan, full output if cycling fan)
REAL(r64) :: OutletAirHumRat         ! Supply air humidity ratio (average value if constant fan, full output if cycling fan)
REAL(r64) :: OutletAirEnthalpy       ! Supply air enthalpy (average value if constant fan, full output if cycling fan)
REAL(r64) :: Adiff                   ! Used for exponential

! If Performance mode not present, then set to 1.  Used only by Multimode/Multispeed DX coil (otherwise mode = 1)
IF (PRESENT(PerfMode)) THEN
  Mode = PerfMode
ELSE
  Mode = 1
END IF

! If AirFlowRatio not present, then set to 1. Used only by DX coils with different air flow
! during cooling and when no cooling is required (constant fan, fan speed changes)
IF (PRESENT(OnOffAirFlowRatio)) THEN
  AirFlowRatio = OnOffAirFlowRatio
ELSE
  AirFlowRatio = 1.0d0
END IF

MaxIter         = 30
RF              = 0.4d0
Counter         = 0
Tolerance       = 0.01d0
CondInletTemp   = 0.0d0
CondInletHumrat = 0.0d0
AirMassFlow         = DXCoil(DXCoilNum)%InletAirMassFlowRate
InletAirDryBulbTemp = DXCoil(DXCoilNum)%InletAirTemp
InletAirEnthalpy    = DXCoil(DXCoilNum)%InletAirEnthalpy
InletAirHumRat      = DXCoil(DXCoilNum)%InletAirHumRat
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!InletAirPressure    = DXCoil(DXCoilNum)%InletAirPressure
HeatReclaimDXCoil(DXCoilNum)%AvailCapacity   = 0.0d0
DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction = 0.0d0
DXCoil(DXCoilNum)%PartLoadRatio              = 0.0d0
DXCoil(DXCoilNum)%BasinHeaterPower           = 0.0d0

IF (DXCoil(DXCoilNum)%CondenserInletNodeNum(Mode) /= 0) THEN
  OutdoorDryBulb  = Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(Mode))%Temp
  OutdoorHumRat   = Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(Mode))%HumRat
  OutdoorPressure = Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(Mode))%Press
  OutdoorWetBulb  = Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(Mode))%OutAirWetBulb
ELSE
  OutdoorDryBulb  = OutDryBulbTemp
  OutdoorHumRat   = OutHumRat
  OutdoorPressure = OutBaroPress
  OutdoorWetBulb  = OutWetBulbTemp
ENDIF


IF (DXCoil(DXCoilNum)%CondenserType(Mode) == AirCooled) THEN
  CondInletTemp   = OutdoorDryBulb ! Outdoor dry-bulb temp
  CompAmbTemp     = OutdoorDryBulb
ELSEIF (DXCoil(DXCoilNum)%CondenserType(Mode) == EvapCooled) THEN
  RhoAir          = PsyRhoAirFnPbTdbW(OutdoorPressure,OutdoorDryBulb,OutdoorHumRat)
  CondAirMassFlow =  RhoAir * DXCoil(DXCoilNum)%EvapCondAirFlow(Mode)
 ! (Outdoor wet-bulb temp from DataEnvironment) + (1.0-EvapCondEffectiveness) * (drybulb - wetbulb)
  CondInletTemp   = OutdoorWetBulb + (OutdoorDryBulb-OutdoorWetBulb)*(1.0d0 - DXCoil(DXCoilNum)%EvapCondEffect(Mode))
  CondInletHumrat = PsyWFnTdbTwbPb(CondInletTemp,OutdoorWetBulb,OutdoorPressure)
  CompAmbTemp     = CondInletTemp
END IF

! Initialize crankcase heater, operates below OAT defined in input deck for HP DX cooling coil
! If used in a heat pump, the value of MaxOAT in the heating coil overrides that in the cooling coil (in GetInput)
IF (CompAmbTemp .LT. DXCoil(DXCoilNum)%MaxOATCrankcaseHeater)THEN
  CrankcaseHeatingPower = DXCoil(DXCoilNum)%CrankcaseHeaterCapacity
ELSE
  CrankcaseHeatingPower = 0.0
END IF

! calculate end time of current time step to determine if error messages should be printed
CurrentEndTime = CurrentTime + SysTimeElapsed

!   Print warning messages only when valid and only for the first ocurrance. Let summary provide statistics.
!   Wait for next time step to print warnings. If simulation iterates, print out
!   the warning for the last iteration only. Must wait for next time step to accomplish this.
!   If a warning occurs and the simulation down shifts, the warning is not valid.
IF(DXCoil(DXCoilNum)%PrintLowAmbMessage)THEN ! .AND. &
  IF(CurrentEndTime .GT. DXCoil(DXCoilNum)%CurrentEndTimeLast .AND. &
     TimeStepSys .GE. DXCoil(DXCoilNum)%TimeStepSysLast)THEN
    IF (DXCoil(DXCoilNum)%LowAmbErrIndex == 0) THEN
      CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%LowAmbBuffer1))
      CALL ShowContinueError(TRIM(DXCoil(DXCoilNum)%LowAmbBuffer2))
      CALL ShowContinueError('... Operation at low ambient temperatures may require special performance curves.')
    ENDIF
    IF (DXCoil(DXCoilNum)%CondenserType(Mode) .EQ. AirCooled) THEN
        CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'&
          //TRIM(DXCoil(DXCoilNum)%Name)//'" - Low condenser dry-bulb temperature error continues...' &
          ,DXCoil(DXCoilNum)%LowAmbErrIndex,DXCoil(DXCoilNum)%LowTempLast,DXCoil(DXCoilNum)%LowTempLast,  &
            ReportMinUnits='[C]',ReportMaxUnits='[C]')
    ELSE
        CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'&
          //TRIM(DXCoil(DXCoilNum)%Name)//'" - Low condenser wet-bulb temperature error continues...' &
          ,DXCoil(DXCoilNum)%LowAmbErrIndex,DXCoil(DXCoilNum)%LowTempLast,DXCoil(DXCoilNum)%LowTempLast,  &
            ReportMinUnits='[C]',ReportMaxUnits='[C]')
    END IF
  END IF
END IF

IF(DXCoil(DXCoilNum)%PrintHighAmbMessage)THEN ! .AND. &
  IF(CurrentEndTime .GT. DXCoil(DXCoilNum)%CurrentEndTimeLast .AND. &
     TimeStepSys .GE. DXCoil(DXCoilNum)%TimeStepSysLast)THEN
    IF (DXCoil(DXCoilNum)%HighAmbErrIndex == 0) THEN
      CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%HighAmbBuffer1))
      CALL ShowContinueError(TRIM(DXCoil(DXCoilNum)%HighAmbBuffer2))
      CALL ShowContinueError('... Operation at high ambient temperatures may require special performance curves.')
    ENDIF
    IF (DXCoil(DXCoilNum)%CondenserType(Mode) .EQ. AirCooled) THEN
        CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'&
          //TRIM(DXCoil(DXCoilNum)%Name)//'" - High condenser dry-bulb temperature error continues...' &
          ,DXCoil(DXCoilNum)%HighAmbErrIndex,DXCoil(DXCoilNum)%HighTempLast,DXCoil(DXCoilNum)%HighTempLast,  &
            ReportMinUnits='[C]',ReportMaxUnits='[C]')
    ELSE
        CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'&
          //TRIM(DXCoil(DXCoilNum)%Name)//'" - High condenser wet-bulb temperature error continues...' &
          ,DXCoil(DXCoilNum)%HighAmbErrIndex,DXCoil(DXCoilNum)%HighTempLast,DXCoil(DXCoilNum)%HighTempLast,  &
            ReportMinUnits='[C]',ReportMaxUnits='[C]')
    END IF
  END IF
END IF

IF(DXCoil(DXCoilNum)%PrintLowOutTempMessage)THEN
  IF(CurrentEndTime .GT. DXCoil(DXCoilNum)%CurrentEndTimeLast .AND. &
       TimeStepSys .GE. DXCoil(DXCoilNum)%TimeStepSysLast)THEN
    IF(DXCoil(DXCoilNum)%LowOutletTempIndex == 0)THEN
      CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%LowOutTempBuffer1))
      CALL ShowContinueError(TRIM(DXCoil(DXCoilNum)%LowOutTempBuffer2))
      CALL ShowContinueError('... Possible reasons for low outlet air dry-bulb temperatures are: This DX coil')
      CALL ShowContinueError('   1) may have a low inlet air dry-bulb temperature. Inlet air temperature = '// &
                                    TRIM(TrimSigDigits(DXCoil(DXCoilNum)%FullLoadInletAirTempLast,3))//' C.')
      CALL ShowContinueError('   2) may have a low air flow rate per watt of cooling capacity. Check inputs.')
      CALL ShowContinueError('   3) is used as part of a HX assisted cooling coil which uses a high sensible'// &
                                   ' effectiveness. Check inputs.')
    END IF
    CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'&
               //TRIM(DXCoil(DXCoilNum)%Name)//'" - Full load outlet temperature'// &
          ' indicates a possibility of frost/freeze error continues. Outlet air temperature statistics follow:', &
          DXCoil(DXCoilNum)%LowOutletTempIndex, DXCoil(DXCoilNum)%FullLoadOutAirTempLast, &
             DXCoil(DXCoilNum)%FullLoadOutAirTempLast)
  END IF
END IF

! save last system time step and last end time of current time step (used to determine if warning is valid)
DXCoil(DXCoilNum)%TimeStepSysLast    = TimeStepSys
DXCoil(DXCoilNum)%CurrentEndTimeLast = CurrentEndTime
DXCoil(DXCoilNum)%PrintLowAmbMessage = .FALSE.
DXCoil(DXCoilNum)%PrintLowOutTempMessage = .FALSE.

IF((AirMassFlow .GT. 0.0) .AND. &
   (GetCurrentScheduleValue(DXCoil(DXCoilNum)%SchedPtr) .GT. 0.0) .AND. &
   (PartLoadRatio .GT. 0.0) .AND. (CompOp == On)) THEN      ! for cycling fan, reset mass flow to full on rate
  IF (FanOpMode .EQ. CycFanCycCoil) THEN
    AirMassFlow = AirMassFlow / PartLoadRatio
  ELSE IF (FanOpMode .EQ. ContFanCycCoil) THEN
    AirMassFlow = AirMassFlow * AirFlowRatio
  ELSE
    AirMassFlow = DXCoil(DXCoilNum)%RatedAirMassFlowRate(Mode)
  END IF

! Check for valid air volume flow per rated total cooling capacity (200 - 500 cfm/ton)

! for some reason there are diff's when using coil inlet air pressure
! these lines (more to follow) are commented out for the time being

  InletAirWetbulbC = PsyTwbFnTdbWPb(InletAirDryBulbTemp,InletAirHumRat,OutdoorPressure)
  AirVolumeFlowRate = AirMassFlow/ PsyRhoAirFnPbTdbW(OutdoorPressure,InletAirDryBulbTemp, InletAirHumRat)
  VolFlowperRatedTotCap = AirVolumeFlowRate/DXCoil(DXCoilNum)%RatedTotCap(Mode)

  IF (DXCoil(DXCoilNum)%RatedTotCap(Mode) .LE. 0.0) THEN
      CALL ShowFatalError(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//&
        '" - Rated total cooling capacity is zero or less.')
  END IF

  IF (.NOT. FirstHVACIteration .AND. .NOT. WarmupFlag .AND. &
      ((VolFlowperRatedTotCap .LT. MinOperVolFlowPerRatedTotCap) .OR. &
       (VolFlowperRatedTotCap .GT. MaxCoolVolFlowPerRatedTotCap))) THEN
    IF (DXCoil(DXCoilNum)%ErrIndex1 == 0) THEN
      CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//&
        '" - Air volume flow rate per watt of rated total cooling capacity is out of range at '//  &
          TRIM(RoundSigDigits(VolFlowperRatedTotCap,3))//' m3/s/W.')
      CALL ShowContinueErrorTimeStamp(' ')
      CALL ShowContinueError('...Expected range for VolumeFlowPerRatedTotalCapacity=['//  &
          TRIM(RoundSigDigits(MinOperVolFlowPerRatedTotCap,3))//'--'//  &
          TRIM(RoundSigDigits(MaxCoolVolFlowPerRatedTotCap,3))//']')
      CALL ShowContinueError('...Possible causes include inconsistent air flow rates in system components,')
      CALL ShowContinueError('...or mixing manual inputs with autosize inputs.'// &
                             ' Also check the following values and calculations.')
      CALL ShowContinueError('...Volume Flow Rate per Rated Total Capacity = Volume Flow Rate / Rated Total Capacity')
      CALL ShowContinueError('...Volume Flow Rate = Air Mass Flow Rate / Air Density')
      CALL ShowContinueError('...Data used for calculations:')
      CALL ShowContinueError('...Rated Total Capacity = '//TRIM(RoundSigDigits(DXCoil(DXCoilNum)%RatedTotCap(Mode),2))//' W.')
      CALL ShowContinueError('...Volume Flow Rate = Air Mass Flow Rate / Air Density')
      CALL ShowContinueError('...Volume Flow Rate   = '//TRIM(RoundSigDigits(AirVolumeFlowRate,8))//' m3/s.')
      CALL ShowContinueError('...Air Mass Flow Rate = '//TRIM(RoundSigDigits(AirMassFlow,8))//' kg/s.')
      CALL ShowContinueError('...Air Density        = '// &
             TRIM(RoundSigDigits(PsyRhoAirFnPbTdbW(OutdoorPressure,InletAirDryBulbTemp, InletAirHumRat),8))//' kg/m3.')
      CALL ShowContinueError('...Data used for air density calculation:')
      CALL ShowContinueError('...Outdoor Air Pressure     = '//TRIM(RoundSigDigits(OutdoorPressure,3))//' Pa.')
      CALL ShowContinueError('...Inlet Air Dry-Bulb Temp  = '//TRIM(RoundSigDigits(InletAirDryBulbTemp,3))//' C.')
      CALL ShowContinueError('...Inlet Air Humidity Ratio = '//TRIM(RoundSigDigits(InletAirHumRat,8))//' kgWater/kgDryAir.')

    ENDIF
    CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//&
          '" - Air volume flow rate per watt of rated total cooling capacity is out ' //&
          'of range error continues...',DXCoil(DXCoilNum)%ErrIndex1,VolFlowperRatedTotCap,VolFlowperRatedTotCap)
  END IF
!
!    Adjust coil bypass factor for actual air flow rate. Use relation CBF = exp(-NTU) where
!    NTU = A0/(m*cp). Relationship models the cooling coil as a heat exchanger with Cmin/Cmax = 0.

  RatedCBF = DXCoil(DXCoilNum)%RatedCBF(Mode)
  IF (RatedCBF .gt. 0.0) THEN
     A0 = -log(RatedCBF)*DXCoil(DXCoilNum)%RatedAirMassFlowRate(Mode)
  ELSE
     A0 = 0.
  END IF
  ADiff=-A0/AirMassFlow
  IF (ADiff >= EXP_LowerLimit) THEN
     CBF = exp(ADiff)
  ELSE
     CBF = 0.0
  END IF

! check boundary for low ambient temperature and post warnings to individual DX coil buffers to print at end of time step
  IF(OutdoorDryBulb .LT. DXCoil(DXCoilNum)%MinOATCompressor .AND. .NOT. WarmupFlag) THEN
    DXCoil(DXCoilNum)%PrintLowAmbMessage = .TRUE.
    DXCoil(DXCoilNum)%LowTempLast = OutdoorDryBulb
    IF(DXCoil(DXCoilNum)%LowAmbErrIndex == 0)THEN
      DXCoil(DXCoilNum)%LowAmbBuffer1 = TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)// &
       '" - Air-cooled condenser inlet dry-bulb temperature below '// &
       TRIM(RoundSigDigits(DXCoil(DXCoilNum)%MinOATCompressor,2))//' C. Outdoor dry-bulb temperature = '//  &
        TRIM(RoundSigDigits(OutdoorDryBulb,2))
      DXCoil(DXCoilNum)%LowAmbBuffer2 = ' '//'... Occurrence info = '//TRIM(EnvironmentName)//', '//Trim(CurMnDy)//' '&
                 //TRIM(CreateSysTimeIntervalString())
    END IF
  END IF

! check boundary for high ambient temperature and post warnings to individual DX coil buffers to print at end of time step
  IF(OutdoorDryBulb .GT. DXCoil(DXCoilNum)%MaxOATCompressor .AND. .NOT. WarmupFlag) THEN
    DXCoil(DXCoilNum)%PrintHighAmbMessage = .TRUE.
    DXCoil(DXCoilNum)%HighTempLast = OutdoorDryBulb
    IF(DXCoil(DXCoilNum)%HighAmbErrIndex == 0)THEN
      DXCoil(DXCoilNum)%HighAmbBuffer1 = TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)// &
       '" - Air-cooled condenser inlet dry-bulb temperature above '// &
       TRIM(RoundSigDigits(DXCoil(DXCoilNum)%MaxOATCompressor,2))//' C. Outdoor dry-bulb temperature = '//  &
        TRIM(RoundSigDigits(OutdoorDryBulb,2))
      DXCoil(DXCoilNum)%HighAmbBuffer2 = ' '//'... Occurrence info = '//TRIM(EnvironmentName)//', '//Trim(CurMnDy)//' '&
                 //TRIM(CreateSysTimeIntervalString())
    END IF
  END IF

!  Get total capacity modifying factor (function of temperature) for off-rated conditions
!  InletAirHumRat may be modified in this ADP/BF loop, use temporary varible for calculations
   InletAirHumRatTemp = InletAirHumRat
50 IF(DXCoil(DXCoilNum)%TotCapTempModFacCurveType(1) .EQ. BiQuadratic)THEN
     TotCapTempModFac = CurveValue(DXCoil(DXCoilNum)%CCapFTemp(Mode),InletAirWetbulbC,CondInletTemp)
   ELSE
     TotCapTempModFac = CurveValue(DXCoil(DXCoilNum)%CCapFTemp(Mode),InletAirWetbulbC)
   END IF

!  Warn user if curve output goes negative
   IF(TotCapTempModFac .LT. 0.0)THEN
     IF(DXCoil(DXCoilNum)%CCapFTempErrorIndex == 0)THEN
       CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//'":')
       CALL ShowContinueError(' Total Cooling Capacity Modifier curve (function of temperature) output is negative (' &
                         //TRIM(TrimSigDigits(TotCapTempModFac,3))//').')
         CALL ShowContinueError(' Negative value occurs using a condenser inlet air temperature of ' &
                          //TRIM(TrimSigDigits(CondInletTemp,1))// &
                             ' and an inlet air wet-bulb temperature of '//TRIM(TrimSigDigits(InletAirWetbulbC,1))//'.')
       IF(Mode .GT. 1)THEN
         Call ShowContinueError(' Negative output results from stage '//TRIM(TrimSigDigits(Mode))// &
                                ' compressor operation.')
       END IF
       CALL ShowContinueErrorTimeStamp(' Resetting curve output to zero and continuing simulation.')
     END IF
     CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//'":'//&
         ' Total Cooling Capacity Modifier curve (function of temperature) output is negative warning continues...' &
         , DXCoil(DXCoilNum)%CCapFTempErrorIndex, TotCapTempModFac, TotCapTempModFac)
     TotCapTempModFac = 0.0
   END IF

!  Get total capacity modifying factor (function of mass flow) for off-rated conditions
   AirMassFlowRatio = AirMassFlow/DXCoil(DXCoilNum)%RatedAirMassFlowRate(Mode)
   TotCapFlowModFac = CurveValue(DXCoil(DXCoilNum)%CCapFFlow(Mode),AirMassFlowRatio)

!  Warn user if curve output goes negative
   IF(TotCapFlowModFac .LT. 0.0)THEN
     IF(DXCoil(DXCoilNum)%CCapFFlowErrorIndex == 0)THEN
       CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//'":')
       CALL ShowContinueError(' Total Cooling Capacity Modifier curve (function of flow fraction) output is negative (' &
                         //TRIM(TrimSigDigits(TotCapFlowModFac,3))//').')
       CALL ShowContinueError(' Negative value occurs using an air flow fraction of ' &
                          //TRIM(TrimSigDigits(AirMassFlowRatio,3))//'.')
       CALL ShowContinueErrorTimeStamp(' Resetting curve output to zero and continuing simulation.')
       IF(Mode .GT. 1)THEN
         Call ShowContinueError(' Negative output results from stage '//TRIM(TrimSigDigits(Mode))// &
                                ' compressor operation.')
       END IF
     END IF
     CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//'":'//&
        ' Total Cooling Capacity Modifier curve (function of flow fraction) output is negative warning continues...' &
        , DXCoil(DXCoilNum)%CCapFFlowErrorIndex, TotCapFlowModFac, TotCapFlowModFac)
     TotCapFlowModFac = 0.0
  END IF

  IF(PRESENT(MaxCoolCap))THEN
    TotCap = MIN(MaxCoolCap,DXCoil(DXCoilNum)%RatedTotCap(Mode) * TotCapFlowModFac * TotCapTempModFac)
  ELSE
    TotCap = DXCoil(DXCoilNum)%RatedTotCap(Mode) * TotCapFlowModFac * TotCapTempModFac
  END IF

  TotCap = TotCap * PartLoadRatio

! Calculate apparatus dew point conditions using TotCap and CBF
  hDelta = TotCap/AirMassFlow
! there is an issue here with using CBF to calculate the ADP enthalpy.
! at low loads the bypass factor increases significantly.
  hADP = InletAirEnthalpy - hDelta/(1.d0-CBF)
  tADP = PsyTsatFnHPb(hADP,OutdoorPressure,'CalcVRFCoolingCoil')
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!  tADP = PsyTsatFnHPb(hADP,InletAirPressure)
  wADP = MIN(InletAirHumRat,PsyWFnTdbH(tADP,hADP,'CalcVRFCoolingCoil'))
  hTinwADP = PsyHFnTdbW(InletAirDryBulbTemp,wADP,'CalcVRFCoolingCoil')
  IF((InletAirEnthalpy-hADP) .NE. 0)THEN
    SHR = MIN((hTinwADP-hADP)/(InletAirEnthalpy-hADP),1.d0)
  ELSE
    SHR = 1.0d0
  END IF
!
! Check for dry evaporator conditions (win < wadp)
!
  IF (wADP .gt. InletAirHumRatTemp .or. (Counter .ge. 1 .and. Counter .lt. MaxIter)) THEN
     If(InletAirHumRatTemp == 0.0)InletAirHumRatTemp=0.00001d0
     werror = (InletAirHumRatTemp - wADP)/InletAirHumRatTemp
!
! Increase InletAirHumRatTemp at constant InletAirTemp to find coil dry-out point. Then use the
! capacity at the dry-out point to determine exiting conditions from coil. This is required
! since the TotCapTempModFac doesn't work properly with dry-coil conditions.
!
      InletAirHumRatTemp = RF*wADP + (1.d0-RF)*InletAirHumRatTemp
      InletAirWetbulbC = PsyTwbFnTdbWPb(InletAirDryBulbTemp,InletAirHumRatTemp,OutdoorPressure)
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!     InletAirWetbulbC = PsyTwbFnTdbWPb(InletAirDryBulbTemp,InletAirHumRatTemp,InletAirPressure)
      Counter = Counter + 1
      IF (ABS(werror) .gt. Tolerance) go to 50   ! Recalculate with modified inlet conditions

  END IF

  IF(DXCoil(DXCoilNum)%PLFFPLR(mode) .GT. 0 .AND. CompCycRatio .LT. 1.d0)THEN
    PLF = CurveValue(DXCoil(DXCoilNum)%PLFFPLR(mode),CompCycRatio) ! Calculate part-load factor
  ELSE
    PLF = 1.0d0
  END IF

  IF (PLF < 0.7d0) THEN
    IF (DXCoil(DXCoilNum)%ErrIndex2 == 0) THEN
      CALL ShowWarningMessage('The PLF curve value for the DX cooling coil '//TRIM(DXCoil(DXCoilNum)%Name)//&
                       ' ='//TRIM(RoundSigDigits(PLF,3))//  &
                       ' for part-load ratio ='//TRIM(RoundSigDigits(PartLoadRatio,3)))
      CALL ShowContinueErrorTimeStamp('PLF curve values must be >= 0.7. PLF has been reset to 0.7 and simulation is continuing.')
      CALL ShowContinueError('Check the IO reference manual for PLF curve guidance [Coil:Cooling:DX:SingleSpeed].')
    END IF
    CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%Name)// &
                  ', DX cooling coil PLF curve < 0.7 warning continues...', &
                   DXCoil(DXCoilNum)%ErrIndex2,PLF,PLF)
    PLF = 0.7d0
  END IF


  DXCoil(DXCoilNum)%PartLoadRatio = PartLoadRatio
  DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction = CompCycRatio / PLF
  IF (DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction > 1.0d0 .and.   &
      ABS(DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction-1.0d0) > .001d0 ) THEN
    IF (DXCoil(DXCoilNum)%ErrIndex3 == 0) THEN
      CALL ShowWarningMessage('The runtime fraction for DX cooling coil '//TRIM(DXCoil(DXCoilNum)%Name)//&
                           ' exceeded 1.0. ['//TRIM(RoundSigDigits(DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction,4))//'].')
      CALL ShowContinueError('Runtime fraction reset to 1 and the simulation will continue.')
      CALL ShowContinueError('Check the IO reference manual for PLF curve guidance [Coil:Cooling:DX:SingleSpeed].')
      CALL ShowContinueErrorTimeStamp(' ')
    END IF
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%Name)//              &
                     ', DX cooling coil runtime fraction > 1.0 warning continues...', &
        DXCoil(DXCoilNum)%ErrIndex3,DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction,DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction)
    DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction = 1.0d0 ! Reset coil runtime fraction to 1.0
  ELSEIF (DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction > 1.0d0) THEN
    DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction = 1.0d0 ! Reset coil runtime fraction to 1.0
  END IF

  ! If cycling fan, send coil part-load fraction to on/off fan via HVACDataGlobals
  IF (FanOpMode .EQ. CycFanCycCoil) OnOffFanPartLoadFraction = PLF

  !  Calculate full load output conditions
  IF (SHR .gt. 1. .OR. Counter .gt. 0) SHR = 1.d0
  FullLoadOutAirEnth = InletAirEnthalpy - TotCap/AirMassFlow
  hTinwout = InletAirEnthalpy - (1.0d0-SHR)*hDelta
  IF (SHR < 1.0d0) THEN
    FullLoadOutAirHumRat = PsyWFnTdbH(InletAirDryBulbTemp,hTinwout)
  ELSE
    FullLoadOutAirHumRat = InletAirHumRat
  END IF
  FullLoadOutAirTemp = PsyTdbFnHW(FullLoadOutAirEnth,FullLoadOutAirHumRat)

! Check for saturation error and modify temperature at constant enthalpy
   IF(FullLoadOutAirTemp .LT. PsyTsatFnHPb(FullLoadOutAirEnth,OutdoorPressure)) THEN
    FullLoadOutAirTemp = PsyTsatFnHPb(FullLoadOutAirEnth,OutdoorPressure)
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!   IF(FullLoadOutAirTemp .LT. PsyTsatFnHPb(FullLoadOutAirEnth,InletAirPressure)) THEN
!    FullLoadOutAirTemp = PsyTsatFnHPb(FullLoadOutAirEnth,InletAirPressure)
    FullLoadOutAirHumRat  = PsyWFnTdbH(FullLoadOutAirTemp,FullLoadOutAirEnth)
   END IF

  ! Store actual outlet conditions when DX coil is ON for use in heat recovery module
  DXCoilFullLoadOutAirTemp(DXCoilNum) = FullLoadOutAirTemp
  DXCoilFullLoadOutAirHumRat(DXCoilNum) = FullLoadOutAirHumRat

! Add warning message for cold cooling coil (FullLoadOutAirTemp < 2 C)
  IF(FullLoadOutAirTemp .LT. 2.0d0 .AND. .NOT. FirstHVACIteration .AND. .NOT. WarmupFlag)THEN
    DXCoil(DXCoilNum)%PrintLowOutTempMessage = .TRUE.
    DXCoil(DXCoilNum)%FullLoadOutAirTempLast = FullLoadOutAirTemp
    IF(DXCoil(DXCoilNum)%LowOutletTempIndex == 0)THEN
      DXCoil(DXCoilNum)%FullLoadInletAirTempLast = InletAirDryBulbTemp
      DXCoil(DXCoilNum)%LowOutTempBuffer1= TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)// &
       '" - Full load outlet air dry-bulb temperature < 2C. This indicates the possibility of coil frost/freeze.'// &
       ' Outlet temperature = '//TRIM(RoundSigDigits(FullLoadOutAirTemp,2))//' C.'
      DXCoil(DXCoilNum)%LowOutTempBuffer2 = ' '//'...Occurrence info = '//TRIM(EnvironmentName)//', '//Trim(CurMnDy)//' '&
                   //TRIM(CreateSysTimeIntervalString())
    END IF
  END IF

  !  If constant fan with cycling compressor, call function to determine "effective SHR"
  !  which includes the part-load degradation on latent capacity
  IF (FanOpMode .EQ. ContFanCycCoil .AND. CompCycRatio .LT. 1.d0) THEN
     QLatRated = DXCoil(DXCoilNum)%RatedTotCap(Mode) * (1.d0 - DXCoil(DXCoilNum)%RatedSHR(Mode))
     QLatActual = TotCap * (1.d0 - SHR)
     SHRUnadjusted = SHR
     SHR = CalcEffectiveSHR(DXCoilNum, SHR, DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction, &
                          QLatRated, QLatActual, InletAirDryBulbTemp, InletAirWetbulbC, Mode)

  !  Calculate full load output conditions
    IF (SHR .gt. 1.d0 .OR. Counter .gt. 0) SHR = 1.d0
    FullLoadOutAirEnth = InletAirEnthalpy - TotCap/AirMassFlow
    hTinwout = InletAirEnthalpy - (1.0d0-SHR)*hDelta
    IF (SHR < 1.0d0) THEN
      FullLoadOutAirHumRat = PsyWFnTdbH(InletAirDryBulbTemp,hTinwout)
    ELSE
      FullLoadOutAirHumRat = InletAirHumRat
    END IF
    FullLoadOutAirTemp = PsyTdbFnHW(FullLoadOutAirEnth,FullLoadOutAirHumRat)

  END IF

!  Calculate actual outlet conditions for the input part load ratio
!  Actual outlet conditions are "average" for time step when compressor cycles

  IF (FanOpMode .EQ. ContFanCycCoil .AND. CompCycRatio .LT. 1.d0) THEN
    ! Continuous fan, cycling compressor
    OutletAirEnthalpy = ((PartLoadRatio * AirFlowRatio)*FullLoadOutAirEnth + &
                                            (1.d0-(PartLoadRatio * AirFlowRatio))*InletAirEnthalpy)
    OutletAirHumRat = ((PartLoadRatio * AirFlowRatio)*FullLoadOutAirHumRat + &
                                            (1.d0-(PartLoadRatio * AirFlowRatio))*InletAirHumRat)
    OutletAirTemp = PsyTdbFnHW(OutletAirEnthalpy,OutletAirHumRat)
  ELSE
    ! Default to cycling fan, cycling compressor
    OutletAirEnthalpy = FullLoadOutAirEnth
    OutletAirHumRat = FullLoadOutAirHumRat
    OutletAirTemp = FullLoadOutAirTemp
  END IF

! Check for saturation error and modify temperature at constant enthalpy
   IF(OutletAirTemp .LT. PsyTsatFnHPb(OutletAirEnthalpy,OutdoorPressure,'CalcVRFCoolingCoil')) THEN
    OutletAirTemp = PsyTsatFnHPb(OutletAirEnthalpy,OutdoorPressure)
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!   IF(OutletAirTemp .LT. PsyTsatFnHPb(OutletAirEnthalpy,InletAirPressure)) THEN
!    OutletAirTemp = PsyTsatFnHPb(OutletAirEnthalpy,InletAirPressure)
    OutletAirHumRat  = PsyWFnTdbH(OutletAirTemp,OutletAirEnthalpy)
   END IF

! Reset AirMassFlow to inlet node air mass flow for final total, sensible and latent calculations
! since AirMassFlow might have been modified above (in this subroutine):
!     IF (FanOpMode .EQ. CycFanCycCoil) AirMassFlow = AirMassFlow / PartLoadRatio
!
! For multimode coil, this should be full flow including bypassed fraction
  AirMassFlow = DXCoil(DXCoilNum)%InletAirMassFlowRate
  DXCoil(DXCoilNum)%TotalCoolingEnergyRate = AirMassFlow * (InletAirEnthalpy - OutletAirEnthalpy)

!! Set DataHeatGlobal heat reclaim variable for use by heat reclaim coil (part load ratio is accounted for)
!! Calculation for heat reclaim needs to be corrected to use compressor power (not including condenser fan power)
!  HeatReclaimDXCoil(DXCoilNum)%AvailCapacity = DXCoil(DXCoilNum)%TotalCoolingEnergyRate + DXCoil(DXCoilNum)%ElecCoolingPower

  MinAirHumRat = MIN(InletAirHumRat,OutletAirHumRat)
  DXCoil(DXCoilNum)%SensCoolingEnergyRate = AirMassFlow * &
                                         (PsyHFnTdbW(InletAirDryBulbTemp,MinAirHumRat) - &
                                          PsyHFnTdbW(OutletAirTemp,MinAirHumRat))
!  Don't let sensible capacity be greater than total capacity
  IF (DXCoil(DXCoilNum)%SensCoolingEnergyRate .GT. DXCoil(DXCoilNum)%TotalCoolingEnergyRate) THEN
     DXCoil(DXCoilNum)%SensCoolingEnergyRate = DXCoil(DXCoilNum)%TotalCoolingEnergyRate
  END IF
!
  DXCoil(DXCoilNum)%LatCoolingEnergyRate = DXCoil(DXCoilNum)%TotalCoolingEnergyRate - DXCoil(DXCoilNum)%SensCoolingEnergyRate
  DXCoil(DXCoilNum)%OutletAirTemp     = OutletAirTemp
  DXCoil(DXCoilNum)%OutletAirHumRat   = OutletAirHumRat
  DXCoil(DXCoilNum)%OutletAirEnthalpy = OutletAirEnthalpy

ELSE

  ! DX coil is off; just pass through conditions
  DXCoil(DXCoilNum)%OutletAirEnthalpy = DXCoil(DXCoilNum)%InletAirEnthalpy
  DXCoil(DXCoilNum)%OutletAirHumRat   = DXCoil(DXCoilNum)%InletAirHumRat
  DXCoil(DXCoilNum)%OutletAirTemp     = DXCoil(DXCoilNum)%InletAirTemp

  DXCoil(DXCoilNum)%ElecCoolingPower = 0.0
  DXCoil(DXCoilNum)%TotalCoolingEnergyRate = 0.0
  DXCoil(DXCoilNum)%SensCoolingEnergyRate = 0.0
  DXCoil(DXCoilNum)%LatCoolingEnergyRate = 0.0
  DXCoil(DXCoilNum)%EvapCondPumpElecPower = 0.0
  DXCoil(DXCoilNum)%EvapWaterConsumpRate = 0.0

! Reset globals when DX coil is OFF for use in heat recovery module
  DXCoilFullLoadOutAirTemp(DXCoilNum) = 0.0
  DXCoilFullLoadOutAirHumRat(DXCoilNum) = 0.0

END IF ! end of on/off if - else

!set water system demand request (if needed)
IF ( DXCoil(DxCoilNum)%EvapWaterSupplyMode == WaterSupplyFromTank) THEN
   WaterStorage(DXCoil(DXCoilNum)%EvapWaterSupTankID)%VdotRequestDemand(DXCoil(DXCoilNum)%EvapWaterTankDemandARRID) &
       = DXCoil(DXCoilNum)%EvapWaterConsumpRate
ENDIF

DXCoilOutletTemp(DXCoilNum)         = DXCoil(DXCoilNum)%OutletAirTemp
DXCoilOutletHumRat(DXCoilNum)       = DXCoil(DXCoilNum)%OutletAirHumRat
DXCoilPartLoadRatio(DXCoilNum)      = DXCoil(DXCoilNum)%PartLoadRatio
DXCoilFanOpMode(DXCoilNum)          = FanOpMode
DXCoil(DXCoilNum)%CondInletTemp     = CondInletTemp
DXCoilTotalCooling(DXCoilNum)       = DXCoil(DXCoilNum)%TotalCoolingEnergyRate
DXCoilCoolInletAirWBTemp(DXCoilNum) = PsyTwbFnTdbWPb(InletAirDryBulbTemp,InletAirHumRat,OutdoorPressure)

RETURN
END SUBROUTINE CalcVRFCoolingCoil

SUBROUTINE CalcDXHeatingCoil(DXCoilNum,PartLoadRatio, FanOpMode, OnOffAirFlowRatio, MaxHeatCap)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Richard Raustad
          !       DATE WRITTEN   October 2001
          !       MODIFIED       Raustad/Shirey Mar 2004
          !                      Kenneth Tang 2004 (Sensitivity of TotCapTempModFac & EIRTempModFac  to indoor dry bulb temp)
          !                      Feb 2005 M. J. Witte, GARD Analytics, Inc.
          !                        Add new coil type COIL:DX:MultiMode:CoolingEmpirical:
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! Calculates the air-side heating performance and electrical heating energy
          ! use of a direct-expansion, air-cooled heat pump unit.

          ! METHODOLOGY EMPLOYED:
          ! This routine simulates the performance of air-cooled DX heating equipment.
          ! The routine requires the user to enter the total heating capacity
          ! and COP for the unit at ARI 210/240 rating conditions (21.11C [70F] dry-bulb,
          ! 15.55C [60F] wet-bulb air entering the heating coil, 8.33C [47F] dry-bulb,
          ! 6.11C [43F] wet-bulb air entering the outdoor condenser. Since different
          ! manufacturer's rate their equipment at different air flow rates, the supply
          ! air flow rate corresponding to the rated capacities and rated COP must also
          ! be entered (should be between 300 cfm/ton and 450 cfm/ton). The rated information
          ! entered by the user should NOT include the thermal or electrical impacts of the
          ! supply air fan, as this is addressed by another module.

          ! With the rated performance data entered by the user, the model employs some of the
          ! DOE-2.1E curve fits to adjust the capacity and efficiency of the unit as a function
          ! of outdoor air temperatures and supply air flow rate (actual vs rated flow). The
          ! model does NOT employ the exact same methodology to calculate performance as DOE-2,
          ! although some of the DOE-2 curve fits are employed by this model.

          ! REFERENCES:
          !
          ! Winkelmann, F.C., Birdsall, B.E., Buhl W.F., Ellington, K.L., Erdem, A.E. 1993.
          ! DOE-2 Supplement Version 2.1E.  Energy and Environment Division, Larwence Berkely
          ! Laboratory.
          !
          ! Henderson, H.I. Jr., Y.J. Huang and Danny Parker. 1999. Residential Equipment Part
          ! Load Curves for Use in DOE-2.  Environmental Energy Technologies Division, Ernest
          ! Orlando Lawrence Berkeley National Laboratory.


          ! USE STATEMENTS:
USE CurveManager, ONLY: CurveValue
USE General,      ONLY: RoundSigDigits

IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
INTEGER,   INTENT(IN) :: DXCoilNum        ! the number of the DX heating coil to be simulated
REAL(r64), INTENT(IN) :: PartLoadRatio    ! sensible cooling load / full load sensible cooling capacity
INTEGER,   INTENT(IN) :: FanOpMode        ! Allows parent object to control fan mode
REAL(r64), INTENT(IN), OPTIONAL :: OnOffAirFlowRatio ! ratio of compressor on airflow to compressor off airflow
REAL(r64), INTENT(IN), OPTIONAL :: MaxHeatCap        ! maximum allowed heating capacity

          ! SUBROUTINE PARAMETER DEFINITIONS:
CHARACTER(len=*), PARAMETER ::  RoutineName='CalcDXHeatingCoil'

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
REAL(r64)     :: AirMassFlow                 ! dry air mass flow rate through coil [kg/s]
REAL(r64)     :: AirMassFlowRatio            ! Ratio of actual air mass flow to rated air mass flow
REAL(r64)     :: AirVolumeFlowRate           ! Air volume flow rate across the cooling coil [m3/s]
REAL(r64)     :: VolFlowperRatedTotCap       ! Air volume flow rate divided by rated total cooling capacity [m3/s-W]
REAL(r64)     :: TotCap                      ! gross total cooling capacity at off-rated conditions [W]
REAL(r64)     :: TotCapTempModFac            ! Total capacity modifier (function of entering drybulb, outside drybulb) depending
                                        ! on the type of curve
REAL(r64)     :: TotCapFlowModFac            ! Total capacity modifier (function of actual supply air flow vs rated flow)
REAL(r64)     :: InletAirDryBulbTemp         ! inlet air dry bulb temperature [C]
REAL(r64)     :: InletAirWetBulbC            ! wetbulb temperature of inlet air [C]
REAL(r64)     :: InletAirEnthalpy            ! inlet air enthalpy [J/kg]
REAL(r64)     :: InletAirHumRat              ! inlet air humidity ratio [kg/kg]
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!REAL(r64)     :: InletAirPressure            ! inlet air pressure [Pa]
REAL(r64)     :: FullLoadOutAirEnth          ! outlet full load enthalpy [J/kg]
REAL(r64)     :: FullLoadOutAirHumRat        ! outlet humidity ratio at full load
REAL(r64)     :: FullLoadOutAirTemp          ! outlet air temperature at full load [C]
REAL(r64)     :: FullLoadOutAirRH            ! outlet air relative humidity at full load
REAL(r64)     :: EIRTempModFac               ! EIR modifier (function of entering drybulb, outside drybulb) depending on the
                                        ! type of curve
REAL(r64)     :: DefrostEIRTempModFac        ! EIR modifier for defrost (function of entering wetbulb, outside drybulb)
REAL(r64)     :: EIRFlowModFac               ! EIR modifier (function of actual supply air flow vs rated flow)
REAL(r64)     :: EIR                         ! EIR at part load and off rated conditions
REAL(r64)     :: PLF                         ! Part load factor, accounts for thermal lag at compressor startup
REAL(r64)     :: PLRHeating                  ! PartLoadRatio in heating
REAL(r64)     :: OutdoorCoilT                ! Outdoor coil temperature (C)
REAL(r64)     :: OutdoorCoildw               ! Outdoor coil delta w assuming coil temp of OutdoorCoilT (kg/kg)
REAL(r64)     :: FractionalDefrostTime       ! Fraction of time step system is in defrost
REAL(r64)     :: HeatingCapacityMultiplier   ! Multiplier for heating capacity when system is in defrost
REAL(r64)     :: InputPowerMultiplier        ! Multiplier for power when system is in defrost
REAL(r64)     :: LoadDueToDefrost            ! Additional load due to defrost
REAL(r64)     :: CrankcaseHeatingPower       ! power due to crankcase heater
REAL(r64)     :: OutdoorDryBulb              ! Outdoor dry-bulb temperature at condenser (C)
REAL(r64)     :: OutdoorWetBulb              ! Outdoor wet-bulb temperature at condenser (C)
REAL(r64)     :: OutdoorHumRat               ! Outdoor humidity ratio at condenser (kg/kg)
REAL(r64)     :: OutdoorPressure             ! Outdoor barometric pressure at condenser (Pa)
INTEGER           :: Mode=1             ! Performance mode for MultiMode DX coil; Always 1 for other coil types
REAL(r64)         :: AirFlowRatio       ! Ratio of compressor on airflow to average timestep airflow
REAL(r64)  :: OutletAirTemp      ! Supply air temperature (average value if constant fan, full output if cycling fan)
REAL(r64)  :: OutletAirHumRat    ! Supply air humidity ratio (average value if constant fan, full output if cycling fan)
REAL(r64)  :: OutletAirEnthalpy  ! Supply air enthalpy (average value if constant fan, full output if cycling fan)

IF (PRESENT(OnOffAirFlowRatio)) THEN
  AirFlowRatio = OnOffAirFlowRatio
ELSE
  AirFlowRatio = 1.0d0
END IF

! Get condenser outdoor node info from DX Heating Coil
IF (DXCoil(DXCoilNum)%CondenserInletNodeNum(1) /= 0) THEN
  OutdoorDryBulb  = Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(1))%Temp
  OutdoorHumRat   = Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(1))%HumRat
  OutdoorPressure = Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(1))%Press
  OutdoorWetBulb  = Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(1))%OutAirWetBulb
ELSE
  OutdoorDryBulb  = OutDryBulbTemp
  OutdoorHumRat   = OutHumRat
  OutdoorPressure = OutBaroPress
  OutdoorWetBulb  = OutWetBulbTemp
ENDIF

AirMassFlow = DXCoil(DXCoilNum)%InletAirMassFlowRate
InletAirDryBulbTemp = DXCoil(DXCoilNum)%InletAirTemp
InletAirEnthalpy = DXCoil(DXCoilNum)%InletAirEnthalpy
InletAirHumRat = DXCoil(DXCoilNum)%InletAirHumRat
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!InletAirPressure = DXCoil(DXCoilNum)%InletAirPressure
!InletAirWetbulbC = PsyTwbFnTdbWPb(InletAirDryBulbTemp,InletAirHumRat,InletAirPressure)
InletAirWetbulbC = PsyTwbFnTdbWPb(InletAirDryBulbTemp,InletAirHumRat,OutdoorPressure)
PLRHeating = 0.0
DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction = 0.0
! Initialize crankcase heater, operates below OAT defined in input deck for HP DX heating coil
IF (OutdoorDryBulb .LT. DXCoil(DXCoilNum)%MaxOATCrankcaseHeater)THEN
  CrankcaseHeatingPower = DXCoil(DXCoilNum)%CrankcaseHeaterCapacity
ELSE
  CrankcaseHeatingPower = 0.0
END IF

IF((AirMassFlow .GT. 0.0) .AND. &
   (GetCurrentScheduleValue(DXCoil(DXCoilNum)%SchedPtr) .GT. 0.0) .AND. &
   (PartLoadRatio .GT. 0.0) .AND. OutdoorDryBulb .GT. DXCoil(DXCoilNum)%MinOATCompressor) THEN
! for cycling fan, reset mass flow to full on rate
  IF (FanOpMode .EQ. CycFanCycCoil) AirMassFlow = AirMassFlow / PartLoadRatio
  IF (FanOpMode .EQ. ContFanCycCoil) AirMassFlow = AirMassFlow * AirFlowRatio
!
! Check for valid air volume flow per rated total cooling capacity (200 - 600 cfm/ton)
!
  AirVolumeFlowRate = AirMassFlow/PsyRhoAirFnPbTdbW(OutdoorPressure,InletAirDryBulbTemp, InletAirHumRat)
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!  AirVolumeFlowRate = AirMassFlow/PsyRhoAirFnPbTdbW(InletAirPressure,InletAirDryBulbTemp, InletAirHumRat)
  VolFlowperRatedTotCap = AirVolumeFlowRate/DXCoil(DXCoilNum)%RatedTotCap(Mode)

  IF ((VolFlowperRatedTotCap.LT.MinOperVolFlowPerRatedTotCap).OR.(VolFlowperRatedTotCap.GT.MaxHeatVolFlowPerRatedTotCap)) THEN
    IF (DXCoil(DXCoilNum)%ErrIndex1 == 0) THEN
      CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//&
        '" - Air volume flow rate per watt of rated total heating capacity is out of range at '// &
          TRIM(RoundSigDigits(VolFlowperRatedTotCap,3))//' m3/s/W.')
      CALL ShowContinueErrorTimeStamp(' ')
      CALL ShowContinueError('Expected range for VolumeFlowPerRatedTotalCapacity=['//  &
          TRIM(RoundSigDigits(MinOperVolFlowPerRatedTotCap,3))//'--'//  &
          TRIM(RoundSigDigits(MaxHeatVolFlowPerRatedTotCap,3))//']')
      CALL ShowContinueError('Possible causes include inconsistent air flow rates in system components or')
      CALL ShowContinueError('inconsistent supply air fan operation modes in coil and unitary system objects.')
    ENDIF
    CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//&
        '" - Air volume flow rate per watt of rated total heating capacity is out ' //&
        'of range error continues...',DXCoil(DXCoilNum)%ErrIndex1,  &
        ReportMinOf=VolFlowperRatedTotCap,ReportMaxOf=VolFlowperRatedTotCap)
  END IF

! Get total capacity modifying factor (function of temperature) for off-rated conditions
! Model was extended to accept bi-quadratic curves. This allows sensitivity of the heating capacity
! to the entering dry-bulb temperature as well as the outside dry-bulb temperature. User is
! advised to use the bi-quaratic curve if sufficient manufacturer data is available.
 IF (DXCoil(DXCoilNum)%TotCapTempModFacCurveType(Mode) == Biquadratic) THEN
   SELECT CASE(DXCoil(DXCoilNum)%HeatingPerformanceOATType)
     CASE(DryBulbIndicator)
       TotCapTempModFac = CurveValue(DXCoil(DXCoilNum)%CCapFTemp(Mode),InletAirDryBulbTemp,OutdoorDryBulb)
     CASE(WetBulbIndicator)
       TotCapTempModFac = CurveValue(DXCoil(DXCoilNum)%CCapFTemp(Mode),InletAirDryBulbTemp,OutdoorWetBulb)
     CASE DEFAULT
       TotCapTempModFac = 1.d0
   END SELECT
 ELSE
   SELECT CASE(DXCoil(DXCoilNum)%HeatingPerformanceOATType)
     CASE(DryBulbIndicator)
       IF(DXCoil(DXCoilNum)%DXCoilType_Num /= CoilVRF_Heating)THEN
         TotCapTempModFac = CurveValue(DXCoil(DXCoilNum)%CCapFTemp(Mode),OutdoorDryBulb)
       ELSE
         TotCapTempModFac = CurveValue(DXCoil(DXCoilNum)%CCapFTemp(Mode),InletAirDryBulbTemp)
       END IF
     CASE(WetBulbIndicator)
       IF(DXCoil(DXCoilNum)%DXCoilType_Num /= CoilVRF_Heating)THEN
         TotCapTempModFac = CurveValue(DXCoil(DXCoilNum)%CCapFTemp(Mode),OutdoorWetBulb)
       ELSE
         TotCapTempModFac = CurveValue(DXCoil(DXCoilNum)%CCapFTemp(Mode),InletAirDryBulbTemp)
       END IF
     CASE DEFAULT
       TotCapTempModFac = 1.d0
   END SELECT
 END IF

!  Get total capacity modifying factor (function of mass flow) for off-rated conditions
  AirMassFlowRatio = AirMassFlow/DXCoil(DXCoilNum)%RatedAirMassFlowRate(Mode)
  TotCapFlowModFac = CurveValue(DXCoil(DXCoilNum)%CCapFFlow(Mode),AirMassFlowRatio)

! Calculate total heating capacity for off-rated conditions
  TotCap = DXCoil(DXCoilNum)%RatedTotCap(Mode) * TotCapFlowModFac * TotCapTempModFac

! Calculating adjustment factors for defrost
! Calculate delta w through outdoor coil by assuming a coil temp of 0.82*DBT-9.7(F) per DOE2.1E
  OutdoorCoilT = 0.82d0 * OutdoorDryBulb - 8.589d0
  OutdoorCoildw = MAX(1.0d-6,(OutdoorHumRat - PsyWFnTdpPb(OutdoorCoilT,OutdoorPressure)))

! Initializing defrost adjustment factors
  LoadDueToDefrost = 0.0
  HeatingCapacityMultiplier = 1.0d0
  FractionalDefrostTime = 0.0
  InputPowerMultiplier = 1.0d0

! Check outdoor temperature to determine of defrost is active
  IF (OutdoorDryBulb .LE. DXCoil(DXCoilNum)%MaxOATDefrost) THEN
! Calculate defrost adjustment factors depending on defrost control type
     IF (DXCoil(DXCoilNum)%DefrostControl .EQ. Timed) THEN
       FractionalDefrostTime = DXCoil(DXCoilNum)%DefrostTime
       HeatingCapacityMultiplier = 0.909d0 - 107.33d0 * OutdoorCoildw
       InputPowerMultiplier = 0.90d0 - 36.45d0*OutdoorCoildw
     ELSE !else defrost control is on-demand
       FractionalDefrostTime = 1.0d0 / (1.0d0 + 0.01446d0 / OutdoorCoildw)
       HeatingCapacityMultiplier = 0.875d0 * ( 1.0d0 - FractionalDefrostTime)
       InputPowerMultiplier = 0.954d0 * ( 1.0d0 - FractionalDefrostTime)
     END IF


     IF (FractionalDefrostTime .GT. 0.0) THEN
! Calculate defrost adjustment factors depending on defrost control strategy
       IF (DXCoil(DXCoilNum)%DefrostStrategy .EQ. ReverseCycle) THEN
         LoadDueToDefrost = (0.01d0 * FractionalDefrostTime) * &
                            (7.222d0 - OutdoorDryBulb) * &
                            (DXCoil(DXCoilNum)%RatedTotCap(Mode)/1.01667d0)
         DefrostEIRTempModFac = CurveValue(DXCoil(DXCoilNum)%DefrostEIRFT,&
                                MAX(15.555d0,InletAirWetbulbC),MAX(15.555d0,OutdoorDryBulb))
         DXCoil(DXCoilNum)%DefrostPower =  DefrostEIRTempModFac * &
                                           (DXCoil(DXCoilNum)%RatedTotCap(Mode)/1.01667d0) &
                                           * FractionalDefrostTime
       ELSE ! Defrost strategy is resistive
         DXCoil(DXCoilNum)%DefrostPower = DXCoil(DXCoilNum)%DefrostCapacity &
                                          * FractionalDefrostTime
       END IF
     ELSE ! Defrost is not active because (OutDryBulbTemp .GT. DXCoil(DXCoilNum)%MaxOATDefrost)
       DXCoil(DXCoilNum)%DefrostPower =  0.0
     END IF
  END IF

! Modify total heating capacity based on defrost heating capacity multiplier
! MaxHeatCap passed from parent object VRF Condenser and is used to limit capacity of TU's to that available from condenser
  IF(PRESENT(MaxHeatCap))THEN
    TotCap = MIN(MaxHeatCap,TotCap * HeatingCapacityMultiplier)
  ELSE
    TotCap = TotCap * HeatingCapacityMultiplier
  END IF

! Calculate full load outlet conditions
  FullLoadOutAirEnth = InletAirEnthalpy + TotCap/AirMassFlow
  FullLoadOutAirHumRat = InletAirHumRat
  FullLoadOutAirTemp = PsyTdbFnHW(FullLoadOutAirEnth,FullLoadOutAirHumRat)
  FullLoadOutAirRH = PsyRhFnTdbWPb(FullLoadOutAirTemp,FullLoadOutAirHumRat,OutdoorPressure,'CalcDXHeatingCoil:fullload')
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!  FullLoadOutAirRH = PsyRhFnTdbWPb(FullLoadOutAirTemp,FullLoadOutAirHumRat,InletAirPressure)
  IF (FullLoadOutAirRH .gt. 1.) THEN  ! Limit to saturated conditions at FullLoadOutAirEnth
    FullLoadOutAirTemp = PsyTsatFnHPb(FullLoadOutAirEnth,OutdoorPressure)
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!    FullLoadOutAirTemp = PsyTsatFnHPb(FullLoadOutAirEnth,InletAirPressure)
    FullLoadOutAirHumRat  = PsyWFnTdbH(FullLoadOutAirTemp,FullLoadOutAirEnth)
  END IF

  ! Calculate actual outlet conditions for the input part load ratio
  ! Actual outlet conditions are "average" for time step
  IF (FanOpMode .EQ. ContFanCycCoil) THEN
    ! continuous fan, cycling compressor
    OutletAirEnthalpy = ((PartLoadRatio * AirFlowRatio)*FullLoadOutAirEnth + (1.-(PartLoadRatio * AirFlowRatio))*InletAirEnthalpy)
    OutletAirHumRat   = (PartLoadRatio*FullLoadOutAirHumRat + (1.-PartLoadRatio)*InletAirHumRat)
    OutletAirTemp     = PsyTdbFnHW(OutletAirEnthalpy,OutletAirHumRat)
  ELSE
    ! default to cycling fan, cycling compressor
    OutletAirEnthalpy = FullLoadOutAirEnth
    OutletAirHumRat   = FullLoadOutAirHumRat
    OutletAirTemp     = FullLoadOutAirTemp
  END IF
  ! Calculate electricity consumed. First, get EIR modifying factors for off-rated conditions
  ! Model was extended to accept bi-quadratic curves. This allows sensitivity of the EIR
  ! to the entering dry-bulb temperature as well as the outside dry-bulb temperature. User is
  ! advised to use the bi-quaratic curve if sufficient manufacturer data is available.
  IF(DXCoil(DXCoilNum)%DXCoilType_Num /= CoilVRF_Heating)THEN
    IF ((DXCoil(DXCoilNum)%EIRTempModFacCurveType(1) == Quadratic).OR.(DXCoil(DXCoilNum)%EIRTempModFacCurveType(1) == Cubic)) THEN
      EIRTempModFac = CurveValue(DXCoil(DXCoilNum)%EIRFTemp(Mode),OutdoorDryBulb)
    ELSEIF (DXCoil(DXCoilNum)%EIRTempModFacCurveType(1) == Biquadratic) THEN
      EIRTempModFac = CurveValue(DXCoil(DXCoilNum)%EIRFTemp(Mode),InletAirDryBulbTemp,OutdoorDryBulb)
    END IF
    EIRFlowModFac = CurveValue(DXCoil(DXCoilNum)%EIRFFlow(Mode), AirMassFlowRatio)
  ELSE
    EIRTempModFac = 1.0d0
    EIRFlowModFac = 1.0d0
  END IF
  EIR = DXCoil(DXCoilNum)%RatedEIR(Mode) * EIRTempModFac * EIRFlowModFac
  ! Calculate modified PartLoadRatio due to defrost (reverse-cycle defrost only)
  PLRHeating = MIN(1.0d0,(PartLoadRatio + LoadDueToDefrost/TotCap))
  IF(DXCoil(DXCoilNum)%DXCoilType_Num /= CoilVRF_Heating)THEN
    PLF = CurveValue(DXCoil(DXCoilNum)%PLFFPLR(Mode),PLRHeating) ! Calculate part-load factor
  ELSE
    PLF = 1.0d0
  END IF

  IF (PLF < 0.7d0) THEN
    IF (DXCoil(DXCoilNum)%PLRErrIndex == 0) THEN
      CALL ShowWarningMessage('The PLF curve value for DX heating coil '//TRIM(DXCoil(DXCoilNum)%Name)//&
                         ' ='//TRIM(RoundSigDigits(PLF,2))//  &
                         ' for part-load ratio ='//TRIM(RoundSigDigits(PLRHeating,2)))
      CALL ShowContinueError('PLF curve values must be >= 0.7. PLF has been reset to 0.7 and simulation is continuing.')
      CALL ShowContinueError('Check the IO reference manual for PLF curve guidance [Coil:Heating:DX:SingleSpeed].')
      CALL ShowContinueErrorTimeStamp(' ')
    ENDIF
    CALL ShowRecurringWarningErrorAtEnd('DX heating coil PLF curve < 0.7 warning continues... ',  &
        DXCoil(DXCoilNum)%PLRErrIndex,ReportMinOf=PLF,ReportMaxOf=PLF)
    PLF = 0.7d0
  END IF

  DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction = (PLRHeating / PLF)
  IF (DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction > 1.0d0 .and. ABS(DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction-1.0) > .001 ) THEN
    IF (DXCoil(DXCoilNum)%ErrIndex4 == 0) THEN
      CALL ShowWarningMessage('The runtime fraction for DX heating coil '//TRIM(DXCoil(DXCoilNum)%Name)//&
                         ' exceeded 1.0. ['//TRIM(RoundSigDigits(DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction,4))//'].')
      CALL ShowContinueError('Runtime fraction is set to 1.0 and the simulation continues...')
      CALL ShowContinueError('Check the IO reference manual for PLF curve guidance [Coil:Heating:DX:SingleSpeed].')
      CALL ShowContinueErrorTimeStamp(' ')
    ENDIF
    CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%Name)//              &
                   ', DX heating coil runtime fraction > 1.0 warning continues...', &
           DXCoil(DXCoilNum)%ErrIndex4,DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction,DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction)
    DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction = 1.0d0 ! Reset coil runtime fraction to 1.0
  ELSEIF (DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction > 1.0) THEN
    DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction = 1.0d0 ! Reset coil runtime fraction to 1.0
  END IF
  ! if cycling fan, send coil part-load fraction to on/off fan via HVACDataGlobals
  IF (FanOpMode .EQ. CycFanCycCoil) OnOffFanPartLoadFraction = PLF
  DXCoil(DXCoilNum)%ElecHeatingPower = TotCap * EIR * DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction * InputPowerMultiplier

! Calculate crankcase heater power using the runtime fraction for this DX heating coil only if there is no companion DX coil.
! Else use the largest runtime fraction of this DX heating coil and the companion DX cooling coil.

  IF(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil .EQ. 0) THEN
    DXCoil(DXCoilNum)%CrankcaseHeaterPower = CrankcaseHeatingPower * &
                                           (1.0d0 - DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction)
  ELSE
    DXCoil(DXCoilNum)%CrankcaseHeaterPower = CrankcaseHeatingPower * &
                                           (1.0d0 - MAX(DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction, &
                                           DXCoil(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil)%CoolingCoilRuntimeFraction))
  END IF

  AirMassFlow = DXCoil(DXCoilNum)%InletAirMassFlowRate
  DXCoil(DXCoilNum)%TotalHeatingEnergyRate = AirMassFlow * (OutletAirEnthalpy - InletAirEnthalpy)
! Adjust defrost power to correct for DOE-2 bug where defrost power is constant regardless of compressor runtime fraction
! Defrosts happen based on compressor run time (frost buildup on outdoor coil), not total elapsed time.
  DXCoil(DXCoilNum)%DefrostPower = DXCoil(DXCoilNum)%DefrostPower * DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction

  DXCoil(DXCoilNum)%OutletAirTemp     = OutletAirTemp
  DXCoil(DXCoilNum)%OutletAirHumRat   = OutletAirHumRat
  DXCoil(DXCoilNum)%OutletAirEnthalpy = OutletAirEnthalpy

ELSE

  ! DX coil is off; just pass through conditions
  DXCoil(DXCoilNum)%OutletAirEnthalpy = DXCoil(DXCoilNum)%InletAirEnthalpy
  DXCoil(DXCoilNum)%OutletAirHumRat   = DXCoil(DXCoilNum)%InletAirHumRat
  DXCoil(DXCoilNum)%OutletAirTemp     = DXCoil(DXCoilNum)%InletAirTemp

  DXCoil(DXCoilNum)%ElecHeatingPower       = 0.0
  DXCoil(DXCoilNum)%TotalHeatingEnergyRate = 0.0
  DXCoil(DXCoilNum)%DefrostPower           = 0.0

! Calculate crankcase heater power using the runtime fraction for this DX heating coil (here DXHeatingCoilRTF=0) if
! there is no companion DX coil, or the runtime fraction of the companion DX cooling coil (here DXCoolingCoilRTF>=0).
  IF(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil .EQ. 0) THEN
    DXCoil(DXCoilNum)%CrankcaseHeaterPower = CrankcaseHeatingPower
  ELSE
    DXCoil(DXCoilNum)%CrankcaseHeaterPower = CrankcaseHeatingPower * &
                                            (1.d0-DXCoil(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil)%CoolingCoilRuntimeFraction)
  END IF

END IF ! end of on/off if - else

DXCoilOutletTemp(DXCoilNum)    = DXCoil(DXCoilNum)%OutletAirTemp
DXCoilOutletHumRat(DXCoilNum)  = DXCoil(DXCoilNum)%OutletAirHumRat
DXCoilFanOpMode(DXCoilNum)     = FanOpMode
DXCoilPartLoadRatio(DXCoilNum) = PLRHeating
DXCoilTotalHeating(DXCoilNum)  = DXCoil(DXCoilNum)%TotalHeatingEnergyRate
DXCoilHeatInletAirDBTemp(DXCoilNum) = InletAirDryBulbTemp
DXCoilHeatInletAirWBTemp(DXCoilNum) = InletAirWetbulbC

RETURN
END SUBROUTINE CalcDXHeatingCoil

SUBROUTINE CalcMultiSpeedDXCoil(DXCoilNum,SpeedRatio, CycRatio)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Fred Buhl
          !       DATE WRITTEN   September 2002
          !       MODIFIED       Raustad/Shirey, Feb 2004
          !                      Feb 2005 M. J. Witte, GARD Analytics, Inc.
          !                        Add new coil type COIL:DX:MultiMode:CoolingEmpirical:
          !                      April 2010, Chandan sharma, FSEC, added basin heater
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! Calculates the air-side performance and electrical energy use of a direct-
          ! expansion, air-cooled cooling unit with a 2 speed or variable speed compressor.

          ! METHODOLOGY EMPLOYED:
          ! Uses the same methodology as the single speed DX unit model (SUBROUTINE CalcDoe2DXCoil).
          ! In addition it assumes that the unit performance is obtained by interpolating between
          ! the performance at high speed and that at low speed. If the output needed is below
          ! that produced at low speed, the compressor cycles between off and low speed.

          ! USE STATEMENTS:
USE CurveManager, ONLY: CurveValue
USE DataWater,    ONLY: WaterStorage
!USE ScheduleManager, ONLY: GetCurrentScheduleValue

IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
INTEGER :: DXCoilNum        ! the number of the DX heating coil to be simulated
REAL(r64) :: SpeedRatio       ! = (CompressorSpeed - CompressorSpeedMin) / (CompressorSpeedMax - CompressorSpeedMin)
                            ! SpeedRatio varies between 1.0 (maximum speed) and 0.0 (minimum speed)
REAL(r64) :: CycRatio         ! cycling part load ratio

          ! SUBROUTINE PARAMETER DEFINITIONS:
CHARACTER(len=*), PARAMETER ::  RoutineName='CalcMultiSpeedDXCoil'

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na


          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
REAL(r64) :: AirMassFlow         ! dry air mass flow rate through coil [kg/s]
REAL(r64) :: AirMassFlowRatio    ! Ratio of max air mass flow to rated air mass flow
REAL(r64) :: InletAirWetBulbC    ! wetbulb temperature of inlet air [C]
REAL(r64) :: InletAirDryBulbTemp ! inlet air dry bulb temperature [C]
REAL(r64) :: InletAirEnthalpy    ! inlet air enthalpy [J/kg]
REAL(r64) :: InletAirHumRat      ! inlet air humidity ratio [kg/kg]
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!REAL(r64) :: InletAirPressure    ! inlet air pressure [Pa]
REAL(r64) :: OutletAirDryBulbTemp    ! outlet air dry bulb temperature [C]
REAL(r64) :: OutletAirEnthalpy   ! outlet air enthalpy [J/kg]
REAL(r64) :: OutletAirHumRat     ! outlet air humidity ratio [kg/kg]
! REAL(r64) :: OutletAirRH         ! outlet air relative humudity [fraction]
REAL(r64) :: OutletAirDryBulbTempSat ! outlet air dry bulb temp at saturation at the outlet enthalpy [C]
REAL(r64) :: LSOutletAirDryBulbTemp ! low speed outlet air dry bulb temperature [C]
REAL(r64) :: LSOutletAirEnthalpy    ! low speed outlet air enthalpy [J/kg]
REAL(r64) :: LSOutletAirHumRat      ! low speed outlet air humidity ratio [kg/kg]
REAL(r64) :: LSOutletAirRH          ! low speed outlet air relative humudity [fraction]
REAL(r64) :: hDelta              ! Change in air enthalpy across the cooling coil [J/kg]
REAL(r64) :: hTinwout            ! Enthalpy at inlet dry-bulb and outlet humidity ratio [J/kg]
REAL(r64) :: hADP                ! Apparatus dew point enthalpy [J/kg]
REAL(r64) :: tADP                ! Apparatus dew point temperature [C]
REAL(r64) :: wADP                ! Apparatus dew point humidity ratio [kg/kg]
REAL(r64) :: hTinwADP            ! Enthalpy at inlet dry-bulb and wADP [J/kg]
REAL(r64) :: RatedCBFHS          ! coil bypass factor at rated conditions (high speed)
REAL(r64) :: CBFHS               ! coil bypass factor at max flow (high speed)
REAL(r64) :: TotCapHS            ! total capacity at high speed [W]
REAL(r64) :: SHRHS               ! sensible heat ratio at high speed
REAL(r64) :: TotCapLS            ! total capacity at low speed [W]
REAL(r64) :: SHRLS               ! sensible heat ratio at low speed
REAL(r64) :: EIRTempModFacHS     ! EIR modifier (function of entering wetbulb, outside drybulb) (high speed)
REAL(r64) :: EIRFlowModFacHS     ! EIR modifier (function of actual supply air flow vs rated flow) (high speed)
REAL(r64) :: EIRHS               ! EIR at off rated conditions (high speed)
REAL(r64) :: EIRTempModFacLS     ! EIR modifier (function of entering wetbulb, outside drybulb) (low speed)
REAL(r64) :: EIRLS               ! EIR at off rated conditions (low speed)
REAL(r64) :: TotCap              ! total capacity at current speed [W]
REAL(r64) :: SHR                 ! sensible heat ratio at current speed
REAL(r64) :: EIR                 ! EIR at current speed
REAL(r64) :: AirMassFlowNom      ! speed ratio weighted average of high and low speed air mass flow rates [kg/s]
REAL(r64) :: CBFNom              ! coil bypass factor corresponding to AirMassFlowNom and SpeedRatio
REAL(r64) :: CBF                 ! CBFNom adjusted for actual air mass flow rate
REAL(r64) :: PLF                 ! Part load factor, accounts for thermal lag at compressor startup, used in
                               ! power calculation
REAL(r64) :: CondInletTemp       ! Condenser inlet temperature (C). Outdoor dry-bulb temp for air-cooled condenser.
                               ! Outdoor Wetbulb +(1 - effectiveness)*(outdoor drybulb - outdoor wetbulb) for evap condenser.
REAL(r64) :: CondInletHumrat     ! Condenser inlet humidity ratio (kg/kg). Zero for air-cooled condenser.
                               ! For evap condenser, its the humidity ratio of the air leaving the evap cooling pads.
REAL(r64) :: RhoAir              ! Density of air [kg/m3]
REAL(r64) :: RhoWater            ! Density of water [kg/m3]
REAL(r64) :: CondAirMassFlow     ! Condenser air mass flow rate [kg/s]
REAL(r64) :: EvapCondPumpElecPower ! Evaporative condenser electric pump power [W]
REAL(r64) :: MinAirHumRat = 0.0    ! minimum of the inlet air humidity ratio and the outlet air humidity ratio
INTEGER,SAVE :: Mode=1  ! Performance mode for MultiMode DX coil; Always 1 for other coil types
REAL(r64) :: OutdoorDryBulb       ! Outdoor dry-bulb temperature at condenser (C)
REAL(r64) :: OutdoorWetBulb       ! Outdoor wet-bulb temperature at condenser (C)
REAL(r64) :: OutdoorHumRat        ! Outdoor humidity ratio at condenser (kg/kg)
REAL(r64) :: OutdoorPressure      ! Outdoor barometric pressure at condenser (Pa)

IF (DXCoil(DXCoilNum)%CondenserInletNodeNum(Mode) /= 0) THEN
  OutdoorDryBulb  = Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(Mode))%Temp
  OutdoorHumRat   = Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(Mode))%HumRat
  OutdoorPressure = Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(Mode))%Press
  OutdoorWetBulb  = PsyTwbFnTdbWPb(OutdoorDryBulb,OutdoorHumRat,OutdoorPressure)
ELSE
  OutdoorDryBulb  = OutDryBulbTemp
  OutdoorHumRat   = OutHumRat
  OutdoorPressure = OutBaroPress
  OutdoorWetBulb  = OutWetBulbTemp
ENDIF

AirMassFlow = DXCoil(DXCoilNum)%InletAirMassFlowRate
AirMassFlowRatio = DXCoil(DXCoilNum)%InletAirMassFlowRateMax / DXCoil(DXCoilNum)%RatedAirMassFlowRate(Mode)
DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction = 0.0
InletAirDryBulbTemp = DXCoil(DXCoilNum)%InletAirTemp
InletAirEnthalpy = DXCoil(DXCoilNum)%InletAirEnthalpy
InletAirHumRat = DXCoil(DXCoilNum)%InletAirHumRat
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!InletAirPressure = DXCoil(DXCoilNum)%InletAirPressure
!InletAirWetbulbC = PsyTwbFnTdbWPb(InletAirDryBulbTemp,InletAirHumRat,InletAirPressure)
InletAirWetbulbC = PsyTwbFnTdbWPb(InletAirDryBulbTemp,InletAirHumRat,OutdoorPressure)
IF (DXCoil(DXCoilNum)%CondenserType(Mode) == AirCooled) THEN
  CondInletTemp = OutdoorDryBulb ! Outdoor dry-bulb temp
ELSEIF (DXCoil(DXCoilNum)%CondenserType(Mode) == EvapCooled) THEN
 ! Outdoor wet-bulb temp from DataEnvironment + (1.0-EvapCondEffectiveness) * (drybulb - wetbulb)
  CondInletTemp = OutdoorWetBulb + (OutdoorDryBulb-OutdoorWetBulb)*(1.0d0 - DXCoil(DXCoilNum)%EvapCondEffect(Mode))
  CondInletHumrat = PsyWFnTdbTwbPb(CondInletTemp,OutdoorWetBulb,OutdoorPressure)
END IF

IF((AirMassFlow .GT. 0.0) .AND. &
   (GetCurrentScheduleValue(DXCoil(DXCoilNum)%SchedPtr) .GT. 0.0) &
   .AND. (SpeedRatio > 0.0 .OR. CycRatio > 0.0) ) THEN

  RhoAir = PsyRhoAirFnPbTdbW(OutdoorPressure,OutdoorDryBulb,OutdoorHumRat)
  IF (SpeedRatio > 0.0) THEN
    ! Adjust high speed coil bypass factor for actual maximum air flow rate.
    RatedCBFHS = DXCoil(DXCoilNum)%RatedCBF(Mode)
    CBFHS = AdjustCBF(RatedCBFHS,DXCoil(DXCoilNum)%RatedAirMassFlowRate(Mode),DXCoil(DXCoilNum)%InletAirMassFlowRateMax)
    ! get high speed total capacity and SHR at current conditions
    CALL CalcTotCapSHR(InletAirDryBulbTemp,InletAirHumRat,InletAirEnthalpy,InletAirWetbulbC,AirMassFlowRatio, &
                       DXCoil(DXCoilNum)%InletAirMassFlowRateMax,DXCoil(DXCoilNum)%RatedTotCap(Mode), &
                       CBFHS,DXCoil(DXCoilNum)%CCapFTemp(Mode),DXCoil(DXCoilNum)%CCapFFlow(Mode),TotCapHS,SHRHS, &
                       CondInletTemp, OutdoorPressure)
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!                       CondInletTemp, Node(DXCoil(DXCoilNum)%AirInNode)%Press)
    ! get low speed total capacity and SHR at current conditions
    CALL CalcTotCapSHR(InletAirDryBulbTemp,InletAirHumRat,InletAirEnthalpy,InletAirWetbulbC,1.0d0, &
                       DXCoil(DXCoilNum)%RatedAirMassFlowRate2,DXCoil(DXCoilNum)%RatedTotCap2, &
                       DXCoil(DXCoilNum)%RatedCBF2,DXCoil(DXCoilNum)%CCapFTemp2, &
                       DXCoil(DXCoilNum)%CCapFFlow(Mode),TotCapLS,SHRLS,CondInletTemp, &
                       OutdoorPressure)
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!                       Node(DXCoil(DXCoilNum)%AirInNode)%Press)
    ! get high speed EIR at current conditions
    EIRTempModFacHS = CurveValue(DXCoil(DXCoilNum)%EIRFTemp(Mode),InletAirWetbulbC,CondInletTemp)
    EIRFlowModFacHS = CurveValue(DXCoil(DXCoilNum)%EIRFFlow(Mode),AirMassFlowRatio)
    EIRHS = DXCoil(DXCoilNum)%RatedEIR(Mode) * EIRFlowModFacHS * EIRTempModFacHS
    ! get low speed EIR at current conditions
!    EIRTempModFacLS = CurveValue(DXCoil(DXCoilNum)%EIRFTemp(Mode),InletAirWetbulbC,CondInletTemp)
!    CR7307 changed EIRTempModFacLS calculation to that shown below.
    EIRTempModFacLS = CurveValue(DXCoil(DXCoilNum)%EIRFTemp2,InletAirWetbulbC,CondInletTemp)
    EIRLS = DXCoil(DXCoilNum)%RatedEIR2 * EIRTempModFacLS

    ! get current total capacity, SHR, EIR
    IF (SpeedRatio >= 1.0d0) THEN
      TotCap = TotCapHS
      SHR = SHRHS
      EIR = EIRHS
      CBFNom = CBFHS
      AirMassFlowNom = DXCoil(DXCoilNum)%InletAirMassFlowRateMax
      CondAirMassFlow =  RhoAir * DXCoil(DXCoilNum)%EvapCondAirFlow(Mode)
      EvapCondPumpElecPower = DXCoil(DXCoilNum)%EvapCondPumpElecNomPower(Mode)
    ELSE
      TotCap = SpeedRatio*TotCapHS + (1.0d0-SpeedRatio)*TotCapLS
      EIR = SpeedRatio*EIRHS + (1.0d0-SpeedRatio)*EIRLS
      CBFNom = SpeedRatio*CBFHS + (1.0d0-SpeedRatio)*DXCoil(DXCoilNum)%RatedCBF2
      AirMassFlowNom = SpeedRatio*DXCoil(DXCoilNum)%InletAirMassFlowRateMax + (1.0d0-SpeedRatio)* &
                       DXCoil(DXCoilNum)%RatedAirMassFlowRate2
      CondAirMassFlow =  RhoAir * (SpeedRatio * DXCoil(DXCoilNum)%EvapCondAirFlow(Mode) + (1.0d0-SpeedRatio)* &
                       DXCoil(DXCoilNum)%EvapCondAirFlow2)
      EvapCondPumpElecPower = SpeedRatio * DXCoil(DXCoilNum)%EvapCondPumpElecNomPower(Mode) + (1.0d0-SpeedRatio)* &
                              DXCoil(DXCoilNum)%EvapCondPumpElecNomPower2
    END IF
    hDelta = TotCap / AirMassFlow
  ! Adjust CBF for off-nominal flow
    CBF = AdjustCBF(CBFNom,AirMassFlowNom,AirMassFlow)
  ! Calculate new apparatus dew point conditions
    hADP = InletAirEnthalpy - hDelta/(1.d0-CBF)
    tADP = PsyTsatFnHPb(hADP,OutdoorPressure)
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!    tADP = PsyTsatFnHPb(hADP,InletAirPressure)
    wADP = PsyWFnTdbH(tADP,hADP)
    hTinwADP = PsyHFnTdbW(InletAirDryBulbTemp,wADP)
  ! get corresponding SHR
    SHR = MIN((hTinwADP-hADP)/(InletAirEnthalpy-hADP),1.d0)
    OutletAirEnthalpy = InletAirEnthalpy - hDelta
  ! get outlet conditions
    hTinwout = InletAirEnthalpy - (1.0d0-SHR)*hDelta
    OutletAirHumRat = PsyWFnTdbH(InletAirDryBulbTemp,hTinwout)
    OutletAirDryBulbTemp = PsyTdbFnHW(OutletAirEnthalpy,OutletAirHumRat)
    ! OutletAirRH = PsyRhFnTdbWPb(OutletAirDryBulbTemp,OutletAirHumRat,OutBaroPress)
    ! IF (OutletAirRH >= 1.) THEN  ! Limit to saturated conditions at OutletAirEnthalpy
    !   OutletAirDryBulbTemp = PsyTsatFnHPb(OutletAirEnthalpy,OutBaroPress)
    !    OutletAirHumRat  = PsyWFnTdbH(OutletAirDryBulbTemp,OutletAirEnthalpy)
    !  END IF
    OutletAirDryBulbTempSat = PsyTsatFnHPb(OutletAirEnthalpy,OutdoorPressure)
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!    OutletAirDryBulbTempSat = PsyTsatFnHPb(OutletAirEnthalpy,InletAirPressure)
    IF (OutletAirDryBulbTemp < OutletAirDryBulbTempSat) THEN ! Limit to saturated conditions at OutletAirEnthalpy
      OutletAirDryBulbTemp = OutletAirDryBulbTempSat
      OutletAirHumRat  = PsyWFnTdbH(OutletAirDryBulbTemp,OutletAirEnthalpy)
    END IF
  ! calculate cooling rate and electrical power
    DXCoil(DXCoilNum)%TotalCoolingEnergyRate = AirMassFlow * (InletAirEnthalpy - OutletAirEnthalpy)
    MinAirHumRat = MIN(InletAirHumRat,OutletAirHumRat)
    DXCoil(DXCoilNum)%SensCoolingEnergyRate = AirMassFlow * (PsyHFnTdbW(InletAirDryBulbTemp,MinAirHumRat) - &
                                                             PsyHFnTdbW(OutletAirDryBulbTemp,MinAirHumRat))
  ! Don't let sensible capacity be greater than total capacity
    IF (DXCoil(DXCoilNum)%SensCoolingEnergyRate > DXCoil(DXCoilNum)%TotalCoolingEnergyRate) THEN
       DXCoil(DXCoilNum)%SensCoolingEnergyRate = DXCoil(DXCoilNum)%TotalCoolingEnergyRate
    END IF
    DXCoil(DXCoilNum)%LatCoolingEnergyRate = DXCoil(DXCoilNum)%TotalCoolingEnergyRate - &
                                             DXCoil(DXCoilNum)%SensCoolingEnergyRate
    DXCoil(DXCoilNum)%ElecCoolingPower = TotCap * EIR
!   Calculation for heat reclaim needs to be corrected to use compressor power (not including condenser fan power)
    HeatReclaimDXCoil(DXCoilNum)%AvailCapacity = DXCoil(DXCoilNum)%TotalCoolingEnergyRate + DXCoil(DXCoilNum)%ElecCoolingPower
    DXCoil(DXCoilNum)%PartLoadRatio    = 1.0d0
    DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction = 1.0d0

    DXCoil(DXCoilNum)%OutletAirEnthalpy = OutletAirEnthalpy
    DXCoil(DXCoilNum)%OutletAirHumRat = OutletAirHumRat
    DXCoil(DXCoilNum)%OutletAirTemp = OutletAirDryBulbTemp

  ELSE IF (CycRatio > 0.0) THEN

    IF (DXCoil(DXCoilNum)%CondenserType(Mode) == EvapCooled) THEN
     ! Outdoor wet-bulb temp from DataEnvironment + (1.0-EvapCondEffectiveness) * (drybulb - wetbulb)
      CondInletTemp = OutdoorWetBulb + (OutdoorDryBulb-OutdoorWetBulb)*(1.0d0 - DXCoil(DXCoilNum)%EvapCondEffect2)
      CondInletHumrat = PsyWFnTdbTwbPb(CondInletTemp,OutdoorWetBulb,OutdoorPressure)
    END IF

    ! Adjust low speed coil bypass factor for actual flow rate.
    ! CBF = AdjustCBF(DXCoil(DXCoilNum)%RatedCBF2,DXCoil(DXCoilNum)%RatedAirMassFlowRate2,AirMassFlow)
    ! get low speed total capacity and SHR at current conditions
    CALL CalcTotCapSHR(InletAirDryBulbTemp,InletAirHumRat,InletAirEnthalpy,InletAirWetbulbC,1.0d0, &
                       DXCoil(DXCoilNum)%RatedAirMassFlowRate2,DXCoil(DXCoilNum)%RatedTotCap2, &
                       DXCoil(DXCoilNum)%RatedCBF2,DXCoil(DXCoilNum)%CCapFTemp2, &
                       DXCoil(DXCoilNum)%CCapFFlow(Mode),TotCapLS,SHRLS,CondInletTemp, &
                       OutdoorPressure)
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!                       Node(DXCoil(DXCoilNum)%AirInNode)%Press)
    hDelta = TotCapLS / AirMassFlow
  ! Adjust CBF for off-nominal flow
    CBF = AdjustCBF(DXCoil(DXCoilNum)%RatedCBF2,DXCoil(DXCoilNum)%RatedAirMassFlowRate2,AirMassFlow)
  ! Calculate new apparatus dew point conditions
    hADP = InletAirEnthalpy - hDelta/(1.d0-CBF)
    tADP = PsyTsatFnHPb(hADP,OutdoorPressure,'CalcMultiSpeedDXCoil:newdewpointconditions')
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!    tADP = PsyTsatFnHPb(hADP,InletAirPressure)
    wADP = PsyWFnTdbH(tADP,hADP,'CalcMultiSpeedDXCoil:newdewpointconditions')
    hTinwADP = PsyHFnTdbW(InletAirDryBulbTemp,wADP,'CalcMultiSpeedDXCoil:newdewpointconditions')
  ! get corresponding SHR
    SHR = MIN((hTinwADP-hADP)/(InletAirEnthalpy-hADP),1.d0)
    ! get low speed outlet conditions
    LSOutletAirEnthalpy = InletAirEnthalpy - hDelta
    hTinwout = InletAirEnthalpy - (1.0d0-SHR)*hDelta
    LSOutletAirHumRat = PsyWFnTdbH(InletAirDryBulbTemp,hTinwout)
    LSOutletAirDryBulbTemp = PsyTdbFnHW(LSOutletAirEnthalpy,LSOutletAirHumRat,'CalcMultiSpeedDXCoil:lowspeedoutlet')
    LSOutletAirRH = PsyRhFnTdbWPb(LSOutletAirDryBulbTemp,LSOutletAirHumRat,OutdoorPressure,'CalcMultiSpeedDXCoil:lowspeedoutlet')
    OutletAirDryBulbTempSat = PsyTsatFnHPb(LSOutletAirEnthalpy,OutdoorPressure,'CalcMultiSpeedDXCoil:lowspeedoutlet')
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!    LSOutletAirRH = PsyRhFnTdbWPb(LSOutletAirDryBulbTemp,LSOutletAirHumRat,InletAirPressure)
!    OutletAirDryBulbTempSat = PsyTsatFnHPb(LSOutletAirEnthalpy,InletAirPressure)
    IF (LSOutletAirDryBulbTemp < OutletAirDryBulbTempSat) THEN ! Limit to saturated conditions at OutletAirEnthalpy
      LSOutletAirDryBulbTemp = OutletAirDryBulbTempSat
      LSOutletAirHumRat  = PsyWFnTdbH(LSOutletAirDryBulbTemp,LSOutletAirEnthalpy,'CalcMultiSpeedDXCoil:lowspeedoutlet')
    END IF
    ! outlet conditions are average of inlet and low speed weighted by CycRatio
    OutletAirEnthalpy = CycRatio*LSOutletAirEnthalpy + (1.d0-CycRatio)*InletAirEnthalpy
    OutletAirHumRat = CycRatio*LSOutletAirHumRat + (1.d0-CycRatio)*InletAirHumRat
    OutletAirDryBulbTemp = PsyTdbFnHW(OutletAirEnthalpy,OutletAirHumRat)
    ! get low speed EIR at current conditions
!    EIRTempModFacLS = CurveValue(DXCoil(DXCoilNum)%EIRFTemp(Mode),InletAirWetbulbC,CondInletTemp)
!    CR7307 changed EIRTempModFacLS calculation to that shown below.
    EIRTempModFacLS = CurveValue(DXCoil(DXCoilNum)%EIRFTemp2,InletAirWetbulbC,CondInletTemp)
    EIRLS = DXCoil(DXCoilNum)%RatedEIR2 * EIRTempModFacLS
    ! get the part load factor that will account for cycling losses
    PLF = CurveValue(DXCoil(DXCoilNum)%PLFFPLR(Mode),CycRatio)
    IF (PLF < 0.7d0) THEN
     PLF = 0.7d0
    END IF
    ! calculate the run time fraction
    DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction = CycRatio / PLF
    DXCoil(DXCoilNum)%PartLoadRatio    = CycRatio
    IF ( DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction > 1.d0 ) THEN
      DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction = 1.0d0 ! Reset coil runtime fraction to 1.0
    END IF
    ! get the eletrical power consumption
    DXCoil(DXCoilNum)%ElecCoolingPower = TotCapLS * EIRLS * DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction
    ! calculate cooling output power
    DXCoil(DXCoilNum)%TotalCoolingEnergyRate = AirMassFlow * (InletAirEnthalpy - OutletAirEnthalpy)
!   Calculation for heat reclaim needs to be corrected to use compressor power (not including condenser fan power)
    HeatReclaimDXCoil(DXCoilNum)%AvailCapacity = DXCoil(DXCoilNum)%TotalCoolingEnergyRate + DXCoil(DXCoilNum)%ElecCoolingPower
    MinAirHumRat = MIN(InletAirHumRat,OutletAirHumRat)
    DXCoil(DXCoilNum)%SensCoolingEnergyRate = AirMassFlow * (PsyHFnTdbW(InletAirDryBulbTemp,MinAirHumRat) - &
                                                             PsyHFnTdbW(OutletAirDryBulbTemp,MinAirHumRat))
  ! Don't let sensible capacity be greater than total capacity
    IF (DXCoil(DXCoilNum)%SensCoolingEnergyRate > DXCoil(DXCoilNum)%TotalCoolingEnergyRate) THEN
       DXCoil(DXCoilNum)%SensCoolingEnergyRate = DXCoil(DXCoilNum)%TotalCoolingEnergyRate
    END IF
    DXCoil(DXCoilNum)%LatCoolingEnergyRate = DXCoil(DXCoilNum)%TotalCoolingEnergyRate - &
                                             DXCoil(DXCoilNum)%SensCoolingEnergyRate
    DXCoil(DXCoilNum)%OutletAirEnthalpy = OutletAirEnthalpy
    DXCoil(DXCoilNum)%OutletAirHumRat = OutletAirHumRat
    DXCoil(DXCoilNum)%OutletAirTemp = OutletAirDryBulbTemp
    CondAirMassFlow =  RhoAir * DXCoil(DXCoilNum)%EvapCondAirFlow2 * DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction
    EvapCondPumpElecPower = DXCoil(DXCoilNum)%EvapCondPumpElecNomPower2 * DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction

  END IF

  IF (DXCoil(DXCoilNum)%CondenserType(Mode) == EvapCooled) THEN
  !******************
  !             WATER CONSUMPTION IN m3 OF WATER FOR DIRECT
  !             H2O [m3/sec] = Delta W[KgH2O/Kg air]*Mass Flow Air[Kg air]
  !                                /RhoWater [kg H2O/m3 H2O]
  !******************
     RhoWater = RhoH2O(OutdoorDryBulb)
     DXCoil(DXCoilNum)%EvapWaterConsumpRate =  (CondInletHumrat - OutdoorHumRat) *  CondAirMassFlow/RhoWater
     DXCoil(DXCoilNum)%EvapCondPumpElecPower = EvapCondPumpElecPower
     !set water system demand request (if needed)
     IF ( DXCoil(DxCoilNum)%EvapWaterSupplyMode == WaterSupplyFromTank) THEN

       WaterStorage(DXCoil(DXCoilNum)%EvapWaterSupTankID)%VdotRequestDemand(DXCoil(DXCoilNum)%EvapWaterTankDemandARRID) &
       = DXCoil(DXCoilNum)%EvapWaterConsumpRate
     ENDIF

     ! Calculate basin heater power
     CALL CalcBasinHeaterPower(DXCoil(DXCoilNum)%BasinHeaterPowerFTempDiff,&
                              DXCoil(DXCoilNum)%BasinHeaterSchedulePtr,&
                              DXCoil(DXCoilNum)%BasinHeaterSetPointTemp,DXCoil(DXCoilNum)%BasinHeaterPower)
     DXCoil(DXCoilNum)%BasinHeaterPower = DXCoil(DXCoilNum)%BasinHeaterPower * &
                                          (1.d0 - DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction)
  ENDIF

ELSE

  ! DX coil is off; just pass through conditions
  DXCoil(DXCoilNum)%OutletAirEnthalpy = DXCoil(DXCoilNum)%InletAirEnthalpy
  DXCoil(DXCoilNum)%OutletAirHumRat = DXCoil(DXCoilNum)%InletAirHumRat
  DXCoil(DXCoilNum)%OutletAirTemp = DXCoil(DXCoilNum)%InletAirTemp

  DXCoil(DXCoilNum)%ElecCoolingPower = 0.0
  DXCoil(DXCoilNum)%TotalCoolingEnergyRate = 0.0
  DXCoil(DXCoilNum)%SensCoolingEnergyRate = 0.0
  DXCoil(DXCoilNum)%LatCoolingEnergyRate = 0.0
  DXCoil(DXCoilNum)%EvapCondPumpElecPower = 0.0
  DXCoil(DXCoilNum)%EvapWaterConsumpRate = 0.0

  ! Calculate basin heater power
  IF (DXCoil(DXCoilNum)%CondenserType(Mode) == EvapCooled) THEN
    CALL CalcBasinHeaterPower(DXCoil(DXCoilNum)%BasinHeaterPowerFTempDiff,&
                              DXCoil(DXCoilNum)%BasinHeaterSchedulePtr,&
                              DXCoil(DXCoilNum)%BasinHeaterSetPointTemp,DXCoil(DXCoilNum)%BasinHeaterPower)
  ENDIF
END IF

DXCoilOutletTemp(DXCoilNum) = DXCoil(DXCoilNum)%OutletAirTemp
DXCoilOutletHumRat(DXCoilNum) = DXCoil(DXCoilNum)%OutletAirHumRat
DXCoil(DXCoilNum)%CondInletTemp = CondInletTemp ! Save condenser inlet temp in the data structure

RETURN

END SUBROUTINE CalcMultiSpeedDXCoil

SUBROUTINE CalcBasinHeaterPowerForMultiModeDXCoil(DXCoilNum, DehumidMode)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Chandan Sharma, FSEC
          !       DATE WRITTEN   May 2010
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! To calculate the basin heater power for multi mode DX cooling coil

          ! METHODOLOGY EMPLOYED:
          ! The methodology employed is as follows:
          ! 1) If the number of capacity stages is equal to 1 and the CondenserType for stage 1
          !    is EvapCooled, then the basin heater power is calculated for (1-runtimefractionstage1) of DX coil
          ! 2) If the number of capacity stages is greater than 1, then
          !    a) If the CondenserType for stage 1 is EvapCooled, then the basin heater power is calculated for
          !       (1-runtimefractionofstage1) of DX coil
          !    b) Elseif the CondenserType for stage 2 is EvapCooled, then the basin heater power is calculated for
          !       (1-runtimefractionofstage2) of DX coil

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
          ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  INTEGER,  INTENT(IN)    :: DehumidMode         ! Dehumidification mode (0=normal, 1=enhanced)
  INTEGER,  INTENT(IN)    :: DXCoilNum           ! Index of coil being simulated

          ! SUBROUTINE PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER :: PerfMode                  ! Performance mode for MultiMode DX coil; Always 1 for other coil types
                                       ! 1-2=normal mode: 1=stage 1 only, 2=stage 1&2
                                       ! 3-4=enhanced dehumidification mode: 3=stage 1 only, 4=stage 1&2

  IF (DXCoil(DXCoilNum)%NumCapacityStages .EQ. 1) THEN
    DXCoil(DXCoilNum)%BasinHeaterPower = DXCoil(DXCoilNum)%BasinHeaterPower * &
                                         (1.d0 - DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction)
  ELSE
    PerfMode = DehumidMode*2 + 1
    IF (DXCoil(DXCoilNum)%CondenserType(PerfMode) .EQ. EvapCooled) THEN
      DXCoil(DXCoilNum)%BasinHeaterPower = DXCoil(DXCoilNum)%BasinHeaterPower * &
                                         (1.d0 - DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction)
    ELSEIF(DXCoil(DXCoilNum)%CondenserType(PerfMode+1) .EQ. EvapCooled) THEN
      CALL CalcBasinHeaterPower(DXCoil(DXCoilNum)%BasinHeaterPowerFTempDiff,&
                          DXCoil(DXCoilNum)%BasinHeaterSchedulePtr,&
                          DXCoil(DXCoilNum)%BasinHeaterSetPointTemp,&
                          DXCoil(DXCoilNum)%BasinHeaterPower)
      DXCoil(DXCoilNum)%BasinHeaterPower = DXCoil(DXCoilNum)%BasinHeaterPower * &
                                       (1.d0 - DXCoil(DXCoilNum)%CoolingCoilStg2RuntimeFrac)
    ENDIF
  ENDIF

RETURN

END SUBROUTINE CalcBasinHeaterPowerForMultiModeDXCoil

FUNCTION AdjustCBF(CBFNom,AirMassFlowRateNom,AirMassFlowRate) RESULT(CBFAdj)

          ! FUNCTION INFORMATION:
          !       AUTHOR         Fred Buhl using Don Shirey's code
          !       DATE WRITTEN   September 2002
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS FUNCTION:
          !    Adjust coil bypass factor for actual air flow rate.

          ! METHODOLOGY EMPLOYED:
          ! Uses relation CBF = exp(-NTU) whereNTU = A0/(m*cp). Relationship models the cooling coil
          ! as a heat exchanger with Cmin/Cmax = 0.
          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
          ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! FUNCTION ARGUMENT DEFINITIONS:
  REAL(r64), INTENT (IN) :: CBFNom                ! nominal coil bypass factor
  REAL(r64), INTENT (IN) :: AirMassFlowRateNom    ! nominal air mass flow rate [kg/s]
  REAL(r64), INTENT (IN) :: AirMassFlowRate       ! actual air mass flow rate [kg/s]
  REAL(r64)         :: CBFAdj                ! the result - the adjusted coil bypass factor


          ! FUNCTION PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  REAL(r64) :: A0  ! intermediate variable
  REAL(r64) :: ADiff  ! intermediate variable

  IF (CBFNom .gt. 0.0) THEN
     A0 = -log(CBFNom)*AirMassFlowRateNom
  ELSE
     A0 = 0.
  END IF
  ADiff=-A0/AirMassFlowRate
  IF (ADiff >= EXP_LowerLimit) THEN
     CBFAdj = exp(ADiff)
  ELSE
     CBFAdj = 0.0
  END IF

  RETURN
END FUNCTION AdjustCBF

FUNCTION CalcCBF(UnitType,UnitName,InletAirTemp,InletAirHumRat,TotCap,AirMassFlowRate,SHR) RESULT(CBF)

          ! FUNCTION INFORMATION:
          !       AUTHOR         Fred Buhl using Don Shirey's code
          !       DATE WRITTEN   September 2002
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS FUNCTION:
          ! Calculate the coil bypass factor for a coil given the total capacity at the entering conditions,
          ! air mass flow rate at the entering conditions, and the sensible heat ratio (SHR) at the
          ! entering conditions.

          ! METHODOLOGY EMPLOYED:
          ! calculate SlopeRated (deltahumrat/deltaT) using rated unit information provided by
          ! user. Then hunt along saturation curve of psychrometric chart until the slope of the line
          ! between the saturation point and rated inlet air humidity ratio and T is the same as SlopeRated.
          ! When the slopes are equal, then we have located the apparatus dewpoint of the coil at rated
          ! conditions. From this information, coil bypass factor is calculated.

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE General, ONLY: RoundSigDigits
  USE DataEnvironment, ONLY: StdRhoAir

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT (IN) :: UnitType
  CHARACTER(len=*), INTENT (IN) :: UnitName
  REAL(r64), INTENT (IN) :: InletAirTemp          ! inlet air temperature [C]
  REAL(r64), INTENT (IN) :: InletAirHumRat        ! inlet air humidity ratio [kg water / kg dry air]
  REAL(r64), INTENT (IN) :: TotCap                ! total cooling  capacity [Watts]
  REAL(r64), INTENT (IN) :: AirMassFlowRate       ! the air mass flow rate at the given capacity [kg/s]
  REAL(r64), INTENT (IN) :: SHR                   ! sensible heat ratio at the given capacity and flow rate
  REAL(r64)         :: CBF                   ! the result - the coil bypass factor


          ! FUNCTION PARAMETER DEFINITIONS:
  CHARACTER(len=*), PARAMETER ::  RoutineName='CalcCBF'
  REAL(r64) :: SmallDifferenceTest=0.00000001d0

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  REAL(r64) :: InletAirEnthalpy  ! Enthalpy of inlet air to evaporator at given conditions [J/kg]
  REAL(r64) :: DeltaH            ! Enthalpy drop across evaporator at given conditions [J/kg]
  REAL(r64) :: DeltaT            ! Temperature drop across evaporator at given conditions [C]
  REAL(r64) :: DeltaHumRat       ! Humidity ratio drop across evaporator at given conditions [kg/kg]
  REAL(r64) :: OutletAirTemp     ! Outlet dry-bulb temperature from evaporator at given conditions [C]
  REAL(r64) :: OutletAirEnthalpy ! Enthalpy of outlet air at given conditions [J/kg]
  REAL(r64) :: OutletAirHumRat   ! Outlet humidity ratio from evaporator at given conditions [kg/kg]
  REAL(r64) :: OutletAirRH       ! relative humidity of the outlet air
  REAL(r64) :: Error                ! Error term used in given coil bypass factor (CBF) calculations
  REAL(r64) :: ErrorLast            ! Error term, from previous iteration
  INTEGER :: Iter                 ! Iteration loop counter in CBF calculations
  INTEGER :: IterMax              ! Maximum number of iterations in CBF calculations
  REAL(r64) :: ADPTemp              ! Apparatus dewpoint temperature used in CBF calculations [C]
  REAL(r64) :: ADPHumRat            ! Apparatus dewpoint humidity used in CBF calculations [kg/kg]
  REAL(r64) :: ADPEnthalpy          ! Air enthalpy at apparatus dew point [J/kg]
  REAL(r64) :: DeltaADPTemp         ! Change in Apparatus Dew Point used in CBF calculations [C]
  REAL(r64) :: SlopeAtConds          ! Slope (DeltaHumRat/DeltaT) at given conditions
  REAL(r64) :: Slope                ! Calculated Slope used while hunting for Tadp
  REAL(r64) :: Tolerance            ! Convergence tolerance for CBF calculations
  REAL(r64) :: HTinHumRatOut        ! Air enthalpy at inlet air temp and outlet air humidity ratio [J/kg]
  LOGICAL :: CBFErrors=.false.    ! Set to true if errors in CBF calculation, fatal at end of routine

  DeltaH = 0.0
  DeltaT = 0.0
  DeltaHumRat = 0.0
  OutletAirTemp =  InletAirTemp
  OutletAirHumRat = InletAirHumRat
  SlopeAtConds = 0.0
  Slope = 0.0
  IterMax = 50
  CBFErrors=.false.

  DeltaH = TotCap/AirMassFlowRate
  InletAirEnthalpy = PsyHFnTdbW(InletAirTemp,InletAirHumRat)
  HTinHumRatOut = InletAirEnthalpy - (1.0d0-SHR)*DeltaH
  OutletAirHumRat = PsyWFnTdbH(InletAirTemp,HTinHumRatOut)
  DeltaHumRat = InletAirHumRat - OutletAirHumRat
  OutletAirEnthalpy = InletAirEnthalpy - DeltaH
  OutletAirTemp = PsyTdbFnHW(OutletAirEnthalpy,OutletAirHumRat)
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!  Pressure will have to be pass into this subroutine to fix this one
  OutletAirRH = PsyRhFnTdbWPb(OutletAirTemp,OutletAirHumRat,StdBaroPress,'CalcCBF')
  IF (OutletAirRH .ge. 1.0d0) THEN
     CALL ShowSevereError ('For object = '//TRIM(UnitType)// ', name = "'//TRIM(UnitName)// '"')
     CALL ShowContinueError ('Calculated outlet air relative humidity greater than 1. The combination of')
     CALL ShowContinueError ('rated air volume flow rate, total cooling capacity and sensible heat ratio yields coil exiting')
     CALL ShowContinueError ('air conditions above the saturation curve. Possible fixes are to reduce the rated total cooling')
     CALL ShowContinueError ('capacity, increase the rated air volume flow rate, or reduce the rated sensible heat'// &
                             ' ratio for this coil.')
     CALL ShowContinueError ('If autosizing, it is recommended that all three of these values be autosized.')
     CALL ShowContinueError('...Inputs used for calculating cooling coil bypass factor.')
     CALL ShowContinueError('...Inlet Air Temperature     = '//TRIM(RoundSigDigits(InletAirTemp,2))//' C')
     CALL ShowContinueError('...Outlet Air Temperature    = '//TRIM(RoundSigDigits(OutletAirTemp,2))//' C')
     CALL ShowContinueError('...Inlet Air Humidity Ratio  = '//TRIM(RoundSigDigits(InletAirHumRat,6))//' kgWater/kgDryAir')
     CALL ShowContinueError('...Outlet Air Humidity Ratio = '//TRIM(RoundSigDigits(OutletAirHumRat,6))//' kgWater/kgDryAir')
     CALL ShowContinueError('...Total Cooling Capacity used in calculation = '//TRIM(RoundSigDigits(TotCap,2))//' W')
     CALL ShowContinueError('...Air Mass Flow Rate used in calculation     = '//TRIM(RoundSigDigits(AirMassFlowRate,6))//' kg/s')
     CALL ShowContinueError('...Air Volume Flow Rate used in calculation   = '// &
       TRIM(RoundSigDigits(AirMassFlowRate/PsyRhoAirFnPbTdbW(StdBaroPress,InletAirTemp,InletAirHumRat,RoutineName),6))//' m3/s')
     IF(TotCap .GT. 0.d0)THEN
       IF (((MinRatedVolFlowPerRatedTotCap - AirMassFlowRate/ &
            PsyRhoAirFnPbTdbW(StdBaroPress,InletAirTemp,InletAirHumRat,RoutineName)/TotCap) > SmallDifferenceTest).OR. &
           ((AirMassFlowRate/PsyRhoAirFnPbTdbW(StdBaroPress,InletAirTemp,InletAirHumRat,RoutineName)/TotCap &
             - MaxRatedVolFlowPerRatedTotCap) > SmallDifferenceTest)) THEN
         CALL ShowContinueError('...Air Volume Flow Rate per Watt of Rated Cooling Capacity is also out of bounds at = '// &
                                TRIM(RoundSigDigits(AirMassFlowRate/ &
                                PsyRhoAirFnPbTdbW(StdBaroPress,InletAirTemp,InletAirHumRat,RoutineName)/TotCap,7))//' m3/s/W')
       END IF
     END IF
     CALL ShowContinueErrorTimeStamp(' ')
     CALL ShowFatalError ('Check and revise the input data for this coil before rerunning the simulation.')
  END IF
  DeltaT = InletAirTemp - OutletAirTemp
  IF (DeltaT .LE. 0.0d0) THEN
     CALL ShowSevereError ('For object = '//TRIM(UnitType)// ', name = "'//TRIM(UnitName)// '"')
     CALL ShowContinueError ('Calculated coil delta T is less than or equal to 0. The combination of')
     CALL ShowContinueError ('rated air volume flow rate, total cooling capacity and sensible heat ratio yields coil exiting')
     CALL ShowContinueError ('air conditions that are not reasonable. Possible fixes are to adjust the rated total cooling')
     CALL ShowContinueError ('capacity, rated air volume flow rate, or rated sensible heat'// &
                             ' ratio for this coil.')
     CALL ShowContinueError ('If autosizing, it is recommended that all three of these values be autosized.')
     CALL ShowContinueError('...Inputs used for calculating cooling coil bypass factor.')
     CALL ShowContinueError('...Inlet Air Temperature     = '//TRIM(RoundSigDigits(InletAirTemp,2))//' C')
     CALL ShowContinueError('...Outlet Air Temperature    = '//TRIM(RoundSigDigits(OutletAirTemp,2))//' C')
     CALL ShowContinueError('...Inlet Air Humidity Ratio  = '//TRIM(RoundSigDigits(InletAirHumRat,6))//' kgWater/kgDryAir')
     CALL ShowContinueError('...Outlet Air Humidity Ratio = '//TRIM(RoundSigDigits(OutletAirHumRat,6))//' kgWater/kgDryAir')
     CALL ShowContinueError('...Total Cooling Capacity used in calculation = '//TRIM(RoundSigDigits(TotCap,2))//' W')
     CALL ShowContinueError('...Air Mass Flow Rate used in calculation     = '//TRIM(RoundSigDigits(AirMassFlowRate,6))//' kg/s')
     CALL ShowContinueError('...Air Volume Flow Rate used in calculation   = '// &
       TRIM(RoundSigDigits(AirMassFlowRate/PsyRhoAirFnPbTdbW(StdBaroPress,InletAirTemp,InletAirHumRat,RoutineName),6))//' m3/s')
     IF(TotCap .GT. 0.d0)THEN
       IF (((MinRatedVolFlowPerRatedTotCap - AirMassFlowRate/ &
            PsyRhoAirFnPbTdbW(StdBaroPress,InletAirTemp,InletAirHumRat,RoutineName)/TotCap) > SmallDifferenceTest).OR. &
           ((AirMassFlowRate/PsyRhoAirFnPbTdbW(StdBaroPress,InletAirTemp,InletAirHumRat,RoutineName)/TotCap &
             - MaxRatedVolFlowPerRatedTotCap) > SmallDifferenceTest)) THEN
         CALL ShowContinueError('...Air Volume Flow Rate per Watt of Rated Cooling Capacity is also out of bounds at = '// &
                                TRIM(RoundSigDigits(AirMassFlowRate/ &
                                PsyRhoAirFnPbTdbW(StdBaroPress,InletAirTemp,InletAirHumRat,RoutineName)/TotCap,7))//' m3/s/W')
       END IF
     END IF
     CALL ShowContinueErrorTimeStamp(' ')
     CALL ShowFatalError ('Check and revise the input data for this coil before rerunning the simulation.')
  END IF
  ! Calculate slope at given conditions
  IF (DeltaT .gt. 0.0) SlopeAtConds = DeltaHumRat/DeltaT

!  IF (SlopeAtConds .le. .0000001d0 .or. OutletAirHumRat .le. 0.) THEN
  IF (SlopeAtConds .lt. 0.0d0 .or. OutletAirHumRat .le. 0.) THEN
!   Invalid conditions, slope can't be less than zero (SHR > 1) or
!   outlet air humidity ratio can't be less than zero.
    CALL ShowSevereError(TRIM(UnitType)//' "'//TRIM(UnitName)//'"')
    CALL ShowContinueError('...Invalid slope or outlet air condition when calculating cooling coil bypass factor.')
    CALL ShowContinueError('...Slope = '//TRIM(RoundSigDigits(SlopeAtConds,8)))
    CALL ShowContinueError('...Inlet Air Temperature     = '//TRIM(RoundSigDigits(InletAirTemp,2))//' C')
    CALL ShowContinueError('...Outlet Air Temperature    = '//TRIM(RoundSigDigits(OutletAirTemp,2))//' C')
    CALL ShowContinueError('...Inlet Air Humidity Ratio  = '//TRIM(RoundSigDigits(InletAirHumRat,6))//' kgWater/kgDryAir')
    CALL ShowContinueError('...Outlet Air Humidity Ratio = '//TRIM(RoundSigDigits(OutletAirHumRat,6))//' kgWater/kgDryAir')
    CALL ShowContinueError('...Total Cooling Capacity used in calculation = '//TRIM(RoundSigDigits(TotCap,2))//' W')
    CALL ShowContinueError('...Air Mass Flow Rate used in calculation     = '//TRIM(RoundSigDigits(AirMassFlowRate,6))//' kg/s')
     CALL ShowContinueError('...Air Volume Flow Rate used in calculation   = '// &
       TRIM(RoundSigDigits(AirMassFlowRate/PsyRhoAirFnPbTdbW(StdBaroPress,InletAirTemp,InletAirHumRat,RoutineName),6))//' m3/s')
    IF(TotCap .GT. 0.d0)THEN
       IF (((MinRatedVolFlowPerRatedTotCap - AirMassFlowRate/ &
            PsyRhoAirFnPbTdbW(StdBaroPress,InletAirTemp,InletAirHumRat,RoutineName)/TotCap) > SmallDifferenceTest).OR. &
           ((AirMassFlowRate/PsyRhoAirFnPbTdbW(StdBaroPress,InletAirTemp,InletAirHumRat,RoutineName)/TotCap &
             - MaxRatedVolFlowPerRatedTotCap) > SmallDifferenceTest)) THEN
         CALL ShowContinueError('...Air Volume Flow Rate per Watt of Rated Cooling Capacity is also out of bounds at = '// &
                                TRIM(RoundSigDigits(AirMassFlowRate/ &
                                PsyRhoAirFnPbTdbW(StdBaroPress,InletAirTemp,InletAirHumRat,RoutineName)/TotCap,7))//' m3/s/W')
       END IF
    END IF
    CALL ShowContinueErrorTimeStamp(' ')
    CBFErrors=.true.
  ELSE

!   First guess for Tadp is outlet air dew point
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!  Pressure will have to be pass into this subroutine to fix this one
    ADPTemp = PsyTdpFnWPb(OutletAirHumRat,StdBaroPress)

    Tolerance = 1.         ! initial conditions for iteration
    ErrorLast = 100.d0
    Iter = 0
    DeltaADPTemp = 5.0d0
    DO WHILE ((Iter .le. IterMax).and.(Tolerance .gt. .001d0))
!     Do for IterMax iterations or until the error gets below .1%
      IF (Iter .gt. 0) ADPTemp = ADPTemp + DeltaADPTemp
      Iter = Iter + 1

!     Find new slope using guessed Tadp

!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!  Pressure will have to be pass into this subroutine to fix this one
      ADPHumRat = PsyWFnTdpPb(ADPTemp,StdBaroPress)
      Slope     = (InletAirHumRat-ADPHumRat)/(InletAirTemp-ADPTemp)

!     check for convergence (slopes are equal to within error tolerance)

      Error     = (Slope-SlopeAtConds)/SlopeAtConds
      IF ((Error .gt. 0.).and.(ErrorLast .lt. 0.)) DeltaADPTemp = -DeltaADPTemp/2.d0
      IF ((Error .lt. 0.).and.(ErrorLast .gt. 0.)) DeltaADPTemp = -DeltaADPTemp/2.d0
      ErrorLast = Error

      Tolerance = ABS(Error)

    END DO

!   Calculate Bypass Factor from Enthalpies

    InletAirEnthalpy=PsyHFnTdbW(InletAirTemp,InletAirHumRat)
    OutletAirEnthalpy=PsyHFnTdbW(OutletAirTemp,OutletAirHumRat)
    ADPEnthalpy=PsyHFnTdbW(ADPTemp,ADPHumRat)
    CBF = (OutletAirEnthalpy-ADPEnthalpy)/(InletAirEnthalpy-ADPEnthalpy)
    IF (Iter .gt. IterMax) THEN
      CALL ShowSevereError(TRIM(UnitType)//' "'//TRIM(UnitName)//&
                          '" -- coil bypass factor calculation did not converge after max iterations.')
      CALL ShowContinueError('The RatedSHR of ['//TRIM(RoundSigDigits(SHR,3))//  &
         '], entered by the user or autosized (see *.eio file),')
      CALL ShowContinueError('may be causing this. The line defined by the coil rated inlet air conditions')
      CALL ShowContinueError('(26.7C drybulb and 19.4C wetbulb) and the RatedSHR (i.e., slope of the line) must intersect')
      CALL ShowContinueError('the saturation curve of the psychrometric chart. If the RatedSHR is too low, then this')
      CALL ShowContinueError('intersection may not occur and the coil bypass factor calculation will not converge.')
      CALL ShowContinueError('If autosizing the SHR, recheck the design supply air humidity ratio and design supply air')
      CALL ShowContinueError('temperature values in the Sizing:System and Sizing:Zone objects. In general, the temperatures')
      CALL ShowContinueError('and humidity ratios specified in these two objects should be the same for each system')
      CALL ShowContinueError('and the zones that it serves.')
      CALL ShowContinueErrorTimeStamp(' ')
      CBFErrors=.true.  ! Didn't converge within MaxIter iterations
    ENDIF
    IF (CBF .lt. 0.) THEN
      CALL ShowSevereError(TRIM(UnitType)//' "'//TRIM(UnitName)//'" -- negative coil bypass factor calculated.')
      CALL ShowContinueErrorTimeStamp(' ')
      CBFErrors=.true. ! Negative CBF not valid
    ENDIF
  END IF

! Show fatal error for specific coil that caused a CBF error
  IF (CBFErrors) THEN
    CALL ShowFatalError(TRIM(UnitType)//' "'//TRIM(UnitName)//&
                        '" Errors found in calculating coil bypass factors')
  END IF

  RETURN
END FUNCTION CalcCBF

FUNCTION CalcEffectiveSHR(DXCoilNum, SHRss, RTF, QLatRated, QLatActual, EnteringDB, EnteringWB, Mode, HeatingRTF) RESULT(SHReff)

        ! FUNCTION INFORMATION:
        !    AUTHOR         Richard Raustad, FSEC
        !    DATE WRITTEN   September 2003
        !                   Feb 2005 M. J. Witte, GARD Analytics, Inc.
        !                    Add new coil type COIL:DX:MultiMode:CoolingEmpirical:
        !                   Nov 2008 R. Raustad, FSEC
        !                    Modified to allow latent degradation with cycling fan
        !    RE-ENGINEERED  na

        ! PURPOSE OF THIS FUNCTION:
        !    Adjust sensible heat ratio to account for degradation of DX coil latent
        !    capacity at part-load (cycling) conditions.

        ! METHODOLOGY EMPLOYED:
        !    With model parameters entered by the user, the part-load latent performance
        !    of a DX cooling coil is determined for a constant air flow system with
        !    a cooling coil that cycles on/off. The model calculates the time
        !    required for condensate to begin falling from the cooling coil.
        !    Runtimes greater than this are integrated to a "part-load" latent
        !    capacity which is used to determine the "part-load" sensible heat ratio.
        !    See reference below for additional details (linear decay model, Eq. 8b).
        !
        ! REFERENCES:
        !   "A Model to Predict the Latent Capacity of Air Conditioners and
        !    Heat Pumps at Part-Load Conditions with Constant Fan Operation"
        !    1996 ASHRAE Transactions, Volume 102, Part 1, Pp. 266 - 274,
        !    Hugh I. Henderson, Jr., P.E., Kannan Rengarajan, P.E.

        ! USE STATEMENTS:
        !    na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! FUNCTION ARGUMENT DEFINITIONS:
  INTEGER, INTENT (IN) :: DXCoilNum   ! Index number for cooling coil
  REAL(r64), INTENT (IN) :: SHRss       ! Steady-state sensible heat ratio
  REAL(r64), INTENT (IN) :: RTF         ! Compressor run-time fraction
  REAL(r64), INTENT (IN) :: QLatRated   ! Rated latent capacity
  REAL(r64), INTENT (IN) :: QLatActual  ! Actual latent capacity
  REAL(r64), INTENT (IN) :: EnteringDB  ! Entering air dry-bulb temperature
  REAL(r64), INTENT (IN) :: EnteringWB  ! Entering air wet-bulb temperature
  INTEGER, INTENT(IN), OPTIONAL :: Mode  ! Performance mode for MultiMode DX coil; Always 1 for other coil types
  REAL(r64), INTENT (IN), OPTIONAL :: HeatingRTF ! Used to recalculate Toff for cycling fan systems
  REAL(r64)            :: SHReff      ! Effective sensible heat ratio, includes degradation due to cycling effects

          ! FUNCTION PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  REAL(r64) :: Twet        ! Nominal time for condensate to begin leaving the coil's condensate drain line
                         !   at the current operating conditions (sec)
  REAL(r64) :: Gamma       ! Initial moisture evaporation rate divided by steady-state AC latent capacity
                         !   at the current operating conditions
  REAL(r64) :: Twet_rated  ! Twet at rated conditions (coil air flow rate and air temperatures), sec
  REAL(r64) :: Gamma_rated ! Gamma at rated conditions (coil air flow rate and air temperatures)
  REAL(r64) :: Twet_max    ! Maximum allowed value for Twet
  REAL(r64) :: Nmax        ! Maximum ON/OFF cycles per hour for the compressor (cycles/hr)
  REAL(r64) :: Tcl         ! Time constant for latent capacity to reach steady state after startup (sec)
  REAL(r64) :: Ton         ! Coil on time (sec)
  REAL(r64) :: Toff        ! Coil off time (sec)
  REAL(r64) :: Toffa       ! Actual coil off time (sec). Equations valid for Toff <= (2.0 * Twet/Gamma)
  REAL(r64) :: aa          ! Intermediate variable
  REAL(r64) :: To1         ! Intermediate variable (first guess at To). To = time to the start of moisture removal
  REAL(r64) :: To2         ! Intermediate variable (second guess at To). To = time to the start of moisture removal
  REAL(r64) :: Error       ! Error for iteration (DO) loop
  REAL(r64) :: LHRmult     ! Latent Heat Ratio (LHR) multiplier. The effective latent heat ratio LHR = (1-SHRss)*LHRmult
  REAL(r64) :: Ton_heating
  REAL(r64) :: Toff_heating

  If (DXCoil(DXCoilNum)%DXCoilType_Num .NE. CoilDX_MultiSpeedCooling) Then
   Twet_rated  = DXCoil(DXCoilNum)%Twet_Rated(Mode)
   Gamma_rated = DXCoil(DXCoilNum)%Gamma_Rated(Mode)
   Nmax        = DXCoil(DXCoilNum)%MaxONOFFCyclesperHour(Mode)
   Tcl         = DXCoil(DXCoilNum)%LatentCapacityTimeConstant(Mode)
  Else
   Twet_rated  = DXCoil(DXCoilNum)%MSTwet_Rated(Mode)
   Gamma_rated = DXCoil(DXCoilNum)%MSGamma_Rated(Mode)
   Nmax        = DXCoil(DXCoilNum)%MSMaxONOFFCyclesperHour(Mode)
   Tcl         = DXCoil(DXCoilNum)%MSLatentCapacityTimeConstant(Mode)
  End If

!  No moisture evaporation (latent degradation) occurs for runtime fraction of 1.0
!  All latent degradation model parameters cause divide by 0.0 if not greater than 0.0
!  Latent degradation model parameters initialize to 0.0 meaning no evaporation model used.
   IF(RTF .GE. 1.0d0 .OR. Twet_rated .LE. 0.0 .OR. &
      Gamma_rated .LE. 0.0 .OR. Nmax .LE. 0.0 .OR. Tcl .LE. 0.0) THEN
     SHReff = SHRss
     RETURN
   ENDIF

   Twet_max   = 9999.0d0 ! high limit for Twet

!  Calculate the model parameters at the actual operating conditions
   Twet    = MIN(Twet_rated*QLatRated /(QLatActual+1.d-10),Twet_max)
   Gamma   = Gamma_rated*QLatRated*(EnteringDB-EnteringWB)/((26.7d0-19.4d0)*QLatActual+1.d-10)

!  Calculate the compressor on and off times using a converntional thermostat curve
   Ton  = 3600.d0/(4.d0*Nmax*(1.-RTF))   ! duration of cooling coil on-cycle (sec)
   Toff = 3600.d0/(4.d0*Nmax*RTF)        ! duration of cooling coil off-cycle (sec)

!  Cap Toff to meet the equation restriction
   IF(Gamma .GT. 0.0)THEN
     Toffa = MIN(Toff, 2.d0*Twet/Gamma)
   ELSE
     Toffa = Toff
   END IF

!  Need to include the reheat coil operation to account for actual fan run time. E+ uses a
!  separate heating coil for heating and reheat (to separate the heating and reheat loads)
!  and real world applications would use a single heating coil for both purposes, the actual
!  fan operation is based on HeatingPLR + ReheatPLR. For cycling fan RH control, latent
!  degradation only occurs when a heating load exists, in this case the reheat load is
!  equal to and oposite in magnitude to the cooling coil sensible output but the reheat
!  coil is not always active. This additional fan run time has not been accounted for at this time.
!
!  Recalculate Toff for cycling fan systems when heating is active
   IF (PRESENT(HeatingRTF)) THEN
     IF (HeatingRTF .LT. 1.0d0 .AND. HeatingRTF .GT. RTF)THEN
       Ton_heating = 3600.d0/(4.d0*Nmax*(1.d0-HeatingRTF))
       Toff_heating = 3600.d0/(4.d0*Nmax*HeatingRTF)
!    add additional heating coil operation during cooling coil off cycle (due to cycling rate difference of coils)
       Ton_heating = Ton_heating + MAX(0.0d0,MIN(Ton_heating, (Ton+Toffa)-(Ton_heating+Toff_heating)))
       Toffa = MIN(Toffa,Ton_heating - Ton)
     END IF
   END IF

!  Use sucessive substitution to solve for To
   aa = (Gamma*Toffa) - (0.25d0/Twet)*(Gamma**2)*(Toffa**2)
   To1 = aa+Tcl
   Error = 1.0d0
   DO WHILE (Error .gt. 0.001d0)
       To2 = aa-Tcl*(EXP(-To1/Tcl)-1.d0)
       Error = ABS((To2-To1)/To1)
       To1 = To2
   END DO

!  Adjust Sensible Heat Ratio (SHR) using Latent Heat Ratio (LHR) multiplier
!  Floating underflow errors occur when -Ton/Tcl is a large negative number.
!  Cap lower limit at -700 to avoid the underflow errors.
   aa = EXP(MAX(-700.0d0,-Ton/Tcl))
!  Calculate latent heat ratio multiplier
   LHRmult = MAX(((Ton-To2)/(Ton+Tcl*(aa-1.0d0))),0.0d0)

!  Calculate part-load or "effective" sensible heat ratio
   SHReff = 1.0d0-(1.0d0-SHRss)*LHRmult

   IF (SHReff .LT. SHRss) SHReff = SHRss ! Effective SHR can be less than the steady-state SHR
   IF (SHReff .GT. 1.0d0) SHReff=1.0d0 ! Effective sensible heat ratio can't be greater than 1.0

 RETURN

END FUNCTION CalcEffectiveSHR

SUBROUTINE CalcTotCapSHR(InletDryBulb,InletHumRat,InletEnthalpy,InletWetBulb,AirMassFlowRatio,&
                         AirMassFlow,TotCapNom,CBF,CCapFTemp,CCapFFlow,TotCap,SHR,CondInletTemp, Pressure)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Fred Buhl using Don Shirey's code
          !       DATE WRITTEN   September 2002
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! Calculates total capacity and sensible heat ratio of a DX coil at the specified conditions

          ! METHODOLOGY EMPLOYED:
          ! With the rated performance data entered by the user, the model employs some of the
          ! DOE-2.1E curve fits to adjust the capacity and SHR of the unit as a function
          ! of entering air temperatures and supply air flow rate (actual vs rated flow). The model
          ! does NOT employ the exact same methodology to calculate performance as DOE-2, although
          ! some of the DOE-2 curve fits are employed by this model.

          ! The model checks for coil dryout conditions, and adjusts the calculated performance
          ! appropriately.

          ! REFERENCES:
          ! ASHRAE HVAC 2 Toolkit page 4-81.
          !
          ! Henderson, H.I. Jr., K. Rengarajan and D.B. Shirey, III. 1992.The impact of comfort
          ! control on air conditioner energy use in humid climates. ASHRAE Transactions 98(2):
          ! 104-113.
          !
          ! Henderson, H.I. Jr., Danny Parker and Y.J. Huang. 2000.Improving DOE-2's RESYS routine:
          ! User Defined Functions to Provide More Accurate Part Load Energy Use and Humidity
          ! Predictions. Proceedings of ACEEE Conference.


          ! USE STATEMENTS:
  USE CurveManager, ONLY: CurveValue

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  REAL(r64), INTENT (IN) :: InletDryBulb       ! inlet air dry bulb temperature [C]
  REAL(r64), INTENT (IN) :: InletHumRat        ! inlet air humidity ratio [kg water / kg dry air]
  REAL(r64), INTENT (IN) :: InletEnthalpy      ! inlet air specific enthalpy [J/kg]
  REAL(r64), INTENT (IN) :: InletWetBulb       ! inlet air wet bulb temperature [C]
  REAL(r64), INTENT (IN) :: AirMassFlowRatio   ! Ratio of actual air mass flow to nominal air mass flow
  REAL(r64), INTENT (IN) :: AirMassFlow        ! actual mass flow for capacity and SHR calculation
  REAL(r64), INTENT (IN) :: TotCapNom          ! nominal total capacity [W]
  REAL(r64), INTENT (IN) :: CBF                ! coil bypass factor
  INTEGER, INTENT (IN) :: CCapFTemp          ! capacity modifier curve index, function of entering wetbulb
                                             ! and outside drybulb
  INTEGER, INTENT (IN) :: CCapFFlow          ! capacity modifier curve, function of actual flow vs rated flow
  REAL(r64), INTENT (OUT)   :: TotCap             ! total capacity at the given conditions [W]
  REAL(r64), INTENT (OUT)   :: SHR                ! sensible heat ratio at the given conditions
  REAL(r64), INTENT (IN) :: CondInletTemp      ! Condenser inlet temperature [C]
  REAL(r64), INTENT (IN) :: Pressure           ! air pressure [Pa]


          ! SUBROUTINE PARAMETER DEFINITIONS:
  CHARACTER(len=*), PARAMETER :: RoutineName='CalcTotCapSHR'

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  REAL(r64) :: InletWetBulbCalc    ! calculated inlet wetbulb temperature used for finding dry coil point [C]
  REAL(r64) :: InletHumRatCalc     ! calculated inlet humidity ratio used for finding dry coil point [kg water / kg dry air]
  REAL(r64) :: TotCapTempModFac    ! Total capacity modifier (function of entering wetbulb, outside drybulb)
  REAL(r64) :: TotCapFlowModFac    ! Total capacity modifier (function of actual supply air flow vs nominal flow)
  REAL(r64) :: hDelta              ! Change in air enthalpy across the cooling coil [J/kg]
  REAL(r64) :: hADP                ! Apparatus dew point enthalpy [J/kg]
  REAL(r64) :: tADP                ! Apparatus dew point temperature [C]
  REAL(r64) :: wADP                ! Apparatus dew point humidity ratio [kg/kg]
  REAL(r64) :: hTinwADP            ! Enthalpy at inlet dry-bulb and wADP [J/kg]
  REAL(r64) :: SHRCalc             ! temporary calculated value of SHR
  REAL(r64) :: TotCapCalc          ! temporary calculated value of total capacity [W]
  INTEGER :: Counter             ! Counter for dry evaporator iterations
  INTEGER :: MaxIter             ! Maximum number of iterations for dry evaporator calculations
  REAL(r64) :: RF                  ! Relaxation factor for dry evaporator iterations
  REAL(r64) :: Tolerance           ! Error tolerance for dry evaporator iterations
  REAL(r64) :: werror              ! Deviation of humidity ratio in dry evaporator iteration loop

  MaxIter = 30
  RF = 0.4d0
  Counter = 0
  Tolerance = 0.01d0
  werror = 0.0

  InletWetBulbCalc = InletWetBulb
  InletHumRatCalc = InletHumRat

!  DO WHILE (ABS(werror) .gt. Tolerance .OR. Counter == 0)
!   Get capacity modifying factor (function of inlet wetbulb & outside drybulb) for off-rated conditions
50  TotCapTempModFac = CurveValue(CCapFTemp,InletWetBulbCalc,CondInletTemp)
!   Get capacity modifying factor (function of mass flow) for off-rated conditions
    TotCapFlowModFac = CurveValue(CCapFFlow,AirMassFlowRatio)
!   Get total capacity
    TotCapCalc = TotCapNom * TotCapFlowModFac * TotCapTempModFac

!   Calculate apparatus dew point conditions using TotCap and CBF
    hDelta = TotCapCalc/AirMassFlow
    hADP = InletEnthalpy - hDelta/(1.d0-CBF)
    tADP = PsyTsatFnHPb(hADP,Pressure)
    wADP = PsyWFnTdbH(tADP,hADP)
    hTinwADP = PsyHFnTdbW(InletDryBulb,wADP)
    SHRCalc = MIN((hTinwADP-hADP)/(InletEnthalpy-hADP),1.d0)
!
!   Check for dry evaporator conditions (win < wadp)
!
    IF (wADP .gt. InletHumRatCalc .or. (Counter .ge. 1 .and. Counter .lt. MaxIter)) THEN
      If(InletHumRatCalc == 0.0)InletHumRatCalc=0.00001d0
      werror = (InletHumRatCalc - wADP)/InletHumRatCalc
!
!     Increase InletHumRatCalc at constant inlet air temp to find coil dry-out point. Then use the
!     capacity at the dry-out point to determine exiting conditions from coil. This is required
!     since the TotCapTempModFac doesn't work properly with dry-coil conditions.
!
      InletHumRatCalc = RF*wADP + (1.d0-RF)*InletHumRatCalc
      InletWetBulbCalc = PsyTwbFnTdbWPb(InletDryBulb,InletHumRatCalc,Pressure)
      Counter = Counter + 1
      IF (ABS(werror) .gt. Tolerance) go to 50   ! Recalculate with modified inlet conditions

    END IF

! END DO

!  Calculate full load output conditions
  IF (SHRCalc .gt. 1.d0 .OR. Counter .gt. 0) SHRCalc = 1.d0

  SHR = SHRCalc
  TotCap = TotCapCalc

  RETURN
END SUBROUTINE CalcTotCapSHR

SUBROUTINE CalcMultiSpeedDXCoilCooling(DXCoilNum,SpeedRatio, CycRatio, SpeedNum, FanOpMode, CompOp)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Lixing Gu, FSEC
          !       DATE WRITTEN   June 2007
          !       MODIFIED       April 2010, Chandan sharma, FSEC, added basin heater
          !       RE-ENGINEERED  Revised based on CalcMultiSpeedDXCoil

          ! PURPOSE OF THIS SUBROUTINE:
          ! Calculates the air-side performance and electrical energy use of a direct-
          ! expansion, air-cooled cooling unit with a multispeed compressor.

          ! METHODOLOGY EMPLOYED:
          ! Uses the same methodology as the single speed DX unit model (SUBROUTINE CalcDoe2DXCoil).
          ! In addition it assumes that the unit performance is obtained by interpolating between
          ! the performance at high speed and that at low speed. If the output needed is below
          ! that produced at low speed, the compressor cycles between off and low speed.

          ! USE STATEMENTS:
USE CurveManager, ONLY: CurveValue
USE DataWater,    ONLY: WaterStorage
USE DataHVACGlobals,     ONLY: MSHPMassFlowRateLow, MSHPMassFlowRateHigh, MSHPWasteHeat
USE General,      ONLY: TrimSigDigits, RoundSigDigits

IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
INTEGER   :: DXCoilNum     ! the number of the DX heating coil to be simulated
REAL(r64) :: SpeedRatio    ! = (CompressorSpeed - CompressorSpeedMin) / (CompressorSpeedMax - CompressorSpeedMin)
                         ! SpeedRatio varies between 1.0 (maximum speed) and 0.0 (minimum speed)
REAL(r64) :: CycRatio      ! cycling part load ratio
INTEGER   :: SpeedNum      ! Speed number
INTEGER   :: FanOpMode     ! Sets fan control to CycFanCycCoil or ContFanCycCoil
INTEGER   :: CompOp        ! Compressor on/off; 1=on, 0=off

          ! SUBROUTINE PARAMETER DEFINITIONS:
CHARACTER(len=*), PARAMETER :: RoutineName='CalcMultiSpeedDXCoilCooling'

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
REAL(r64) :: AirMassFlow         ! dry air mass flow rate through coil [kg/s]
REAL(r64) :: InletAirWetBulbC    ! wetbulb temperature of inlet air [C]
REAL(r64) :: InletAirDryBulbTemp ! inlet air dry bulb temperature [C]
REAL(r64) :: InletAirEnthalpy    ! inlet air enthalpy [J/kg]
REAL(r64) :: InletAirHumRat      ! inlet air humidity ratio [kg/kg]
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!REAL(r64)   :: InletAirPressure    ! inlet air pressure [Pa]
REAL(r64) :: OutletAirDryBulbTemp    ! outlet air dry bulb temperature [C]
REAL(r64) :: OutletAirEnthalpy   ! outlet air enthalpy [J/kg]
REAL(r64) :: OutletAirHumRat     ! outlet air humidity ratio [kg/kg]
!REAL(r64)   :: OutletAirRH         ! outlet air relative humudity [fraction]
REAL(r64) :: OutletAirDryBulbTempSat ! outlet air dry bulb temp at saturation at the outlet enthalpy [C]
REAL(r64) :: LSOutletAirDryBulbTemp ! low speed outlet air dry bulb temperature [C]
REAL(r64) :: LSOutletAirEnthalpy    ! low speed outlet air enthalpy [J/kg]
REAL(r64) :: LSOutletAirHumRat      ! low speed outlet air humidity ratio [kg/kg]
REAL(r64) :: LSOutletAirRH          ! low speed outlet air relative humudity [fraction]
REAL(r64) :: HSOutletAirDryBulbTemp ! hihg speed outlet air dry bulb temperature [C]
REAL(r64) :: HSOutletAirEnthalpy    ! high speed outlet air enthalpy [J/kg]
REAL(r64) :: HSOutletAirHumRat      ! high speed outlet air humidity ratio [kg/kg]
REAL(r64) :: HSOutletAirRH          ! high speed outlet air relative humudity [fraction]
REAL(r64) :: hDelta              ! Change in air enthalpy across the cooling coil [J/kg]
REAL(r64) :: hTinwout            ! Enthalpy at inlet dry-bulb and outlet humidity ratio [J/kg]
REAL(r64) :: hADP                ! Apparatus dew point enthalpy [J/kg]
REAL(r64) :: tADP                ! Apparatus dew point temperature [C]
REAL(r64) :: wADP                ! Apparatus dew point humidity ratio [kg/kg]
REAL(r64) :: hTinwADP            ! Enthalpy at inlet dry-bulb and wADP [J/kg]
REAL(r64) :: RatedCBFHS          ! coil bypass factor at rated conditions (high speed)
REAL(r64) :: CBFHS               ! coil bypass factor at max flow (high speed)
REAL(r64) :: RatedCBFLS          ! coil bypass factor at rated conditions (low speed)
REAL(r64) :: CBFLS               ! coil bypass factor at max flow (low speed)
REAL(r64) :: TotCapHS            ! total capacity at high speed [W]
REAL(r64) :: SHRHS               ! sensible heat ratio at high speed
REAL(r64) :: TotCapLS            ! total capacity at low speed [W]
REAL(r64) :: SHRLS               ! sensible heat ratio at low speed
REAL(r64) :: EIRTempModFacHS     ! EIR modifier (function of entering wetbulb, outside drybulb) (high speed)
REAL(r64) :: EIRFlowModFacHS     ! EIR modifier (function of actual supply air flow vs rated flow) (high speed)
REAL(r64) :: EIRHS               ! EIR at off rated conditions (high speed)
REAL(r64) :: EIRTempModFacLS     ! EIR modifier (function of entering wetbulb, outside drybulb) (low speed)
REAL(r64) :: EIRFlowModFacLS     ! EIR modifier (function of actual supply air flow vs rated flow) (low speed)
REAL(r64) :: EIRLS               ! EIR at off rated conditions (low speed)
REAL(r64) :: SHR                 ! sensible heat ratio at current speed
REAL(r64) :: EIR                 ! EIR at current speed
REAL(r64) :: CBF                 ! CBFNom adjusted for actual air mass flow rate
REAL(r64) :: PLF                 ! Part load factor, accounts for thermal lag at compressor startup, used in
                               ! power calculation
REAL(r64) :: CondInletTemp       ! Condenser inlet temperature (C). Outdoor dry-bulb temp for air-cooled condenser.
                               ! Outdoor Wetbulb +(1 - effectiveness)*(outdoor drybulb - outdoor wetbulb) for evap condenser.
REAL(r64) :: CondInletHumrat     ! Condenser inlet humidity ratio (kg/kg). Zero for air-cooled condenser.
                               ! For evap condenser, its the humidity ratio of the air leaving the evap cooling pads.
REAL(r64) :: RhoAir              ! Density of air [kg/m3]
REAL(r64) :: RhoWater            ! Density of water [kg/m3]
REAL(r64) :: CondAirMassFlow     ! Condenser air mass flow rate [kg/s]
REAL(r64) :: EvapCondPumpElecPower ! Evaporative condenser electric pump power [W]
REAL(r64) :: MinAirHumRat = 0.0    ! minimum of the inlet air humidity ratio and the outlet air humidity ratio
INTEGER,SAVE :: DXMode=1        ! Performance mode for MultiMode DX coil; Always 1 for other coil types
REAL(r64) :: OutdoorDryBulb       ! Outdoor dry-bulb temperature at condenser (C)
REAL(r64) :: OutdoorWetBulb       ! Outdoor wet-bulb temperature at condenser (C)
REAL(r64) :: OutdoorHumRat        ! Outdoor humidity ratio at condenser (kg/kg)
REAL(r64) :: OutdoorPressure      ! Outdoor barometric pressure at condenser (Pa)
INTEGER :: SpeedNumHS           ! High speed number
INTEGER :: SpeedNumLS           ! Low speed number
REAL(r64) :: SHRUnadjusted        ! Temp SHR
REAL(r64)              :: QLatRated                ! Qlatent at rated conditions of indoor(TDB,TWB)=(26.7C,19.4C)
REAL(r64)              :: QLatActual               ! Qlatent at actual operating conditions
REAL(r64) :: AirMassFlowRatioLS   ! airflow ratio at low speed
REAL(r64) :: AirMassFlowRatioHS   ! airflow ratio at high speed
REAL(r64) :: WasteHeatLS          ! Waste heat at low speed
REAL(r64) :: WasteHeatHS          ! Waste heat at high speed
REAL(r64) :: LSElecCoolingPower   ! low speed power [W]
REAL(r64) :: HSElecCoolingPower   ! high speed power [W]
REAL(r64) :: CrankcaseHeatingPower ! Power due to crank case heater
REAL(r64) :: Hfg
REAL(r64) :: AirVolumeFlowRate    ! Air volume flow rate across the heating coil
REAL(r64) :: VolFlowperRatedTotCap ! Air volume flow rate divided by rated total heating capacity

IF (DXCoil(DXCoilNum)%CondenserInletNodeNum(DXMode) /= 0) THEN
  OutdoorDryBulb  = Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(DXMode))%Temp
  OutdoorHumRat   = Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(DXMode))%HumRat
  OutdoorPressure = Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(DXMode))%Press
  OutdoorWetBulb  = PsyTwbFnTdbWPb(OutdoorDryBulb,OutdoorHumRat,OutdoorPressure,RoutineName)
ELSE
  OutdoorDryBulb  = OutDryBulbTemp
  OutdoorHumRat   = OutHumRat
  OutdoorPressure = OutBaroPress
  OutdoorWetBulb  = OutWetBulbTemp
ENDIF

If (SpeedNum > 1) Then
   SpeedNumLS = SpeedNum-1
   SpeedNumHS = SpeedNum
  If (SpeedNum .GT. DXCoil(DXCoilNum)%NumOfSpeeds) Then
    SpeedNumLS = DXCoil(DXCoilNum)%NumOfSpeeds-1
    SpeedNumHS = DXCoil(DXCoilNum)%NumOfSpeeds
  End If
Else
  SpeedNumLS = 1
  SpeedNumHS = 1
End If

MSHPWasteHeat = 0.0
AirMassFlow = DXCoil(DXCoilNum)%InletAirMassFlowRate
AirMassFlowRatioLS = MSHPMassFlowRateLow/DXCoil(DXCoilNum)%MSRatedAirMassFlowRate(SpeedNumLS)
AirMassFlowRatioHS = MSHPMassFlowRateHigh/DXCoil(DXCoilNum)%MSRatedAirMassFlowRate(SpeedNumHS)

DXCoil(DXCoilNum)%PartLoadRatio              = 0.0
HeatReclaimDXCoil(DXCoilNum)%AvailCapacity   = 0.0
DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction = 0.0
InletAirDryBulbTemp = DXCoil(DXCoilNum)%InletAirTemp
InletAirEnthalpy = DXCoil(DXCoilNum)%InletAirEnthalpy
InletAirHumRat = DXCoil(DXCoilNum)%InletAirHumRat
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!InletAirPressure = DXCoil(DXCoilNum)%InletAirPressure
!InletAirWetbulbC = PsyTwbFnTdbWPb(InletAirDryBulbTemp,InletAirHumRat,InletAirPressure)
InletAirWetbulbC = PsyTwbFnTdbWPb(InletAirDryBulbTemp,InletAirHumRat,OutdoorPressure,RoutineName)
IF (DXCoil(DXCoilNum)%CondenserType(DXMode) == AirCooled) THEN
  CondInletTemp = OutdoorDryBulb ! Outdoor dry-bulb temp
ELSEIF (DXCoil(DXCoilNum)%CondenserType(DXMode) == EvapCooled) THEN
 ! Outdoor wet-bulb temp from DataEnvironment + (1.0-EvapCondEffectiveness) * (drybulb - wetbulb)
  CondInletTemp = OutdoorWetBulb + (OutdoorDryBulb-OutdoorWetBulb)*(1.0d0 - DXCoil(DXCoilNum)%MSEvapCondEffect(SpeedNumHS))
  CondInletHumrat = PsyWFnTdbTwbPb(CondInletTemp,OutdoorWetBulb,OutdoorPressure,RoutineName)
END IF
IF (OutdoorDryBulb .LT. DXCoil(DXCoilNum)%MaxOATCrankcaseHeater)THEN
  CrankcaseHeatingPower = DXCoil(DXCoilNum)%CrankcaseHeaterCapacity
ELSE
  CrankcaseHeatingPower = 0.0
END IF

IF((AirMassFlow .GT. 0.0) .AND. &
   (GetCurrentScheduleValue(DXCoil(DXCoilNum)%SchedPtr) .GT. 0.0) &
   .AND. (SpeedRatio > 0.0 .OR. CycRatio > 0.0) .AND. (CompOp == On)) THEN

  RhoAir = PsyRhoAirFnPbTdbW(OutdoorPressure,OutdoorDryBulb,OutdoorHumRat,RoutineName)
  IF (SpeedNum > 1) THEN

    ! Check for valid air volume flow per rated total cooling capacity (200 - 500 cfm/ton) at low speed
    AirVolumeFlowRate = MSHPMassFlowRateLow/PsyRhoAirFnPbTdbW(OutdoorPressure,InletAirDryBulbTemp, InletAirHumRat,RoutineName)
    !  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
    !  AirVolumeFlowRate = AirMassFlow/PsyRhoAirFnPbTdbW(InletAirPressure,InletAirDryBulbTemp, InletAirHumRat)
    VolFlowperRatedTotCap = AirVolumeFlowRate/DXCoil(DXCoilNum)%MSRatedTotCap(SpeedNumLS)
    IF ((VolFlowperRatedTotCap.LT.MinOperVolFlowPerRatedTotCap).OR.(VolFlowperRatedTotCap.GT.MaxCoolVolFlowPerRatedTotCap)) THEN
      IF (DXCoil(DXCoilNum)%MSErrIndex(SpeedNumLS) == 0) THEN
        CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//&
          '" - Air volume flow rate per watt of rated total cooling capacity is out of range at speed ' &
          //TRIM(TrimSigDigits(SpeedNumLS))//'.')
        CALL ShowContinueErrorTimeStamp(' ')
        CALL ShowContinueError('Expected range for VolumeFlowPerRatedTotalCapacity=['//  &
          TRIM(RoundSigDigits(MinOperVolFlowPerRatedTotCap,3))//'--'//  &
          TRIM(RoundSigDigits(MaxCoolVolFlowPerRatedTotCap,3))//']')
        CALL ShowContinueError('Possible causes include inconsistent air flow rates in system components or')
        CALL ShowContinueError('inconsistent supply air fan operation modes in coil and unitary system objects.')
      END IF
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//&
          '" - Air volume flow rate per watt of rated total cooling capacity is out ' //&
         'of range at speed '//TRIM(TrimSigDigits(SpeedNumLS))//'error continues...',   &
         DXCoil(DXCoilNum)%MSErrIndex(SpeedNumLS),ReportMinOf=VolFlowperRatedTotCap,ReportMaxOf=VolFlowperRatedTotCap)
    END IF

    ! Check for valid air volume flow per rated total cooling capacity (200 - 500 cfm/ton) at high speed
    AirVolumeFlowRate = MSHPMassFlowRateHigh/PsyRhoAirFnPbTdbW(OutdoorPressure,InletAirDryBulbTemp, InletAirHumRat,RoutineName)
    !  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
    !  AirVolumeFlowRate = AirMassFlow/PsyRhoAirFnPbTdbW(InletAirPressure,InletAirDryBulbTemp, InletAirHumRat)
    VolFlowperRatedTotCap = AirVolumeFlowRate/DXCoil(DXCoilNum)%MSRatedTotCap(SpeedNumHS)
    IF ((VolFlowperRatedTotCap.LT.MinOperVolFlowPerRatedTotCap).OR.(VolFlowperRatedTotCap.GT.MaxCoolVolFlowPerRatedTotCap)) THEN
      IF (DXCoil(DXCoilNum)%MSErrIndex(SpeedNumHS) == 0) THEN
        CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//&
          '" - Air volume flow rate per watt of rated total cooling capacity is out of range at speed ' &
          //Trim(TrimSigDigits(SpeedNumHS))//'.')
        CALL ShowContinueErrorTimeStamp(' ')
        CALL ShowContinueError('Expected range for VolumeFlowPerRatedTotalCapacity=['//  &
          TRIM(RoundSigDigits(MinOperVolFlowPerRatedTotCap,3))//'--'//  &
          TRIM(RoundSigDigits(MaxCoolVolFlowPerRatedTotCap,3))//']')
        CALL ShowContinueError('Possible causes include inconsistent air flow rates in system components or')
        CALL ShowContinueError('inconsistent supply air fan operation modes in coil and unitary system objects.')
      END IF
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//&
          '" - Air volume flow rate per watt of rated total cooling capacity is out ' //&
         'of range at speed '//TRIM(TrimSigDigits(SpeedNumHS))//'error continues...',   &
         DXCoil(DXCoilNum)%MSErrIndex(SpeedNumHS),ReportMinOf=VolFlowperRatedTotCap,ReportMaxOf=VolFlowperRatedTotCap)
    END IF

    ! Adjust high speed coil bypass factor for actual maximum air flow rate.
    RatedCBFHS = DXCoil(DXCoilNum)%MSRatedCBF(SpeedNumHS)
    CBFHS = AdjustCBF(RatedCBFHS,DXCoil(DXCoilNum)%MSRatedAirMassFlowRate(SpeedNumHS),MSHPMassFlowRateHigh)
    RatedCBFLS = DXCoil(DXCoilNum)%MSRatedCBF(SpeedNumLS)
    CBFLS = AdjustCBF(RatedCBFLS,DXCoil(DXCoilNum)%MSRatedAirMassFlowRate(SpeedNumLS),MSHPMassFlowRateLow)
    ! get low speed total capacity and SHR at current conditions
    CALL CalcTotCapSHR(InletAirDryBulbTemp,InletAirHumRat,InletAirEnthalpy,InletAirWetbulbC,AirMassFlowRatioLS, &
                       MSHPMassFlowRateLow,DXCoil(DXCoilNum)%MSRatedTotCap(SpeedNumLS), &
                       CBFLS,DXCoil(DXCoilNum)%MSCCapFTemp(SpeedNumLS), &
                       DXCoil(DXCoilNum)%MSCCapFFlow(SpeedNumLS),TotCapLS,SHRLS,CondInletTemp, &
                       OutdoorPressure)
    ! get low speed outlet conditions
    hDelta = TotCapLS / MSHPMassFlowRateLow
    ! Calculate new apparatus dew point conditions
    hADP = InletAirEnthalpy - hDelta/(1.d0-CBFLS)
    tADP = PsyTsatFnHPb(hADP,OutdoorPressure,'CalcMultiSpeedDXCoilCooling highspeed')
    !  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
    !  tADP = PsyTsatFnHPb(hADP,InletAirPressure)
    wADP = PsyWFnTdbH(tADP,hADP,'CalcMultiSpeedDXCoilCooling')
    hTinwADP = PsyHFnTdbW(InletAirDryBulbTemp,wADP,'CalcMultiSpeedDXCoilCooling highspeed')
    ! get corresponding SHR
    SHRLS = MIN((hTinwADP-hADP)/(InletAirEnthalpy-hADP),1.d0)
    ! get low speed outlet conditions
    LSOutletAirEnthalpy = InletAirEnthalpy - hDelta
    hTinwout = InletAirEnthalpy - (1.0d0-SHRLS)*hDelta
    LSOutletAirHumRat = PsyWFnTdbH(InletAirDryBulbTemp,hTinwout,RoutineName)
    LSOutletAirDryBulbTemp = PsyTdbFnHW(LSOutletAirEnthalpy,LSOutletAirHumRat,RoutineName)
    LSOutletAirRH = PsyRhFnTdbWPb(LSOutletAirDryBulbTemp,LSOutletAirHumRat,OutdoorPressure, &
                    'CalcMultiSpeedDXCoilCooling:highspeed')
    OutletAirDryBulbTempSat = PsyTsatFnHPb(LSOutletAirEnthalpy,OutdoorPressure,RoutineName)
    !  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
    !    LSOutletAirRH = PsyRhFnTdbWPb(LSOutletAirDryBulbTemp,LSOutletAirHumRat,InletAirPressure)
    !    OutletAirDryBulbTempSat = PsyTsatFnHPb(LSOutletAirEnthalpy,InletAirPressure)
    IF (LSOutletAirDryBulbTemp < OutletAirDryBulbTempSat) THEN ! Limit to saturated conditions at OutletAirEnthalpy
      LSOutletAirDryBulbTemp = OutletAirDryBulbTempSat
      LSOutletAirHumRat  = PsyWFnTdbH(LSOutletAirDryBulbTemp,LSOutletAirEnthalpy,RoutineName)
    END IF

    ! get high speed total capacity and SHR at current conditions

    CALL CalcTotCapSHR(InletAirDryBulbTemp,InletAirHumRat,InletAirEnthalpy,InletAirWetbulbC,AirMassFlowRatioHS, &
                       MSHPMassFlowRateHigh,DXCoil(DXCoilNum)%MSRatedTotCap(SpeedNumHS), &
                       CBFHS,DXCoil(DXCoilNum)%MSCCapFTemp(SpeedNumHS),DXCoil(DXCoilNum)%MSCCapFFlow(SpeedNumHS),TotCapHS,SHRHS, &
                       CondInletTemp, OutdoorPressure)
    hDelta = TotCapHS / MSHPMassFlowRateHigh
    ! Calculate new apparatus dew point conditions
    hADP = InletAirEnthalpy - hDelta/(1.d0-CBFHS)
    tADP = PsyTsatFnHPb(hADP,OutdoorPressure,RoutineName)
    !  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
    !  tADP = PsyTsatFnHPb(hADP,InletAirPressure)
    wADP = PsyWFnTdbH(tADP,hADP,RoutineName)
    hTinwADP = PsyHFnTdbW(InletAirDryBulbTemp,wADP,RoutineName)
    ! get corresponding SHR
    SHRHS = MIN((hTinwADP-hADP)/(InletAirEnthalpy-hADP),1.d0)
    ! get the part load factor that will account for cycling losses
    PLF = CurveValue(DXCoil(DXCoilNum)%MSPLFFPLR(SpeedNumHS),SpeedRatio)
    IF (PLF < 0.7d0) THEN
     PLF = 0.7d0
    END IF
    ! calculate the run time fraction
    DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction = SpeedRatio / PLF
    DXCoil(DXCoilNum)%PartLoadRatio    = SpeedRatio

    IF ( DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction > 1. ) THEN
      DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction = 1.0d0 ! Reset coil runtime fraction to 1.0
    END IF

    ! get high speed outlet conditions
    HSOutletAirEnthalpy = InletAirEnthalpy - hDelta
    hTinwout = InletAirEnthalpy - (1.0d0-SHRHS)*hDelta
    HSOutletAirHumRat = PsyWFnTdbH(InletAirDryBulbTemp,hTinwout,RoutineName)
    HSOutletAirDryBulbTemp = PsyTdbFnHW(HSOutletAirEnthalpy,HSOutletAirHumRat,RoutineName)
    HSOutletAirRH = PsyRhFnTdbWPb(HSOutletAirDryBulbTemp,HSOutletAirHumRat,OutdoorPressure, &
                    RoutineName//':highspeedoutlet')
    OutletAirDryBulbTempSat = PsyTsatFnHPb(HSOutletAirEnthalpy,OutdoorPressure,RoutineName)
    !  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
    !    LSOutletAirRH = PsyRhFnTdbWPb(LSOutletAirDryBulbTemp,LSOutletAirHumRat,InletAirPressure)
    !    OutletAirDryBulbTempSat = PsyTsatFnHPb(LSOutletAirEnthalpy,InletAirPressure)
    IF (HSOutletAirDryBulbTemp < OutletAirDryBulbTempSat) THEN ! Limit to saturated conditions at OutletAirEnthalpy
      HSOutletAirDryBulbTemp = OutletAirDryBulbTempSat
      HSOutletAirHumRat  = PsyWFnTdbH(HSOutletAirDryBulbTemp,HSOutletAirEnthalpy,RoutineName)
    END IF

    !  If constant fan with cycling compressor, call function to determine "effective SHR"
    !  which includes the part-load degradation on latent capacity
    IF (DXCoil(DXCoilNum)%LatentImpact .AND. FanOpMode .EQ. ContFanCycCoil .AND. SpeedRatio .GT. 0.0) THEN
      QLatRated = DXCoil(DXCoilNum)%MSRatedTotCap(SpeedNumHS) * (1.d0 - DXCoil(DXCoilNum)%MSRatedSHR(SpeedNumHS))
      QLatActual = TotCapHS * (1.d0 - SHRHS)
      SHRUnadjusted = SHRHS
      SHR = CalcEffectiveSHR(DXCoilNum, SHRHS, DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction, &
                          QLatRated, QLatActual, InletAirDryBulbTemp, InletAirWetbulbC, SpeedNumHS)
      ! Calculate full load output conditions
      If (SHR .GT. 1.0d0) SHR=1.0d0
      hTinwout = InletAirEnthalpy - (1.0d0-SHR)*hDelta
      If (SHR < 1.0d0) Then
        HSOutletAirHumRat = PsyWFnTdbH(InletAirDryBulbTemp,hTinwout,RoutineName)
      Else
        HSOutletAirHumRat = InletAirHumRat
      End If
      HSOutletAirDryBulbTemp = PsyTdbFnHW(HSOutletAirEnthalpy,HSOutletAirHumRat,RoutineName)
    END IF

    ! get high speed EIR at current conditions
    EIRTempModFacHS = CurveValue(DXCoil(DXCoilNum)%MSEIRFTemp(SpeedNumHS),InletAirWetbulbC,CondInletTemp)
    EIRFlowModFacHS = CurveValue(DXCoil(DXCoilNum)%MSEIRFFlow(SpeedNumHS),AirMassFlowRatioHS)
    EIRHS = 1.0d0/DXCoil(DXCoilNum)%MSRatedCOP(SpeedNumHS) * EIRFlowModFacHS * EIRTempModFacHS
    ! get low speed EIR at current conditions
    EIRTempModFacLS = CurveValue(DXCoil(DXCoilNum)%MSEIRFTemp(SpeedNumLS),InletAirWetbulbC,CondInletTemp)
    EIRFlowModFacLS = CurveValue(DXCoil(DXCoilNum)%MSEIRFFlow(SpeedNumLS),AirMassFlowRatioLS)
    EIRLS = 1.0d0/DXCoil(DXCoilNum)%MSRatedCOP(SpeedNumLS) * EIRTempModFacLS * EIRFlowModFacLS

    ! get current total capacity, SHR, EIR
    IF (SpeedRatio >= 1.0d0) THEN
      SHR = SHRHS
      EIR = EIRHS
      CondAirMassFlow =  RhoAir * DXCoil(DXCoilNum)%MSEvapCondAirFlow(SpeedNumHS)
      EvapCondPumpElecPower = DXCoil(DXCoilNum)%MSEvapCondPumpElecNomPower(SpeedNumHS)
    ELSE
      EIR = SpeedRatio*EIRHS + (1.0d0-SpeedRatio)*EIRLS
      SHR = SpeedRatio*SHRHS + (1.0d0-SpeedRatio)*SHRLS
      CondAirMassFlow =  RhoAir * (SpeedRatio * DXCoil(DXCoilNum)%MSEvapCondAirFlow(SpeedNumHS) + (1.0d0-SpeedRatio)* &
                       DXCoil(DXCoilNum)%MSEvapCondAirFlow(SpeedNumLS))
      EvapCondPumpElecPower = SpeedRatio * DXCoil(DXCoilNum)%MSEvapCondPumpElecNomPower(SpeedNumHS) + (1.0d0-SpeedRatio)* &
                              DXCoil(DXCoilNum)%MSEvapCondPumpElecNomPower(SpeedNumLS)
    END IF
    ! Outlet calculation
    DXCoil(DXCoilNum)%TotalCoolingEnergyRate = MSHPMassFlowRateHigh*(InletAirEnthalpy - HSOutletAirEnthalpy)*SpeedRatio + &
                                               MSHPMassFlowRateLow*(InletAirEnthalpy - LSOutletAirEnthalpy)*(1.0d0-SpeedRatio)
    ! Average outlet enthalpy
    OutletAirEnthalpy = InletAirEnthalpy - DXCoil(DXCoilNum)%TotalCoolingEnergyRate/DXCoil(DXCoilNum)%InletAirMassFlowRate
    MinAirHumRat = MIN(InletAirHumRat,SpeedRatio*HSOutletAirHumRat+(1.0d0-SpeedRatio)*LSOutletAirHumRat)
    DXCoil(DXCoilNum)%SensCoolingEnergyRate = MSHPMassFlowRateHigh*(PsyHFnTdbW(InletAirDryBulbTemp,MinAirHumRat,RoutineName) - &
                                              PsyHFnTdbW(HSOutletAirDryBulbTemp,MinAirHumRat,RoutineName))*SpeedRatio + &
                                              MSHPMassFlowRateLow*(PsyHFnTdbW(InletAirDryBulbTemp,MinAirHumRat,RoutineName) - &
                                              PsyHFnTdbW(LSOutletAirDryBulbTemp,MinAirHumRat,RoutineName))*(1.0d0-SpeedRatio)
    IF (DXCoil(DXCoilNum)%SensCoolingEnergyRate > DXCoil(DXCoilNum)%TotalCoolingEnergyRate) THEN
      DXCoil(DXCoilNum)%SensCoolingEnergyRate = DXCoil(DXCoilNum)%TotalCoolingEnergyRate
    END IF
    DXCoil(DXCoilNum)%LatCoolingEnergyRate = DXCoil(DXCoilNum)%TotalCoolingEnergyRate - &
                                             DXCoil(DXCoilNum)%SensCoolingEnergyRate

    IF (FanOpMode .EQ. CycFanCycCoil) OnOffFanPartLoadFraction = 1.0d0
    ! Update outlet conditions
    If (SpeedRatio .EQ. 0.0 .AND. FanOpMode .EQ. CycFanCycCoil) Then
      OutletAirEnthalpy = LSOutletAirEnthalpy
      OutletAirHumRat = LSOutletAirHumRat
      OutletAirDryBulbTemp = LSOutletAirDryBulbTemp
    Else If (SpeedRatio >= 1.0d0 .AND. FanOpMode .EQ. CycFanCycCoil) Then
      OutletAirEnthalpy = HSOutletAirEnthalpy
      OutletAirHumRat = HSOutletAirHumRat
      OutletAirDryBulbTemp = HSOutletAirDryBulbTemp
    Else
      Hfg = PsyHfgAirFnWTdb(MinAirHumRat,HSOutletAirDryBulbTemp*SpeedRatio+(1.0d0-SpeedRatio)*LSOutletAirDryBulbTemp, &
          RoutineName//':highspeed')
      ! Average outlet HR
      OutletAirHumRat = InletAirHumRat-DXCoil(DXCoilNum)%LatCoolingEnergyRate/Hfg/DXCoil(DXCoilNum)%InletAirMassFlowRate
      OutletAirDryBulbTemp = PsyTdbFnHW(OutletAirEnthalpy,OutletAirHumRat,RoutineName)
      IF (OutletAirDryBulbTemp < OutletAirDryBulbTempSat) THEN ! Limit to saturated conditions at OutletAirEnthalpy
        OutletAirDryBulbTemp = OutletAirDryBulbTempSat
        OutletAirHumRat  = PsyWFnTdbH(OutletAirDryBulbTemp,OutletAirEnthalpy,RoutineName)
        MinAirHumRat = MIN(InletAirHumRat,OutletAirHumRat)
        DXCoil(DXCoilNum)%SensCoolingEnergyRate = AirMassFlow*(PsyHFnTdbW(InletAirDryBulbTemp,MinAirHumRat,RoutineName) - &
                                                            PsyHFnTdbW(OutletAirDryBulbTemp,MinAirHumRat,RoutineName))
        IF (DXCoil(DXCoilNum)%SensCoolingEnergyRate > DXCoil(DXCoilNum)%TotalCoolingEnergyRate) THEN
          DXCoil(DXCoilNum)%SensCoolingEnergyRate = DXCoil(DXCoilNum)%TotalCoolingEnergyRate
        END IF
        DXCoil(DXCoilNum)%LatCoolingEnergyRate = DXCoil(DXCoilNum)%TotalCoolingEnergyRate - &
                                             DXCoil(DXCoilNum)%SensCoolingEnergyRate
      END IF
    End If

    LSElecCoolingPower = TotCapLS*EIRLS
    HSElecCoolingPower = TotCapHS*EIRHS

    ! Power calculation
    If (.NOT. DXCoil(DXCoilNum)%PLRImpact) Then
      DXCoil(DXCoilNum)%ElecCoolingPower = SpeedRatio*HSElecCoolingPower+(1.0d0-SpeedRatio)*LSElecCoolingPower
    Else
      DXCoil(DXCoilNum)%ElecCoolingPower = DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction*HSElecCoolingPower + &
                   (1.0d0-DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction)*LSElecCoolingPower
    End If
!   Calculation for heat reclaim needs to be corrected to use compressor power (not including condenser fan power)
    HeatReclaimDXCoil(DXCoilNum)%AvailCapacity = DXCoil(DXCoilNum)%TotalCoolingEnergyRate + DXCoil(DXCoilNum)%ElecCoolingPower

    ! Waste heat calculation
    WasteHeatLS = CurveValue(DXCoil(DXCoilNum)%MSWasteHeat(SpeedNumLS),OutdoorDryBulb,InletAirDryBulbTemp) * &
                  DXCoil(DXCoilNum)%MSWasteHeatFrac(SpeedNumLS)
    WasteHeatHS = CurveValue(DXCoil(DXCoilNum)%MSWasteHeat(SpeedNumHS),OutdoorDryBulb,InletAirDryBulbTemp) * &
                  DXCoil(DXCoilNum)%MSWasteHeatFrac(SpeedNumHS)
    MSHPWasteHeat = (SpeedRatio*WasteHeatHS + (1.0d0-SpeedRatio)*WasteHeatLS)*DXCoil(DXCoilNum)%ElecCoolingPower

    ! Energy use for other fuel types
    If (DXCoil(DXCoilNum)%FuelType .NE. FuelTypeElectricity) Then
      DXCoil(DXCoilNum)%FuelUsed = DXCoil(DXCoilNum)%ElecCoolingPower
      DXCoil(DXCoilNum)%ElecCoolingPower = 0.0
    End If

    DXCoil(DXCoilNum)%OutletAirEnthalpy = OutletAirEnthalpy
    DXCoil(DXCoilNum)%OutletAirHumRat = OutletAirHumRat
    DXCoil(DXCoilNum)%OutletAirTemp = OutletAirDryBulbTemp
    DXCoil(DXCoilNum)%CrankcaseHeaterPower = 0.0

  ELSE IF (CycRatio > 0.0) THEN

    IF (FanOpMode .EQ. CycFanCycCoil) AirMassFlow = AirMassFlow/CycRatio
    IF (FanOpMode .EQ. ContFanCycCoil) AirMassFlow = MSHPMassFlowRateLow

    ! Check for valid air volume flow per rated total cooling capacity (200 - 500 cfm/ton) at low speed
    AirVolumeFlowRate = MSHPMassFlowRateLow/PsyRhoAirFnPbTdbW(OutdoorPressure,InletAirDryBulbTemp, InletAirHumRat,RoutineName)
    !  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
    !  AirVolumeFlowRate = AirMassFlow/PsyRhoAirFnPbTdbW(InletAirPressure,InletAirDryBulbTemp, InletAirHumRat)
    VolFlowperRatedTotCap = AirVolumeFlowRate/DXCoil(DXCoilNum)%MSRatedTotCap(SpeedNumLS)
    IF ((VolFlowperRatedTotCap.LT.MinOperVolFlowPerRatedTotCap).OR.(VolFlowperRatedTotCap.GT.MaxCoolVolFlowPerRatedTotCap)) THEN
      IF (DXCoil(DXCoilNum)%MSErrIndex(SpeedNumLS) == 0) THEN
        CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//&
          '" - Air volume flow rate per watt of rated total cooling capacity is out of range at speed ' &
          //TRIM(TrimSigDigits(SpeedNumLS))//'.')
        CALL ShowContinueErrorTimeStamp(' ')
        CALL ShowContinueError('Expected range for VolumeFlowPerRatedTotalCapacity=['//  &
          TRIM(RoundSigDigits(MinOperVolFlowPerRatedTotCap,3))//'--'//  &
          TRIM(RoundSigDigits(MaxCoolVolFlowPerRatedTotCap,3))//']')
        CALL ShowContinueError('Possible causes include inconsistent air flow rates in system components or')
        CALL ShowContinueError('inconsistent supply air fan operation modes in coil and unitary system objects.')
      END IF
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//&
          '" - Air volume flow rate per watt of rated total cooling capacity is out ' //&
         'of range at speed '//TRIM(TrimSigDigits(SpeedNumHS))//'error continues...',   &
         DXCoil(DXCoilNum)%MSErrIndex(SpeedNumHS),ReportMinOf=VolFlowperRatedTotCap,ReportMaxOf=VolFlowperRatedTotCap)
    END IF

    IF (DXCoil(DXCoilNum)%CondenserType(SpeedNumLS) == EvapCooled) THEN
     ! Outdoor wet-bulb temp from DataEnvironment + (1.0-EvapCondEffectiveness) * (drybulb - wetbulb)
      CondInletTemp = OutdoorWetBulb + (OutdoorDryBulb-OutdoorWetBulb)*(1.0d0 - DXCoil(DXCoilNum)%MSEvapCondEffect(SpeedNumLS))
      CondInletHumrat = PsyWFnTdbTwbPb(CondInletTemp,OutdoorWetBulb,OutdoorPressure,RoutineName)
    END IF

    RatedCBFLS = DXCoil(DXCoilNum)%MSRatedCBF(SpeedNumLS)
    CBFLS = AdjustCBF(RatedCBFLS,DXCoil(DXCoilNum)%MSRatedAirMassFlowRate(SpeedNumLS),MSHPMassFlowRateLow)

    ! Adjust low speed coil bypass factor for actual flow rate.
    ! CBF = AdjustCBF(DXCoil(DXCoilNum)%RatedCBF2,DXCoil(DXCoilNum)%RatedAirMassFlowRate2,AirMassFlow)
    ! get low speed total capacity and SHR at current conditions
    CALL CalcTotCapSHR(InletAirDryBulbTemp,InletAirHumRat,InletAirEnthalpy,InletAirWetbulbC,AirMassFlowRatioLS, &
                       MSHPMassFlowRateLow,DXCoil(DXCoilNum)%MSRatedTotCap(SpeedNumLS), &
                       CBFLS,DXCoil(DXCoilNum)%MSCCapFTemp(SpeedNumLS), &
                       DXCoil(DXCoilNum)%MSCCapFFlow(SpeedNumLS),TotCapLS,SHRLS,CondInletTemp, &
                       OutdoorPressure)
    !  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
    !  Node(DXCoil(DXCoilNum)%AirInNode)%Press)
    hDelta = TotCapLS / AirMassFlow
    ! Adjust CBF for off-nominal flow
    CBF = AdjustCBF(DXCoil(DXCoilNum)%MSRatedCBF(SpeedNumLS),DXCoil(DXCoilNum)%MSRatedAirMassFlowRate(SpeedNumLS),AirMassFlow)
    ! Calculate new apparatus dew point conditions
    hADP = InletAirEnthalpy - hDelta/(1.d0-CBF)
    tADP = PsyTsatFnHPb(hADP,OutdoorPressure,RoutineName)
    !  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
    !  tADP = PsyTsatFnHPb(hADP,InletAirPressure)
    wADP = PsyWFnTdbH(tADP,hADP,RoutineName)
    hTinwADP = PsyHFnTdbW(InletAirDryBulbTemp,wADP,RoutineName)
    ! get corresponding SHR
    SHR = MIN((hTinwADP-hADP)/(InletAirEnthalpy-hADP),1.d0)

    ! get the part load factor that will account for cycling losses
    PLF = CurveValue(DXCoil(DXCoilNum)%MSPLFFPLR(SpeedNumLS),CycRatio)
    IF (FanOpMode .EQ. CycFanCycCoil .AND. CycRatio .EQ. 1.0d0 .AND. PLF .NE. 1.0d0) Then
      IF (DXCoil(DXCoilNum)%PLFErrIndex == 0) THEN
        CALL ShowWarningMessage('The PLF curve value for DX cooling coil '//TRIM(DXCoil(DXCoilNum)%Name)//&
                         ' ='//TRIM(RoundSigDigits(PLF,2))//' for part-load ratio = 1')
        CALL ShowContinueError('PLF curve value must be = 1.0 and has been reset to 1.0. Simulation is continuing.')
        CALL ShowContinueErrorTimeStamp(' ')
      ENDIF
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%Name)//'":'//&
          ' DX cooling coil PLF curve value <> 1.0 warning continues...' &
          , DXCoil(DXCoilNum)%PLFErrIndex, PLF, PLF)
      PLF = 1.0d0
    END IF

    IF (PLF < 0.7d0) THEN
     PLF = 0.7d0
    END IF
    ! calculate the run time fraction
    DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction = CycRatio / PLF
    DXCoil(DXCoilNum)%PartLoadRatio    = CycRatio

    IF ( DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction > 1.d0 ) THEN
      DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction = 1.0d0 ! Reset coil runtime fraction to 1.0
    END IF

    ! get low speed outlet conditions
    LSOutletAirEnthalpy = InletAirEnthalpy - hDelta
    hTinwout = InletAirEnthalpy - (1.0d0-SHR)*hDelta
    LSOutletAirHumRat = PsyWFnTdbH(InletAirDryBulbTemp,hTinwout,RoutineName)
    LSOutletAirDryBulbTemp = PsyTdbFnHW(LSOutletAirEnthalpy,LSOutletAirHumRat,RoutineName)
    LSOutletAirRH = PsyRhFnTdbWPb(LSOutletAirDryBulbTemp,LSOutletAirHumRat,OutdoorPressure, &
                    RoutineName//':lowspeedoutlet')
    OutletAirDryBulbTempSat = PsyTsatFnHPb(LSOutletAirEnthalpy,OutdoorPressure,RoutineName)
    !  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
    !    LSOutletAirRH = PsyRhFnTdbWPb(LSOutletAirDryBulbTemp,LSOutletAirHumRat,InletAirPressure)
    !    OutletAirDryBulbTempSat = PsyTsatFnHPb(LSOutletAirEnthalpy,InletAirPressure)
    IF (LSOutletAirDryBulbTemp < OutletAirDryBulbTempSat) THEN ! Limit to saturated conditions at OutletAirEnthalpy
      LSOutletAirDryBulbTemp = OutletAirDryBulbTempSat
      LSOutletAirHumRat  = PsyWFnTdbH(LSOutletAirDryBulbTemp,LSOutletAirEnthalpy,RoutineName)
    END IF

    !  If constant fan with cycling compressor, call function to determine "effective SHR"
    !  which includes the part-load degradation on latent capacity
    IF (FanOpMode .EQ. ContFanCycCoil) THEN
      QLatRated = DXCoil(DXCoilNum)%MSRatedTotCap(SpeedNumLS) * (1.d0 - DXCoil(DXCoilNum)%MSRatedSHR(SpeedNumLS))
      QLatActual = TotCapLS * (1.d0 - SHR)
      SHRUnadjusted = SHR
      SHR = CalcEffectiveSHR(DXCoilNum, SHR, DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction, &
                          QLatRated, QLatActual, InletAirDryBulbTemp, InletAirWetbulbC, SpeedNumLS)
      ! Calculate full load output conditions
      If (SHR .GT. 1.0d0) SHR=1.0d0
      hTinwout = InletAirEnthalpy - (1.0d0-SHR)*hDelta
      If (SHR < 1.0d0) Then
        LSOutletAirHumRat = PsyWFnTdbH(InletAirDryBulbTemp,hTinwout,RoutineName)
      Else
        LSOutletAirHumRat = InletAirHumRat
      End If
      LSOutletAirDryBulbTemp = PsyTdbFnHW(LSOutletAirEnthalpy,LSOutletAirHumRat,RoutineName)
    END IF

    IF (FanOpMode .EQ. CycFanCycCoil) OnOffFanPartLoadFraction = PLF
    ! outlet conditions are average of inlet and low speed weighted by CycRatio
    OutletAirEnthalpy = LSOutletAirEnthalpy
    OutletAirHumRat = LSOutletAirHumRat
    OutletAirDryBulbTemp = LSOutletAirDryBulbTemp
    ! get low speed EIR at current conditions
    EIRTempModFacLS = CurveValue(DXCoil(DXCoilNum)%MSEIRFTemp(SpeedNumLS),InletAirWetbulbC,CondInletTemp)
    EIRFlowModFacLS = CurveValue(DXCoil(DXCoilNum)%MSEIRFFlow(SpeedNumLS),AirMassFlowRatioLS)
    EIRLS = 1.0d0/DXCoil(DXCoilNum)%MSRatedCOP(SpeedNumLS) * EIRTempModFacLS * EIRFlowModFacLS

    ! get the eletrical power consumption
    DXCoil(DXCoilNum)%ElecCoolingPower = TotCapLS * EIRLS * DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction
    ! calculate cooling output power
!    AirMassFlow = DXCoil(DXCoilNum)%InletAirMassFlowRate
    DXCoil(DXCoilNum)%TotalCoolingEnergyRate = AirMassFlow * (InletAirEnthalpy - LSOutletAirEnthalpy)*CycRatio
    IF (FanOpMode .EQ. ContFanCycCoil) THEN
      OutletAirEnthalpy = InletAirEnthalpy - DXCoil(DXCoilNum)%TotalCoolingEnergyRate/DXCoil(DXCoilNum)%InletAirMassFlowRate
      MinAirHumRat = MIN(InletAirHumRat,LSOutletAirHumRat)
      DXCoil(DXCoilNum)%SensCoolingEnergyRate = AirMassFlow * (PsyHFnTdbW(InletAirDryBulbTemp,MinAirHumRat,RoutineName) - &
                        PsyHFnTdbW(LSOutletAirDryBulbTemp,MinAirHumRat,RoutineName))*CycRatio
      IF (DXCoil(DXCoilNum)%SensCoolingEnergyRate > DXCoil(DXCoilNum)%TotalCoolingEnergyRate) THEN
        DXCoil(DXCoilNum)%SensCoolingEnergyRate = DXCoil(DXCoilNum)%TotalCoolingEnergyRate
      END IF
      DXCoil(DXCoilNum)%LatCoolingEnergyRate = DXCoil(DXCoilNum)%TotalCoolingEnergyRate - &
                                             DXCoil(DXCoilNum)%SensCoolingEnergyRate
      ! Calculate avarage outlet conditions
      Hfg = PsyHfgAirFnWTdb(MinAirHumRat,OutletAirDryBulbTemp*CycRatio+(1.0d0-CycRatio)*InletAirDryBulbTemp,'MultiSpeedCooling ')
      OutletAirHumRat = InletAirHumRat-DXCoil(DXCoilNum)%LatCoolingEnergyRate/Hfg/DXCoil(DXCoilNum)%InletAirMassFlowRate
      OutletAirDryBulbTemp = PsyTdbFnHW(OutletAirEnthalpy,OutletAirHumRat,RoutineName)
      IF (OutletAirDryBulbTemp < OutletAirDryBulbTempSat) THEN ! Limit to saturated conditions at OutletAirEnthalpy
        OutletAirDryBulbTemp = OutletAirDryBulbTempSat
        OutletAirHumRat  = PsyWFnTdbH(OutletAirDryBulbTemp,OutletAirEnthalpy,RoutineName)
        MinAirHumRat = MIN(InletAirHumRat,OutletAirHumRat)
        DXCoil(DXCoilNum)%SensCoolingEnergyRate = DXCoil(DXCoilNum)%InletAirMassFlowRate* &
          (PsyHFnTdbW(InletAirDryBulbTemp,MinAirHumRat,RoutineName) - PsyHFnTdbW(OutletAirDryBulbTemp,MinAirHumRat,RoutineName))
        IF (DXCoil(DXCoilNum)%SensCoolingEnergyRate > DXCoil(DXCoilNum)%TotalCoolingEnergyRate) THEN
          DXCoil(DXCoilNum)%SensCoolingEnergyRate = DXCoil(DXCoilNum)%TotalCoolingEnergyRate
        END IF
        DXCoil(DXCoilNum)%LatCoolingEnergyRate = DXCoil(DXCoilNum)%TotalCoolingEnergyRate - &
                                             DXCoil(DXCoilNum)%SensCoolingEnergyRate
      END IF
    Else
      MinAirHumRat = MIN(InletAirHumRat,OutletAirHumRat)
      DXCoil(DXCoilNum)%SensCoolingEnergyRate = AirMassFlow * (PsyHFnTdbW(InletAirDryBulbTemp,MinAirHumRat,RoutineName) - &
                                                             PsyHFnTdbW(LSOutletAirDryBulbTemp,MinAirHumRat,RoutineName))*CycRatio
      ! Don't let sensible capacity be greater than total capacity
      IF (DXCoil(DXCoilNum)%SensCoolingEnergyRate > DXCoil(DXCoilNum)%TotalCoolingEnergyRate) THEN
         DXCoil(DXCoilNum)%SensCoolingEnergyRate = DXCoil(DXCoilNum)%TotalCoolingEnergyRate
      END IF
      DXCoil(DXCoilNum)%LatCoolingEnergyRate = DXCoil(DXCoilNum)%TotalCoolingEnergyRate - &
                                             DXCoil(DXCoilNum)%SensCoolingEnergyRate
    End If
    !   Calculation for heat reclaim needs to be corrected to use compressor power (not including condenser fan power)
    HeatReclaimDXCoil(DXCoilNum)%AvailCapacity = DXCoil(DXCoilNum)%TotalCoolingEnergyRate + DXCoil(DXCoilNum)%ElecCoolingPower
    DXCoil(DXCoilNum)%OutletAirEnthalpy = OutletAirEnthalpy
    DXCoil(DXCoilNum)%OutletAirHumRat = OutletAirHumRat
    DXCoil(DXCoilNum)%OutletAirTemp = OutletAirDryBulbTemp
    CondAirMassFlow =  RhoAir * DXCoil(DXCoilNum)%MSEvapCondAirFlow(SpeedNumLS) * DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction
    EvapCondPumpElecPower = DXCoil(DXCoilNum)%MSEvapCondPumpElecNomPower(SpeedNumLS) * DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction

    ! Waste heat
    MSHPWasteHeat = CurveValue(DXCoil(DXCoilNum)%MSWasteHeat(SpeedNumLS),OutdoorDryBulb,InletAirDryBulbTemp) * &
                  DXCoil(DXCoilNum)%MSWasteHeatFrac(SpeedNumLS)*DXCoil(DXCoilNum)%ElecCoolingPower

    ! Energy use for other fuel types
    If (DXCoil(DXCoilNum)%FuelType .NE. FuelTypeElectricity) Then
      DXCoil(DXCoilNum)%FuelUsed = DXCoil(DXCoilNum)%ElecCoolingPower
      DXCoil(DXCoilNum)%ElecCoolingPower = 0.0
    End If

    IF(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil .EQ. 0) THEN
      DXCoil(DXCoilNum)%CrankcaseHeaterPower = CrankcaseHeatingPower * &
                                           (1.0d0 - DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction)
    ELSE
      DXCoil(DXCoilNum)%CrankcaseHeaterPower = CrankcaseHeatingPower * &
                                           (1.0d0 - MAX(DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction, &
                                           DXCoil(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil)%HeatingCoilRuntimeFraction))
    END IF

  END IF

  IF (DXCoil(DXCoilNum)%CondenserType(DXMode) == EvapCooled) THEN
  !******************
  !             WATER CONSUMPTION IN m3 OF WATER FOR DIRECT
  !             H2O [m3/sec] = Delta W[KgH2O/Kg air]*Mass Flow Air[Kg air]
  !                                /RhoWater [kg H2O/m3 H2O]
  !******************
     RhoWater = RhoH2O(OutdoorDryBulb)
     DXCoil(DXCoilNum)%EvapWaterConsumpRate =  (CondInletHumrat - OutdoorHumRat) *  CondAirMassFlow/RhoWater
     DXCoil(DXCoilNum)%EvapCondPumpElecPower = EvapCondPumpElecPower
     !set water system demand request (if needed)
     IF ( DXCoil(DxCoilNum)%EvapWaterSupplyMode == WaterSupplyFromTank) THEN
       WaterStorage(DXCoil(DXCoilNum)%EvapWaterSupTankID)%VdotRequestDemand(DXCoil(DXCoilNum)%EvapWaterTankDemandARRID) &
       = DXCoil(DXCoilNum)%EvapWaterConsumpRate
     ENDIF

     ! Calculate basin heater power
     CALL CalcBasinHeaterPower(DXCoil(DXCoilNum)%BasinHeaterPowerFTempDiff,&
                              DXCoil(DXCoilNum)%BasinHeaterSchedulePtr,&
                              DXCoil(DXCoilNum)%BasinHeaterSetPointTemp,DXCoil(DXCoilNum)%BasinHeaterPower)
     DXCoil(DXCoilNum)%BasinHeaterPower = DXCoil(DXCoilNum)%BasinHeaterPower * &
                                      (1.d0 - DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction)
  ENDIF

ELSE

  ! DX coil is off; just pass through conditions
  DXCoil(DXCoilNum)%OutletAirEnthalpy = DXCoil(DXCoilNum)%InletAirEnthalpy
  DXCoil(DXCoilNum)%OutletAirHumRat = DXCoil(DXCoilNum)%InletAirHumRat
  DXCoil(DXCoilNum)%OutletAirTemp = DXCoil(DXCoilNum)%InletAirTemp

  DXCoil(DXCoilNum)%FuelUsed = 0.0
  DXCoil(DXCoilNum)%ElecCoolingPower = 0.0
  DXCoil(DXCoilNum)%TotalCoolingEnergyRate = 0.0
  DXCoil(DXCoilNum)%SensCoolingEnergyRate = 0.0
  DXCoil(DXCoilNum)%LatCoolingEnergyRate = 0.0
  DXCoil(DXCoilNum)%EvapCondPumpElecPower = 0.0
  DXCoil(DXCoilNum)%EvapWaterConsumpRate = 0.0

  IF(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil .EQ. 0) THEN
    DXCoil(DXCoilNum)%CrankcaseHeaterPower = CrankcaseHeatingPower
  ELSE
    DXCoil(DXCoilNum)%CrankcaseHeaterPower = CrankcaseHeatingPower * &
                                         (1.0d0 - DXCoil(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil)%HeatingCoilRuntimeFraction)
  END IF

  ! Calculate basin heater power
  IF (DXCoil(DXCoilNum)%CondenserType(DXMode) == EvapCooled) THEN
    CALL CalcBasinHeaterPower(DXCoil(DXCoilNum)%BasinHeaterPowerFTempDiff,&
                              DXCoil(DXCoilNum)%BasinHeaterSchedulePtr,&
                              DXCoil(DXCoilNum)%BasinHeaterSetPointTemp,DXCoil(DXCoilNum)%BasinHeaterPower)
  ENDIF
END IF

DXCoilOutletTemp(DXCoilNum)     = DXCoil(DXCoilNum)%OutletAirTemp
DXCoilOutletHumRat(DXCoilNum)   = DXCoil(DXCoilNum)%OutletAirHumRat
DXCoilPartLoadRatio(DXCoilNum)  = DXCoil(DXCoilNum)%PartLoadRatio
DXCoilFanOpMode(DXCoilNum)      = FanOpMode
DXCoil(DXCoilNum)%CondInletTemp = CondInletTemp ! Save condenser inlet temp in the data structure

RETURN

END SUBROUTINE CalcMultiSpeedDXCoilCooling

SUBROUTINE CalcMultiSpeedDXCoilHeating(DXCoilNum,SpeedRatio, CycRatio, SpeedNum, FanOpMode)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Lixing Gu, FSEC
          !       DATE WRITTEN   June 2007
          !       MODIFIED       na
          !       RE-ENGINEERED  Revised based on CalcDXHeatingCoil

          ! PURPOSE OF THIS SUBROUTINE:
          ! Calculates the air-side performance and electrical energy use of a direct-
          ! expansion, air-cooled cooling unit with a multispeed compressor.

          ! METHODOLOGY EMPLOYED:
          ! Uses the same methodology as the single speed DX heating unit model (SUBROUTINE CalcDXHeatingCoil).
          ! In addition it assumes that the unit performance is obtained by interpolating between
          ! the performance at high speed and that at low speed. If the output needed is below
          ! that produced at low speed, the compressor cycles between off and low speed.

          ! USE STATEMENTS:
USE CurveManager, ONLY: CurveValue
USE General,      ONLY: TrimSigDigits, RoundSigDigits
USE DataWater,    ONLY: WaterStorage
USE DataHVACGlobals,     ONLY: MSHPMassFlowRateLow, MSHPMassFlowRateHigh, MSHPWasteHeat

IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
INTEGER   :: DXCoilNum      ! the number of the DX heating coil to be simulated
REAL(r64) :: SpeedRatio     ! = (CompressorSpeed - CompressorSpeedMin) / (CompressorSpeedMax - CompressorSpeedMin)
                          ! SpeedRatio varies between 1.0 (maximum speed) and 0.0 (minimum speed)
REAL(r64) :: CycRatio       ! cycling part load ratio
INTEGER   :: SpeedNum       ! Speed number
INTEGER   :: FanOpMode      ! Fan operation mode

          ! SUBROUTINE PARAMETER DEFINITIONS:
CHARACTER(len=*), PARAMETER :: RoutineName='CalcMultiSpeedDXCoilHeating'

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
REAL(r64) :: AirMassFlow         ! dry air mass flow rate through coil [kg/s]
REAL(r64) :: InletAirWetBulbC    ! wetbulb temperature of inlet air [C]
REAL(r64) :: InletAirDryBulbTemp ! inlet air dry bulb temperature [C]
REAL(r64) :: InletAirEnthalpy    ! inlet air enthalpy [J/kg]
REAL(r64) :: InletAirHumRat      ! inlet air humidity ratio [kg/kg]
REAL(r64) :: OutletAirEnthalpy   ! outlet air enthalpy [J/kg]
REAL(r64) :: OutletAirHumRat     ! outlet air humidity ratio [kg/kg]
REAL(r64) :: TotCapHS            ! total capacity at high speed [W]
REAL(r64) :: TotCapLS            ! total capacity at low speed [W]
REAL(r64) :: EIRHS               ! EIR at off rated conditions (high speed)
REAL(r64) :: EIRLS               ! EIR at off rated conditions (low speed)
REAL(r64) :: TotCap              ! total capacity at current speed [W]
REAL(r64) :: EIR                 ! EIR at current speed
REAL(r64) :: PLF                 ! Part load factor, accounts for thermal lag at compressor startup, used in
                               ! power calculation
REAL(r64) :: OutdoorDryBulb       ! Outdoor dry-bulb temperature at condenser (C)
REAL(r64) :: OutdoorHumRat        ! Outdoor humidity ratio at condenser (kg/kg)
REAL(r64) :: OutdoorPressure      ! Outdoor barometric pressure at condenser (Pa)
INTEGER :: SpeedNumHS           ! High speed number
INTEGER :: SpeedNumLS           ! Low speed number
REAL(r64) :: AirMassFlowRatioLS   ! airflow ratio at low speed
REAL(r64) :: AirMassFlowRatioHS   ! airflow ratio at high speed
REAL(r64) :: AirFlowRatio         ! Airflow ratio
REAL(r64) :: PLRHeating           ! Part load ratio in heating
REAL(r64) :: CrankcaseHeatingPower ! Power due to crank case heater
REAL(r64) :: AirVolumeFlowRate    ! Air volume flow rate across the heating coil
REAL(r64) :: VolFlowperRatedTotCap ! Air volume flow rate divided by rated total heating capacity
REAL(r64) :: TotCapTempModFac     ! Total capacity modifier as a function ot temperature
REAL(r64) :: TotCapFlowModFac     ! Total capacity modifier as a function of flow ratio
REAL(r64) :: OutdoorCoilT         ! Outdoor coil temperature
REAL(r64) :: OutdoorCoildw        ! Outdoor coil delta w assuming coil temperature of OutdoorCoilT
REAL(r64) :: LoadDueToDefrost     ! Additonal load due to defrost
REAL(r64) :: LoadDueToDefrostLS   ! Additonal load due to defrost at low speed
REAL(r64) :: LoadDueToDefrostHS   ! Additonal load due to defrost at high speed
REAL(r64) :: HeatingCapacityMultiplier ! Multiplier for heating capacity when system is in defrost
REAL(r64) :: FractionalDefrostTime ! Fraction of time step when system is in defrost
REAL(r64) :: InputPowerMultiplier  ! Multiplier for poer when system is in defrost
REAL(r64) :: DefrostEIRTempModFac ! EIR modifier for defrost
REAL(r64) :: FullLoadOutAirEnth   ! Outlet full load enthalpy
REAL(r64) :: FullLoadOutAirHumRat ! Outlet humidity ratio at full load
REAL(r64) :: FullLoadOutAirTemp   ! Outlet temperature at full load
REAL(r64) :: FullLoadOutAirRH     ! Outler relative humidity at full load
REAL(r64) :: OutletAirTemp        ! Supply ari temperature
REAL(r64) :: EIRTempModFac        ! EIR modifier as a function of temperature
REAL(r64) :: EIRFlowModFac        ! EIR modifier as a function of airflow ratio
REAL(r64) :: WasteHeatLS          ! Waste heat at low speed
REAL(r64) :: WasteHeatHS          ! Waste heat at high speed
REAL(r64) :: LSFullLoadOutAirEnth ! Outlet full load enthalpy at low speed
REAL(r64) :: HSFullLoadOutAirEnth ! Outlet full load enthalpy at high speed
REAL(r64) :: LSElecHeatingPower   ! Full load power at low speed
REAL(r64) :: HSElecHeatingPower   ! Full load power at high speed
REAL(r64) :: DefrostPowerLS       ! Defrost power at low speed [W]
REAL(r64) :: DefrostPowerHS       ! Defrost power at high speed [W]

! FLOW
If (SpeedNum > 1) Then
   SpeedNumLS = SpeedNum-1
   SpeedNumHS = SpeedNum
  If (SpeedNum .GT. DXCoil(DXCoilNum)%NumOfSpeeds) Then
    SpeedNumLS = DXCoil(DXCoilNum)%NumOfSpeeds-1
    SpeedNumHS = DXCoil(DXCoilNum)%NumOfSpeeds
  End If
Else
  SpeedNumLS = 1
  SpeedNumHS = 1
End If

AirMassFlow = DXCoil(DXCoilNum)%InletAirMassFlowRate
AirMassFlowRatioLS = MSHPMassFlowRateLow/DXCoil(DXCoilNum)%MSRatedAirMassFlowRate(SpeedNumLS)
AirMassFlowRatioHS = MSHPMassFlowRateHigh/DXCoil(DXCoilNum)%MSRatedAirMassFlowRate(SpeedNumHS)

AirFlowRatio = 1.0d0
IF(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil .EQ. 0) MSHPWasteHeat = 0.0

! Get condenser outdoor node info from DX Heating Coil
IF (DXCoil(DXCoilNum)%CondenserInletNodeNum(1) /= 0) THEN
  OutdoorDryBulb  = Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(1))%Temp
  OutdoorHumRat   = Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(1))%HumRat
  OutdoorPressure = Node(DXCoil(DXCoilNum)%CondenserInletNodeNum(1))%Press
ELSE
  OutdoorDryBulb  = OutDryBulbTemp
  OutdoorHumRat   = OutHumRat
  OutdoorPressure = OutBaroPress
ENDIF

InletAirDryBulbTemp = DXCoil(DXCoilNum)%InletAirTemp
InletAirEnthalpy = DXCoil(DXCoilNum)%InletAirEnthalpy
InletAirHumRat = DXCoil(DXCoilNum)%InletAirHumRat
!  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
!InletAirPressure = DXCoil(DXCoilNum)%InletAirPressure
!InletAirWetbulbC = PsyTwbFnTdbWPb(InletAirDryBulbTemp,InletAirHumRat,InletAirPressure)
InletAirWetbulbC = PsyTwbFnTdbWPb(InletAirDryBulbTemp,InletAirHumRat,OutdoorPressure,RoutineName)
PLRHeating = 0.0
DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction = 0.0
! Initialize crankcase heater, operates below OAT defined in input deck for HP DX heating coil
IF (OutdoorDryBulb .LT. DXCoil(DXCoilNum)%MaxOATCrankcaseHeater)THEN
  CrankcaseHeatingPower = DXCoil(DXCoilNum)%CrankcaseHeaterCapacity
ELSE
  CrankcaseHeatingPower = 0.0
END IF
DXCoil(DXCoilNum)%PartLoadRatio              = 0.0
HeatReclaimDXCoil(DXCoilNum)%AvailCapacity   = 0.0

IF((AirMassFlow .GT. 0.0) .AND. &
   (GetCurrentScheduleValue(DXCoil(DXCoilNum)%SchedPtr) .GT. 0.0) .AND. &
   ((CycRatio .GT. 0.0) .OR. (SpeedRatio .GT. 0.0)) .AND. OutdoorDryBulb .GT. DXCoil(DXCoilNum)%MinOATCompressor) THEN

  If (SpeedNum > 1) Then

    ! Check for valid air volume flow per rated total cooling capacity (200 - 600 cfm/ton) at low speed
    AirVolumeFlowRate = MSHPMassFlowRateLow/PsyRhoAirFnPbTdbW(OutdoorPressure,InletAirDryBulbTemp, InletAirHumRat,RoutineName)
    !  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
    !  AirVolumeFlowRate = AirMassFlow/PsyRhoAirFnPbTdbW(InletAirPressure,InletAirDryBulbTemp, InletAirHumRat)
    VolFlowperRatedTotCap = AirVolumeFlowRate/DXCoil(DXCoilNum)%MSRatedTotCap(SpeedNumLS)
    IF ((VolFlowperRatedTotCap.LT.MinOperVolFlowPerRatedTotCap).OR.(VolFlowperRatedTotCap.GT.MaxHeatVolFlowPerRatedTotCap)) THEN
      IF (DXCoil(DXCoilNum)%MSErrIndex(SpeedNumLS) == 0) THEN
        CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//&
          '" - Air volume flow rate per watt of rated total heating capacity is out of range at speed ' &
          //TRIM(TrimSigDigits(SpeedNumLS))//'.')
        CALL ShowContinueErrorTimeStamp(' ')
        CALL ShowContinueError('Expected range for VolumeFlowPerRatedTotalCapacity=['//  &
          TRIM(RoundSigDigits(MinOperVolFlowPerRatedTotCap,3))//'--'//  &
          TRIM(RoundSigDigits(MaxHeatVolFlowPerRatedTotCap,3))//']')
        CALL ShowContinueError('Possible causes include inconsistent air flow rates in system components or')
        CALL ShowContinueError('inconsistent supply air fan operation modes in coil and unitary system objects.')
      END IF
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//&
          '" - Air volume flow rate per watt of rated total heating capacity is out ' //&
         'of range at speed '//TRIM(TrimSigDigits(SpeedNumLS))//'error continues...',   &
         DXCoil(DXCoilNum)%MSErrIndex(SpeedNumLS),ReportMinOf=VolFlowperRatedTotCap,ReportMaxOf=VolFlowperRatedTotCap)
    END IF

    ! Check for valid air volume flow per rated total cooling capacity (200 - 600 cfm/ton) at high speed
    AirVolumeFlowRate = MSHPMassFlowRateHigh/PsyRhoAirFnPbTdbW(OutdoorPressure,InletAirDryBulbTemp, InletAirHumRat,RoutineName)
    !  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
    !  AirVolumeFlowRate = AirMassFlow/PsyRhoAirFnPbTdbW(InletAirPressure,InletAirDryBulbTemp, InletAirHumRat)
    VolFlowperRatedTotCap = AirVolumeFlowRate/DXCoil(DXCoilNum)%MSRatedTotCap(SpeedNumHS)
    IF ((VolFlowperRatedTotCap.LT.MinOperVolFlowPerRatedTotCap).OR.(VolFlowperRatedTotCap.GT.MaxHeatVolFlowPerRatedTotCap)) THEN
      IF (DXCoil(DXCoilNum)%MSErrIndex(SpeedNumHS) == 0) THEN
        CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//&
          '" - Air volume flow rate per watt of rated total heating capacity is out of range at speed ' &
          //TRIM(TrimSigDigits(SpeedNumHS))//'.')
        CALL ShowContinueErrorTimeStamp(' ')
        CALL ShowContinueError('Expected range for VolumeFlowPerRatedTotalCapacity=['//  &
          TRIM(RoundSigDigits(MinOperVolFlowPerRatedTotCap,3))//'--'//  &
          TRIM(RoundSigDigits(MaxHeatVolFlowPerRatedTotCap,3))//']')
        CALL ShowContinueError('Possible causes include inconsistent air flow rates in system components or')
        CALL ShowContinueError('inconsistent supply air fan operation modes in coil and unitary system objects.')
      END IF
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//&
          '" - Air volume flow rate per watt of rated total heating capacity is out ' //&
         'of range at speed '//TRIM(TrimSigDigits(SpeedNumHS))//'error continues...',   &
         DXCoil(DXCoilNum)%MSErrIndex(SpeedNumHS),ReportMinOf=VolFlowperRatedTotCap,ReportMaxOf=VolFlowperRatedTotCap)
    END IF

    ! Get total capacity modifying factor (function of temperature) for off-rated conditions
    ! Model was extended to accept bi-quadratic curves. This allows sensitivity of the heating capacity
    ! to the entering dry-bulb temperature as well as the outside dry-bulb temperature. User is
    ! advised to use the bi-quaratic curve if sufficient manufacturer data is available.
    ! Low speed
    IF ((DXCoil(DXCoilNum)%MSTotCapTempModFacCurveType(SpeedNumLS) == Quadratic).OR. &
      (DXCoil(DXCoilNum)%MSTotCapTempModFacCurveType(SpeedNumLS) == Cubic)) THEN
        TotCapTempModFac = CurveValue(DXCoil(DXCoilNum)%MSCCapFTemp(SpeedNumLS),OutdoorDryBulb)
    ELSEIF (DXCoil(DXCoilNum)%MSTotCapTempModFacCurveType(SpeedNumLS) == Biquadratic) THEN
      TotCapTempModFac = CurveValue(DXCoil(DXCoilNum)%MSCCapFTemp(SpeedNumLS),InletAirDryBulbTemp,OutdoorDryBulb)
    END IF
    !  Get total capacity modifying factor (function of mass flow) for off-rated conditions
    TotCapFlowModFac = CurveValue(DXCoil(DXCoilNum)%MSCCapFFlow(SpeedNumLS),AirMassFlowRatioLS)
    ! Calculate total heating capacity for off-rated conditions
    TotCapLS = DXCoil(DXCoilNum)%MSRatedTotCap(SpeedNumLS) * TotCapFlowModFac * TotCapTempModFac
    ! High speed
    IF ((DXCoil(DXCoilNum)%MSTotCapTempModFacCurveType(SpeedNumHS) == Quadratic).OR. &
      (DXCoil(DXCoilNum)%MSTotCapTempModFacCurveType(SpeedNumHS) == Cubic)) THEN
        TotCapTempModFac = CurveValue(DXCoil(DXCoilNum)%MSCCapFTemp(SpeedNumHS),OutdoorDryBulb)
    ELSEIF (DXCoil(DXCoilNum)%MSTotCapTempModFacCurveType(SpeedNumHS) == Biquadratic) THEN
      TotCapTempModFac = CurveValue(DXCoil(DXCoilNum)%MSCCapFTemp(SpeedNumHS),InletAirDryBulbTemp,OutdoorDryBulb)
    END IF
    !  Get total capacity modifying factor (function of mass flow) for off-rated conditions
    TotCapFlowModFac = CurveValue(DXCoil(DXCoilNum)%MSCCapFFlow(SpeedNumHS),AirMassFlowRatioHS)
    ! Calculate total heating capacity for off-rated conditions
    TotCapHS = DXCoil(DXCoilNum)%MSRatedTotCap(SpeedNumHS) * TotCapFlowModFac * TotCapTempModFac
    ! Calculate electricity consumed. First, get EIR modifying factors for off-rated conditions
    ! Model was extended to accept bi-quadratic curves. This allows sensitivity of the EIR
    ! to the entering dry-bulb temperature as well as the outside dry-bulb temperature. User is
    ! advised to use the bi-quaratic curve if sufficient manufacturer data is available.
    ! Low Speed
    IF ((DXCoil(DXCoilNum)%MSEIRTempModFacCurveType(SpeedNumLS) == Quadratic).OR. &
       (DXCoil(DXCoilNum)%MSEIRTempModFacCurveType(SpeedNumLS) == Cubic)) THEN
      EIRTempModFac = CurveValue(DXCoil(DXCoilNum)%MSEIRFTemp(SpeedNumLS),OutdoorDryBulb)
    ELSEIF (DXCoil(DXCoilNum)%MSEIRTempModFacCurveType(SpeedNumLS) == Biquadratic) THEN
      EIRTempModFac = CurveValue(DXCoil(DXCoilNum)%MSEIRFTemp(SpeedNumLS),InletAirDryBulbTemp,OutdoorDryBulb)
    END IF
    EIRFlowModFac = CurveValue(DXCoil(DXCoilNum)%MSEIRFFlow(SpeedNumLS), AirMassFlowRatioLS)
    EIRLS = 1.0d0/DXCoil(DXCoilNum)%MSRatedCOP(SpeedNumLS) * EIRTempModFac * EIRFlowModFac
    ! High Speed
    IF ((DXCoil(DXCoilNum)%MSEIRTempModFacCurveType(SpeedNumHS) == Quadratic).OR. &
       (DXCoil(DXCoilNum)%MSEIRTempModFacCurveType(SpeedNumHS) == Cubic)) THEN
      EIRTempModFac = CurveValue(DXCoil(DXCoilNum)%MSEIRFTemp(SpeedNumHS),OutdoorDryBulb)
    ELSEIF (DXCoil(DXCoilNum)%MSEIRTempModFacCurveType(SpeedNumHS) == Biquadratic) THEN
      EIRTempModFac = CurveValue(DXCoil(DXCoilNum)%MSEIRFTemp(SpeedNumHS),InletAirDryBulbTemp,OutdoorDryBulb)
    END IF
    EIRFlowModFac = CurveValue(DXCoil(DXCoilNum)%MSEIRFFlow(SpeedNumHS), AirMassFlowRatioHS)
    EIRHS = 1.0d0/DXCoil(DXCoilNum)%MSRatedCOP(SpeedNumHS) * EIRTempModFac * EIRFlowModFac

    ! Calculating adjustment factors for defrost
    ! Calculate delta w through outdoor coil by assuming a coil temp of 0.82*DBT-9.7(F) per DOE2.1E
    OutdoorCoilT = 0.82d0 * OutdoorDryBulb - 8.589d0
    OutdoorCoildw = MAX(1.0d-6,(OutdoorHumRat - PsyWFnTdpPb(OutdoorCoilT,OutdoorPressure,RoutineName)))

    ! Initializing defrost adjustment factors
    LoadDueToDefrostLS = 0.0
    LoadDueToDefrostHS = 0.0
    HeatingCapacityMultiplier = 1.0d0
    FractionalDefrostTime = 0.0
    InputPowerMultiplier = 1.0d0
    DefrostPowerLS = 0.0
    DefrostPowerHS = 0.0

    ! Check outdoor temperature to determine of defrost is active
    IF (OutdoorDryBulb .LE. DXCoil(DXCoilNum)%MaxOATDefrost) THEN
      ! Calculate defrost adjustment factors depending on defrost control type
      IF (DXCoil(DXCoilNum)%DefrostControl .EQ. Timed) THEN
        FractionalDefrostTime = DXCoil(DXCoilNum)%DefrostTime
        HeatingCapacityMultiplier = 0.909d0 - 107.33d0 * OutdoorCoildw
        InputPowerMultiplier = 0.90d0 - 36.45d0*OutdoorCoildw
      ELSE !else defrost control is on-demand
        FractionalDefrostTime = 1.0d0 / (1.0d0 + 0.01446d0 / OutdoorCoildw)
        HeatingCapacityMultiplier = 0.875d0 * ( 1.0d0 - FractionalDefrostTime)
        InputPowerMultiplier = 0.954d0 * ( 1.0d0 - FractionalDefrostTime)
      END IF

      IF (FractionalDefrostTime .GT. 0.0) THEN
        ! Calculate defrost adjustment factors depending on defrost control strategy
        IF (DXCoil(DXCoilNum)%DefrostStrategy .EQ. ReverseCycle) THEN
          DefrostEIRTempModFac = CurveValue(DXCoil(DXCoilNum)%DefrostEIRFT,&
                                MAX(15.555d0,InletAirWetbulbC),MAX(15.555d0,OutdoorDryBulb))
          LoadDueToDefrostLS = (0.01d0 * FractionalDefrostTime) * (7.222d0 - OutdoorDryBulb) * &
                            (DXCoil(DXCoilNum)%MSRatedTotCap(SpeedNumLS)/1.01667d0)
          DefrostPowerLS =  DefrostEIRTempModFac *(DXCoil(DXCoilNum)%MSRatedTotCap(SpeedNumLS)/1.01667d0)* FractionalDefrostTime
          LoadDueToDefrostHS = (0.01d0 * FractionalDefrostTime) * (7.222d0 - OutdoorDryBulb) * &
                            (DXCoil(DXCoilNum)%MSRatedTotCap(SpeedNumHS)/1.01667d0)
          DefrostPowerHS =  DefrostEIRTempModFac *(DXCoil(DXCoilNum)%MSRatedTotCap(SpeedNumHS)/1.01667d0)* FractionalDefrostTime
        ELSE ! Defrost strategy is resistive
          DXCoil(DXCoilNum)%DefrostPower = DXCoil(DXCoilNum)%DefrostCapacity &
                                          * FractionalDefrostTime
        END IF
      ELSE ! Defrost is not active because (OutDryBulbTemp .GT. DXCoil(DXCoilNum)%MaxOATDefrost)
        DXCoil(DXCoilNum)%DefrostPower =  0.0
      END IF
    END IF

    TotCapLS = TotCapLS*HeatingCapacityMultiplier
    TotCapHS = TotCapHS*HeatingCapacityMultiplier

    ! Calculate modified PartLoadRatio due to defrost (reverse-cycle defrost only)
    PLRHeating = MIN(1.0d0,(SpeedRatio + LoadDueToDefrostHS/TotCapHS))
    PLF = CurveValue(DXCoil(DXCoilNum)%MSPLFFPLR(SpeedNumHS),PLRHeating) ! Calculate part-load factor

    IF (PLF < 0.7d0) THEN
      IF (DXCoil(DXCoilNum)%PlrErrIndex == 0) THEN
        CALL ShowWarningMessage('The PLF curve value at high speed for DX multispeed heating coil '//  &
           TRIM(DXCoil(DXCoilNum)%Name)//' ='//TRIM(RoundSigDigits(PLF,2))//  &
           ' for part-load ratio ='//TRIM(RoundSigDigits(PLRHeating,2)))
        CALL ShowContinueError('PLF curve values must be >= 0.7. PLF has been reset to 0.7 and simulation is continuing.')
        CALL ShowContinueError('Check the IO reference manual for PLF curve guidance [Coil:Heating:DX:MultiSpeed].')
        CALL ShowContinueErrorTimeStamp(' ')
      END IF
      CALL ShowRecurringWarningErrorAtEnd('DX heating coil PLF curve < 0.7 warning continues... ',   &
          DXCoil(DXCoilNum)%PlrErrIndex,ReportMinOf=PLF,ReportMaxOf=PLF)
      PLF = 0.7d0
    END IF

    DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction = (PLRHeating / PLF)
    IF (DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction > 1.0d0 .and.   &
        ABS(DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction-1.0d0) > .001d0)THEN
      IF (DXCoil(DXCoilNum)%ErrIndex4 == 0) THEN
        CALL ShowWarningMessage('The runtime fraction at high speed for DX multispeed heating coil '//  &
           TRIM(DXCoil(DXCoilNum)%Name)//&
           ' exceeded 1.0. ['//TRIM(RoundSigDigits(DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction,4))//'].')
        CALL ShowContinueError('Runtime fraction is set to 1.0 and the simulation continues...')
        CALL ShowContinueError('Check the IO reference manual for PLF curve guidance [Coil:Heating:DX:SingleSpeed].')
        CALL ShowContinueErrorTimeStamp(' ')
      END IF
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%Name)//              &
                   ', DX heating coil runtime fraction > 1.0 warning continues...', &
           DXCoil(DXCoilNum)%ErrIndex4,DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction,DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction)
      DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction = 1.0d0 ! Reset coil runtime fraction to 1.0
    ELSEIF (DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction > 1.0d0) THEN
      DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction = 1.0d0 ! Reset coil runtime fraction to 1.0
    END IF

    ! Get full load output and power
    LSFullLoadOutAirEnth = InletAirEnthalpy + TotCapLS/MSHPMassFlowRateLow
    HSFullLoadOutAirEnth = InletAirEnthalpy + TotCapHS/MSHPMassFlowRateHigh
    LSElecHeatingPower = TotCapLS*EIRLS*InputPowerMultiplier
    HSElecHeatingPower = TotCapHS*EIRHS*InputPowerMultiplier
    OutletAirHumRat   = InletAirHumRat

    ! if cycling fan, send coil part-load fraction to on/off fan via HVACDataGlobals
    IF (FanOpMode .EQ. CycFanCycCoil) OnOffFanPartLoadFraction = 1.0d0

    ! Power calculation
    If (.NOT. DXCoil(DXCoilNum)%PLRImpact) Then
      DXCoil(DXCoilNum)%ElecHeatingPower = SpeedRatio*HSElecHeatingPower+(1.0d0-SpeedRatio)*LSElecHeatingPower
    Else
      DXCoil(DXCoilNum)%ElecHeatingPower = DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction*HSElecHeatingPower + &
                   (1.0d0-DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction)*LSElecHeatingPower
    End If

    DXCoil(DXCoilNum)%TotalHeatingEnergyRate = MSHPMassFlowRateHigh*(HSFullLoadOutAirEnth-InletAirEnthalpy)*SpeedRatio + &
                                               MSHPMassFlowRateLow*(LSFullLoadOutAirEnth-InletAirEnthalpy)*(1.0d0-SpeedRatio)
    OutletAirEnthalpy = InletAirEnthalpy + DXCoil(DXCoilNum)%TotalHeatingEnergyRate/DXCoil(DXCoilNum)%InletAirMassFlowRate
    OutletAirTemp = PsyTdbFnHW(OutletAirEnthalpy,OutletAirHumRat,RoutineName)
    FullLoadOutAirRH = PsyRhFnTdbWPb(OutletAirTemp,OutletAirHumRat,OutdoorPressure,RoutineName//':Averageload')
    IF (FullLoadOutAirRH .gt. 1.d0) THEN  ! Limit to saturated conditions at FullLoadOutAirEnth
      OutletAirTemp = PsyTsatFnHPb(FullLoadOutAirEnth,OutdoorPressure,RoutineName)
      OutletAirHumRat = PsyWFnTdbH(OutletAirTemp,FullLoadOutAirEnth,RoutineName)
    END IF

    ! Waste heat calculation
    WasteHeatLS = CurveValue(DXCoil(DXCoilNum)%MSWasteHeat(SpeedNumLS),OutdoorDryBulb,InletAirDryBulbTemp) * &
                  DXCoil(DXCoilNum)%MSWasteHeatFrac(SpeedNumLS)
    WasteHeatHS = CurveValue(DXCoil(DXCoilNum)%MSWasteHeat(SpeedNumHS),OutdoorDryBulb,InletAirDryBulbTemp) * &
                  DXCoil(DXCoilNum)%MSWasteHeatFrac(SpeedNumHS)
    MSHPWasteHeat = (SpeedRatio*WasteHeatHS + (1.0d0-SpeedRatio)*WasteHeatLS)*DXCoil(DXCoilNum)%ElecHeatingPower
    If (DXCoil(DXCoilNum)%FuelType .NE. FuelTypeElectricity) Then
      DXCoil(DXCoilNum)%FuelUsed = DXCoil(DXCoilNum)%ElecHeatingPower
      DXCoil(DXCoilNum)%ElecHeatingPower = 0.0
    End If

    ! Adjust defrost power to correct for DOE-2 bug where defrost power is constant regardless of compressor runtime fraction
    ! Defrosts happen based on compressor run time (frost buildup on outdoor coil), not total elapsed time.
    If (DXCoil(DXCoilNum)%DefrostStrategy .EQ. ReverseCycle) Then
      If (.NOT. DXCoil(DXCoilNum)%PLRImpact) Then
        DXCoil(DXCoilNum)%DefrostPower = DefrostPowerHS * SpeedRatio + DefrostPowerLS * (1.0d0-SpeedRatio)
      Else
        DXCoil(DXCoilNum)%DefrostPower = DefrostPowerHS * DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction + &
                                         DefrostPowerLS * (1.0d0-DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction)
      End If
    End If
    DXCoil(DXCoilNum)%OutletAirTemp     = OutletAirTemp
    DXCoil(DXCoilNum)%OutletAirHumRat   = OutletAirHumRat
    DXCoil(DXCoilNum)%OutletAirEnthalpy = OutletAirEnthalpy
    DXCoil(DXCoilNum)%CrankcaseHeaterPower = 0.0

  ! Stage 1
  Else If (CycRatio > 0.0) Then

    ! for cycling fan, reset mass flow to full on rate
    IF (FanOpMode .EQ. CycFanCycCoil) AirMassFlow = AirMassFlow / CycRatio
    IF (FanOpMode .EQ. ContFanCycCoil) AirMassFlow = MSHPMassFlowRateLow
    !
    ! Check for valid air volume flow per rated total cooling capacity (200 - 600 cfm/ton)
    !
    AirVolumeFlowRate = AirMassFlow/PsyRhoAirFnPbTdbW(OutdoorPressure,InletAirDryBulbTemp, InletAirHumRat,RoutineName)
    !  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
    !  AirVolumeFlowRate = AirMassFlow/PsyRhoAirFnPbTdbW(InletAirPressure,InletAirDryBulbTemp, InletAirHumRat)
    VolFlowperRatedTotCap = AirVolumeFlowRate/DXCoil(DXCoilNum)%MSRatedTotCap(SpeedNumLS)

    IF ((VolFlowperRatedTotCap.LT.MinOperVolFlowPerRatedTotCap).OR.(VolFlowperRatedTotCap.GT.MaxHeatVolFlowPerRatedTotCap)) THEN
      IF (DXCoil(DXCoilNum)%ErrIndex1 == 0) THEN
        CALL ShowWarningMessage(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//&
          '" - Air volume flow rate per watt of rated total heating capacity is out of range at speed 1.')
        CALL ShowContinueErrorTimeStamp(' ')
        CALL ShowContinueError('Expected range for VolumeFlowPerRatedTotalCapacity=['//  &
          TRIM(RoundSigDigits(MinOperVolFlowPerRatedTotCap,3))//'--'//  &
          TRIM(RoundSigDigits(MaxHeatVolFlowPerRatedTotCap,3))//']')
        CALL ShowContinueError('Possible causes include inconsistent air flow rates in system components or')
        CALL ShowContinueError('inconsistent supply air fan operation modes in coil and unitary system objects.')
      END IF
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%DXCoilType)//' "'//TRIM(DXCoil(DXCoilNum)%Name)//&
          '" - Air volume flow rate per watt of rated total heating capacity is out ' //&
          'of range error continues at speed 1...',  &
          DXCoil(DXCoilNum)%ErrIndex1,ReportMinOf=VolFlowperRatedTotCap,ReportMaxOf=VolFlowperRatedTotCap)
    END IF

    ! Get total capacity modifying factor (function of temperature) for off-rated conditions
    ! Model was extended to accept bi-quadratic curves. This allows sensitivity of the heating capacity
    ! to the entering dry-bulb temperature as well as the outside dry-bulb temperature. User is
    ! advised to use the bi-quaratic curve if sufficient manufacturer data is available.
    IF ((DXCoil(DXCoilNum)%MSTotCapTempModFacCurveType(SpeedNumLS) == Quadratic).OR. &
       (DXCoil(DXCoilNum)%MSTotCapTempModFacCurveType(SpeedNumLS) == Cubic)) THEN
      TotCapTempModFac = CurveValue(DXCoil(DXCoilNum)%MSCCapFTemp(SpeedNumLS),OutdoorDryBulb)
    ELSEIF (DXCoil(DXCoilNum)%MSTotCapTempModFacCurveType(SpeedNumLS) == Biquadratic) THEN
      TotCapTempModFac = CurveValue(DXCoil(DXCoilNum)%MSCCapFTemp(SpeedNumLS),InletAirDryBulbTemp,OutdoorDryBulb)
    END IF

    !  Get total capacity modifying factor (function of mass flow) for off-rated conditions
!    AirMassFlowRatio = AirMassFlow/DXCoil(DXCoilNum)%MSRatedAirMassFlowRate(SpeedNumLS)
!    TotCapFlowModFac = CurveValue(DXCoil(DXCoilNum)%MSCCapFFlow(SpeedNumLS),AirMassFlowRatio)
    TotCapFlowModFac = CurveValue(DXCoil(DXCoilNum)%MSCCapFFlow(SpeedNumLS),AirMassFlowRatioLS)
    ! Calculate total heating capacity for off-rated conditions
    TotCap = DXCoil(DXCoilNum)%MSRatedTotCap(SpeedNumLS) * TotCapFlowModFac * TotCapTempModFac

    ! Calculating adjustment factors for defrost
    ! Calculate delta w through outdoor coil by assuming a coil temp of 0.82*DBT-9.7(F) per DOE2.1E
    OutdoorCoilT = 0.82d0 * OutdoorDryBulb - 8.589d0
    OutdoorCoildw = MAX(1.0d-6,(OutdoorHumRat - PsyWFnTdpPb(OutdoorCoilT,OutdoorPressure,RoutineName)))

    ! Initializing defrost adjustment factors
    LoadDueToDefrost = 0.0
    HeatingCapacityMultiplier = 1.0d0
    FractionalDefrostTime = 0.0
    InputPowerMultiplier = 1.0d0

    ! Check outdoor temperature to determine of defrost is active
    IF (OutdoorDryBulb .LE. DXCoil(DXCoilNum)%MaxOATDefrost) THEN
      ! Calculate defrost adjustment factors depending on defrost control type
      IF (DXCoil(DXCoilNum)%DefrostControl .EQ. Timed) THEN
        FractionalDefrostTime = DXCoil(DXCoilNum)%DefrostTime
        HeatingCapacityMultiplier = 0.909d0 - 107.33d0 * OutdoorCoildw
        InputPowerMultiplier = 0.90d0 - 36.45d0*OutdoorCoildw
      ELSE !else defrost control is on-demand
        FractionalDefrostTime = 1.0d0 / (1.0d0 + 0.01446d0 / OutdoorCoildw)
        HeatingCapacityMultiplier = 0.875d0 * ( 1.0d0 - FractionalDefrostTime)
        InputPowerMultiplier = 0.954d0 * ( 1.0d0 - FractionalDefrostTime)
      END IF

      IF (FractionalDefrostTime .GT. 0.0) THEN
        ! Calculate defrost adjustment factors depending on defrost control strategy
        IF (DXCoil(DXCoilNum)%DefrostStrategy .EQ. ReverseCycle) THEN
          LoadDueToDefrost = (0.01d0 * FractionalDefrostTime) * &
                            (7.222d0 - OutdoorDryBulb) * &
                            (DXCoil(DXCoilNum)%MSRatedTotCap(1)/1.01667d0)
          DefrostEIRTempModFac = CurveValue(DXCoil(DXCoilNum)%DefrostEIRFT,&
                                MAX(15.555d0,InletAirWetbulbC),MAX(15.555d0,OutdoorDryBulb))
          DXCoil(DXCoilNum)%DefrostPower =  DefrostEIRTempModFac * &
                                           (DXCoil(DXCoilNum)%MSRatedTotCap(1)/1.01667d0) &
                                           * FractionalDefrostTime
        ELSE ! Defrost strategy is resistive
          DXCoil(DXCoilNum)%DefrostPower = DXCoil(DXCoilNum)%DefrostCapacity &
                                          * FractionalDefrostTime
        END IF
      ELSE ! Defrost is not active because (OutDryBulbTemp .GT. DXCoil(DXCoilNum)%MaxOATDefrost)
        DXCoil(DXCoilNum)%DefrostPower =  0.0
      END IF
    END IF

    ! Modify total heating capacity based on defrost heating capacity multiplier
    TotCap = TotCap * HeatingCapacityMultiplier

    ! Calculate full load outlet conditions
    FullLoadOutAirEnth = InletAirEnthalpy + TotCap/AirMassFlow
    FullLoadOutAirHumRat = InletAirHumRat
    FullLoadOutAirTemp = PsyTdbFnHW(FullLoadOutAirEnth,FullLoadOutAirHumRat,RoutineName)
    FullLoadOutAirRH = PsyRhFnTdbWPb(FullLoadOutAirTemp,FullLoadOutAirHumRat,OutdoorPressure,RoutineName//':fullload')
    !  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
    !  FullLoadOutAirRH = PsyRhFnTdbWPb(FullLoadOutAirTemp,FullLoadOutAirHumRat,InletAirPressure)
    IF (FullLoadOutAirRH .gt. 1.d0) THEN  ! Limit to saturated conditions at FullLoadOutAirEnth
      FullLoadOutAirTemp = PsyTsatFnHPb(FullLoadOutAirEnth,OutdoorPressure,RoutineName)
      !  Eventually inlet air conditions will be used in DX Coil, these lines are commented out and marked with this comment line
      !  FullLoadOutAirTemp = PsyTsatFnHPb(FullLoadOutAirEnth,InletAirPressure)
      FullLoadOutAirHumRat  = PsyWFnTdbH(FullLoadOutAirTemp,FullLoadOutAirEnth,RoutineName)
    END IF

    ! Set outlet conditions from the full load calculation
    OutletAirEnthalpy = FullLoadOutAirEnth
    OutletAirHumRat   = FullLoadOutAirHumRat
    OutletAirTemp     = FullLoadOutAirTemp
    ! Calculate electricity consumed. First, get EIR modifying factors for off-rated conditions
    ! Model was extended to accept bi-quadratic curves. This allows sensitivity of the EIR
    ! to the entering dry-bulb temperature as well as the outside dry-bulb temperature. User is
    ! advised to use the bi-quaratic curve if sufficient manufacturer data is available.
    IF ((DXCoil(DXCoilNum)%MSEIRTempModFacCurveType(1) == Quadratic).OR.  &
        (DXCoil(DXCoilNum)%MSEIRTempModFacCurveType(1) == Cubic)) THEN
      EIRTempModFac = CurveValue(DXCoil(DXCoilNum)%MSEIRFTemp(1),OutdoorDryBulb)
    ELSEIF (DXCoil(DXCoilNum)%MSEIRTempModFacCurveType(1) == Biquadratic) THEN
      EIRTempModFac = CurveValue(DXCoil(DXCoilNum)%MSEIRFTemp(1),InletAirDryBulbTemp,OutdoorDryBulb)
    END IF
    EIRFlowModFac = CurveValue(DXCoil(DXCoilNum)%MSEIRFFlow(1), AirMassFlowRatioLS)
    EIR = 1.0d0/DXCoil(DXCoilNum)%MSRatedCOP(1) * EIRTempModFac * EIRFlowModFac
    ! Calculate modified PartLoadRatio due to defrost (reverse-cycle defrost only)
    PLRHeating = MIN(1.0d0,(CycRatio + LoadDueToDefrost/TotCap))
    PLF = CurveValue(DXCoil(DXCoilNum)%MSPLFFPLR(1),PLRHeating) ! Calculate part-load factor
    IF (FanOpMode .EQ. CycFanCycCoil .AND. CycRatio .EQ. 1.0d0 .AND. PLF .NE. 1.0d0) Then
      IF (DXCoil(DXCoilNum)%PLFErrIndex == 0) THEN
        CALL ShowWarningMessage('The PLF curve value for DX heating coil '//TRIM(DXCoil(DXCoilNum)%Name)//&
                         ' ='//TRIM(RoundSigDigits(PLF,2))//' for part-load ratio = 1')
        CALL ShowContinueError('PLF curve value must be = 1.0 and has been reset to 1.0. Simulation is continuing.')
        CALL ShowContinueErrorTimeStamp(' ')
      END IF
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%Name)//'":'//&
          ' DX heating coil PLF curve value <> 1.0 warning continues...' &
          , DXCoil(DXCoilNum)%PLFErrIndex, PLF, PLF)
      PLF = 1.0d0
    END IF

    IF (PLF < 0.7d0) THEN
      IF (DXCoil(DXCoilNum)%PLRErrIndex == 0) THEN
        CALL ShowWarningMessage('The PLF curve value for DX heating coil '//TRIM(DXCoil(DXCoilNum)%Name)//&
                         ' ='//TRIM(RoundSigDigits(PLF,2))//  &
                         ' for part-load ratio ='//TRIM(RoundSigDigits(PLRHeating,2)))
        CALL ShowContinueError('PLF curve values must be >= 0.7. PLF has been reset to 0.7 and simulation is continuing.')
        CALL ShowContinueError('Check the IO reference manual for PLF curve guidance [Coil:Heating:DX:SingleSpeed].')
        CALL ShowContinueErrorTimeStamp(' ')
      END IF
      CALL ShowRecurringWarningErrorAtEnd('DX heating coil PLF curve < 0.7 warning continues... ',  &
        DXCoil(DXCoilNum)%PLRErrIndex,ReportMinOf=PLF,ReportMaxOf=PLF)
      PLF = 0.7d0
    END IF

    DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction = (PLRHeating / PLF)
    IF (DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction > 1.0d0 .and.   &
        ABS(DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction-1.0d0) > .001 ) THEN
      IF (DXCoil(DXCoilNum)%ErrIndex4 == 0) THEN
        CALL ShowWarningMessage('The runtime fraction for DX heating coil '//TRIM(DXCoil(DXCoilNum)%Name)//&
                         ' exceeded 1.0. ['//TRIM(RoundSigDigits(DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction,4))//'].')
        CALL ShowContinueError('Runtime fraction is set to 1.0 and the simulation continues...')
        CALL ShowContinueError('Check the IO reference manual for PLF curve guidance [Coil:Heating:DX:SingleSpeed].')
        CALL ShowContinueErrorTimeStamp(' ')
      END IF
      CALL ShowRecurringWarningErrorAtEnd(TRIM(DXCoil(DXCoilNum)%Name)//              &
                   ', DX heating coil runtime fraction > 1.0 warning continues...', &
           DXCoil(DXCoilNum)%ErrIndex4,DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction,DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction)
      DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction = 1.0d0 ! Reset coil runtime fraction to 1.0
    ELSEIF (DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction > 1.0d0) THEN
      DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction = 1.0d0 ! Reset coil runtime fraction to 1.0
    END IF
    ! if cycling fan, send coil part-load fraction to on/off fan via HVACDataGlobals
    IF (FanOpMode .EQ. CycFanCycCoil) OnOffFanPartLoadFraction = PLF
    DXCoil(DXCoilNum)%ElecHeatingPower = TotCap * EIR * DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction * InputPowerMultiplier

    ! Calculate crankcase heater power using the runtime fraction for this DX heating coil only if there is no companion DX coil.
    ! Else use the largest runtime fraction of this DX heating coil and the companion DX cooling coil.

    IF(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil .EQ. 0) THEN
      DXCoil(DXCoilNum)%CrankcaseHeaterPower = CrankcaseHeatingPower * &
                                           (1.0d0 - DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction)
    ELSE
      DXCoil(DXCoilNum)%CrankcaseHeaterPower = CrankcaseHeatingPower * &
                                           (1.0d0 - MAX(DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction, &
                                          DXCoil(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil)%CoolingCoilRuntimeFraction))
    END IF

    DXCoil(DXCoilNum)%TotalHeatingEnergyRate = AirMassFlow * (FullLoadOutAirEnth - InletAirEnthalpy)*CycRatio
    IF (FanOpMode .EQ. ContFanCycCoil) THEN
      OutletAirEnthalpy = InletAirEnthalpy+DXCoil(DXCoilNum)%TotalHeatingEnergyRate/DXCoil(DXCoilNum)%InletAirMassFlowRate
      OutletAirTemp = PsyTdbFnHW(OutletAirEnthalpy,OutletAirHumRat,RoutineName)
    END IF
    MSHPWasteHeat = CurveValue(DXCoil(DXCoilNum)%MSWasteHeat(SpeedNumLS),OutdoorDryBulb,InletAirDryBulbTemp) * &
                  DXCoil(DXCoilNum)%MSWasteHeatFrac(SpeedNumLS)*DXCoil(DXCoilNum)%ElecHeatingPower

    If (DXCoil(DXCoilNum)%FuelType .NE. FuelTypeElectricity) Then
      DXCoil(DXCoilNum)%FuelUsed = DXCoil(DXCoilNum)%ElecHeatingPower
      DXCoil(DXCoilNum)%ElecHeatingPower = 0.0
    End If
    ! Adjust defrost power to correct for DOE-2 bug where defrost power is constant regardless of compressor runtime fraction
    ! Defrosts happen based on compressor run time (frost buildup on outdoor coil), not total elapsed time.
    DXCoil(DXCoilNum)%DefrostPower = DXCoil(DXCoilNum)%DefrostPower * DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction

    DXCoil(DXCoilNum)%OutletAirTemp     = OutletAirTemp
    DXCoil(DXCoilNum)%OutletAirHumRat   = OutletAirHumRat
    DXCoil(DXCoilNum)%OutletAirEnthalpy = OutletAirEnthalpy
  End If

ELSE

  ! DX coil is off; just pass through conditions
  DXCoil(DXCoilNum)%OutletAirEnthalpy = DXCoil(DXCoilNum)%InletAirEnthalpy
  DXCoil(DXCoilNum)%OutletAirHumRat   = DXCoil(DXCoilNum)%InletAirHumRat
  DXCoil(DXCoilNum)%OutletAirTemp     = DXCoil(DXCoilNum)%InletAirTemp

  DXCoil(DXCoilNum)%ElecHeatingPower       = 0.0
  DXCoil(DXCoilNum)%FuelUsed               = 0.0
  DXCoil(DXCoilNum)%TotalHeatingEnergyRate = 0.0
  DXCoil(DXCoilNum)%DefrostPower           = 0.0

! Calculate crankcase heater power using the runtime fraction for this DX heating coil (here DXHeatingCoilRTF=0) if
! there is no companion DX coil, or the runtime fraction of the companion DX cooling coil (here DXCoolingCoilRTF>=0).
  IF(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil .EQ. 0) THEN
    DXCoil(DXCoilNum)%CrankcaseHeaterPower = CrankcaseHeatingPower
  ELSE
    DXCoil(DXCoilNum)%CrankcaseHeaterPower = CrankcaseHeatingPower * &
                                            (1.d0-DXCoil(DXCoil(DXCoilNum)%CompanionUpstreamDXCoil)%CoolingCoilRuntimeFraction)
  END IF

END IF ! end of on/off if - else

DXCoilOutletTemp(DXCoilNum)     = DXCoil(DXCoilNum)%OutletAirTemp
DXCoilOutletHumRat(DXCoilNum)   = DXCoil(DXCoilNum)%OutletAirHumRat
DXCoil(DXCoilNum)%PartLoadRatio = PLRHeating
DXCoilFanOpMode(DXCoilNum)      = FanOpMode
DXCoilPartLoadRatio(DXCoilNum)  = PLRHeating

RETURN

END SUBROUTINE CalcMultiSpeedDXCoilHeating

SUBROUTINE UpdateDXCoil(DXCoilNum)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Fred Buhl
          !       DATE WRITTEN   May 2000
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! This subroutine is for passing results to the outlet air node.

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE DataContaminantBalance, ONLY: Contaminant

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  INTEGER, INTENT (IN) :: DXCoilNum ! number of the current fan coil unit being simulated

          ! SUBROUTINE PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER :: AirOutletNode ! air outlet node number
  INTEGER :: AirInletNode ! air inlet node number

AirOutletNode = DXCoil(DXCoilNum)%AirOutNode
AirInletNode = DXCoil(DXCoilNum)%AirInNode
! changed outputs
Node(AirOutletNode)%Enthalpy     = DXCoil(DXCoilNum)%OutletAirEnthalpy
Node(AirOutletNode)%Temp         = DXCoil(DXCoilNum)%OutletAirTemp
Node(AirOutletNode)%HumRat       = DXCoil(DXCoilNum)%OutletAirHumRat
Node(AirOutletNode)%MassFlowRate = DXCoil(DXCoilNum)%InletAirMassFlowRate
! pass through outputs
Node(AirOutletNode)%Quality              = Node(AirInletNode)%Quality
Node(AirOutletNode)%Press                = Node(AirInletNode)%Press
Node(AirOutletNode)%MassFlowRateMin      = Node(AirInletNode)%MassFlowRateMin
Node(AirOutletNode)%MassFlowRateMax      = Node(AirInletNode)%MassFlowRateMax
Node(AirOutletNode)%MassFlowRateMinAvail = Node(AirInletNode)%MassFlowRateMinAvail
Node(AirOutletNode)%MassFlowRateMaxAvail = Node(AirInletNode)%MassFlowRateMaxAvail

IF (Contaminant%CO2Simulation) Then
  Node(AirOutletNode)%CO2 = Node(AirInletNode)%CO2
End If
IF (Contaminant%GenericContamSimulation) Then
  Node(AirOutletNode)%GenContam = Node(AirInletNode)%GenContam
End If

RETURN
END SUBROUTINE UpdateDXCoil

SUBROUTINE ReportDXCoil(DXCoilNum)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Fred Buhl
          !       DATE WRITTEN   May 2000
          !       MODIFIED       Richard Raustad/Don Shirey Oct 2001, Feb 2004
          !                      Feb 2005 M. J. Witte, GARD Analytics, Inc.
          !                        Always update evap value to support new coil type COIL:DX:MultiMode:CoolingEmpirical:
          !                      Lixing Gu. Jan. 5, 2007, pass information to the AirflowNetwork model
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! Fills some of the report variables for the DX coils

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE DataHVACGlobals, ONLY: TimeStepSys, DXElecCoolingPower, DXElecHeatingPower
  USE Psychrometrics,  ONLY: RhoH2O
  USE DataWater,       ONLY: WaterStorage
  USE DataAirLoop,     ONLY: LoopDXCoilRTF

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  INTEGER, INTENT (IN) :: DXCoilNum ! number of the current fan coil unit being simulated

          ! SUBROUTINE PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
REAL(r64) :: rhoWater
REAL(r64) :: Tavg
REAL(r64) :: SpecHumOut
REAL(r64) :: SpecHumIn
REAL(r64) :: ReportingConstant ! Number of seconds per HVAC system time step, to convert from W (J/s) to J

ReportingConstant = TimeStepSys*SecInHour

SELECT CASE(DXCoil(DXCoilNum)%DXCoilType_Num)

CASE (CoilDX_HeatingEmpirical)
  DXCoil(DXCoilNum)%TotalHeatingEnergy = DXCoil(DXCoilNum)%TotalHeatingEnergyRate*ReportingConstant
  DXCoil(DXCoilNum)%ElecHeatingConsumption = DXCoil(DXCoilNum)%ElecHeatingPower*ReportingConstant
  DXCoil(DXCoilNum)%DefrostConsumption = DXCoil(DXCoilNum)%DefrostPower*ReportingConstant
  DXCoil(DXCoilNum)%CrankcaseHeaterConsumption = DXCoil(DXCoilNum)%CrankcaseHeaterPower*ReportingConstant
  DXElecHeatingPower = DXCoil(DXCoilNum)%ElecHeatingPower + DXCoil(DXCoilNum)%CrankcaseHeaterPower
CASE (CoilDX_MultiSpeedHeating)
  DXCoil(DXCoilNum)%TotalHeatingEnergy = DXCoil(DXCoilNum)%TotalHeatingEnergyRate*ReportingConstant
  If (DXCoil(DXCoilNum)%FuelType .EQ. FuelTypeElectricity) Then
    DXCoil(DXCoilNum)%ElecHeatingConsumption = DXCoil(DXCoilNum)%ElecHeatingPower*ReportingConstant
  Else
    DXCoil(DXCoilNum)%FuelConsumed = DXCoil(DXCoilNum)%FuelUsed*ReportingConstant
  End If
  DXCoil(DXCoilNum)%DefrostConsumption = DXCoil(DXCoilNum)%DefrostPower*ReportingConstant
  DXCoil(DXCoilNum)%CrankcaseHeaterConsumption = DXCoil(DXCoilNum)%CrankcaseHeaterPower*ReportingConstant
  DXElecHeatingPower = DXCoil(DXCoilNum)%ElecHeatingPower + DXCoil(DXCoilNum)%CrankcaseHeaterPower
CASE (CoilDX_MultiSpeedCooling)
  DXCoil(DXCoilNum)%TotalCoolingEnergy = DXCoil(DXCoilNum)%TotalCoolingEnergyRate*ReportingConstant
  DXCoil(DXCoilNum)%SensCoolingEnergy = DXCoil(DXCoilNum)%SensCoolingEnergyRate*ReportingConstant
  DXCoil(DXCoilNum)%LatCoolingEnergy = DXCoil(DXCoilNum)%TotalCoolingEnergy - DXCoil(DXCoilNum)%SensCoolingEnergy
  DXCoil(DXCoilNum)%CrankcaseHeaterConsumption = DXCoil(DXCoilNum)%CrankcaseHeaterPower*ReportingConstant
  DXElecCoolingPower = DXCoil(DXCoilNum)%ElecCoolingPower
  DXCoil(DXCoilNum)%EvapCondPumpElecConsumption = DXCoil(DXCoilNum)%EvapCondPumpElecPower*ReportingConstant
  DXCoil(DXCoilNum)%EvapWaterConsump = DXCoil(DXCoilNum)%EvapWaterConsumpRate*ReportingConstant
  If (DXCoil(DXCoilNum)%FuelType .EQ. FuelTypeElectricity) Then
    DXCoil(DXCoilNum)%ElecCoolingConsumption = DXCoil(DXCoilNum)%ElecCoolingPower*ReportingConstant
  Else
    DXCoil(DXCoilNum)%FuelConsumed = DXCoil(DXCoilNum)%FuelUsed*ReportingConstant
  End If
  If (ANY(DXCoil(DXCoilNum)%CondenserType == EvapCooled)) THEN
    DXCoil(DXCoilNum)%BasinHeaterConsumption = DXCoil(DXCoilNum)%BasinHeaterPower*ReportingConstant
  Endif
CASE (CoilDX_HeatPumpWaterHeater)
! water heating energy for HP water heater DX Coil condenser
  DXCoil(DXCoilNum)%TotalHeatingEnergy = DXCoil(DXCoilNum)%TotalHeatingEnergyRate*ReportingConstant
! water heating power for HP water heater
  DXCoil(DXCoilNum)%ElecWaterHeatingConsumption = DXCoil(DXCoilNum)%ElecWaterHeatingPower*ReportingConstant
! other usual DX cooling coil outputs
  DXCoil(DXCoilNum)%TotalCoolingEnergy = DXCoil(DXCoilNum)%TotalCoolingEnergyRate*ReportingConstant
  DXCoil(DXCoilNum)%SensCoolingEnergy = DXCoil(DXCoilNum)%SensCoolingEnergyRate*ReportingConstant
  DXCoil(DXCoilNum)%LatCoolingEnergy = DXCoil(DXCoilNum)%TotalCoolingEnergy - DXCoil(DXCoilNum)%SensCoolingEnergy
  DXCoil(DXCoilNum)%ElecCoolingConsumption = DXCoil(DXCoilNum)%ElecCoolingPower*ReportingConstant
  DXCoil(DXCoilNum)%CrankcaseHeaterConsumption = DXCoil(DXCoilNum)%CrankcaseHeaterPower*ReportingConstant
! DXElecCoolingPower global is only used for air-to-air cooling and heating coils
  DXElecCoolingPower = 0.0
CASE DEFAULT
  DXCoil(DXCoilNum)%TotalCoolingEnergy = DXCoil(DXCoilNum)%TotalCoolingEnergyRate*ReportingConstant
  DXCoil(DXCoilNum)%SensCoolingEnergy = DXCoil(DXCoilNum)%SensCoolingEnergyRate*ReportingConstant
  DXCoil(DXCoilNum)%LatCoolingEnergy = DXCoil(DXCoilNum)%TotalCoolingEnergy - DXCoil(DXCoilNum)%SensCoolingEnergy
  DXCoil(DXCoilNum)%ElecCoolingConsumption = DXCoil(DXCoilNum)%ElecCoolingPower*ReportingConstant
  DXCoil(DXCoilNum)%CrankcaseHeaterConsumption = DXCoil(DXCoilNum)%CrankcaseHeaterPower*ReportingConstant
  DXElecCoolingPower = DXCoil(DXCoilNum)%ElecCoolingPower
  DXCoil(DXCoilNum)%EvapCondPumpElecConsumption = DXCoil(DXCoilNum)%EvapCondPumpElecPower*ReportingConstant
  DXCoil(DXCoilNum)%EvapWaterConsump = DXCoil(DXCoilNum)%EvapWaterConsumpRate*ReportingConstant
  If (ANY(DXCoil(DXCoilNum)%CondenserType == EvapCooled)) THEN
    DXCoil(DXCoilNum)%BasinHeaterConsumption = DXCoil(DXCoilNum)%BasinHeaterPower*ReportingConstant
  Endif
END SELECT

IF (DXCoil(DXCoilNum)%CondensateCollectMode == CondensateToTank) THEN
  ! calculate and report condensation rates  (how much water extracted from the air stream)
  ! water flow of water in m3/s for water system interactions
  !  put here to catch all types of DX coils
  Tavg =( DXCoil(DXCoilNum)%InletAirTemp - DXCoil(DXCoilNum)%OutletAirTemp)/2.0d0
  rhoWater = RhoH2O(Tavg)
  SpecHumIn = DXCoil(DXCoilNum)%InletAirHumRat / ( 1.d0 + DXCoil(DXCoilNum)%InletAirHumRat) !eq. 9b ASHRAE HOF 2001 page 6.8
  SpecHumOut = DXCoil(DXCoilNum)%OutletAirHumRat / ( 1.d0 + DXCoil(DXCoilNum)%OutletAirHumRat )
  !  mdot * del HumRat / rho water
  DXCoil(DXCoilNum)%CondensateVdot = MAX(0.0d0, (DXCoil(DXCoilNum)%InletAirMassFlowRate *   &
            (SpecHumIn - SpecHumOut) / rhoWater) )
  DXCoil(DXCoilNum)%CondensateVol  = DXCoil(DXCoilNum)%CondensateVdot *ReportingConstant

  WaterStorage(DXCoil(DXCoilNum)%CondensateTankID)%VdotAvailSupply(DXCoil(DXCoilNum)%CondensateTankSupplyARRID) &
      =  DXCoil(DXCoilNum)%CondensateVdot
  WaterStorage(DXCoil(DXCoilNum)%CondensateTankID)%TwaterSupply(DXCoil(DXCoilNum)%CondensateTankSupplyARRID) &
      =  DXCoil(DXCoilNum)%OutletAirTemp
ENDIF

LoopDXCoilRTF = MAX(DXCoil(DXCoilNum)%CoolingCoilRuntimeFraction,DXCoil(DXCoilNum)%HeatingCoilRuntimeFraction)

RETURN
END SUBROUTINE ReportDXCoil

! ======================  Utility routines ======================================

SUBROUTINE GetDXCoilIndex(DXCoilName,DXCoilIndex,ErrorsFound,ThisObjectType,SuppressWarning)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Richard Raustad
          !       DATE WRITTEN   March 2005
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! This subroutine sets an index for a given DX Coil -- issues error message if that
          ! DX Coil is not a legal DX Coil.

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor, ONLY: FindItemInList
  !USE DataGlobals, ONLY: ShowSevereError

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: DXCoilName
  INTEGER, INTENT(INOUT)       :: DXCoilIndex
  LOGICAL, INTENT(INOUT)       :: ErrorsFound
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: ThisObjectType
  LOGICAL, INTENT(IN), OPTIONAL :: SuppressWarning

          ! SUBROUTINE PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
          ! na
  IF (GetCoilsInputFlag) THEN
    CALL GetDXCoils
    GetCoilsInputFlag = .FALSE.
  END IF

  DXCoilIndex = FindItemInList(DXCoilName,DXCoil%Name,NumDXCoils)
  IF (DXCoilIndex == 0) THEN
    IF(PRESENT(SuppressWarning)) THEN
!     No warning printed if only searching for the existence of a DX Coil
    ELSE
      IF (PRESENT(ThisObjectType)) THEN
        CALL ShowSevereError(TRIM(ThisObjectType)//', GetDXCoilIndex: DX Coil not found='//TRIM(DXCoilName))
      ELSE
        CALL ShowSevereError('GetDXCoilIndex: DX Coil not found='//TRIM(DXCoilName))
      ENDIF
    END IF
    ErrorsFound = .TRUE.
  ENDIF

  RETURN

END SUBROUTINE GetDXCoilIndex

FUNCTION GetCoilCapacity(CoilType,CoilName,ErrorsFound) RESULT(CoilCapacity)

          ! FUNCTION INFORMATION:
          !       AUTHOR         Linda Lawrie
          !       DATE WRITTEN   February 2006
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS FUNCTION:
          ! This function looks up the coil capacity for the given coil and returns it.  If
          ! incorrect coil type or name is given, errorsfound is returned as true and capacity is returned
          ! as negative.

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor,  ONLY: FindItem,SameString

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

          ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: CoilType     ! must match coil types in this module
  CHARACTER(len=*), INTENT(IN) :: CoilName     ! must match coil names for the coil type
  LOGICAL, INTENT(INOUT)       :: ErrorsFound  ! set to true if problem
  REAL(r64)                    :: CoilCapacity ! returned capacity of matched coil

          ! FUNCTION PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  INTEGER :: WhichCoil

  ! Obtains and Allocates DXCoils
  IF (GetCoilsInputFlag) THEN
    CALL GetDXCoils
    GetCoilsInputFlag = .FALSE.
  END IF

  IF (SameString(CoilType,'Coil:Heating:DX:SingleSpeed') .or.   &
      SameString(CoilType,'Coil:Cooling:DX:SingleSpeed')) THEN
    WhichCoil=FindItem(CoilName,DXCoil%Name,NumDXCoils)
    IF (WhichCoil /= 0) THEN
      CoilCapacity=DXCoil(WhichCoil)%RatedTotCap(1)
    ENDIF
  ELSE IF(SameString(CoilType,'Coil:Cooling:DX:TwoStageWithHumidityControlMode')) THEN
    WhichCoil=FindItem(CoilName,DXCoil%Name,NumDXCoils)
    IF (WhichCoil /= 0) THEN
      CoilCapacity=DXCoil(WhichCoil)%RatedTotCap(DXCoil(WhichCoil)%NumCapacityStages)
    ENDIF
  ELSE
    WhichCoil=0
  ENDIF

  IF (WhichCoil == 0) THEN
    CALL ShowSevereError('GetCoilCapacity: Could not find Coil, Type="'//TRIM(CoilType)//'" Name="'//TRIM(CoilName)//'"')
    CALL ShowContinueError('... returning capacity as -1000.')
    ErrorsFound=.true.
    CoilCapacity=-1000.
  ENDIF

  RETURN

END FUNCTION GetCoilCapacity


FUNCTION GetCoilCapacityByIndexType(CoilIndex,CoilType_Num,ErrorsFound) RESULT(CoilCapacity)

          ! FUNCTION INFORMATION:
          !       AUTHOR         Richard Raustad
          !       DATE WRITTEN   October 2010
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS FUNCTION:
          ! This function looks up the coil capacity for the given coil and returns it.  If
          ! incorrect coil index or type is given, errorsfound is returned as true and capacity is returned
          ! as negative.

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor,  ONLY: FindItem,SameString

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

          ! FUNCTION ARGUMENT DEFINITIONS:
  INTEGER, INTENT(IN)    :: CoilIndex    ! must match coil index for the coil type
  INTEGER, INTENT(IN)    :: CoilType_Num ! must match coil types in this module
  LOGICAL, INTENT(INOUT) :: ErrorsFound  ! set to true if problem
  REAL(r64)              :: CoilCapacity ! returned capacity of matched coil

          ! FUNCTION PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! FUNCTION LOCAL VARIABLE DECLARATIONS:
          ! na

  ! Obtains and Allocates DXCoils
  IF (GetCoilsInputFlag) THEN
    CALL GetDXCoils
    GetCoilsInputFlag = .FALSE.
  END IF

  IF(CoilIndex .EQ. 0)THEN
    CALL ShowSevereError('GetCoilCapacityByIndexType: Invalid index passed = 0')
    CALL ShowContinueError('... returning capacity as -1000.')
    ErrorsFound=.TRUE.
    CoilCapacity=-1000.
    RETURN
  END IF

  IF (CoilType_Num /= DXCoil(CoilIndex)%DXCoilType_Num) THEN
    CALL ShowSevereError('GetCoilCapacityByIndexType: Index passed does not match DX Coil type passed.')
    CALL ShowContinueError('... returning capacity as -1000.')
    ErrorsFound=.TRUE.
    CoilCapacity=-1000.
  ELSE
    CoilCapacity=DXCoil(CoilIndex)%RatedTotCap(DXCoil(CoilIndex)%NumCapacityStages)
  ENDIF

  RETURN

END FUNCTION GetCoilCapacityByIndexType

FUNCTION GetCoilTypeNum(CoilType,CoilName,ErrorsFound,PrintWarning) RESULT(TypeNum)

          ! FUNCTION INFORMATION:
          !       AUTHOR         R. Raustad - FSEC
          !       DATE WRITTEN   August 2008
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS FUNCTION:
          ! This function looks up the integerized coil type for the given coil and returns it.  If
          ! incorrect coil type or name is given, errorsfound is returned as true and capacity is returned
          ! as negative.

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor,  ONLY: FindItemInList

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

          ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: CoilType     ! must match coil types in this module
  CHARACTER(len=*), INTENT(IN) :: CoilName     ! must match coil names for the coil type
  LOGICAL, INTENT(INOUT)       :: ErrorsFound  ! set to true if problem
  LOGICAL, OPTIONAL, INTENT(IN):: PrintWarning ! prints warning when true
  INTEGER                      :: TypeNum      ! returned integerized type of matched coil

          ! FUNCTION PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  INTEGER :: WhichCoil
  LOGICAL :: PrintMessage

  ! Obtains and Allocates DXCoils
  IF (GetCoilsInputFlag) THEN
    CALL GetDXCoils
    GetCoilsInputFlag = .FALSE.
  END IF

  IF (PRESENT(PrintWarning)) THEN
    PrintMessage = PrintWarning
  ELSE
    PrintMessage = .TRUE.
  END IF

  WhichCoil=FindItemInList(CoilName,DXCoil%Name,NumDXCoils)
  IF (WhichCoil /= 0) THEN
    TypeNum=DXCoil(WhichCoil)%DXCoilType_Num
  ELSE
   IF(PrintMessage)THEN
     CALL ShowSevereError('GetCoilTypeNum: Could not find Coil, Type="'//TRIM(CoilType)//'" Name="'//TRIM(CoilName)//'"')
   END IF
     ErrorsFound=.true.
     TypeNum=0
  ENDIF

  RETURN

END FUNCTION GetCoilTypeNum

FUNCTION GetMinOATCompressor(CoilType,CoilName,ErrorsFound) RESULT(MinOAT)

          ! FUNCTION INFORMATION:
          !       AUTHOR         Linda Lawrie
          !       DATE WRITTEN   February 2006
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS FUNCTION:
          ! This function looks up the the min oat for the heating coil compressor and returns it.  If
          ! incorrect coil type or name is given, errorsfound is returned as true and value is returned
          ! as negative.

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor,  ONLY: FindItem,SameString

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

          ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: CoilType     ! must match coil types in this module
  CHARACTER(len=*), INTENT(IN) :: CoilName     ! must match coil names for the coil type
  LOGICAL, INTENT(INOUT)       :: ErrorsFound  ! set to true if problem
  REAL(r64)                    :: MinOAT       ! returned min oa temperature of matched coil

          ! FUNCTION PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  INTEGER :: WhichCoil

  ! Obtains and Allocates DXCoils
  IF (GetCoilsInputFlag) THEN
    CALL GetDXCoils
    GetCoilsInputFlag = .FALSE.
  END IF

  IF (SameString(CoilType,'Coil:Heating:DX:SingleSpeed')) THEN
    WhichCoil=FindItem(CoilName,DXCoil%Name,NumDXCoils)
    IF (WhichCoil /= 0) THEN
      MinOAT=DXCoil(WhichCoil)%MinOATCompressor
    ENDIF
  ELSE
    WhichCoil=0
  ENDIF

  IF (WhichCoil == 0) THEN
    CALL ShowSevereError('Could not find CoilType="'//TRIM(CoilType)//'" with Name="'//TRIM(CoilName)//'"')
    CALL ShowContinueError('... returning Min OAT as -1000.')
    ErrorsFound=.true.
    MinOAT=-1000.
  ENDIF

  RETURN

END FUNCTION GetMinOATCompressor

FUNCTION GetCoilInletNode(CoilType,CoilName,ErrorsFound) RESULT(NodeNumber)

          ! FUNCTION INFORMATION:
          !       AUTHOR         Linda Lawrie
          !       DATE WRITTEN   February 2006
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS FUNCTION:
          ! This function looks up the given coil and returns the inlet node number.  If
          ! incorrect coil type or name is given, errorsfound is returned as true and node number is returned
          ! as zero.

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor,  ONLY: FindItemInList

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

          ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: CoilType     ! must match coil types in this module
  CHARACTER(len=*), INTENT(IN) :: CoilName     ! must match coil names for the coil type
  LOGICAL, INTENT(INOUT)       :: ErrorsFound  ! set to true if problem
  INTEGER                      :: NodeNumber   ! returned node number of matched coil

          ! FUNCTION PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  INTEGER :: WhichCoil

  ! Obtains and Allocates DXCoils
  IF (GetCoilsInputFlag) THEN
    CALL GetDXCoils
    GetCoilsInputFlag = .FALSE.
  END IF

  WhichCoil=FindItemInList(CoilName,DXCoil%Name,NumDXCoils)
  IF (WhichCoil /= 0) THEN
    NodeNumber=DXCoil(WhichCoil)%AirInNode
  ELSE
    CALL ShowSevereError('GetCoilInletNode: Could not find Coil, Type="'//TRIM(CoilType)//'" Name="'//TRIM(CoilName)//'"')
    ErrorsFound=.true.
    NodeNumber=0
  ENDIF

  RETURN

END FUNCTION GetCoilInletNode

FUNCTION GetCoilOutletNode(CoilType,CoilName,ErrorsFound) RESULT(NodeNumber)

          ! FUNCTION INFORMATION:
          !       AUTHOR         Linda Lawrie
          !       DATE WRITTEN   February 2006
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS FUNCTION:
          ! This function looks up the given coil and returns the inlet node number.  If
          ! incorrect coil type or name is given, errorsfound is returned as true and node number is returned
          ! as zero.

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor,  ONLY: FindItemInList

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

          ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: CoilType     ! must match coil types in this module
  CHARACTER(len=*), INTENT(IN) :: CoilName     ! must match coil names for the coil type
  LOGICAL, INTENT(INOUT)       :: ErrorsFound  ! set to true if problem
  INTEGER                      :: NodeNumber   ! returned node number of matched coil

          ! FUNCTION PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  INTEGER :: WhichCoil

  ! Obtains and Allocates DXCoils
  IF (GetCoilsInputFlag) THEN
    CALL GetDXCoils
    GetCoilsInputFlag = .FALSE.
  END IF

  WhichCoil=FindItemInList(CoilName,DXCoil%Name,NumDXCoils)
  IF (WhichCoil /= 0) THEN
    NodeNumber=DXCoil(WhichCoil)%AirOutNode
  ELSE
    CALL ShowSevereError('GetCoilOutletNode: Could not find Coil, Type="'//TRIM(CoilType)//'" Name="'//TRIM(CoilName)//  &
          '" when accessing coil outlet node number.')
    ErrorsFound=.true.
    NodeNumber=0
  ENDIF

  RETURN

END FUNCTION GetCoilOutletNode

FUNCTION GetCoilCondenserInletNode(CoilType,CoilName,ErrorsFound) RESULT(CondNode)

          ! FUNCTION INFORMATION:
          !       AUTHOR         R. Raustad
          !       DATE WRITTEN   January 2007
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS FUNCTION:
          ! This function looks up the given coil and returns the condenser inlet node.  If
          ! incorrect coil type or name is given, errorsfound is returned as true.

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor,  ONLY: FindItemInList

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

          ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: CoilType     ! must match coil types in this module
  CHARACTER(len=*), INTENT(IN) :: CoilName     ! must match coil names for the coil type
  LOGICAL, INTENT(INOUT)       :: ErrorsFound  ! set to true if problem
  INTEGER                      :: CondNode     ! returned condenser node number of matched coil

          ! FUNCTION PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  INTEGER :: WhichCoil

  ! Obtains and Allocates DXCoils
  IF (GetCoilsInputFlag) THEN
    CALL GetDXCoils
    GetCoilsInputFlag = .FALSE.
  END IF

  WhichCoil=FindItemInList(CoilName,DXCoil%Name,NumDXCoils)
  IF (WhichCoil /= 0) THEN
    CondNode=DXCoil(WhichCoil)%CondenserInletNodeNum(1)
   ELSE
    CALL ShowSevereError('GetCoilCondenserInletNode: Invalid DX Coil, Type= "'//TRIM(CoilType)//'" Name="'//TRIM(CoilName)//'"')
    ErrorsFound=.true.
    CondNode=0
  ENDIF

  RETURN

END FUNCTION GetCoilCondenserInletNode

FUNCTION GetDXCoilBypassedFlowFrac(CoilType,CoilName,ErrorsFound) RESULT(BypassFraction)

          ! FUNCTION INFORMATION:
          !       AUTHOR         R. Raustad
          !       DATE WRITTEN   June 2007
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS FUNCTION:
          ! This function looks up the given coil and returns the bypassed air flow fraction.
          ! Bypassed air flow fraction can only be greater than 0 for multimode DX cooling coils and is typical for 1st stage
          ! If incorrect coil type or name is given, errorsfound is returned as true.

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor,  ONLY: FindItemInList

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

          ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: CoilType       ! must match coil types in this module
  CHARACTER(len=*), INTENT(IN) :: CoilName       ! must match coil names for the coil type
  LOGICAL, INTENT(INOUT)       :: ErrorsFound    ! set to true if problem
  REAL(r64)                    :: BypassFraction ! returned bypass air fraction of matched coil

          ! FUNCTION PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  INTEGER :: WhichCoil

  ! Obtains and Allocates DXCoils
  IF (GetCoilsInputFlag) THEN
    CALL GetDXCoils
    GetCoilsInputFlag = .FALSE.
  END IF

  WhichCoil=FindItemInList(CoilName,DXCoil%Name,NumDXCoils)
  IF (WhichCoil /= 0) THEN
    BypassFraction=DXCoil(WhichCoil)%BypassedFlowFrac(1)
   ELSE
    CALL ShowSevereError('GetDXCoilBypassedFlowFrac: Invalid DX Coil Type="'//TRIM(CoilType)//'" Name="'//TRIM(CoilName)//'"')
    ErrorsFound=.true.
    BypassFraction=0.0
  ENDIF

  RETURN

END FUNCTION GetDXCoilBypassedFlowFrac

FUNCTION GetHPCoolingCoilIndex(HeatingCoilType, HeatingCoilName, HeatingCoilIndex) RESULT(DXCoolingCoilIndex)

          ! FUNCTION INFORMATION:
          !       AUTHOR         R. Raustad
          !       DATE WRITTEN   February 2007
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS FUNCTION:
          ! This function looks up the given DX heating coil and returns the companion DX cooling coil.

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor,            ONLY: FindItemInList, SameString
  USE DataBranchNodeConnections, ONLY: CompSets, NumCompSets

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

          ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: HeatingCoilType  ! Type of DX heating coil used in HP
  CHARACTER(len=*), INTENT(IN) :: HeatingCoilName  ! Name of DX heating coil used in HP
  INTEGER,          INTENT(IN) :: HeatingCoilIndex ! Index of DX heating coil used in HP

          ! FUNCTION PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  INTEGER                      :: DXCoolingCoilIndex   ! Index of HP DX cooling coil returned from this function
  INTEGER                      :: WhichComp            ! DO loop counter to find correct comp set
  INTEGER                      :: WhichCompanionComp   ! DO loop counter to find companion coil comp set
  INTEGER                      :: WhichHXAssistedComp  ! DO loop counter when DX coil is used in a HX assisted cooling coil
  CHARACTER(len=MaxNameLength) :: CompSetsParentType   ! Parent object type which uses DX heating coil pass into this function
  CHARACTER(len=MaxNameLength) :: CompSetsParentName   ! Parent object name which uses DX heating coil pass into this function
  CHARACTER(len=MaxNameLength) :: HXCompSetsParentType ! Used when DX cooling coil is a child of a HX assisted cooling coil
  CHARACTER(len=MaxNameLength) :: HXCompSetsParentName ! Used when DX cooling coil is a child of a HX assisted cooling coil

  DXCoolingCoilIndex = 0

  DO WhichComp = 1, NumCompSets
    IF(.NOT. SameString(HeatingCoilType, CompSets(WhichComp)%CTYPE) .OR. &
       .NOT. SameString(HeatingCoilName, CompSets(WhichComp)%CNAME))CYCLE
      CompSetsParentType = CompSets(WhichComp)%PARENTCTYPE
      CompSetsParentName = CompSets(WhichComp)%PARENTCNAME
      IF(SameString(CompSetsParentType,'AirLoopHVAC:UnitaryHeatPump:AirToAir') .OR. &
         SameString(CompSetsParentType,'ZoneHVAC:PackagedTerminalHeatPump') .OR. &
         SameString(CompSetsParentType,'AirLoopHVAC:UnitaryHeatPump:AirToAir:MultiSpeed') .OR. &
         SameString(CompSetsParentType,'AirLoopHVAC:UnitaryHeatCool:VAVChangeoverBypass'))THEN
!       Search for DX cooling coils
        DO WhichCompanionComp = 1, NumCompSets
          IF(.NOT. SameString(CompSets(WhichCompanionComp)%PARENTCNAME,CompSetsParentName) .OR. &
             .NOT. SameString(CompSets(WhichCompanionComp)%CTYPE,'Coil:Cooling:DX:SingleSpeed'))CYCLE
            DXCoolingCoilIndex = FindItemInList(CompSets(WhichCompanionComp)%CNAME,DXCoil%Name,NumDXCoils)
            EXIT
        END DO
        DO WhichCompanionComp = 1, NumCompSets
          IF(.NOT. SameString(CompSets(WhichCompanionComp)%PARENTCNAME,CompSetsParentName) .OR. &
             .NOT. SameString(CompSets(WhichCompanionComp)%CTYPE,'Coil:Cooling:DX:MultiSpeed'))CYCLE
            DXCoolingCoilIndex = FindItemInList(CompSets(WhichCompanionComp)%CNAME,DXCoil%Name,NumDXCoils)
            EXIT
        END DO
!       Search for Heat Exchanger Assisted DX cooling coils
        IF(DXCoolingCoilIndex .EQ. 0)THEN
          DO WhichHXAssistedComp = 1, NumCompSets
            IF(.NOT. SameString(CompSets(WhichHXAssistedComp)%PARENTCNAME,CompSetsParentName) .OR. &
               .NOT. SameString(CompSets(WhichHXAssistedComp)%CTYPE,'CoilSystem:Cooling:DX:HeatExchangerAssisted'))CYCLE
            HXCompSetsParentType = CompSets(WhichHXAssistedComp)%CTYPE
            HXCompSetsParentName = CompSets(WhichHXAssistedComp)%CNAME
            DO WhichCompanionComp = 1, NumCompSets
              IF(.NOT. SameString(CompSets(WhichCompanionComp)%PARENTCNAME,HXCompSetsParentName) .OR. &
                 .NOT. SameString(CompSets(WhichCompanionComp)%CTYPE,'Coil:Cooling:DX:SingleSpeed'))CYCLE
                DXCoolingCoilIndex = FindItemInList(CompSets(WhichCompanionComp)%CNAME,DXCoil%Name,NumDXCoils)
                EXIT
            END DO
            EXIT
          END DO
        END IF
      ELSE
!     ErrorFound, Coil:Heating:DX:SingleSpeed is used in wrong type of parent object (should never get here)
        CALL ShowSevereError('Configuration error in '//TRIM(CompSetsParentType)//' "'//TRIM(CompSetsParentName)//'"')
        CALL ShowContinueError('DX heating coil not allowed in this configuration.')
        CALL ShowFatalError('Preceding condition(s) causes termination.')
      END IF
    EXIT
  END DO

! Check and warn user is crankcase heater power or max OAT for crankcase heater differs in DX cooling and heating coils
  IF(DXCoolingCoilIndex .GT. 0)THEN
    IF(DXCoil(DXCoolingCoilIndex)%CrankcaseHeaterCapacity /= 0.) THEN
      IF(DXCoil(DXCoolingCoilIndex)%CrankcaseHeaterCapacity /= DXCoil(HeatingCoilIndex)%CrankcaseHeaterCapacity .OR. &
         DXCoil(DXCoolingCoilIndex)%MaxOATCrankcaseHeater /= DXCoil(HeatingCoilIndex)%MaxOATCrankcaseHeater) THEN
        CALL ShowWarningError('Crankcase heater capacity or max outdoor temp for crankcase heater operation specified in')
        CALL ShowContinueError('Coil:Cooling:DX:SingleSpeed = '//TRIM(DXCoil(DXCoolingCoilIndex)%Name))
        CALL ShowContinueError('is different than that specified in Coil:Heating:DX:SingleSpeed = ' &
                               //TRIM(HeatingCoilName)//'.')
        CALL ShowContinueError('Both of these DX coils are part of '//TRIM(CompSetsParentType)//' = '//  &
                               TRIM(CompSetsParentName)//'.')
        CALL ShowContinueError('The value specified in the DX heating coil will be used and the simulation continues...')
      END IF
    END IF
  END IF

  RETURN

END FUNCTION GetHPCoolingCoilIndex

FUNCTION GetDXCoilNumberOfSpeeds(CoilType,CoilName,ErrorsFound) RESULT(NumberOfSpeeds)

          ! FUNCTION INFORMATION:
          !       AUTHOR         L. Gu
          !       DATE WRITTEN   July 2007
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS FUNCTION:
          ! This function looks up the given coil and returns the number of speeds for multispeed coils.
          ! If incorrect coil type or name is given, errorsfound is returned as true.

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor,  ONLY: FindItemInList

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

          ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: CoilType       ! must match coil types in this module
  CHARACTER(len=*), INTENT(IN) :: CoilName       ! must match coil names for the coil type
  LOGICAL, INTENT(INOUT)       :: ErrorsFound    ! set to true if problem
  INTEGER                      :: NumberOfSpeeds ! returned the number of speed of matched coil

          ! FUNCTION PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  INTEGER :: WhichCoil

  ! Obtains and Allocates DXCoils
  IF (GetCoilsInputFlag) THEN
    CALL GetDXCoils
    GetCoilsInputFlag = .FALSE.
  END IF

  WhichCoil=FindItemInList(CoilName,DXCoil%Name,NumDXCoils)
  IF (WhichCoil /= 0) THEN
    NumberOfSpeeds=DXCoil(WhichCoil)%NumOfSpeeds
   ELSE
    CALL ShowSevereError('GetDXCoilNumberOfSpeeds: Invalid DX Coil Type="'//TRIM(CoilType)//'" Name="'//TRIM(CoilName)//'"')
    ErrorsFound=.true.
    NumberOfSpeeds=0
  ENDIF

  RETURN

END FUNCTION GetDXCoilNumberOfSpeeds

SUBROUTINE SetDXCoolingCoilData(DXCoilNum,ErrorsFound,HeatingCoilPLFCurvePTR,CondenserType,CondenserInletNodeNum, &
                                MaxOATCrankcaseHeater,MinOATCooling,MaxOATCooling,MinOATHeating,MaxOATHeating, &
                                HeatingPerformanceOATType,DefrostStrategy,DefrostControl,DefrostEIRPtr, &
                                DefrostFraction,DefrostCapacity,MaxOATDefrost,CoolingCoilPresent,HeatingCoilPresent)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Richard Raustad, FSEC
          !       DATE WRITTEN   December 2008
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! This routine was designed to allow the DX coil to access information from a gas or
          ! electric heating coil when these coils are each used in a parent object.
          ! Also, this is an illustration of setting Data from an outside source.

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE General, ONLY: TrimSigDigits

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  INTEGER, INTENT(IN)       :: DXCoilNum   ! Number of DX Cooling Coil
  LOGICAL, INTENT(INOUT)    :: ErrorsFound ! Set to true if certain errors found
  INTEGER,OPTIONAL          :: HeatingCoilPLFCurvePTR   ! Parameter equivalent of heating coil PLR curve index
  INTEGER,OPTIONAL          :: CondenserType            ! Parameter equivalent of condenser type parameter
  INTEGER,OPTIONAL          :: CondenserInletNodeNum    ! Parameter equivalent of condenser inlet node number
  REAL(r64),OPTIONAL        :: MaxOATCrankcaseHeater    ! Parameter equivalent of condenser Max OAT for Crank Case Heater temp
  REAL(r64),OPTIONAL        :: MaxOATCooling        ! Parameter equivalent of condenser Max OAT for compressor cooling operation
  REAL(r64),OPTIONAL        :: MinOATCooling        ! Parameter equivalent of condenser Min OAT for compressor cooling operation
  REAL(r64),OPTIONAL        :: MaxOATHeating        ! Parameter equivalent of condenser Max OAT for compressor heating operation
  REAL(r64),OPTIONAL        :: MinOATHeating        ! Parameter equivalent of condenser Min OAT for compressor heating operation
  INTEGER,OPTIONAL          :: HeatingPerformanceOATType ! Parameter equivalent to condenser entering air temp type (1-db, 2=wb)
  INTEGER,OPTIONAL          :: DefrostStrategy
  INTEGER,OPTIONAL          :: DefrostControl
  INTEGER,OPTIONAL          :: DefrostEIRPtr
  REAL(r64),OPTIONAL        :: DefrostFraction
  REAL(r64),OPTIONAL        :: DefrostCapacity
  REAL(r64),OPTIONAL        :: MaxOATDefrost
  LOGICAL,OPTIONAL          :: CoolingCoilPresent
  LOGICAL,OPTIONAL          :: HeatingCoilPresent
          ! SUBROUTINE PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
          ! na

  ! Obtains and Allocates DXCoils
  IF (GetCoilsInputFlag) THEN
    CALL GetDXCoils
    GetCoilsInputFlag = .FALSE.
  END IF

  IF (DXCoilNum <= 0 .or. DXCoilNum > NumDXCoils) THEN
    CALL ShowSevereError('SetDXCoolingCoilData: called with DX Cooling Coil Number out of range='//  &
         TRIM(TrimSigDigits(DXCoilNum))//' should be >0 and <'//TRIM(TrimSigDigits(NumDXCoils)))
    ErrorsFound=.true.
    RETURN
  ENDIF

  IF (PRESENT(HeatingCoilPLFCurvePTR)) THEN
    DXCoil(DXCoilNum)%HeatingCoilPLFCurvePTR=HeatingCoilPLFCurvePTR
  ENDIF

  IF (PRESENT(CondenserType)) THEN
    DXCoil(DXCoilNum)%CondenserType=CondenserType
  ENDIF

  IF (PRESENT(CondenserInletNodeNum)) THEN
    DXCoil(DXCoilNum)%CondenserInletNodeNum(1)=CondenserInletNodeNum
  ENDIF

  IF (PRESENT(MaxOATCrankcaseHeater)) THEN
    DXCoil(DXCoilNum)%MaxOATCrankcaseHeater=MaxOATCrankcaseHeater
  ENDIF

  IF (PRESENT(MaxOATCooling)) THEN
    DXCoil(DXCoilNum)%MaxOATCompressor=MaxOATCooling
  ENDIF

  IF (PRESENT(MaxOATHeating)) THEN
    DXCoil(DXCoilNum)%MaxOATCompressor=MaxOATHeating
  ENDIF

  IF (PRESENT(MinOATCooling)) THEN
    DXCoil(DXCoilNum)%MinOATCompressor=MinOATCooling
  ENDIF

  IF (PRESENT(MinOATHeating)) THEN
    DXCoil(DXCoilNum)%MinOATCompressor=MinOATHeating
  ENDIF

  IF(PRESENT(HeatingPerformanceOATType))THEN
    DXCoil(DXCoilNum)%HeatingPerformanceOATType=HeatingPerformanceOATType
  END IF

  IF(PRESENT(DefrostStrategy))THEN
    DXCoil(DXCoilNum)%DefrostStrategy=DefrostStrategy
  END IF

  IF(PRESENT(DefrostControl))THEN
    DXCoil(DXCoilNum)%DefrostControl=DefrostControl
  END IF

  IF(PRESENT(DefrostEIRPtr))THEN
    DXCoil(DXCoilNum)%DefrostEIRFT=DefrostEIRPtr
  END IF

  IF(PRESENT(DefrostFraction))THEN
    DXCoil(DXCoilNum)%DefrostTime=DefrostFraction
  END IF

  IF(PRESENT(DefrostCapacity))THEN
    DXCoil(DXCoilNum)%DefrostCapacity=DefrostCapacity
  END IF

  IF(PRESENT(MaxOATDefrost))THEN
    DXCoil(DXCoilNum)%MaxOATDefrost=MaxOATDefrost
  END IF

  IF(PRESENT(CoolingCoilPresent))THEN
    DXCoil(DXCoilNum)%CoolingCoilPresent=CoolingCoilPresent
  END IF

  IF(PRESENT(HeatingCoilPresent))THEN
    DXCoil(DXCoilNum)%HeatingCoilPresent=HeatingCoilPresent
  END IF

  RETURN

END SUBROUTINE SetDXCoolingCoilData

SUBROUTINE SetCoilSystemHeatingDXFlag(CoilType, CoilName)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         B. Griffith
          !       DATE WRITTEN   Jan. 2012
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! inform DX heating coil that is is part of a CoilSystem:Heating:DX
          ! and therefore it need not find its companion cooling coil

          ! METHODOLOGY EMPLOYED:
          ! set value of logical flag FindCompanionUpStreamCoil to true

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor,  ONLY: FindItemInList

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: CoilType     ! must match coil types in this module
  CHARACTER(len=*), INTENT(IN) :: CoilName     ! must match coil names for the coil type

          ! SUBROUTINE PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS:
          ! na

          ! DERIVED TYPE DEFINITIONS:
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER :: WhichCoil

  ! Obtains and Allocates DXCoils
  IF (GetCoilsInputFlag) THEN
    CALL GetDXCoils
    GetCoilsInputFlag = .FALSE.
  END IF

  WhichCoil=FindItemInList(CoilName,DXCoil%Name,NumDXCoils)
  IF (WhichCoil /= 0) THEN
    DXCoil(WhichCoil)%FindCompanionUpStreamCoil = .FALSE.
  ELSE
    CALL ShowSevereError('SetCoilSystemHeatingDXFlag: Could not find Coil, Type="'//TRIM(CoilType)//'"Name="'//TRIM(CoilName)//'"')

  ENDIF

  RETURN

END SUBROUTINE SetCoilSystemHeatingDXFlag

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

END MODULE DXCoils
