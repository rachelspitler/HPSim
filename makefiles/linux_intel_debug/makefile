# Makefile created by mkmf.rb $Id: mkmf,v 14.0 2007/03/20 22:13:27 fms Exp $ 

SRCROOT = /media/Store_Shared/workspace/YorkHP_IVF/src/


.DEFAULT:
	-touch $@
all: /media/Store_Shared/workspace/YorkHP_IVF/YorkHP_IVF_intel_debug
AIRPR.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./AIRPR.FOR
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./AIRPR.FOR
Accumulator.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./Accumulator.f90 DataModuleSim.o Refprops.o
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./Accumulator.f90
AirTempLoop.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./AirTempLoop.f90 Refprops.o Evaporator.o Accumulator.o DataModuleSim.o
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./AirTempLoop.f90
CapillaryTube.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./CapillaryTube.F90 Refprops.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./CapillaryTube.F90
ChargeLoop.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./ChargeLoop.f90 DataModuleSim.o
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./ChargeLoop.f90
CoilCalc.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./CoilCalc.f90 DataModuleSim.o OilMixtureCalc.o AIRPR.o
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./CoilCalc.f90
Compressor.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./Compressor.f90 DataModuleSim.o Refprops.o
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./Compressor.f90
Condenser.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./Condenser.f90 CoilCalc.o Refprops.o AIRPR.o OilMixtureCalc.o UnitConvert.o
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./Condenser.f90
DRIVER.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./DRIVER.F90 Refprops.o RefName.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./DRIVER.F90
DataErrorTracking.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./DataErrorTracking.f90
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./DataErrorTracking.f90
DataGlobals.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./DataGlobals.f90
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./DataGlobals.f90
DataHVACGlobals.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./DataHVACGlobals.f90 DataGlobals.o
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./DataHVACGlobals.f90
DataModuleSim.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./DataModuleSim.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./DataModuleSim.F90
DataSizing.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./DataSizing.f90 DataGlobals.o
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./DataSizing.f90
DataStringGlobals.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./DataStringGlobals.f90
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./DataStringGlobals.f90
Distributor.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./Distributor.f90 Refprops.o
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./Distributor.f90
DumpHPsimOutputs.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./DumpHPsimOutputs.f90 Refprops.o AIRPR.o DataModuleSim.o
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./DumpHPsimOutputs.f90
EnergyPlus.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./EnergyPlus.f90 InputProcessor.o DataStringGlobals.o DataGlobals.o
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./EnergyPlus.f90
Evaporator.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./Evaporator.f90 CoilCalc.o DataModuleSim.o Refprops.o AIRPR.o OilMixtureCalc.o ReversingValve.o UnitConvert.o
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./Evaporator.f90
FlowRateLoop.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./FlowRateLoop.f90 Refprops.o Condenser.o Compressor.o ShortTube.o CapillaryTube.o DataModuleSim.o
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./FlowRateLoop.f90
FrostModel.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./FrostModel.F90 DataModuleSim.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./FrostModel.F90
GUESS1.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./GUESS1.FOR
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./GUESS1.FOR
GUESS3.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./GUESS3.FOR
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./GUESS3.FOR
General.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./General.f90 DataGlobals.o InputProcessor.o DataStringGlobals.o DataHVACGlobals.o
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./General.f90
GeneralRoutines.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./GeneralRoutines.f90
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./GeneralRoutines.f90
GetHPsimInputs.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./GetHPsimInputs.f90 DataModuleSim.o RefName.o
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./GetHPsimInputs.f90
HPdesignMod.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./HPdesignMod.f90 Refprops.o Condenser.o Evaporator.o ShortTube.o CapillaryTube.o TXV.o Accumulator.o DataModuleSim.o GeneralRoutines.o
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./HPdesignMod.f90
InputProcessor.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./InputProcessor.f90 DataStringGlobals.o DataGlobals.o DataSizing.o RefName.o
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./InputProcessor.f90
LiquidLineFilterDrier.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./LiquidLineFilterDrier.f90
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./LiquidLineFilterDrier.f90
MatrixCalc.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./MatrixCalc.f90
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./MatrixCalc.f90
ORNLsolver.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./ORNLsolver.f90 Refprops.o GetHPsimInputs.o Compressor.o Condenser.o Evaporator.o Accumulator.o UnitConvert.o DataModuleSim.o FrostModel.o GeneralRoutines.o
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./ORNLsolver.f90
OilMixtureCalc.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./OilMixtureCalc.f90 Refprops.o
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./OilMixtureCalc.f90
PsyChart.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./PsyChart.f90
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./PsyChart.f90
Psychometric.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./Psychometric.f90
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./Psychometric.f90
RefName.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./RefName.f90
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./RefName.f90
Refprops.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./Refprops.f90 DataGlobals.o InputProcessor.o General.o
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./Refprops.f90
ReversingValve.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./ReversingValve.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./ReversingValve.F90
ShortTube.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./ShortTube.f90 Refprops.o
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./ShortTube.f90
TXV.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./TXV.F90 Refprops.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./TXV.F90
UnitConvert.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./UnitConvert.f90
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./UnitConvert.f90
UtilityRoutines.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./UtilityRoutines.f90 DataStringGlobals.o DataErrorTracking.o InputProcessor.o DataGlobals.o
	$(FC) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./UtilityRoutines.f90
WATPR.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./WATPR.FOR
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./WATPR.FOR
ZERO3.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./ZERO3.FOR
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./ZERO3.FOR
ZEROCH.o: /media/Store_Shared/workspace/YorkHP_IVF/src/./ZEROCH.FOR
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	/media/Store_Shared/workspace/YorkHP_IVF/src/./ZEROCH.FOR
./ReversingValve.F90: /media/Store_Shared/workspace/YorkHP_IVF/src/./ReversingValve.F90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./ReversingValve.F90 .
./ZEROCH.FOR: /media/Store_Shared/workspace/YorkHP_IVF/src/./ZEROCH.FOR
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./ZEROCH.FOR .
./GUESS1.FOR: /media/Store_Shared/workspace/YorkHP_IVF/src/./GUESS1.FOR
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./GUESS1.FOR .
./Accumulator.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./Accumulator.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./Accumulator.f90 .
./AIRPR.FOR: /media/Store_Shared/workspace/YorkHP_IVF/src/./AIRPR.FOR
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./AIRPR.FOR .
./ChargeLoop.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./ChargeLoop.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./ChargeLoop.f90 .
./DataGlobals.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./DataGlobals.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./DataGlobals.f90 .
./WATPR.FOR: /media/Store_Shared/workspace/YorkHP_IVF/src/./WATPR.FOR
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./WATPR.FOR .
./CoilCalc.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./CoilCalc.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./CoilCalc.f90 .
./UnitConvert.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./UnitConvert.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./UnitConvert.f90 .
./UtilityRoutines.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./UtilityRoutines.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./UtilityRoutines.f90 .
./DataStringGlobals.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./DataStringGlobals.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./DataStringGlobals.f90 .
./FlowRateLoop.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./FlowRateLoop.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./FlowRateLoop.f90 .
./OilMixtureCalc.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./OilMixtureCalc.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./OilMixtureCalc.f90 .
./DataErrorTracking.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./DataErrorTracking.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./DataErrorTracking.f90 .
./LiquidLineFilterDrier.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./LiquidLineFilterDrier.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./LiquidLineFilterDrier.f90 .
./DataHVACGlobals.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./DataHVACGlobals.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./DataHVACGlobals.f90 .
./PsyChart.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./PsyChart.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./PsyChart.f90 .
./General.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./General.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./General.f90 .
./Distributor.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./Distributor.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./Distributor.f90 .
./Evaporator.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./Evaporator.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./Evaporator.f90 .
./GUESS3.FOR: /media/Store_Shared/workspace/YorkHP_IVF/src/./GUESS3.FOR
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./GUESS3.FOR .
./DataModuleSim.F90: /media/Store_Shared/workspace/YorkHP_IVF/src/./DataModuleSim.F90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./DataModuleSim.F90 .
./DataSizing.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./DataSizing.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./DataSizing.f90 .
./MatrixCalc.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./MatrixCalc.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./MatrixCalc.f90 .
./CapillaryTube.F90: /media/Store_Shared/workspace/YorkHP_IVF/src/./CapillaryTube.F90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./CapillaryTube.F90 .
./AirTempLoop.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./AirTempLoop.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./AirTempLoop.f90 .
./FrostModel.F90: /media/Store_Shared/workspace/YorkHP_IVF/src/./FrostModel.F90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./FrostModel.F90 .
./GetHPsimInputs.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./GetHPsimInputs.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./GetHPsimInputs.f90 .
./Compressor.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./Compressor.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./Compressor.f90 .
./HPdesignMod.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./HPdesignMod.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./HPdesignMod.f90 .
./TXV.F90: /media/Store_Shared/workspace/YorkHP_IVF/src/./TXV.F90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./TXV.F90 .
./Refprops.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./Refprops.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./Refprops.f90 .
./InputProcessor.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./InputProcessor.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./InputProcessor.f90 .
./DRIVER.F90: /media/Store_Shared/workspace/YorkHP_IVF/src/./DRIVER.F90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./DRIVER.F90 .
./Psychometric.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./Psychometric.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./Psychometric.f90 .
./DumpHPsimOutputs.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./DumpHPsimOutputs.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./DumpHPsimOutputs.f90 .
./GeneralRoutines.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./GeneralRoutines.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./GeneralRoutines.f90 .
./RefName.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./RefName.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./RefName.f90 .
./EnergyPlus.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./EnergyPlus.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./EnergyPlus.f90 .
./ZERO3.FOR: /media/Store_Shared/workspace/YorkHP_IVF/src/./ZERO3.FOR
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./ZERO3.FOR .
./ORNLsolver.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./ORNLsolver.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./ORNLsolver.f90 .
./ShortTube.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./ShortTube.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./ShortTube.f90 .
./Condenser.f90: /media/Store_Shared/workspace/YorkHP_IVF/src/./Condenser.f90
	cp /media/Store_Shared/workspace/YorkHP_IVF/src/./Condenser.f90 .
SRC = /media/Store_Shared/workspace/YorkHP_IVF/src/./Evaporator.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./InputProcessor.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./ZERO3.FOR /media/Store_Shared/workspace/YorkHP_IVF/src/./EnergyPlus.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./Distributor.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./Refprops.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./CapillaryTube.F90 /media/Store_Shared/workspace/YorkHP_IVF/src/./DataSizing.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./ORNLsolver.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./GeneralRoutines.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./UnitConvert.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./PsyChart.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./CoilCalc.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./MatrixCalc.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./ReversingValve.F90 /media/Store_Shared/workspace/YorkHP_IVF/src/./HPdesignMod.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./GetHPsimInputs.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./TXV.F90 /media/Store_Shared/workspace/YorkHP_IVF/src/./Psychometric.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./DataGlobals.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./Condenser.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./Accumulator.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./DataHVACGlobals.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./ChargeLoop.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./UtilityRoutines.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./LiquidLineFilterDrier.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./DumpHPsimOutputs.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./General.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./GUESS1.FOR /media/Store_Shared/workspace/YorkHP_IVF/src/./FlowRateLoop.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./DRIVER.F90 /media/Store_Shared/workspace/YorkHP_IVF/src/./Compressor.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./AirTempLoop.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./DataErrorTracking.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./ShortTube.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./RefName.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./OilMixtureCalc.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./DataModuleSim.F90 /media/Store_Shared/workspace/YorkHP_IVF/src/./AIRPR.FOR /media/Store_Shared/workspace/YorkHP_IVF/src/./FrostModel.F90 /media/Store_Shared/workspace/YorkHP_IVF/src/./WATPR.FOR /media/Store_Shared/workspace/YorkHP_IVF/src/./DataStringGlobals.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./GUESS3.FOR /media/Store_Shared/workspace/YorkHP_IVF/src/./ZEROCH.FOR
OBJ = Evaporator.o InputProcessor.o ZERO3.o EnergyPlus.o Distributor.o Refprops.o CapillaryTube.o DataSizing.o ORNLsolver.o GeneralRoutines.o UnitConvert.o PsyChart.o CoilCalc.o MatrixCalc.o ReversingValve.o HPdesignMod.o GetHPsimInputs.o TXV.o Psychometric.o DataGlobals.o Condenser.o Accumulator.o DataHVACGlobals.o ChargeLoop.o UtilityRoutines.o LiquidLineFilterDrier.o DumpHPsimOutputs.o General.o GUESS1.o FlowRateLoop.o DRIVER.o Compressor.o AirTempLoop.o DataErrorTracking.o ShortTube.o RefName.o OilMixtureCalc.o DataModuleSim.o AIRPR.o FrostModel.o WATPR.o DataStringGlobals.o GUESS3.o ZEROCH.o
OFF = /media/Store_Shared/workspace/YorkHP_IVF/src/./ReversingValve.F90 /media/Store_Shared/workspace/YorkHP_IVF/src/./ZEROCH.FOR /media/Store_Shared/workspace/YorkHP_IVF/src/./GUESS1.FOR /media/Store_Shared/workspace/YorkHP_IVF/src/./Accumulator.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./AIRPR.FOR /media/Store_Shared/workspace/YorkHP_IVF/src/./ChargeLoop.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./DataGlobals.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./WATPR.FOR /media/Store_Shared/workspace/YorkHP_IVF/src/./CoilCalc.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./UnitConvert.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./UtilityRoutines.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./DataStringGlobals.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./FlowRateLoop.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./OilMixtureCalc.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./DataErrorTracking.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./LiquidLineFilterDrier.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./DataHVACGlobals.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./PsyChart.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./General.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./Distributor.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./Evaporator.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./GUESS3.FOR /media/Store_Shared/workspace/YorkHP_IVF/src/./DataModuleSim.F90 /media/Store_Shared/workspace/YorkHP_IVF/src/./DataSizing.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./MatrixCalc.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./CapillaryTube.F90 /media/Store_Shared/workspace/YorkHP_IVF/src/./AirTempLoop.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./FrostModel.F90 /media/Store_Shared/workspace/YorkHP_IVF/src/./GetHPsimInputs.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./Compressor.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./HPdesignMod.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./TXV.F90 /media/Store_Shared/workspace/YorkHP_IVF/src/./Refprops.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./InputProcessor.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./DRIVER.F90 /media/Store_Shared/workspace/YorkHP_IVF/src/./Psychometric.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./DumpHPsimOutputs.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./GeneralRoutines.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./RefName.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./EnergyPlus.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./ZERO3.FOR /media/Store_Shared/workspace/YorkHP_IVF/src/./ORNLsolver.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./ShortTube.f90 /media/Store_Shared/workspace/YorkHP_IVF/src/./Condenser.f90
clean: neat
	-rm -f .cppdefs $(OBJ) /media/Store_Shared/workspace/YorkHP_IVF/YorkHP_IVF_intel_debug
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
/media/Store_Shared/workspace/YorkHP_IVF/YorkHP_IVF_intel_debug: $(OBJ) 
	$(LD) $(OBJ) -o /media/Store_Shared/workspace/YorkHP_IVF/YorkHP_IVF_intel_debug  $(LDFLAGS)
