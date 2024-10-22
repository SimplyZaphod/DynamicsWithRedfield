ifort -c InputOutput.f90 -traceback -qopenmp
ifort -c MathFunc.f90 -traceback -qopenmp
ifort -c conversion.f90 -traceback -qopenmp
ifort -c Redfieldmodule.f90 -traceback -qopenmpEvolutionModules.f90 -traceback -qopenmp
ifort RedfieldDynamicsGeneral.f90 -o RedfieldDynamicsGeneral.e -llapack MathFunc.o conversion.o InputOutput.o EvolutionModules.o Redfieldmodule.o -traceback -qopenmp
ifort CorrelationFunctionGeneral.f90 -o CorrelationFunctionGeneral.e -llapack MathFunc.o conversion.o InputOutput.o EvolutionModules.o Redfieldmodule.o -traceback -qopenmp
