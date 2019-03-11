MPICC = mpic++
CXX = g++
CFLAGS = -O3 -I./include
LDFLAGS = -lumfpack
DEBUG = -g -ggdb

# inputs for the executables
PROCS = 4
ARGS = input.inp

OBJS = main.o 									\
	   src/utilities.o 							\
	   src/sparseDiracMatrixClass.o				\
	   src/importanceSamplingIntegrator.o 		\
	   src/gaussianRandomGenerator.o 			\
	   src/ranlxd.o 							\
	   src/determinantWrapper.o 				\
	   src/suiteSparseDet.o

EXE = gn

default: $(EXE) 

$(EXE): $(OBJS) 
	$(MPICC) $^ -o $@  $(LDFLAGS)

%.o: %.cc
	$(MPICC) $< -c -o $@  $(CFLAGS) 

main.o: include/utilities.h 					\
		include/sparseDiracMatrixClass.h 		\
		include/importanceSamplingIntegrator.h	\
		include/gaussianRandomGenerator.h		\
		include/mpiCommManager.h				\
		include/conditionalOStreamClass.h


src/utilities.o: include/utilities.h
src/sparseDiracMatrixClass.o: include/sparseDiracMatrixClass.h
src/ranlxd.o: include/ranlxd.h
src/suiteSparseDet.o: include/suiteSparseDet.h	
src/gaussianRandomGenerator.o: include/gaussianRandomGenerator.h 			\
							   include/ranlxd.h
src/importanceSamplingIntegrator.o: include/importanceSamplingIntegrator.h 	\
									include/gaussianRandomGenerator.h		\
									include/observables.h
src/determinantWrapper.o: include/determinantWrapper.h						\
						  include/suiteSparseDet.h

run: clean default
	mpirun -np $(PROCS) $(EXE) $(ARGS)

tests:

clean:
	@rm -f *~ $(EXE) $(OBJS) *.dat *.png 

.PHONY: clean default test

