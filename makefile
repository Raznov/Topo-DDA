
CC = g++
cuFFT_INCLUDE_PATH = "/opt/apps/software/Compiler/GCC/8.3.0/CUDA/10.1.168/include"
cuFFT_LIB_DIR = "/opt/apps/software/Compiler/GCC/8.3.0/CUDA/10.1.168/lib64"
FLAG = -fopenmp -lpthread -O3 -Wall -Wno-maybe-uninitialized -fPIC -m64 -I${cuFFT_INCLUDE_PATH} -L${cuFFT_LIB_DIR} -Wl,--no-as-needed -lcuda -lcudart -lcufft -lm -ldl
OBJ = test2d_periodic_TiO2.o Structure.o Space.o Model.o EvoModel.o tools.o AProductCore.o DDAModel.o EvoDDAModel.o Objective.o ObjectiveDDAModel.o SiCi.o kernel.o

te : $(OBJ)
	$(CC) $(FLAG) -o te $(OBJ)
test2d_periodic_TiO2.o : test2d_periodic_TiO2.cpp
	$(CC) $(FLAG) -c test2d_periodic_TiO2.cpp
Structure.o : Structure.cpp
	$(CC) $(FLAG) -c Structure.cpp
Space.o : Space.cpp
	$(CC) $(FLAG) -c Space.cpp
Model.o : Model.cpp
	$(CC) $(FLAG) -c Model.cpp
EvoModel.o : EvoModel.cpp
	$(CC) $(FLAG) -c EvoModel.cpp
tools.o : tools.cpp
	$(CC) $(FLAG) -c tools.cpp
Objective.o : Objective.cpp
	$(CC) $(FLAG) -c Objective.cpp
AProductCore.o : AProductCore.cpp
	$(CC) $(FLAG) -c AProductCore.cpp
DDAModel.o : DDAModel.cpp
	$(CC) $(FLAG) -c DDAModel.cpp
EvoDDAModel.o : EvoDDAModel.cpp
	$(CC) $(FLAG) -c EvoDDAModel.cpp
ObjectiveDDAModel.o : ObjectiveDDAModel.cpp
	$(CC) $(FLAG) -c ObjectiveDDAModel.cpp
SiCi.o : SiCi.cpp
	$(CC) $(FLAG) -c SiCi.cpp
kernel.o : kernel.cu
	nvcc -w -c kernel.cu


.PHONY : clean
clean:
	rm $(OBJ) te

.PHONY : cleandata
cleandata:
	rm EXT_mine.txt time_mine_bicgstab.txt wavelength.txt iteration.txt


