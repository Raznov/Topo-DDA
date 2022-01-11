
CC = g++
##cuFFT_INCLUDE_PATH = "/opt/apps/software/Compiler/GCC/8.3.0/CUDA/10.1.168/include"
cuFFT_INCLUDE_PATH = "/usr/local/cuda/include"
##cuFFT_LIB_DIR = "/opt/apps/software/Compiler/GCC/8.3.0/CUDA/10.1.168/lib64"
cuFFT_LIB_DIR = "/usr/local/cuda/lib64"
Eigen_PATH="/home/yage/Topo-DDA/eigen"

FLAG = -fopenmp -lpthread -O3 -Wall -Wno-maybe-uninitialized -fPIC -m64 -I${Eigen_PATH} -I${cuFFT_INCLUDE_PATH} -L${cuFFT_LIB_DIR} -Wl,--no-as-needed -lcuda -lcudart -lcufft -lm -ldl
OBJ = test.o Structure.o Space.o SpacePara.o CoreStructure.o tools.o AProductCore.o DDAModel.o EvoDDAModel.o ObjectiveDDAModel.o SiCi.o kernel.o FOM.o FilterOption.o

te : $(OBJ)
	$(CC) $(FLAG) -o te $(OBJ)
test.o : test.cpp
	$(CC) $(FLAG) -c test.cpp
Structure.o : Structure.cpp
	$(CC) $(FLAG) -c Structure.cpp
Space.o : Space.cpp
	$(CC) $(FLAG) -c Space.cpp
SpacePara.o : SpacePara.cpp
	$(CC) $(FLAG) -c SpacePara.cpp
CoreStructure.o : CoreStructure.cpp
	$(CC) $(FLAG) -c CoreStructure.cpp		
tools.o : tools.cpp
	$(CC) $(FLAG) -c tools.cpp
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
FOM.o : FOM.cpp
	$(CC) $(FLAG) -c FOM.cpp
FilterOption.o : FilterOption.cpp
	$(CC) $(FLAG) -c FilterOption.cpp
kernel.o : kernel.cu
	/usr/local/cuda/bin/nvcc -w -c kernel.cu


.PHONY : clean
clean:
	rm $(OBJ) te

.PHONY : cleandata
cleandata:
	rm EXT_mine.txt time_mine_bicgstab.txt wavelength.txt iteration.txt


