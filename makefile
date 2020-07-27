
CC = g++
cuFFT_INCLUDE_PATH = "/opt/apps/software/Compiler/GCC/9.3.0/CUDA/11.0.182/include"
cuFFT_LIB_DIR = "/opt/apps/software/Compiler/GCC/9.3.0/CUDA/11.0.182/lib64"
FLAG = -fopenmp -lpthread -O3 -Wall -Wno-maybe-uninitialized -fPIC -m64 -I${cuFFT_INCLUDE_PATH} -L${cuFFT_LIB_DIR} -Wl,--no-as-needed -lcuda -lcudart -lcufft -lm -ldl
OBJ = test.o Structure.o Space.o Model.o EvoModel.o tools.o Objective.o kernel.o

te : $(OBJ)
	$(CC) $(FLAG) -o te $(OBJ)
test.o : test.cpp
	$(CC) $(FLAG) -c test.cpp
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
kernel.o : kernel.cu
	nvcc -w -c kernel.cu


.PHONY : clean
clean:
	rm $(OBJ) te

.PHONY : cleandata
cleandata:
	rm EXT_mine.txt time_mine_bicgstab.txt wavelength.txt iteration.txt


