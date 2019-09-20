GENCODE_SM35     := #-gencode arch=compute_35,code=sm_35
GENCODE_FLAGS    := #$(#GENCODE_SM35)

LDFLAGS   := -L/usr/local/cuda/lib64 -lcudart -lcudadevrt -L./lib -lm -lcxsparse -lpthread -lcusolver -lcusparse  -lcublas_static -lcusparse_static -lculibos -ldl
CCFLAGS   := -m64 -std=c++17

#NVCCFLAGS := -m64 -dc -L/usr/local/cuda/lib64 -lcudart -lcudadevrt -L./lib
NVCCFLAGS := -arch=compute_61 -m64 -std=c++11 -L/usr/local/cuda/lib64 -lcudart -lcudadevrt -lcusolver -lcusparse -lcublas_static -lcusparse_static -lculibos

NVCC = nvcc
CXX  = g++
RM   = rm -rf
CP   = cp -r

################################################################################
# Setup debug and release variants
 RELDIR = Release
 DEBUGDIR = Debug

ifeq "$(findstring release, $(MAKECMDGOALS))" ""
 BUILD=$(DEBUGDIR)
 CCFLAGS += -g -O0 -D_DEBUG='1'  -Wall
 NVCCFLAGS += -g -O0 -D_DEBUG='1' 
 else
 BUILD=$(RELDIR)
 CCFLAGS += -O3 -D_RELEASE='1'
 NVCCFLAGS += -O3 -D_RELEASE='1'
 endif


# Common includes and paths for CUDA
INCLUDES      :=  -I. -I.. -I./src -I./include -I/usr/local/cuda/include/ 

# Additional parameters
#MAXRREGCOUNT  :=  -po maxrregcount=16

OBJ = obj
SRC = src


## List of objs
# nvidia OBJS
OBJ_CU_1  = Sparse/Sparse.cu SparseMatrix/Convert2CSR.cu\
            StiffnessMatrix/StiffnessMatrix.cu\
            StiffnessMatrix/StiffnessMatrixFirstOrder/StiffnessMatrixFirstOrder.cu\
            StiffnessMatrix/StiffnessMatrixGPU/StiffnessMatrixGPU.cu\
			#Assembly/Assembly.cu\

# CPU OBJS
OBJ_CPP_1 = main.cpp Timer/Timer.cpp Log/Log.cpp Application/Application.cpp\
            Geometry/Geometry.cpp Material/Material.cpp Material/ElasticPlainStrain.cpp Material/ElasticPlainStress.cpp Material/Concrete01.cpp\
            StiffnessMatrix/StiffnessMatrixSingleCPU/StiffnessMatrixSingleCPU.cpp\
            StiffnessMatrix/StiffnessMatrixParallelCPU/StiffnessMatrixParallelCPU.cpp\
            Recorder/Recorder.cpp \
			SparseMatrix/SparseMatrix.cpp\
			Solver/Cholesky.cpp Solver/ConjugateGradient.cpp Solver/ConjugateGradientPrecond.cpp Solver/Solver.cpp \
			#Assembly/AssemblyParCpu.cpp Assembly/AssemblySingleCpu.cpp\

OBJ_CPP_2 = $(patsubst %,$(SRC)/%,$(OBJ_CPP_1))
OBJ_CU_2  = $(patsubst %,$(SRC)/%,$(OBJ_CU_1))
OBJ_CPP   = $(patsubst %,%.o,$(OBJ_CPP_2))
OBJ_CU    = $(patsubst %,%.o,$(OBJ_CU_2))

DEL = $(patsubst %,%~,$(OBJ_CPP_2)) $(patsubst %,%~,$(OBJ_CU_2)) $(patsubst %.cpp,%.h~,$(OBJ_CPP_2)) $(patsubst %.cu,%.h~,$(OBJ_CU_2))

APP_NAME = finite
# Target rules
debug : all
release : all

all: $(APP_NAME)

$(APP_NAME): $(OBJ_CU) $(OBJ)/link.o $(OBJ_CPP) 
	$(CXX) $+ -o $@ $(LDFLAGS) $(CCFLAGS) $(EXTRA_LDFLAGS)
	$(RM) obj/link.o

$(OBJ)/link.o: 
	$(NVCC) $(NVCCFLAGS) -dlink  $(OBJ_CU) -o $@ 

$(SRC)/%.cu.o: $(SRC)/%.cu $(SRC)/%.h
	$(NVCC) $(NVCCFLAGS) --device-c $< -o $@ $(INCLUDES)    

$(SRC)/%.cpp.o: $(SRC)/%.cpp $(SRC)/%.h
	$(CXX) $(CCFLAGS) $(INCLUDES) -o $@ -c $<

clean:
	$(RM) $(OBJ_CPP) $(OBJ_CU) $(OBJ)/link.o 
	$(RM) $(DEL)
