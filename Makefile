OPT= 
CCFLAGS = ${OPT} -Wall 
LDFLAGS = ${OPT} 

DEPS_DS=dspaces/get_regions.h dspaces/put_regions.h src/divide.h cluster/cluster.h src/get_divs.h src/read_file.h
RM= rm -rf
BIN = bin

H5_ROOT=/opt/hdf5/intel/mvapich2_ib
## dataspaces configurations
DS_ROOT=/home/rlu/Dataspacesroot
DS_INC=-I ${DS_ROOT}/include -I${H5_ROOT}/include -Icluster -Isrc -Idspaces
DS_LIB=-L ${DS_ROOT}/lib -L $(H5_ROOT)/lib -ldspaces -ldscommon -ldart -lrdmacm -libverbs -lm -lpthread

CC=mpicc                       
.PHONY: clean get_regions put_regions all analysis

%.o : %.c $(DEPS_DS)
	    $(CC) -c -o $@ $< $(DS_INC) $(CCFLAGS)


get_regions: dspaces/get_regions.o  src/get_divs.o src/common_utility.o
	    $(CC) -o $(BIN)/get_regions $^ $(DS_LIB) $(LDFLAGS)

    
put_regions: $(OBJ_PUT) dspaces/put_regions.o dspaces/run_with_dspaces.o src/read_file.o src/divide.o  src/common_utility.o
	    h5pcc -o $(BIN)/put_regions $^ $(DS_LIB)  $(LDFLAGS)

analysis: dspaces/analysis.o cluster/cluster.o
	    $(CC) -o $(BIN)/analysis  $^ $(DS_LIB) $(LDFLAGS)

all: get_regions put_regions analysis

clean:
	$(RM) cluster/*.o
	$(RM) src/*.o
	$(RM) dspaces/*.o
	$(RM) bin/*
