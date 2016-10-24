OPT= -g
CCFLAGS = ${OPT} -Wall 
LDFLAGS = ${OPT} 

DEPS_DS=dspaces/get_regions.h dspaces/put_regions.h src/divide.h cluster/cluster.h src/get_divs.h src/read_file.h
ODIR=obj
_OBJ_GET = get_regions.o cluster.o
OBJ_GET = $(patsubst %,$(ODIR)/%,$(_OBJ_GET))
_OBJ_PUT = put_regions.o read_file.o divide.o
OBJ_PUT = $(patsubst %,$(ODIR)/%,$(_OBJ_PUT))
RM= rm -rf

H5_ROOT=/home/lifen/tools/hdf5-1.8.17/hdf5
## dataspaces configurations
DS_ROOT=/home/user/software/dataspaces-1.5.0/bin
DS_INC=-I ${DS_ROOT}/include ${H5_ROOT}/include
DS_LIB=-L ${DS_ROOT}/lib -ldspaces -ldscommon -ldart -lrdmacm -libverbs -lm

DS_CC=mpicc                       

$(ODIR)/%.o: %.c $(DEPS_DS)
	    $(DS_CC) -c -o $@ $< $(DS_INC) $(CCFLAGS)


get_regions: $(OBJ_GET) 
	    $(DS_CC) -o get_regions $^ $(DS_LIB) $(LDFLAGS)

    
put_regions: $(OBJ_PUT)
	    $(DS_CC) -o put_regions $^ $(DS_LIB)  $(LDFLAGS)

.PHONY: clean
clean:
	$(RM) cluster/*.o
	$(RM) src/*.o
	$(RM) dspaces/*.o
	$(RM) bin/*
