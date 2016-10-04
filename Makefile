HDFPATH=/home/lifen/tools/hdf5-1.8.17/hdf5
all: read_file get_divs main
	h5cc -o bin/run  bin/main.o bin/read_file.o bin/get_divs.o -lm

withlg: read_file get_divs main
	h5cc -o bin/run_withlg  bin/main.o bin/read_file.o bin/get_divs.o -lm

main: src/main.c
	gcc src/main.c -c -g -o bin/main.o 

read_file: src/read_file.c
	gcc src/read_file.c -c -g -o bin/read_file.o -I ${HDFPATH}/include
get_divs: src/get_divs.c
	gcc src/get_divs.c -c -g -o bin/get_divs.o 
clustering: src/do_clustering.c
	gcc src/do_clustering.c cluster/cluster.c -g -W -lm -o bin/do_clustering

test_file: src/test_hdf.c
	gcc src/test_hdf.c -c -g  -I ${HDFPATH}/include -o bin/test_hdf.o
test_div: src/test_div.c
	gcc src/test_div.c -o bin/test_div -lm -g

clean:
	rm bin/*.o -rf
