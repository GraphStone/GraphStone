CC = g++
CFLAGS = -O3 -g -fopenmp -std=c++11 -fpermissive
#CFLAGS = -g -std=c++11 -fpermissive
LDFLAGS = -lconfig++  -ljemalloc -lboost_system -lboost_program_options

app: patch_pagerank6 patch_rw patch_cc patch_dij patch_bfs

patch_% : src/patch_%.cpp inc/*pp inc/*h
	@echo $@.cpp
	$(CC) -o bin/$@ $(CFLAGS) src/$@.cpp   $(LDFLAGS)  -I./inc

clean:
	rm -f bin/*
