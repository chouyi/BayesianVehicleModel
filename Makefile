CXX = g++
HOME= /usr/local/include 

LIBS = -lmpfr -lgmp -lgsl -lgslcblas -lm -lglpk 
CFLAGS = -I . -I $(HOME) -g -O3 -std=c++11
LINK_FLAGS = -g -L/usr/local/lib
OBJS = dynamic.o inference.o linpack_d.o blas0.o blas1_d.o

all: UAV_CTmodel 


UAV_CTmodel: $(OBJS) UAV_CTmodel.o 
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)



%.o: %.cc
	$(CXX) -O3 -c $(CFLAGS) -o $@ $<
%.o: %.cpp
	$(CXX) -O3 -c $(CFLAGS) -o $@ $<
%.o: %.c
	$(CXX) -O3 -c $(CFLAGS) -o $@ $<


clean: 
	rm -f *.o *.exe UAV_CTmodel 
