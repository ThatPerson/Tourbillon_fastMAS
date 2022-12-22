# choice of compiler
CXX = g++
# compiler flags
CXXFLAGS = -fopenmp -O3
# linker flags
LDFLAGS =  -fopenmp -O3
# cleanup command
RM = /bin/rm
# program executable filename
PROG = ssd
# object files needed by main program
OBJS = main.o density_matrix.o crystal.o hamiltonian.o simulation.o tensor.o

all: $(PROG)


$(PROG): $(OBJS)
	$(CXX) $(OBJS) -o $(PROG) $(LDFLAGS) 


main.o: main.cpp crystal.h hamiltonian.h density_matrix.h constants.h simulation.h tensor.h
	$(CXX) $(CXXFLAGS) -c main.cpp

crystal.o: crystal.cpp crystal.h constants.h tensor.h
	$(CXX) $(CXXFLAGS) -c crystal.cpp

density_matrix.o: density_matrix.cpp density_matrix.h hamiltonian.h constants.h
	$(CXX) $(CXXFLAGS) -c density_matrix.cpp

hamiltonian.o: hamiltonian.cpp hamiltonian.h crystal.h constants.h tensor.h
	$(CXX) $(CXXFLAGS) -c hamiltonian.cpp

simulation.o: simulation.cpp simulation.h constants.h
	$(CXX) $(CXXFLAGS) -c simulation.cpp

tensor.o: tensor.cpp tensor.h constants.h
	$(CXX) $(CXXFLAGS) -c tensor.cpp

clean:
	$(RM) main.o density_matrix.o crystal.o hamiltonian.o simulation.o tensor.o
