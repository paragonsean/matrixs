# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -I../include -I. -fmessage-length=0 -Wno-sign-compare -O3 -g0
LDFLAGS = -lpthread

# Target executable
TARGET = main

# Object files
OBJS = main.o matrixSolvers.o matrixGenerator.o IterativeSolvers.o

# Default target
all: $(TARGET)

# Linking the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

# Compiling source files into object files
main.o: main.cpp matrixSolvers.h matrixGenerator.h IterativeSolvers.h
	$(CXX) $(CXXFLAGS) -MMD -o main.o -c main.cpp

matrixSolvers.o: matrixSolvers.cpp matrixSolvers.h
	$(CXX) $(CXXFLAGS) -MMD -o matrixSolvers.o -c matrixSolvers.cpp

matrixGenerator.o: matrixGenerator.cpp matrixGenerator.h
	$(CXX) $(CXXFLAGS) -MMD -o matrixGenerator.o -c matrixGenerator.cpp

IterativeSolvers.o: IterativeSolvers.cpp IterativeSolvers.h
	$(CXX) $(CXXFLAGS) -MMD -o IterativeSolvers.o -c IterativeSolvers.cpp

# Clean up the generated files
clean:
	rm -f $(OBJS) $(TARGET) *.d

# Include dependency files
-include $(OBJS:.o=.d)

