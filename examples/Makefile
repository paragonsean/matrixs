# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -I../include -I. -fmessage-length=0 -Wno-sign-compare -O3 -g0
LDFLAGS = -lpthread

# Target executable
TARGET = main

# Object files
OBJS = main.o iterative_solvers.o

# Default target
all: $(TARGET)

# Linking the executable
$(TARGET): $(OBJS)
	@echo "Linking $(TARGET)..."
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

# Compiling main.cpp, includes template headers directly
main.o: main.cpp matrix_file_handler.h matrix_solvers.h
	@echo "Compiling main.cpp..."
	$(CXX) $(CXXFLAGS) -MMD -c main.cpp -o main.o

# Compiling iterative_solvers.cpp
iterative_solvers.o: iterative_solvers.cpp iterative_solvers.h
	@echo "Compiling iterative_solvers.cpp..."
	$(CXX) $(CXXFLAGS) -MMD -c iterative_solvers.cpp -o iterative_solvers.o

# Clean up the generated files
clean:
	@echo "Cleaning up..."
	rm -f $(OBJS) $(TARGET) *.d

# Include dependency files
-include $(OBJS:.o=.d)
