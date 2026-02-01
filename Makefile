CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra
INCLUDES = -I include

TARGET = raytracer
SOURCES = src/main.cpp

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^

clean:
	rm -f $(TARGET) *.ppm *.png

run: $(TARGET)
	./$(TARGET) --builtin two_spheres -w 800 -h 450 -o two_spheres.ppm

run-all: $(TARGET)
	./$(TARGET) --builtin two_spheres -w 800 -h 450 -o two_spheres.ppm
	./$(TARGET) --builtin three_spheres -w 800 -h 450 -o three_spheres.ppm
	./$(TARGET) --builtin glass -w 800 -h 450 -o glass.ppm
	./$(TARGET) --builtin csg -w 800 -h 450 -o csg.ppm
