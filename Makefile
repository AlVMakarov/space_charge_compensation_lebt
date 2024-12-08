CC = g++
LIBS = `pkg-config --libs ibsimu-1.0.6dev`
LDFLAGS = -Wall -g
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6dev`

all: simu analysis

simu: simu.o
	$(CC) $(LDFLAGS) $^ $(LIBS) -o $@
analysis: analysis.o
	$(CC) $(LDFLAGS) $^ $(LIBS) -o $@

.cpp.o:
	$(CC) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) *.png *.o
