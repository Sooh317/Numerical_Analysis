MAIN = main
OBJS = load.o stencil5.o

CXX = g++
CXXFLAGS = -O3

all: $(MAIN)

$(MAIN): $(MAIN).o $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)


.SUFFIXES: .c .cpp .o
.c.o:
	$(CXX) $(CXXFLAGS) -c $<

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: clean
clean:
	rm -f *.bmp *.o $(MAIN)
