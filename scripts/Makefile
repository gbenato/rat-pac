CFLAGS = -Wall -Wno-write-strings -DSWAP_BYTES -g -fdiagnostics-show-option
LFLAGS = -L/usr/lib -L$(ROOTSYS)/lib -L$(RATROOT)/lib
CC= g++

#LOCAL
INCLUDES= `root-config --cflags --glibs` -lGeom\
-I$(RATROOT)/include \
-I$(RATROOT)/include/RAT \
-I$(RATROOT)/include/RAT/DS \
-lRATEvent\
-lhdf5_cpp\
-lhdf5

.PHONY: all

#Tell make to make one .out file for each .cpp file found in the current directory
all: $(patsubst %.cpp, %.exe, $(wildcard *.cpp))

#Rule how to create arbitary .exe files.
#First state what is needed for them e.g. additional headers, .cpp files in an include folder...
#Then the command to create the .exe file, probably you want to add further options to the g++ call.
%.exe: %.cpp Makefile
	$(CC) $(CFLAGS) $(LFLAGS) $< -o $@ $(INCLUDES) -std=c++0x

clean:
	rm -v *.exe
	rm -vrf *.dSYM
