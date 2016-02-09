CFLAGS = -Wall -Wno-write-strings -DSWAP_BYTES -g -fdiagnostics-show-option -std=c++0x -D__WAVEFORMS_IN_DS__
LFLAGS = -L/usr/lib -L$(ROOTSYS)/lib -L$(RATROOT)/lib
CC= g++

INCLUDES= `root-config --cflags --glibs` -lGeom\
-I$(RATROOT)/include \
-I$(RATROOT)/include/RAT \
-I$(RATROOT)/include/RAT/DS \
-lRATEvent

INCLUDES2= `root-config --cflags` \
-I$(RATROOT)/include \
-I$(RATROOT)/include/RAT \
-I$(RATROOT)/include/RAT/DS \
-lRATEvent

SRCS= EventGeometry.cc EventDisplay.cc EventDisplayMain.cc
OBJS=$(subst .cc,.o,$(SRCS))

all: EventDisplayMain.exe

EventDisplayMain.exe: $(OBJS)
	$(CC) $(CFLAGS) $(LFLAGS) -o $@ $(INCLUDES) $(OBJS)

EventDisplayMain.o: EventDisplayMain.cc
	$(CC) $(CFLAGS) -c EventDisplayMain.cc $(INCLUDES2)

EventDisplay.o: EventDisplay.cc EventDisplay.hh
	$(CC) $(CFLAGS) -c EventDisplay.cc $(INCLUDES2)

EventGeometry.o: EventGeometry.cc EventGeometry.hh
	$(CC) $(CFLAGS) -c EventGeometry.cc $(INCLUDES2)

clean:
	rm $(OBJS)
	rm -v *.exe
