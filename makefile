ifeq ($(LHAPDF),)
   LHAPDF = $(HOME)/local/local
endif 

LHAPDFINC = $(shell $(LHAPDF)/bin/lhapdf-config --incdir)
LHAPDFLIB = $(shell $(LHAPDF)/bin/lhapdf-config --libdir)


HATHORPATH = .

CC  = gcc
CXX = g++
FC  = gfortran
AR  = ar
RANLIB = ranlib

IFLAGS = -I. -I$(LHAPDFINC) -I$(HATHORPATH)/include -I$(shell root-config --incdir)
MYLIBS =  -L $(HATHORPATH)/lib -lHathor -L $(LHAPDFLIB) -lLHAPDF -lff $(shell root-config --libs)

#IFLAGS = -I. -I$(LHAPDFINC) -I$(HATHORPATH)/include
#MYLIBS =  -L $(HATHORPATH)/lib -lHathor -L $(LHAPDFLIB) -lLHAPDF -lff 

LFLAGS := $(MYLIBS) $(LFLAGS) -lgfortranbegin -lgfortran -lm

# default configuration
CFLAGS := $(CFLAGS) -O2 -Wall
FFLAGS := $(FFLAGS)

DEMOS =  Tunegt

all: $(DEMOS)

%: %.cxx
	$(CXX) $(CFLAGS) $(IFLAGS) -o $@ $< $(LFLAGS)

clean:
	rm -f $(DEMOS) 

distclean: clean



