
# CONFIG REQUIRED: SET THIS VARIABLE TO THE PATH OF THE INDRI SOURCE CODE
ATIRESRC=../external/atire
# CONFIG REQUIRED: SET THIS TO THE SDSL INSTALL PATH
SDSLPREFIX=../build/external/sdsl-lite/
DIVSUF=../build/external/sdsl-lite/external/libdivsufsort/
NEWTSRC=../include
FASTPFOR=../external/fastpfor/headers/

BASE_DIR=$(ATIRESRC)
include $(ATIRESRC)/GNUmakefile.defns
SHARED=
INCPATH=-I$(ATIRESRC)/source/ -I$(ATIRESRC)/atire/ -I $(SDSLPREFIX)/include/ -I$(NEWTSRC) -I$(FASTPFOR) -I$(DIVSUF)/include -I..
LIBPATH=-L $(SDSLPREFIX)/lib/ -L$(DIVSUF)/lib
LIBS= $(patsubst %, -l%, $(DEPENDENCIES)) -lsdsl -ldivsufsort ../build/libfastpfor_lib.a

all:
	$(CXX) $(CXXFLAGS) -Wno-write-strings -O3 -msse4.2 -std=c++11 build_index.cpp compress_qmx_d4.cpp compress_qmx.cpp -o build_index $(INCPATH) $(LIBPATH) $(LIBS) $(CPPLDFLAGS) $(SOURCES_OBJECTS) $(EXTRA_OBJS) $(LDFLAGS) $(MINUS_D)

clean:
	rm -f build_index *~
