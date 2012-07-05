# Makefile for building library and C++ test files for Rice Wavelet Toolbox.
# 
# Robert Brockman II
# See README file in this directory.

# Executable files
ECLIPSEBUILD := RwtLib RwtTest
CLEAN_ECLIPSEBUILD := $(addprefix Clean,$(ECLIPSEBUILD))

# The default target: run all examples
all: $(ECLIPSEBUILD)

# Library for RWT
RwtLib: RwtLib/Release/libRwtLib.a 
RwtLib/Release/libRwtLib.a:
	cd RwtLib/Release; make
CleanRwtLib:
	cd RwtLib/Release; make clean

# C++ Tests for RWT
RwtTest: RwtLib/Debug/libRwtLib.a "RwtTest/Debug Gcov/RwtTest" 
RwtLib/Debug/libRwtLib.a:
	cd RwtLib/Debug; make
"RwtTest/Debug Gcov/RwtTest":
	cd RwtTest/Debug\ Gcov; make
CleanRwtTest:
	cd RwtLib/Debug; make clean
	cd "RwtTest/Debug Gcov"; make clean

clean: $(CLEAN_ECLIPSEBUILD)

distclean: clean
	
