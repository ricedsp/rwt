cmake_minimum_required (VERSION 2.6)
project (rwt)
subdirs(lib/src)
subdirs(doc)

# This section based on http://www.cmake.org/pipermail/cmake/2003-June/003953.html
IF (UNIX)
  ADD_CUSTOM_TARGET (distclean @echo cleaning for source distribution)
  SET(DISTCLEANED
   CMakeFiles
   cmake.depends
   cmake.check_depends
   CMakeCache.txt
   cmake.check_cache
   Makefile
   *.cmake
   lib/src/*.cmake
   lib/src/Makefile
   lib/src/CMakeFiles
   lib/src/*.a
   doc/CMakeFiles
   doc/Doxyfile
   doc/Makefile
   doc/*.cmake
   doc/html
   doc/latex
   core core.*
   gmon.out
   bin/*.mexa64
   bin/*.mexmaci64
   src/*.mexa64
   src/*.mexmaci64
   tests/*.mexa64
   tests/*.mexmaci64
   *~
  )
  
  ADD_CUSTOM_COMMAND(
    DEPENDS clean
    COMMENT "distribution clean"
    COMMAND rm
    ARGS    -Rf CMakeTmp ${DISTCLEANED}
    TARGET  distclean
  )
ENDIF(UNIX)


