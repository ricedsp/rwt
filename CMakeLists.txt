cmake_minimum_required (VERSION 2.6)
project (rwt)
subdirs(lib/src)

# This section based on http://www.cmake.org/pipermail/cmake/2003-June/003953.html
IF (UNIX)
  ADD_CUSTOM_TARGET (distclean @echo cleaning for source distribution)
  SET(DISTCLEANED
   CMakeFiles
   lib/src/CMakeFiles
   cmake.depends
   cmake.check_depends
   CMakeCache.txt
   cmake.check_cache
   *.cmake
   lib/src/*.cmake
   Makefile
   lib/src/Makefile
   core core.*
   gmon.out
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
