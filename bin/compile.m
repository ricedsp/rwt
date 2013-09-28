%    COMPILE compiles the c files and generates mex files.
%

mex -v ../mex/mdwt.c ../lib/src/dwt.c ../lib/src/init.c ../lib/src/platform.c -I../lib/inc 
mex -v ../mex/midwt.c ../lib/src/idwt.c ../lib/src/init.c ../lib/src/platform.c -I../lib/inc 
mex -v ../mex/mrdwt.c ../lib/src/rdwt.c ../lib/src/init.c ../lib/src/platform.c -I../lib/inc 
mex -v ../mex/mirdwt.c ../lib/src/irdwt.c ../lib/src/init.c ../lib/src/platform.c -I../lib/inc 
