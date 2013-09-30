%    COMPILE compiles the c files and generates mex files.
%

x = computer();
if (x(length(x)-1:length(x)) == '64')
  mex -v -largeArrayDims ../mex/mdwt.c ../lib/src/dwt.c ../lib/src/init.c ../lib/src/platform.c -I../lib/inc 
  mex -v -largeArrayDims ../mex/midwt.c ../lib/src/idwt.c ../lib/src/init.c ../lib/src/platform.c -I../lib/inc 
  mex -v -largeArrayDims ../mex/mrdwt.c ../lib/src/rdwt.c ../lib/src/init.c ../lib/src/platform.c -I../lib/inc 
  mex -v -largeArrayDims ../mex/mirdwt.c ../lib/src/irdwt.c ../lib/src/init.c ../lib/src/platform.c -I../lib/inc 
else
  mex -v -compatibleArrayDims ../mex/mdwt.c ../lib/src/dwt.c ../lib/src/init.c ../lib/src/platform.c -I../lib/inc 
  mex -v -compatibleArrayDims ../mex/midwt.c ../lib/src/idwt.c ../lib/src/init.c ../lib/src/platform.c -I../lib/inc 
  mex -v -compatibleArrayDims ../mex/mrdwt.c ../lib/src/rdwt.c ../lib/src/init.c ../lib/src/platform.c -I../lib/inc 
  mex -v -compatibleArrayDims ../mex/mirdwt.c ../lib/src/irdwt.c ../lib/src/init.c ../lib/src/platform.c -I../lib/inc 
end


