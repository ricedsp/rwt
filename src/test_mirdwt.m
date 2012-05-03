function test_suite = test_mirdwt
initTestSuite;



function test_mirdwt_1
       x = makesig('LinChirp',8);
       h = daubcqf(4,'min');
       L = 2;
       [y,L] = mdwt(x,h,L);
       [x_new,L] = midwt(y,h,L);
       
       
       xin = makesig('Leopold',8);
       h = daubcqf(4,'min');
       Lin = 1;
       [yl,yh,L] = mrdwt(xin,h,Lin);
       [x,L] = mirdwt(yl,yh,h,L);
       
       disp(x);
       disp(xin);

assertEqual(L,Lin);
assertVectorsAlmostEqual(x, xin,'relative',0.0001);

