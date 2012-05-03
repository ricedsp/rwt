function test_suite = test_midwt
initTestSuite;



function test_midwt_1
       x = makesig('LinChirp',8);
       h = daubcqf(4,'min');
       L = 2;
       [y,L] = mdwt(x,h,L);
       [x_new,L] = midwt(y,h,L);
     %  disp(x);
     %  disp(y);
     %  disp(x_new);

assertVectorsAlmostEqual(x, x_new,'relative',0.0001);

