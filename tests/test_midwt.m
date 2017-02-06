function test_suite = test_midwt
initTestSuite;


function test_midwt_1D
       x = makesig('LinChirp',8);
       h = daubcqf(4,'min');
       L = 2;
       [y,L] = mdwt(x,h,L);
       [x_new,L] = midwt(y,h,L);
assertVectorsAlmostEqual(x, x_new,'relative',0.0001);

function test_midwt_1Dc
       x = randn(8,2)*[1;1j];
       h = daubcqf(4,'min');
       L = 2;
       [y,L] = midwt(x,h,L);
       [yr,L] = midwt(real(x),h,L);
       [yi,L] = midwt(imag(x),h,L);
assertVectorsAlmostEqual(y, yr+yi*1j,'relative',0.0001);

function test_midwt_1Ds
       x = single(makesig('LinChirp',8));
       h = single(daubcqf(4,'min'));
       L = 2;
       [y,L] = mdwt(x,h,L);
       [x_new,L] = midwt(y,h,L);
assertVectorsAlmostEqual(x, x_new,'relative',0.0001);

function test_midwt_2D
       load lena512;
       x = lena512;
       h = daubcqf(6);
       [y,L] = mdwt(x,h);
       [x_new,L] = midwt(y,h);
assertEqual(L,9);
assertVectorsAlmostEqual(x, x_new,'relative',0.0001);

function test_midwt_2Ds
       load lena512;
       x = lena512;
       h = daubcqf(6);
       [y,L] = mdwt(x,h);
       [x_new,L] = midwt(single(y),single(h));
assertEqual(L,9);
assertVectorsAlmostEqual(x, x_new,'relative',0.0001);
