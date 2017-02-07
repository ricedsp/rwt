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

function test_midwt_mat1D
       x = randn(16,4);
       h = daubcqf(4,'min');
       y1 = mdwt(x,h,[],1);
       y2 = nan(size(y1));
       for i2=1:size(x,2)
           y2(:,i2) = mdwt(squeeze(x(:,i2)),h,[],1);
       end
assertVectorsAlmostEqual(y1, y2,'relative',0.0001);

function test_midwt_tensor2D
       x = randn(16,4,3,2);
       h = daubcqf(4,'min');
       y1 = mdwt(x,h);  % this should default to a 2d transform across the first 2 dimensions
       y2 = nan(size(y1));
       for i3=1:size(x,3)
           for i4=1:size(x,4)
               y2(:,:,i3,i4) = mdwt(squeeze(x(:,:,i3,i4)),h,[],2);
           end
       end
assertVectorsAlmostEqual(y1, y2,'relative',0.0001);

function test_midwt_tensor1Dc
       sz = [16,5,3,2];
       x = randn(sz) + 1j*randn(sz);
       h = daubcqf(4,'min');
       L = 2;
       y = mdwt(x,h,L,1);
       [x2,L] = midwt(y,h,L,1);
assertVectorsAlmostEqual(x, x2,'relative',0.0001);
