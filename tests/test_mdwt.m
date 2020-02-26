function test_suite = test_mdwt
initTestSuite;

function test_mdwt_1D
  x = makesig('LinChirp', 8);
  h = daubcqf(4, 'min');
  L = 2;  % For 8 values in x we would normally be L=2 
  [y, L] = mdwt(x, h, L);
  y_corr = [1.1097 0.8767 0.8204 -0.5201 -0.0339 0.1001 0.2201 -0.1401];
  L_corr = 2;
assertVectorsAlmostEqual(y, y_corr, 'relative', 0.001);
assertEqual(L, L_corr);

function test_mdwt_1Dcpx
  x = randn(8,2)*[1;1j];
  h = daubcqf(4, 'min');
  L = 2;  % For 8 values in x we would normally be L=2 
  [y, L] = mdwt(x, h, L);
  yr = mdwt(real(x), h, L);
  yi = mdwt(imag(x), h, L);
  L_corr = 2;
assertVectorsAlmostEqual(y, yr+1j*yi, 'relative', 0.001);
assertEqual(L, L_corr);

function test_mdwt_1Ds
  x = single(makesig('LinChirp', 8));
  h = single(daubcqf(4, 'min'));
  L = 2;  % For 8 values in x we would normally be L=2 
  [y, L] = mdwt(x, h, L);
  y_corr = [1.1097 0.8767 0.8204 -0.5201 -0.0339 0.1001 0.2201 -0.1401];
  L_corr = 2;
assertVectorsAlmostEqual(y, y_corr, 'relative', 0.001);
assertEqual(L, L_corr);

function test_mdwt_2D
  x = [1 2 3 4; 5 6 7 8 ; 9 10 11 12; 13 14 15 16];
  h = daubcqf(4);
  y = mdwt(x, h);
  y_corr = [34.0000 -3.4641 0.0000 -2.0000; -13.8564 0.0000 0.0000 -2.0000; -0.0000 0.0000 -0.0000 -0.0000; -8.0000 -8.0000 0.0000 -0.0000];
assertVectorsAlmostEqual(y, y_corr, 'relative', 0.001);

function test_mdwt_compute_L1
  x = [1 2];
  h = daubcqf(4, 'min');
  [y, L] = mdwt(x, h);
assertEqual(L, 1);

function test_mdwt_compute_L2
  x = [1 2 3 4];
  h = daubcqf(4, 'min');
  [y, L] = mdwt(x, h);
assertEqual(L, 2);

function test_mdwt_compute_L3
  x = [1 2 3 4 5 6 7 8];
  h = daubcqf(4, 'min');
  [y, L] = mdwt(x, h);
assertEqual(L, 3);

function test_tensor_mdwt_1D
  x = randn(8,3,2);
  h = daubcqf(4);
  y1 = mdwt(x, h,[],1);
  y2 = nan(size(x));
  for i3=1:size(x,3)
    for i4=1:size(x,4)
      y2(:,:,i3,i4) = mdwt( squeeze(x(:,:,i3,i4)),h,[],1);
    end
  end
assertVectorsAlmostEqual(y1, y2, 'relative', 0.001);

function test_tensor_mdwt_2D
  x = randn(8,16,3,2);
  h = daubcqf(4);
  y1 = mdwt(x, h);
  y2 = nan(size(x));
  for i3=1:size(x,3)
    for i4=1:size(x,4)
      y2(:,:,i3,i4) = mdwt( squeeze(x(:,:,i3,i4)),h);
    end
  end
assertVectorsAlmostEqual(y1, y2, 'relative', 0.001);


function test_mdwt_compute_bad_L
  L = -1;
  x = [1 2 3 4 5 6 7 8 9];
  h = daubcqf(4, 'min');
  mdwtHandle = @() mdwt(x, h);
assertExceptionThrown(mdwtHandle, '');

function test_mdwt_empty_input
  mdwtHandle = @() mdwt([], [0 0 0 0]);
assertExceptionThrown(mdwtHandle, '');
