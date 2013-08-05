function test_suite = test_mrdwt
initTestSuite;

function test_mrdwt_1
  x = makesig('Leopold',8);
  h = daubcqf(4,'min');
  L = 1;
  [yl, yh, L] = mrdwt(x, h, L);
  yl_corr = [0.8365  0.4830 0 0 0 0 -0.1294 0.2241];
  yh_corr = [-0.2241 -0.1294 0 0 0 0 -0.4830 0.8365];
  L_corr = 1;
assertVectorsAlmostEqual(yl, yl_corr,'relative',0.001);
assertVectorsAlmostEqual(yh, yh_corr,'relative',0.001);
assertEqual(L,L_corr);

function test_mrdwt_2
  x = makesig('Doppler', 8);
  x2 = [x*.5;x*1.5;x*2.7;x*.4];
  h = daubcqf(4,'min');
  L = 2;
  [yl, yh, L] = mrdwt(x2, h, L);
  yl_slice =  yl(1,:);
  yl_slice_corr = [-0.3749 -0.1749 0.3773 0.5764 0.7591 0.5591 0.0069 -0.1922];
  yh_slice = yh(1,1:12);
  yh_slice_corr = [0.0914 0.0587 -0.4528 -0.3521 0.5165 0.6673 0.1996 -0.0192 0.4201 -0.3404 -0.5793 0.5879];
assertVectorsAlmostEqual(yl_slice, yl_slice_corr,'relative',0.001);
assertVectorsAlmostEqual(yh_slice, yh_slice_corr,'relative',0.001);
