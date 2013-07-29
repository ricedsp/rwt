function test_suite = test_mrdwt
initTestSuite;

function test_mrdwt_1
    x = makesig('Leopold',8);
    h = daubcqf(4,'min');
    L = 1;
    [yl,yh,L] = mrdwt(x,h,L);
    yl_corr = [0.8365  0.4830 0 0 0 0 -0.1294 0.2241];
    yh_corr = [-0.2241 -0.1294 0 0 0 0 -0.4830 0.8365];
    L_corr = 1;
       
assertVectorsAlmostEqual(yl, yl_corr,'relative',0.001);
assertVectorsAlmostEqual(yh, yh_corr,'relative',0.001);
assertEqual(L,L_corr);
