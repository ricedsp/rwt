function test_mdwt
       x = makesig('LinChirp',8);
       h = daubcqf(4,'min');
       L = 2;
       [y,L] = mdwt(x,h,L);
%
%    1D Example's  output and explanation:
%
       y_corr = [1.1097 0.8767 0.8204 -0.5201 -0.0339 0.1001 0.2201 -0.1401];
       L_corr = 2;

assertVectorsAlmostEqual(y, y_corr,'relative',0.001);
assertEqual(L,L_corr);
end