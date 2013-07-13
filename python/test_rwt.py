#!/usr/bin/python

import unittest
from numpy import *
from scipy.io import loadmat
from rwt import *

class TestRWT(unittest.TestCase):

  def setUp(self):
      pass

  def test_dwt(self):
    x = makesig('LinChirp', 8)
    h = daubcqf(4, 'min')[0]
    L = 2
    y, L = dwt(x, h, L)
    y_corr = array([1.1097, 0.8767, 0.8204, -0.5201, -0.0339, 0.1001, 0.2201, -0.1401])
    self.assertTrue(allclose(y, y_corr, 0.00082))

  def test_dwt_2d(self):
    x = array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [3, 14, 15, 16.0]])
    h = daubcqf(4)[0]
    L = 2
    y, L = dwt(x, h, L)
    y_corr = array([[34.0000, -3.4641, 0.0000, -2.0000], [-13.8564, 0.0000, 0.0000, -2.0000], [-0.0000, 0.0000, -0.0000, -0.0000], [-8.0000, -8.0000, 0.0000, -0.0000]])
    self.assertTrue(allclose(y, y_corr, 0.1))

  def test_idwt(self):
    x = makesig('LinChirp', 8)
    h = daubcqf(4, 'min')[0]
    L = 2
    y, L = dwt(x, h, L)
    x_new, L = idwt(y, h, L)
    self.assertTrue(allclose(x, x_new, 0.0005))

  def test_idwt_2d(self):
    x = loadmat('../tests/lena512.mat')['lena512'] * 1.0
    h = daubcqf(6)[0]
    L = 9
    y, L = dwt(x, h, L)
    x_new, L = idwt(y, h, L)
    self.assertTrue(allclose(x, x_new, 0.0005))

  def test_makesig_heavisine(self):
    x = makesig('HeaviSine', 8)
    y = array([4.0000, 0.0000, -6.0000, -2.0000, 2.0000, 0.0000, -4.0000, -0.0000])
    self.assertTrue(allclose(x, y, 0.0005))
  
  def test_makesig_bumps(self):
    x = around(makesig('Bumps', 8), 4)
    y = array([0.3206, 5.0527, 0.3727, 0.0129, 0.0295, 0.0489, 0.0004, 0.0000])
    self.assertTrue(allclose(x, y, 0.0005))
  
  def test_makesig_blocks(self):
    x = makesig('Blocks', 8)
    y = array([4.0000, 0.5000, 3.0000, 0.9000, 0.9000, 5.2000, -0.0000, -0.0000])
    self.assertTrue(allclose(x, y, 0.0005))
  
  def test_makesig_doppler(self):
    x = makesig('Doppler', 12)
    y = array([-0.1954, -0.3067, 0.0000, -0.4703, 0.4930, -0.2703, -0.4127, 0.1025, 0.4001, 0.3454, 0.1425, 0])
    self.assertTrue(allclose(x, y, 0.0005))
  
  def test_makesig_ramp(self):
    x = makesig('Ramp', 8)
    y = array([0.1250, 0.2500, -0.6250, -0.5000, -0.3750, -0.2500, -0.1250, 0])
    self.assertTrue(allclose(x, y, 0.0005))
  
  def test_makesig_cusp(self):
    x = makesig('Cusp', 8)
    y = array([0.4950, 0.3464, 0.0707, 0.3606, 0.5050, 0.6164, 0.7106, 0.7937])
    self.assertTrue(allclose(x, y, 0.0005))
  
  def test_makesig_sing(self):
    x = makesig('Sing', 8)
    y = array([5.3333, 16.0000, 16.0000, 5.3333, 3.2000, 2.2857, 1.7778, 1.4545])
    self.assertTrue(allclose(x, y, 0.0005))
  
  def test_makesig_hisine(self):
    x = makesig('HiSine', 8)
    y = array([0.8267, -0.9302, 0.2200, 0.6827, -0.9882, 0.4292, 0.5053, -0.9977])
    self.assertTrue(allclose(x, y, 0.0005))
  
  def test_makesig_losine(self):
    x = makesig('LoSine', 8)
    y = array([0.8660, 0.8661, 0.0003, -0.8658, -0.8663, -0.0006, 0.8657, 0.8664])
    self.assertTrue(allclose(x, y, 0.0472))
  
  def test_makesig_linchirp(self):
    x = makesig('LinChirp', 8)
    y = array([0.0491, 0.1951, 0.4276, 0.7071, 0.9415, 0.9808, 0.6716, 0.0000])
    self.assertTrue(allclose(x, y, 0.0007))
  
  def test_makesig_twochirp(self):
    x = makesig('TwoChirp', 8)
    y = array([0.5132, 1.5000, 0.5412, 0.8660, -0.5132, 0, 0.5132, 0.8660])
    self.assertTrue(allclose(x, y, 0.0005))
  
  def test_makesig_quadchirp(self):
    x = makesig('QuadChirp', 8)
    y = array([0.0164, 0.1305, 0.4276, 0.8660, 0.8895, -0.3827, -0.6217, 0.8660])
    self.assertTrue(allclose(x, y, 0.0024))
  
  def test_makesig_mishmash(self):
    x = makesig('MishMash', 8)
    y = array([0.8922, -0.6046, 1.0751, 2.2558, 0.8429, 1.0273, 0.5551, -0.1317])
    self.assertTrue(allclose(x, y, 0.0005))
  
  def test_makesig_wernersorrows(self):
    x = makesig('WernerSorrows', 8)
    y = array([1.5545, 5.3175, 0.8252, 1.6956, -1.2678, 0.6466, 1.7332, -0.9977])
    self.assertTrue(allclose(x, y, 0.0005))
  
  def test_makesig_leopold(self):
    x = makesig('Leopold', 8)
    y = array([0, 1, 0, 0, 0, 0, 0, 0])
    self.assertTrue(allclose(x, y, 0.0005))
    
if __name__ == '__main__':
    unittest.main()
