"""Unit tests"""
import unittest

import pandas as pd
import numpy as np
import os
from uncertainties import ufloat,unumpy
from uncertainties.umath import *
import sys
sys.path.insert(0, './')
from flash_helper.flash_functions import mean_array_calculator, find_nearest, get_index


class TestFLASH(unittest.TestCase):  # pylint:disable=R0902
    """Unit test class"""

    @classmethod
    def setUpClass(cls):  # pylint: disable=R0915
        cls.x=np.array([1,0,0,0])
        cls.y=np.array([0,0,0,1])
        cls.x_find=np.array([2,3])
        cls.y_find=np.array([1,2,3,4])



    @classmethod
    def tearDownClass(cls):
        return

    def SetUp(self):
        return

    def TearDown(self):
        return


    def test_find_nearest(self):

        idx = find_nearest(self.x,1)
        self.assertTrue(idx==0)

        with self.assertRaises(Exception):
            find_nearest(self.x,"string")


    def test_get_index(self):

        idx = get_index(self.x_find,self.y_find)
        self.assertTrue(idx==[1,2])

    def test_mean_array_calculator(self):

        tester = mean_array_calculator(self.x,self.y,self.y_find)

        tester_mean = unumpy.nominal_values(tester)

        tester_sigma = unumpy.std_devs(tester)

        self.assertTrue(tester_mean[0]==2/3)

        self.assertTrue(tester_mean[1]==2/3)

        self.assertTrue(tester_mean[2]==1)

        self.assertTrue(tester_mean[3]==5/3)


        self.assertTrue(tester_sigma[0]==np.sqrt(2/9))

        self.assertTrue(tester_sigma[1]==np.sqrt(8/9))

        self.assertTrue(tester_sigma[2]==np.sqrt(2))

        self.assertTrue(tester_sigma[3]==np.sqrt(78/27))

        with self.assertRaises(Exception):
            mean_array_calculator(self.x,self.x_find)



if __name__ == "__main__":
    unittest.main()
