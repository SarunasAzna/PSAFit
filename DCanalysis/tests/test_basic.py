import unittest
from ..modules import basic_test as basic
import numpy as np


class TestBasic(unittest.TestCase):

    def setUp(self):
        super(TestBasic, self).setUp()
        self.x = [i for i in range(1, 10)]
        self.data = np.array([
            [i for i in range(1, 10)],
            [i*3.0 for i in range(1, 10)]
        ])

    def test_true(self):
        """ Test True"""
        self.assertTrue(True)

    def test_log_derivative(self):
        """ Test logarithmic derivative output value."""
        res = basic.logDerivative(self.x, self.data)
        self.assertEqual(3.3219280948873631, res[0][0])
