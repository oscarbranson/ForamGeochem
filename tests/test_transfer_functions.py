import unittest
from foramgeochem import transfer as tfr

class test_transfer_functions(unittest.TestCase):

    def test_exponential(self):
        temp = tfr.MgCa.exp_mgca_2_temp(5, 0.3, 0.09)
        mgca = tfr.MgCa.exp_temp_2_mgca(31.26011907511152, 0.3, 0.09)

        self.assertAlmostEqual(temp, 31.26011907511152)
        self.assertAlmostEqual(mgca, 5.0)

    def test_linear(self):
        temp = tfr.MgCa.lin_mgca_2_temp(5, 0.3, 0.09)
        mgca = tfr.MgCa.lin_temp_2_mgca(1.59, 0.3, 0.09)

        self.assertAlmostEqual(temp, 1.59)
        self.assertAlmostEqual(mgca, 5.0)

    def test_evans2012(self):
        mgca_modern = 5.0
        mgca_ancient = 2.12
        A, B, H = 0.09, 0.3, 0.53

        temp = tfr.MgCa.evans2012_mgca_2_temp(5, A, B, H, mgca_modern, mgca_ancient)
        mgca = tfr.MgCa.evans2012_temp_2_mgca(36.31291425941813, A, B, H, mgca_modern, mgca_ancient)

        mgca_sw_calc = tfr.MgCa.evans2012_mgca_2_mgca_sw(5, 36.31291425941813, A, B, H, mgca_modern)

        self.assertAlmostEqual(mgca, 5)
        self.assertAlmostEqual(temp, 36.31291425941813)
        self.assertAlmostEqual(mgca_sw_calc, mgca_ancient)
    
    def test_evans2015(self):
        temp = tfr.MgCa.evans2015_mgca_2_temp(126.2)
        mgca = tfr.MgCa.evans2015_temp_2_mgca(17.02967505036528)

        self.assertAlmostEqual(temp, 17.02967505036528)
        self.assertAlmostEqual(mgca, 126.2)

if __name__ == '__main__':
    unittest.main()