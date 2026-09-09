import unittest
import numpy as np
from audit_pppar_same_epoch_joint_nis import analyse

def case(v):
    return {'x':np.zeros((2,1)),'P':np.array([[1.,.999],[.999,1.]]),
            'H':np.eye(2),'v':np.array(v)[:,None],
            'R':np.eye(2)*1e-8,'rhs':np.array(v)[:,None]}, {
            0:{'type':'AMBIGUITY','receiver':'A','satellite':'G01','signal_number':1},
            1:{'type':'AMBIGUITY','receiver':'B','satellite':'G01','signal_number':1}}

class SamePriorTests(unittest.TestCase):
    def test_opposed_marginally_small_residuals_conflict_jointly(self):
        m,k=case([.1,-.1]); r=analyse(m,k)
        self.assertLess(r['subsets']['A']['nis'],.011)
        self.assertGreater(r['full_nis'],19)
        self.assertLess(r['cross_station_covariance_removed_nis'],.021)
        self.assertAlmostEqual(r['full_nis'],r['conditional_on_other_stations']['A']['conditional_nis']+r['subsets']['B']['nis'],places=8)

    def test_compatible_common_residual_is_not_flagged_by_correlation(self):
        m,k=case([.1,.1]); r=analyse(m,k)
        self.assertLess(r['full_nis'],.011)

    def test_rejects_mismatched_constraint_residual(self):
        m,k=case([.1,-.1]); m['rhs'][0,0]=1
        with self.assertRaises(AssertionError): analyse(m,k)

    def test_rejects_asymmetric_covariance(self):
        m,k=case([.1,-.1]); m['P'][0,1]=.5
        with self.assertRaises(AssertionError): analyse(m,k)

if __name__=='__main__': unittest.main()
