import unittest
import numpy as np
from lattice import Lattice

class TestLattice(unittest.TestCase):
    def setUp(self):
        N = 4
        dim = 4
        self.lattice = Lattice(N=4, lambda_=0.5, N_measurements=100, N_thermalization=0, width=1, HMC=True, epsilon=0, N_steps=0, dim=dim)

    def test_initialization(self):
        self.assertEqual(self.lattice.N, 4)
        self.assertEqual(self.lattice.lambda_, 0.5)
        self.assertEqual(self.lattice.N_measurements, 100)
        self.assertEqual(self.lattice.N_thermalization, 0)
        self.assertEqual(self.lattice.width, 1)
        self.assertTrue(self.lattice.HMCs)
        self.assertEqual(self.lattice.epsilon, 0)
        self.assertEqual(self.lattice.N_steps, 0)
        self.assertEqual(self.lattice.dim, 4)

    def test_action_type(self):
        action = self.lattice.action()
        self.assertIsInstance(action, float)
    
    def test_action_zero(self):
        self.lattice.lattice = np.zeros((4, 4, 4, 4))
        action = self.lattice.action()
        self.assertEqual(action, 0.0)
        random_number = np.random.normal(size=(1))

        self.lattice.lattice = np.ones((4, 4, 4, 4))*random_number
        action = self.lattice.action()
        self.assertAlmostEqual(action, 0.0)
    
    def test_action_inneficient(self):
        N = 4
        self.lattice.lattice = np.random.normal(size=(N, N,N,N))
        self.lattice.dim = 4
        #self.lattice.lattice = np.ones((4, 4, 4, 4))
        action = self.lattice.action()
        action_true = 0
        squarable = np.zeros((N, N, N, N))
        
      
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    for l in range(N):
                        direction_x = self.lattice.lattice[(i+1)%N,j,k,l]+self.lattice.lattice[(i-1)%N,j,k,l]-2*self.lattice.lattice[i, j, k, l] + (self.lattice.lattice[(i+1)%N,j,k,l]-self.lattice.lattice[i,j,k,l])**2
                        direction_y = self.lattice.lattice[i,(j+1)%N,k,l]+self.lattice.lattice[i,(j-1)%N,k,l]-2*self.lattice.lattice[i, j, k, l] + (self.lattice.lattice[i,(j+1)%N,k,l]-self.lattice.lattice[i,j,k,l])**2
                        direction_z = self.lattice.lattice[i,j,(k+1)%N,l]+self.lattice.lattice[i,j,(k-1)%N,l]-2*self.lattice.lattice[i, j, k, l] + (self.lattice.lattice[i,j,(k+1)%N,l]-self.lattice.lattice[i,j,k,l])**2
                        direction_t = self.lattice.lattice[i,j,k,(l+1)%N]+self.lattice.lattice[i,j,k,(l-1)%N]-2*self.lattice.lattice[i, j, k, l] + (self.lattice.lattice[i,j,k,(l+1)%N]-self.lattice.lattice[i,j,k,l])**2
                        squarable[i,j,k,l] = direction_x + direction_y + direction_z + direction_t
                        action_true += 1/(2*np.abs(self.lattice.lambda_))*squarable[i,j,k,l]**2
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    for l in range(N):
                        self.assertAlmostEqual(squarable[i,j,k,l],self.lattice.I(self.lattice.lattice)[i,j,k,l])
    
        self.assertAlmostEqual(action, action_true)
    def test_action_known_values(self):
        self.lattice.dim = 2
        self.lattice.N = 2
        self.lattice.lattice = np.array([[1,2],[3,4]])
        action = self.lattice.action()
        self.assertEqual(action, 180)
        I = self.lattice.I(self.lattice.lattice)
        Result = None 
        self.assertEqual(I, Result)


if __name__ == '__main__':
    unittest.main()