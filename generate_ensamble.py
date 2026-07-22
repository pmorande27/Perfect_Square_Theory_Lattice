from lattice import Lattice
import LatticeRun
import numpy as np
import plotting
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import processing
import scipy.optimize as optimize


def main():
    N = 10
    beta = 5
    N_th = 1000
    N_measure = 200

    #LatticeRun.calibration(N,beta,0.25,True,300,calibration_runs = 100)
    LatticeRun.generate_phis(N,beta,N_measure,True,4,guess = 200)


if __name__ == "__main__":
    main()