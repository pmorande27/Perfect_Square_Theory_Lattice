from lattice import Lattice
import LatticeRun
import numpy as np
import plotting
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import processing
import scipy.optimize as optimize




def plotNmeasurements_dependence_tderiv():
    fig,ax = plt.subplots()
    N_measure = 10000


    for i,N_measure in enumerate([1000,10000]):
        beta =1
        N = 10
        results,errs =LatticeRun.evaluate_observable_1D_dontsave(N,beta,N_measure,lambda l: Lattice.measure_time_derivative_correlator_zero_momentum_operator(4,l))


        
      

        results = results
        print(results)
        forward = np.roll(results,-1)
        popt,pcov = optimize.curve_fit(model, range(len(results)-1), forward[:-1], p0=(1, 0.1), bounds=([0, 0], [np.inf, np.inf]))
        print(f'Optimal parameters for N={N}: {popt}')
        xs = np.linspace(0, len(results)-1, 100)
        ys = model(xs, *popt)
        #ax.plot(xs, ys, label=f'Fit for $N$ = {N}')


        

        
        
        ax.errorbar(np.array(range(len(results)))+i*0.01,results,yerr = errs,fmt = 's',label = f'$N_m$ = {N_measure}')
        #ax.errorbar(np.array(range(len(check_t)))+0.01,check_t,yerr = check_t_err,fmt = 's',label = f'Check $N_m$ = {N_measure}')
        
        
        ax.set_xlabel('t')
        operator_label = '$\langle(\\nabla_t\chi_t)(t)(\\nabla_t\chi_t)(0)\\rangle$'
        ax.set_ylabel(operator_label)
        ax.legend(frameon = False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    plt.show()

def plotbeta_dependence_tderiv():
    fig,ax = plt.subplots()
    N_measure = 1000


    for i,beta in enumerate([1,2,3,4,5]):
        N = 10
        results,errs =LatticeRun.evaluate_observable_1D_dontsave(N,beta,N_measure,lambda l: Lattice.measure_time_derivative_correlator_zero_momentum_operator(4,l))


        
      

        results = results
        forward = np.roll(results,-1)
        popt,pcov = optimize.curve_fit(model, range(len(results)-1), forward[:-1], p0=(1, 0.1), bounds=([0, 0], [np.inf, np.inf]))
        print(f'Optimal parameters for N={N}: {popt}')
        xs = np.linspace(0, len(results)-1, 100)
        ys = model(xs, *popt)
        #ax.plot(xs, ys, label=f'Fit for $N$ = {N}')


        

        
        
        ax.errorbar(np.array(range(len(results)))+i*0.01,results,yerr = errs,fmt = 's',label = f'$\\beta$ = {beta}')
        #ax.errorbar(np.array(range(len(check_t)))+0.01,check_t,yerr = check_t_err,fmt = 's',label = f'Check $N_m$ = {N_measure}')
        
        
        ax.set_xlabel('t')
        operator_label = '$\langle(\\nabla_t\chi_t)(t)(\\nabla_t\chi_t)(0)\\rangle$'
        ax.set_ylabel(operator_label)
        ax.legend(frameon = False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    plt.show()
    
    
def model(x, a, b):
    return a * np.exp(-b * x)

def plotNmeasurements_dependence():
    fig,ax = plt.subplots()
    N_measure = 10000


    for i,N_measure in enumerate([100000]):
        beta =6 
        N = 10
        #LatticeRun.calibration(N,beta,0.25,True,300,calibration_runs = 100)
        #LatticeRun.generate_phis(N,beta,N_th,N_measure,True,4,guess = 200)
        operator = lambda l: Lattice.operator_gradient_sq(l)

        
        result,err =LatticeRun.evaluate_observable_dontsave(N,beta,N_measure,lambda l: Lattice.measure_average_local_operator(4,l,operator))
        results,errs =LatticeRun.evaluate_observable_1D_dontsave(N,beta,N_measure,lambda l: Lattice.measure_correlator_zero_momentum_operator(4,l,operator))


        r = result**2
        r_err = 2*result*err
        
        
        err = np.sqrt(errs**2)
        results = results-r
        forward = np.roll(results,-1)
        popt,pcov = optimize.curve_fit(model, range(len(results)-1), forward[:-1], p0=(1, 0.1), bounds=([0, 0], [np.inf, np.inf]))
        print(f'Optimal parameters for N={N}: {popt}')
        xs = np.linspace(0, len(results)-1, 100)
        ys = model(xs, *popt)
        #ax.plot(xs, ys, label=f'Fit for $N$ = {N}')


        

        
        
        ax.errorbar(np.array(range(len(results)))+i*0.1,results,yerr = err,fmt = 's',label = f'$N_m$ = {N_measure}')
        
        
        ax.set_xlabel('t')
        operator_label = '$\langle(\\nabla\chi)^2(t)(\\nabla\chi)^2(0)\\rangle$'
        ax.set_ylabel(operator_label)
        ax.legend(frameon = False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    plt.show()
def plotN_dependence():
    fig,ax = plt.subplots()
    N_measure = 1000


    for N in [8,10,16]:
        beta =2
        #LatticeRun.calibration(N,beta,0.25,True,300,calibration_runs = 100)
        #LatticeRun.generate_phis(N,beta,N_th,N_measure,True,4,guess = 200)
        operator = lambda l: Lattice.operator_gradient_sq(l)

        
        result,err =LatticeRun.evaluate_observable_dontsave(N,beta,N_measure,lambda l: Lattice.measure_average_local_operator(4,l,operator))
        results,errs =LatticeRun.evaluate_observable_1D_dontsave(N,beta,N_measure,lambda l: Lattice.measure_correlator_zero_momentum_operator(4,l,operator))


        r = result**2
        r_err = 2*result*err
        
        
        err = np.sqrt(errs**2)
        results = results-r
        forward = np.roll(results,-1)
        middle = len(results)//2
        results_mid = results[:middle]
        err_mid = err[:middle]
        popt,pcov = optimize.curve_fit(model, range(1,middle), results_mid[1:], p0=(1, 0.1), bounds=([0, 0], [np.inf, np.inf]))
        print(f'Optimal parameters for N={N}: {popt}')
        xs = np.linspace(1, len(results_mid)-1, 100)
        ys = model(xs, *popt)
        #ax.plot(xs, ys, label=f'Fit for $N$ = {N}')


        

        
        
        ax.errorbar(range(len(results)),results,yerr = err,fmt = 's',label = f'$N$ = {N}')
        
        
        ax.set_xlabel('t')
        operator_label = '$\langle(\\nabla\chi)^2(t)(\\nabla\chi)^2(0)\\rangle$'
        ax.set_ylabel(operator_label)
        ax.legend(frameon = False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    plt.show()


def correlate_operator(N, beta, N_measure, operator):
    results,errs =LatticeRun.evaluate_observable_1D_dontsave(N,beta,N_measure,lambda l: Lattice.measure_correlator_general_t_operator(4,l,operator))
    #results,errs =LatticeRun.evaluate_observable_1D_dontsave(N,beta,N_measure,lambda l: Lattice.measure_time_derivative_correlator_zero_momentum_operator(4,l))

    return results,errs




def vev_operator(N, beta, N_measure, operator):
    #LatticeRun.check_if_ensemble_exists(N, beta, N_measure)
    result,err =LatticeRun.evaluate_observable_dontsave(N,beta,N_measure,lambda l: Lattice.measure_vev_general_t_operator(4,l,operator))
    return result,err

def plot_correlator(ax,N, beta, N_measure, operator, operator_label):
    LatticeRun.check_if_ensemble_exists(N, beta, N_measure)

    results,errs = correlate_operator(N, beta, N_measure, operator)
   
    result,err = vev_operator(N, beta, N_measure, operator)
    disconnected = result**2
    disconnected_err = 2*result*err
    results = results-disconnected
    errs = np.sqrt(errs**2+disconnected_err**2)
    print("Results: ", results)
    print("Errors: ", errs)
    print("Disconnected contribution: ", disconnected)
    ax.errorbar(range(len(results)),results,yerr = errs,fmt = 's',label = f'$N$ = {N}, $\\beta$ = {beta}')
    ax.set_xlabel('t')
    ax.set_ylabel(operator_label)
    ax.legend(frameon = False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    


def projected_chi(lattice):
    dim = 4
    chi_t = np.sum(lattice, axis=tuple(range(dim - 1)))/(lattice.shape[0]**(dim-1)) # O(t) = sum_x O(x,t), shape (N,)
    return chi_t.copy()

def deriv_t(lattice):

    chi_t = projected_chi(lattice)
    dchi_t = np.roll(chi_t, -1) - chi_t
    return dchi_t

def local_gradient_sq(lattice):
    dim = 4
    gradient_sq = np.zeros_like(lattice)
    for i in range(dim):
           
        forward = np.roll(lattice, -1, axis = i)
        backward = np.roll(lattice, 1, axis = i)
        gradient_sq +=  0.5*(forward-lattice)**2+0.5*(backward-lattice)**2
    
    return gradient_sq

def zero_momentum_local_operator(lattice, operator):
    O = operator(lattice)
    O_t = np.sum(O, axis=tuple(range(3)))/(lattice.shape[0]**3) # O(t) = sum_x O(x,t), shape (N,)
    return O_t.copy()
def zero_momentum_gradient_sq(lattice):
    O = local_gradient_sq(lattice)
    O_t = np.sum(O, axis=tuple(range(3)))/(lattice.shape[0]**3) # O(t) = sum_x O(x,t), shape (N
    return O_t.copy()
def measure_F_local(lattice):
    dim = 4
    result = 0
    for i in range(dim):
        forward = np.roll(lattice, -1, axis = i)
        backward = np.roll(lattice, 1, axis = i)
        result += forward + backward - 2*lattice + 0.5*(forward-lattice)**2+0.5*(backward-lattice)**2
    return result
def measure_phi_without_zero_mode(lattice):
    phi = lattice.copy()
    phi -= np.mean(phi)
    return phi

def plot_correlator_gradient_sq_t(N, beta, N_measure):
    operator = lambda l: zero_momentum_gradient_sq(l)
   
    label = '$\langle(\\nabla\chi)^2(t)(\\nabla\chi)^2(0)\\rangle$'
    fig,ax = plt.subplots()
    plot_correlator(ax,N,beta,N_measure,operator,label)
    plt.show()
def plot_correlator_Neil_t(N, beta, N_measure):
    operator = lambda l: deriv_t(l)
    label = '$\langle(\\nabla_t\chi_t)(t)(\\nabla_t\chi_t)(0)\\rangle$'
    fig,ax = plt.subplots()
    plot_correlator(ax,N,beta,N_measure,operator,label)
    plt.show()
    
def plot_correlator_F_t(N, beta, N_measure):
    O_t = lambda l: zero_momentum_local_operator(l, measure_F_local)
    label = '$\langle F(t)F(0)\\rangle$'
   
    fig,ax = plt.subplots()
    plot_correlator(ax,N,beta,N_measure,O_t,label)
    plt.show()
def plot_correlator_phi_t(N, beta, N_measure):
    O_t = lambda l: zero_momentum_local_operator(l, measure_phi_without_zero_mode)
    label = "$\langle \\chi(t)\'\\chi(0)\'\\rangle$"
    
    fig,ax = plt.subplots()
    plot_correlator(ax,N,beta,N_measure,O_t,label)
    plt.show()

def main():
    N = 10
    beta = 6
    N_measure = 10**5

    #plot_correlator_gradient_sq_t(N, beta, N_measure)
    #plot_correlator_Neil_t(N, beta, N_measure)
    #plot_correlator_F_t(N, beta, N_measure)
    #plot_correlator_phi_t(N, beta, N_measure)
    plot_correlator_Neil_t(N, beta, N_measure)


   
if __name__ == "__main__":
    main()