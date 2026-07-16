from lattice import Lattice
import LatticeRun
import numpy as np
import plotting
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import processing
def plot_moving_average():
    N = 8
    lambda_ = -100
    N_measurements = 100000
    N_thermalization = 100000
    processing.analysis_corr_moving_average(N,lambda_,N_measurements,N_thermalization,True)
def plot_multiple_volumes_diff_sq_corr():
    Ns = [12]
    lambdas = [-1,-10,-100]
    N_thermalization = 1000
    N_measurements = 1000
    HMC = True
    plotting.plot_difference_lambdas_volumes_phis_corr(Ns,lambdas,N_measurements,N_thermalization,HMC)
def plot_multiple_volumes_diff_sq():
    Ns = [10]
    lambdas = [-1,-10,-100]
    N_thermalization = 1000
    N_measurements = 1000
    HMC = True
    plotting.plot_difference_lambdas_volumes_phis(Ns,lambdas,N_measurements,N_thermalization,HMC)

def run_dissipartive_separation(N,measurements):
    spacing = [n for n in range(N//2+1)]
    lambda_ = -1
    a  = Lattice(N,1,0,0,0)
    actions = []

    for i in spacing:
        pass
        """a.lattice = a.define_potential(0,i)

        #plt.plot(a.return_slice_lattice_x())
        #plt.show()
        action = a.run_dissipative(10,0.005,measurements)
        actions.append(action)"""
    return actions
def save_final_speration(N,measurements,sep):
    a = Lattice(N,1,0,0,0)
    a.lattice = a.define_potential(3,(3+sep)%N)
    a.run_dissipative2(10,0.005,measurements)
    plt.plot(a.return_slice_lattice_x())
    plt.show()
    return a.lattice
def run_dissipartive_separation2(N,measurements):
    spacing = [n for n in range(N//2+1)]
    lambda_ = -1
    a  = Lattice(N,1,0,0,0)
    actions = []
    #a.lattice = a.define_potential(0,0)
    #actions = a.run_dissipative2(10,0.001,measurements)
    print(actions)
    for i in spacing:
        pass
        

        a.lattice = a.define_potential(3,(i+3)%N)

        plt.plot(a.return_slice_lattice_x())
        
        action = a.run_dissipative2(10,0.001,measurements)
        plt.plot(a.return_slice_lattice_x())
        plt.show()
        
        #actions.append(action)
    return actions
def vev_analysis():
    N = 12
    lambda_ = -10
    N_measurements = [1000,1000,1000,1000]
    N_thermalization = [1000,1000,1000,1000]
    guess = 30

    partition = 100
    lambdas = [np.pi**2/4,np.pi**2/3]
    lambdas = [-1.0]
    N_measurements = [1000 for i in lambdas]
    N_thermalization = [1000 for i in lambdas]
        

    mode = [0,1]

        #plot_multiple_volumes_diff_sq_corr()
        #plot_moving_average()
        
        #print(run_dissipartive_separation2(8,1000))
        

        #print(a.lattice[1][0][0][0])
        #print(a.lattice)
        
    c =4*np.pi
    b = []
    for i,lambda_ in enumerate(lambdas):
        for msq in [0.46*lambda_+0.1,0.46*lambda_+0.2,0.46*lambda_+0.3,0.46*lambda_+0.4]:

            pass
        #guess = LatticeRun.calibration(N,lambda_,0.25,True,100,mode =1)
            guess = LatticeRun.generate_phis(N,lambda_,N_measurements[i],N_thermalization[i],True,4,False,guess = guess,msq=msq)
            LatticeRun.turn_phis_to_measurements_1D(N,lambda_,N_measurements[i],N_thermalization[i],lambda l: Lattice.measure_exponential_f(4,l,c,lambda_), "Exponential",True,4,False,msq = msq)       
        #print(processing.process_corr_deri(N,lambda_,N_measurements[i],N_thermalization[i],True))
        #LatticeRun.turn_phis_to_measurements(N,lambda_,N_measurements[i],N_thermalization[i],lambda l: Lattice.measure_exponential_average_f(4,l,c,lambda_), "Exponential Average",True,4)
        #LatticeRun.turn_phis_to_measurements_1D(N,lambda_,N_measurements[i],N_thermalization[i],lambda l: Lattice.measure_exponential_f(4,l,c,lambda_), "Exponential",True,4,False,msq = msq)
            processing.process_exponential(N,lambda_,N_measurements[i],N_thermalization[i],c,True,msq)
        #b+= [processing.process_exponential_average(N,lambda_,N_measurements[i],N_thermalization[i],c,True)[0]]
        #LatticeRun.measure_func_1D_3(N,lambda_,N_measurements[i],N_thermalization[i],lambda l: Lattice.measure_difference(4,l), "Difference",lambda l: Lattice.measure_correlation(4,l), "Correlation",lambda l: Lattice.measure_sq_array(4,l), "Square",True,4,mode=mode[i])
            #processing.process_exponential(N,lambda_,N_measurements[i],N_thermalization[i],c,True)
            #print(c)
        #processing.process_corr(N,lambda_,N_measurements[i],N_thermalization[i],True)
        #processing.process_Square(N,lambda_,N_measurements[i],N_thermalization[i],True)
            #plotting.plot_exponential(N,lambda_,N_measurements[i],N_thermalization[i],c,True,msq = msq)
        #processing.analysis_corr_two_portions_print(N,lambda_,N_measurements[i],N_thermalization[i],good_indices,bad_indices,True)
"""print(np.array(b)**0.5)
    fig,axes = plt.subplots(3,sharex=True,figsize=(10,10))

    exp = np.exp(-8*np.pi/(5*(np.array(lambdas))))
    axes[0].plot(lambdas,np.array(b)**0.5,'o')
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)
    axes[0].set_ylabel("$\sqrt{\\langle hh^{\dagger}\\rangle}$")
    axes[0].set_xlabel("$\lambda$")
    #plt.plot(lambdas,exp*np.array(b)**0.5,'o')
    axes[1].plot(lambdas,exp*np.array(b)**0.5,'o')
    z = np.linspace(lambdas[0],lambdas[len(lambdas)-1],1000)
    y = np.exp(-8*np.pi/(5*(np.array(z))))
    

    axes[2].plot(z,y,'--')
    axes[2].set_ylim(0,10**3)
    print(exp*np.array(b)**0.5)"""
    #plt.show()

def main():
    N = 10
    beta = 5
    N_th = 1000
    N_measure = 10000
    #LatticeRun.calibration(N,beta,0.25,True,300,calibration_runs = 100)
    #LatticeRun.generate_phis(N,beta,1000,1000,True,4,guess = 200)
    #LatticeRun.turn_phis_to_measurements_1D(N,beta,1000,1000,lambda l: Lattice.measure_difference(4,l), "Difference",True,4)
    
    for beta in [2,3,4,5]:
        pass
        LatticeRun.generate_phis(N,beta,N_measure,True,4,guess = 150)
    
    
    
    
    
def plotN_dependence():
    fig,ax = plt.subplots()
    N_measure = 1000


    for N in [8,10]:
        beta =2
        #LatticeRun.calibration(N,beta,0.25,True,300,calibration_runs = 100)
        #LatticeRun.generate_phis(N,beta,N_th,N_measure,True,4,guess = 200)
        operator = lambda l: Lattice.operator_gradient_sq(l)

        
        result,err =LatticeRun.evaluate_observable_dontsave(N,beta,N_measure,lambda l: Lattice.measure_average_local_operator(4,l,operator))
        results,errs =LatticeRun.evaluate_observable_1D_dontsave(N,beta,N_measure,lambda l: Lattice.measure_correlator_zero_momentum_operator(4,l,operator))


        r = result**2
        r_err = 2*result*err
        print(r,r_err)
        
        
        err = np.sqrt(errs**2)
        results = results-r
        print(results,err)
        forward = np.roll(results,-1)

        

        
        
        ax.errorbar(range(len(results)),results,yerr = err,fmt = 's',label = f'$N$ = {N}')
        ax.set_xlabel('t')
        operator_label = '$\langle(\\nabla\chi)^2(t)(\\nabla\chi)^2(0)\\rangle$'
        ax.set_ylabel(operator_label)
        ax.legend(frameon = False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    plt.show()
    
    
def main2():

    N = 15
    beta = 5
    #LatticeRun.calibration(N,beta,0.25,True,300,calibration_runs = 100)
    #LatticeRun.generate_phis(N,beta,1000,1000,True,4,guess = 200)
    #LatticeRun.turn_phis_to_measurements_1D(N,beta,1000,1000,lambda l: Lattice.measure_difference(4,l), "Difference",True,4)
    fig,ax = plt.subplots()
    for beta in [1,2,3,4,5]:
        #LatticeRun.calibration(N,beta,0.25,True,300,calibration_runs = 100)
        LatticeRun.generate_phis(N,beta,1000,1000,True,4,guess = 200)
        LatticeRun.turn_phis_to_measurements_1D(N,beta,1000,1000,lambda l: Lattice.measure_difference(4,l), "Difference",True,4)
    

        lambda_ = np.sqrt(1/beta)
        results,err = processing.process_difference(N,beta,1000,1000,True)
        phi_r = results/lambda_
        err_r = err/lambda_
        ax.errorbar(range(len(phi_r)),phi_r,yerr = err_r,fmt = 's',label = f'$\\beta$ = {beta}')
    ax.set_xlabel('r')
    ax.set_ylabel('$\langle\phi(r)-\phi(0)\\rangle$')
    ax.legend(frameon = False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.show()
    
# Call the function with the test array
    #result = print(Lattice.measure_difference(lattice))
main()
#run_dissipartive_separation2(8,1000)
#np.save("data.npy",save_final_speration(8,1000,3))
"""b = np.load("data.npy")
n = Lattice(8,1,0,0,0)
n.lattice = b
print(n.action2())
def format_array_for_mathematica(array):
    if isinstance(array, np.ndarray):
        return "{" + ", ".join(format_array_for_mathematica(subarray) for subarray in array) + "}"
    else:
        return str(array)


# Print the array in Mathematica style
a = format_array_for_mathematica(b)
with open("data.txt" , "w") as f:
    f.write(a)
print(b[0,0,0,0])"""

"""
 if myobject_m.updated <= 500:
        myobject_m.actions += [model.run_dissipative2_force_two(10,0.005,1,20,20,5,0,0,0,7,0,0,0)]
    elif myobject_m.updated <= 1000:
        myobject_m.actions += [model.run_dissipative2(10,0.005,1)]
    else:
        myobject_m.actions += [model.run_dissipative2(10,0.001,1)]
"""
def record_configurations_evolving_push(N,diss_runs,force_runs,time_step,x,x_2,c_1,c_2):
    model = Lattice(N,1,0,0,0)
    configs = [model.return_slice_lattice_x().copy()]
    fig,(ax1,ax2) = plt.subplots(1,2)
    actions = [model.action2()]
    ax1.plot(model.return_slice_lattice_x(),label = 'Initial two vortices',color = 'r')
    for i in range(diss_runs):
        if i<= force_runs:
            actions += [model.run_dissipative2_force_two(10,0.005,1,c_1,c_2,x,0,0,0,x_2,0,0,0)]
        elif i == 2000:
            actions += [model.run_dissipative2(10,0.001,1)]
            model.lattice = model.lattice/5
        else:
            actions += [model.run_dissipative2(10,time_step,1)]
        configs += [model.return_slice_lattice_x().copy()]
        print(i)
    path = f"Vortices/Two_vortices_push_{force_runs}_{diss_runs}_time_step_{time_step}_x_{x}_x2_{x_2}_c_{c_1}_c_2_{c_2}_N_{N}"
    np.save(path+"_configs.npy",configs)
    np.save(path+"_actions.npy",actions)
    ax1.plot(model.return_slice_lattice_x(),label = 'final configuration',color = 'b')
    ax2.plot(range(diss_runs+1),actions,label = 'Action',color = 'b')
    np.save("vortex.npy",model.lattice)
    plt.legend()
    plt.show()
def record_configurations_evolving_smoothed(N,diss_runs,time_step,separation,shift,Ls):
    model = Lattice(N,1,0,0,0)
    model.lattice = model.two_vortices_smoothed_initial_configuration(separation,shift,Ls)
    configs = [model.return_slice_lattice_x().copy()]
    fig,(ax1,ax2) = plt.subplots(1,2)
    actions = [model.action2()]
    ax1.plot(model.return_slice_lattice_x(),label = 'Initial two vortices',color = 'r')
    for i in range(diss_runs):
        actions += [model.run_dissipative2(10,0.005,1)]
        configs += [model.return_slice_lattice_x().copy()]
        print(i)
    path = f"Vortices/Two_vortices_smoothed_runs_{diss_runs}_time_step_{time_step}_separation_{separation}_shift_{shift}_Ls_{Ls}_N_{N}"
    np.save(path+"_configs.npy",configs)
    np.save(path+"_actions.npy",actions)
    ax1.plot(model.return_slice_lattice_x(),label = 'final configuration',color = 'b')
    ax2.plot(range(diss_runs+1),actions,label = 'Action',color = 'b')

    plt.legend()
    plt.show()
def animate_2(N,diss_runs,force_runs,time_step,x,x_2,c_1,c_2):
    frames_configs = np.load(f"Vortices/Two_vortices_push_{force_runs}_{diss_runs}_time_step_{time_step}_x_{x}_x2_{x_2}_c_{c_1}_c_2_{c_2}_N_{N}_configs.npy")
    frames_actions = np.load(f"Vortices/Two_vortices_push_{force_runs}_{diss_runs}_time_step_{time_step}_x_{x}_x2_{x_2}_c_{c_1}_c_2_{c_2}_N_{N}_actions.npy")
    x = np.arange(0,N,1)
    fig, (ax,ax2) = plt.subplots(1,2)
    def init_func():
        ax.clear()
        ax2.clear()
        ax.set_xlim((x[0], N))
        
        ax.set_ylim(min(frames_configs[0]), max(frames_configs[0]))
    def update_plot(i):
        ax.clear()
        ax2.clear()
        actions = frames_actions[:i]
        ax.plot(x,frames_configs[i],label = 'running configuration',color = 'b')
        ax.set_xlim(0,N)
        #ax.set_ylim(0,1)
        
        #ax.plot(x,frames_configs[0],label = 'two vortices initital config',color = 'r')

        ax2.plot(range(i),actions,label = 'Action',color = 'b')
        ax.legend()

    anim = FuncAnimation(fig,
                     update_plot,
                     frames=np.arange(0, len(frames_configs), 1),
                     init_func=init_func,
                     interval=0.1)
    plt.show()
def animate(N,diss_runs,time_step,separation,shift,Ls):
    frames_configs = np.load(f"Vortices/Two_vortices_smoothed_runs_{diss_runs}_time_step_{time_step}_separation_{separation}_shift_{shift}_Ls_{Ls}_N_{N}_configs.npy")
    frames_actions = np.load(f"Vortices/Two_vortices_smoothed_runs_{diss_runs}_time_step_{time_step}_separation_{separation}_shift_{shift}_Ls_{Ls}_N_{N}_actions.npy")
    x = np.arange(0,N,1)
    fig, (ax,ax2) = plt.subplots(1,2)
    def init_func():
        ax.clear()
        ax2.clear()
        ax.set_xlim((x[0], N))
        
        ax.set_ylim(min(frames_configs[0]), max(frames_configs[0]))
    def update_plot(i):
        ax.clear()
        ax2.clear()
        actions = frames_actions[:i]
        ax.plot(x,frames_configs[i],label = 'running configuration',color = 'b')
        ax.plot(x,frames_configs[0],label = 'two vortices initital config',color = 'r')

        ax2.plot(range(i),actions,label = 'Action',color = 'b')
        ax.legend()

    anim = FuncAnimation(fig,
                     update_plot,
                     frames=np.arange(0, len(frames_configs), 1),
                     init_func=init_func,
                     interval=20)
    plt.show()
def main2():
    N = 20
    diss_runs = 1000
    time_step = 0.005
    separation = 3
    shift = 5
    Ls = 6
    record_configurations_evolving_smoothed(N,diss_runs,time_step,separation,shift,Ls)
    #
    animate(N,diss_runs,time_step,separation,shift,Ls)
#main()

    #
#record_configurations_evolving_push(N,diss_runs,force_runs,time_step,x,x_2,c_1,c_2)
#animate_2(N,diss_runs,force_runs,time_step,x,x_2,c_1,c_2)
#plt.plot(model.return_slice_lattice_x())
#plt.show()
"""x = np.arange(0, N, 1)
y = model.return_slice_lattice_x()
fig, (ax,ax2,ax3) = plt.subplots(1,3)
actions = []
class myobject(object):
    def __init__(self):
        self.actions = []
        self.distances = []
myobjects = myobject()
#ax = plt.subplot(1, 1, 1)

data_skip = 50


def init_func():
    ax.clear()
    ax3.clear()
    
    ax.set_xlim((x[0], N))
    #plt.ylim((-1, 1))


def update_plot(i):
    #ax.plot(x[i:i+data_skip], y[i:i+data_skip], color='k')
    #ax.plot(x[i], y[i], color='r')


    myobjects.actions+= [model.run_dissipative2(10,0.001,1)]

    #ax.set_ylim(min(model.return_slice_lattice_x().copy()), max(model.return_slice_lattice_x().copy()))

    ax.plot(model.return_slice_lattice_x().copy(), color='g',label ='two vortices')
    z = np.arange(0, len(myobjects.actions),1)
    myobjects.distances += [abs(max(model.return_slice_lattice_x().copy())-min(model.return_slice_lattice_x().copy()))]

    ax2.plot( z,myobjects.actions, color='g')
    ax3.plot(z,myobjects.distances, color='g')



anim = FuncAnimation(fig,
                     update_plot,
                     frames=np.arange(0, len(x), data_skip),
                     init_func=init_func,
                     interval=20)
plt.show()
#anim.save('sine.mp4', dpi=150, fps = 30, writer='ffmpeg')"""


