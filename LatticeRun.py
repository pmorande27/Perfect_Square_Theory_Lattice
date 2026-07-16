import numpy as np
from Stats import Stats
from lattice import Lattice
from alive_progress import alive_bar
import os
import processing

def calibration( N, beta, width_guess,HMC,  N_steps_guess=0,calibration_runs=10**2,msq = 0):


    accel = False
    print('Calibration with beta = ' + str(beta) + " N = " +str(N) )
    up = 0.8
    low = 0.4
    max_count = 10
    results = [0 for i in range(max_count)]
    width = width_guess
    if not HMC:
        for i in range(max_count):
            N_measurements = 0
            N_thermalization = 0
            lat = Lattice( N, beta, N_measurements, N_thermalization,width,False,mode = 1)
            rate = lat.calibration_runs(calibration_runs, calibration_runs)
            
            
            results[i] = (rate-up,width)
            print(rate,width)
            

            if rate <=up and rate >= low:
                if msq == 0:
                    file_name = "Parameters/Calibration parameters beta = " + str(beta) + " N = " + str(N)
                else:
                    file_name = "Parameters/Calibration parameters beta = " + str(beta) + " N = " + str(N) + f" msq = {msq}"
                print("-----------------")
                np.save(file_name, [width])
                return width
            
            else:
                new_width = width
                if rate > up:
                    new_width *= 2
                else:
                    new_width *=    0.5
            
            width = new_width
    else:
        minimum = 10
        N_tau = N_steps_guess
        epsilon = 1/(N_tau)
        width = 0
        print('Calibration with beta = ' + str(beta) + " N = " +str(N)+ " N_tau = " + str(N_tau))
        up = 1.0
        low = 0.75
        max_count = 10
        results = [0 for i in range(max_count)]
        for i in range(max_count):
            if N_tau < minimum and i != 0:
                N_tau = minimum
                epsilon = 1/(N_tau)
                if msq == 0:
                    file_name = "Parameters/Calibration parameters beta = " + str(beta) + " N = " + str(N) + " HMC"
                else:
                    file_name = "Parameters/Calibration parameters beta = " + str(beta) + " N = " + str(N) + f" HMC msq = {msq}"
            
                print(file_name)
                np.save(file_name, [N_tau,1/N_tau])
                print("-----------------")
                print(rate, N_tau)
                return N_tau
            epsilon = 1/(N_tau)
            lat = Lattice(N, beta,0,width, HMC, epsilon,N_tau)
            lat.calibration_runs(calibration_runs, calibration_runs)
            rate = lat.accepted/calibration_runs
            d_rate = 0.85-rate
            results[i] = (rate-up,N_tau)
            print(rate*100,N_tau)
            

            new_N = int(np.rint(N_tau*(1+d_rate)))
            if rate <=up and rate >= low:
                if msq == 0:
                        file_name = "Parameters/Calibration parameters beta = " + str(beta) + " N = " + str(N) + " HMC"
                else:
                        file_name = "Parameters/Calibration parameters beta = " + str(beta) + " N = " + str(N) + f" HMC msq = {msq}"
                    
                    
                np.save(file_name, [N_tau,1/N_tau])
                print("-----------------")
                print(file_name,rate, N_tau)
                return N_tau
            if new_N == N_tau:
                if d_rate <0:
                    new_N -= 1
                else:
                    new_N +=1
            
            N_tau = new_N

       

    print("-----------------")
    print("Calibration Unsucessful, better run:")
    results_abs = [(abs(x),y) for (x,y) in results]
    d_rate, width = min(results_abs)
    d_rate_2 = lookup(d_rate,width,results)
    rate = (d_rate_2+up)*100
    print(rate,width)
    file_name = "ChiralParams/Chiral Calibration parameters beta = " + str(beta) + " N = " + str(N) 
    np.save(file_name, [width])
    return width
def load_calibration( N, lambda_ ,HMC = False,msq = 0):

    file_name = "Parameters/Calibration parameters beta = " + str(lambda_) + " N = " + str(N) + ".npy"
    if HMC:
        if msq == 0:
            file_name = "Parameters/Calibration parameters beta = " + str(lambda_) + " N = " + str(N) + " HMC.npy"
        else:
            file_name = "Parameters/Calibration parameters beta = " + str(lambda_) + " N = " + str(N) + f" HMC msq = {msq}.npy"
    values = np.load(    file_name)
    if HMC:
        return values[0],values[1]
    else:  

        return values[0]


def lookup(d_rate,N_tau,results):
    for (x,y) in results:
        if abs(x) == d_rate and y == N_tau:
            return x
        
def generate_phis( N,lambda_, N_measure,HMC = False,dim=4,guess = 40,mode=1,msq=0):
    msq =0
    if msq == 0:
            file_name = "Parameters/Calibration parameters beta = " + str(lambda_) + " N = " + str(N)+' HMC.npy'
    else:  
        file_name = "Parameters/Calibration parameters beta = " + str(lambda_) + " N = " + str(N)+f' HMC msq = {msq}.npy'
    if os.path.exists(file_name):
        print('Calibration already done',file_name)
        pass
    else:
        print(file_name)
        calibration(N,lambda_,1,HMC,guess)

    observable_name = 'phi'
    count = 0
    while True:
        try:
            if count == 10:    
                count = 0
                print('Recalibration')
                calibration(N,lambda_,1)
            
            if msq == 0:
                    file_name = "Results/"+observable_name+"/"+observable_name+" beta = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measure)+'.npy'
            else:
                    file_name = "Results/"+observable_name+"/"+observable_name+" beta = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measure)+' msq = {msq}.npy'
            if HMC:
                N_tau,epsilon = load_calibration(N,lambda_,HMC)
                N_tau = int(N_tau)
                model = Lattice(N,lambda_,N_measure,1,HMC,epsilon,N_tau,dim,mode=mode)
            
            else:
                width = load_calibration(N,lambda_,)
                model = Lattice(N,lambda_,N_measure,width,False,0,0,dim,msq =msq) 

            results = model.generate_phis()

        except (ValueError) as e:
            print(e)
            count+= 1
            continue
        break
    #print(Stats(vals).estimate())
    if os.path.exists("Results/"+observable_name) == False:
        os.makedirs("Results/"+observable_name)
    np.save(file_name,results)
    return N_tau




def turn_phis_to_measurements( N,lambda_, N_measure, observable, observable_name,HMC = False,thermalization_percent=0.1):
    file_name = "Results/"+"phi"+"/"+"phi"+" beta = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measure)+'.npy'
    configurations = np.load(file_name)
    measurements = []
    n_thermal = int(thermalization_percent * len(configurations))
    configurations = configurations[n_thermal:]
    with alive_bar(len(configurations)) as bar:
        for config in configurations:
            measurements.append(observable(config))
            bar()
    measurements = np.array(measurements)
    if os.path.exists("Results/"+observable_name) == False:
        os.makedirs("Results/"+observable_name)
    file_name = "Results/"+observable_name+"/"+observable_name+" beta = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measure)+'.npy'
    np.save(file_name,measurements)
def evaluate_observable_dontsave(N,lambda_, N_measure, observable,thermalization_percent=0.1):
    file_name = "Results/"+"phi"+"/"+"phi"+" beta = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measure)+'.npy'
    configurations = np.load(file_name)
    measurements = []
    n_thermal = int(thermalization_percent * len(configurations))
    configurations = configurations[n_thermal:]
    with alive_bar(len(configurations)) as bar:
        for config in configurations:
            measurements.append(observable(config))
            bar()
    measurements = np.array(measurements)
    
    vals = measurements
    results =0
    errs = 0
    result ,err,_,_= Stats(vals).estimate()
    
   

    return result,err

def evaluate_observable_1D_dontsave(N,lambda_, N_measure, observable,msq = 0,thermalization_percent=0.1):
   
    if msq == 0:
            
        file_name = "Results/"+"phi"+"/"+"phi"+" beta = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measure)+'.npy'
    else:
        file_name = "Results/"+"phi"+"/"+"phi"+" beta = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measure)+f' msq = {msq}.npy'
    configurations = np.load(file_name)
    print('Loaded configurations from',file_name)
    measurements = []
    n_thermal = int(thermalization_percent * len(configurations))
    configurations = configurations[n_thermal:]
    with alive_bar(len(configurations)) as bar:

        for config in configurations:
                measurements.append(observable(config))
                bar()
    measurements = np.array(measurements)
    vals = measurements.swapaxes(0,1)
    results = np.zeros((N+1))
    errs = np.zeros((N+1))
    for i in range(0,N+1):
        
        result ,err,_,_= Stats(vals[i]).estimate()
        results[i] = result
        errs[i] = err
   
    return results,errs

    np.save(file_name,vals)
def turn_phis_to_measurements_1D( N,lambda_, N_measure, observable, observable_name,HMC = False,dim=4,accel =False, mass = 0.1,msq = 0,thermalization_percent=0.1):
    
    if accel == False:
        if msq == 0:
                
            file_name = "Results/"+"phi"+"/"+"phi"+" beta = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measure)+'.npy'
        else:
            file_name = "Results/"+"phi"+"/"+"phi"+" beta = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measure)+f' msq = {msq}.npy'
    else:
        file_name = "ChiralResults/"+"phi"+"/"+"phi"+" beta = " + str(lambda_) + " N = " + str(N)  + " N measurements = "  + str(N_measure)+" Accel.npy"
    configurations = np.load(file_name)
    measurements = []
    n_thermal = int(thermalization_percent * len(configurations))
    configurations = configurations[n_thermal:]
    with alive_bar(len(configurations)) as bar:

        for config in configurations:
                measurements.append(observable(config))
                bar()
    measurements = np.array(measurements)
    if os.path.exists("Results/"+observable_name) == False:
        os.makedirs("Results/"+observable_name)
   
    file_name = "Results/"+observable_name+"/"+observable_name+" beta = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measure)+'.npy'
    vals = measurements.swapaxes(0,1)
    np.save(file_name,vals)