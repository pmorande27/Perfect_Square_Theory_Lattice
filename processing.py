from Stats import Stats
import numpy as np
import matplotlib.pyplot as plt
import os

def process_exponential(N,lambda_,N_measurements,N_thermalization,c,HMC=False,msq=0):
    observable_name = "Exponential"  
    if msq == 0:
        file_name = "Results/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
    else:
        file_name = "Results/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+ ' msq = ' + str(msq) + '.npy'
    vals = np.load(file_name)
    results = np.zeros(N)
    errs = np.zeros(N)
    for i in range(N):
        
        result ,err,_,_= Stats(vals[i]).estimate()
        results[i] = result
        errs[i] = err
    file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" c = "+str(c)+'.npy'
    if HMC:
        if msq == 0:
            file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" c = "+str(c)+" HMC.npy"
        else:
            file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" c = "+str(c)+" msq = "+str(msq)+" HMC.npy"



    np.save(file_name,(results,errs))
def process_exponential_average(N,lambda_,N_measurements,N_thermalization,c,HMC=False):
    observable_name = "Exponential Average"  
    file_name = "Results/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'

    vals = np.load(file_name)
    results = np.zeros(N)
    errs = np.zeros(N)

    
        
    result ,err,_,_= Stats(vals).estimate()
    
    file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" c = "+str(c)+'.npy'
    if HMC:
        file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" c = "+str(c)+" HMC.npy"


    np.save(file_name,(result,err))

    return result,err
def process_derivative_correlator(N,lambda_,N_measurements,N_thermalization,HMC=False):
    observable_name = "Derivative Correlator"  
    file_name = "Results/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'

    vals = np.load(file_name)
    results = np.zeros(N)
    errs = np.zeros(N)

    
        
    result ,err,_,_= Stats(vals).estimate()
    
    file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
    if HMC:
        file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"


    np.save(file_name,(result,err))

    return result,err
def process_difference(N,lambda_,N_measurements,N_thermalization,HMC=False):
    observable_name = "Difference"  
    file_name = "Results/"+observable_name+"/"+observable_name+" beta = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'

    vals = np.load(file_name)
    results = np.zeros((N+1))
    errs = np.zeros((N+1))
    for i in range(0,N+1):
        
        result ,err,_,_= Stats(vals[i]).estimate()
        results[i] = result
        errs[i] = err
    file_name = "Results/Processed/"+observable_name+"/"+observable_name+" beta = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
    if HMC:
        file_name = "Results/Processed/"+observable_name+"/"+observable_name+" beta = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"


    np.save(file_name,(results,errs))

    return results,errs
def process_zero_momentum_correlator(N,lambda_,N_measurements,N_thermalization,HMC=False):
    observable_name = "Zero Momentum Correlator"
    file_name = "Results/"+observable_name+"/"+observable_name+" beta = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'

    vals = np.load(file_name)
    results = np.zeros((N+1))
    errs = np.zeros((N+1))
    for i in range(0,N+1):
        
        result ,err,_,_= Stats(vals[i]).estimate()
        results[i] = result
        errs[i] = err
    file_name = "Results/Processed/"+observable_name+"/"+observable_name+" beta = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
    if HMC:
        file_name = "Results/Processed/"+observable_name+"/"+observable_name+" beta = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"


    np.save(file_name,(results,errs))

    return results,errs

def process_correlator_I_zero_momentum(N,lambda_,N_measurements,N_thermalization,HMC=False):
    observable_name = "I Correlator Zero Momentum"
    file_name = "Results/"+observable_name+"/"+observable_name+" beta = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'

    vals = np.load(file_name)
    results = np.zeros((N+1))
    errs = np.zeros((N+1))
    for i in range(0,N+1):
        
        result ,err,_,_= Stats(vals[i]).estimate()
        results[i] = result
        errs[i] = err
    file_name = "Results/Processed/"+observable_name+"/"+observable_name+" beta = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
    if HMC:
        file_name = "Results/Processed/"+observable_name+"/"+observable_name+" beta = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"


    np.save(file_name,(results,errs))

    return results,errs

def process_I_zero_momentum(N,lambda_,N_measurements,N_thermalization,HMC=False):
    observable_name = "I Zero Momentum"
    if not os.path.exists("Results/Processed/"+observable_name):
        os.makedirs("Results/Processed/"+observable_name)

    file_name = "Results/"+observable_name+"/"+observable_name+" beta = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'

    vals = np.load(file_name)
    results =0
    errs = 0
    result ,err,_,_= Stats(vals).estimate()
    
    file_name = "Results/Processed/"+observable_name+"/"+observable_name+" beta = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
    if HMC:
        file_name = "Results/Processed/"+observable_name+"/"+observable_name+" beta = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"


    np.save(file_name,(result,err))

    return result,err

def process_corr(N,lambda_,N_measurements,N_thermalization,HMC=False):
    observable_name = "Correlation"  
    file_name = "Results/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'

    vals = np.load(file_name)
    results = np.zeros(N+1)
    errs = np.zeros(N+1)
    for i in range(0,N+1):
        
        result ,err,_,_= Stats(vals[i]).estimate()
        results[i] = result
        errs[i] = err
    file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
    if HMC:
        file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"


    np.save(file_name,(results,errs))


    return results,errs
def process_corr_deri(N,lambda_,N_measurements,N_thermalization,HMC=False):
    observable_name = "Derivative Correlator r"  
    file_name = "Results/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'

    vals = np.load(file_name)
    results = np.zeros(N+1)
    errs = np.zeros(N+1)
    for i in range(0,N+1):
        
        result ,err,_,_= Stats(vals[i]).estimate()
        results[i] = result
        errs[i] = err
    file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
    file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"


    #np.save(file_name,(results,errs))


    return results,errs
def process_Square(N,lambda_,N_measurements,N_thermalization,HMC=False):
    observable_name = "Square"  
    file_name = "Results/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'

    vals = np.load(file_name)
    results = np.zeros(N+1)
    errs = np.zeros(N+1)
    for i in range(0,N+1):
        
        result ,err,_,_= Stats(vals[i]).estimate()
        results[i] = result
        errs[i] = err
    file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
    if HMC:
        file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"


    np.save(file_name,(results,errs))
def process_field(N,lambda_,N_measurements,N_thermalization,HMC=False):
    observable_name = "Field"  
    file_name = "Results/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'

    val = np.load(file_name)
    result ,err,_,_= Stats(val).estimate()
    file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
    if HMC:
        file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"


    np.save(file_name,(result,err))


    return result,err
def process_field(N,lambda_,N_measurements,N_thermalization,HMC=False):
    observable_name = "Exponetial Average"  
    file_name = "Results/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'

    val = np.load(file_name)
    result ,err,_,_= Stats(val).estimate()
    file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
    if HMC:
        file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"


    np.save(file_name,(result,err))

    return result,err

def process_corr_portion(N,lambda_,N_measurements,N_thermalization,indices,HMC=False, ):
    observable_name = "Correlation"  
    file_name = "Results/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'

    vals = np.load(file_name)
    val_r = [[] for i in range(0,N+1)]
    print(vals.shape)
    for j in range(0,N+1):
        for i in indices:
            val_r[j].append(vals[j][i])
    vals = val_r.copy()
    results = np.zeros(N+1)
    errs = np.zeros(N+1)
    for i in range(0,N+1):
        
        result ,err,_,_= Stats(vals[i]).estimate()
        results[i] = result
        errs[i] = err
    file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
    if HMC:
        file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"


    np.save(file_name,(results,errs))

    return results,errs

def analysis_corr_two_portions(N,lambda_,N_measurements,N_thermalization,indices_one,indices_two,HMC=False, ):
    observable_name = "Correlation"  
    file_name = "Results/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'

    vals = np.load(file_name)
    val_r = [[] for i in range(0,N+1)]
    for j in range(0,N+1):
        for i in indices_one:
            val_r[j].append(vals[j][i])
    vals = val_r.copy()
    results = np.zeros(N+1)
    errs = np.zeros(N+1)
    for i in range(0,N+1):
        
        result ,err,_,_= Stats(vals[i]).estimate()
        results[i] = result
        errs[i] = err
    results_one = results.copy()
    print(results)
    vals = np.load(file_name)

    val_r = [[] for i in range(0,N+1)]
    for j in range(0,N+1):
        for i in indices_two:
            val_r[j].append(vals[j][i])
    vals = val_r.copy()
    results = np.zeros(N+1)
    errs = np.zeros(N+1)
    for i in range(0,N+1):
        
        result ,err,_,_= Stats(vals[i]).estimate()
        results[i] = result
        errs[i] = err
    results_two = results.copy()
    print(results)
    print(abs(results_one-results_two))



    return results,errs
def analysis_corr_two_portions_print(N,lambda_,N_measurements,N_thermalization,indices_one,indices_two,HMC=False, ):
    observable_name = "Correlation"  
    file_name = "Results/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'

    vals = np.load(file_name)
    val_r = [[] for i in range(0,N+1)]
    for j in range(0,N+1):
        for i in indices_one:
            val_r[j].append(vals[j][i])
    vals = val_r.copy()
    results = np.zeros(N+1)
    errs = np.zeros(N+1)
    print(vals[0])
    for i in range(0,N+1):
        pass
        
        #result ,err,_,_= Stats(vals[i]).estimate()
        
    print(results)
    vals = np.load(file_name)

    val_r = [[] for i in range(0,N+1)]
    for j in range(0,N+1):
        for i in indices_two:
            val_r[j].append(vals[j][i])
    vals = val_r.copy()
    results = np.zeros(N+1)
    errs = np.zeros(N+1)
    print(vals[0])
        
        
    print(results)

def analysis_corr_moving_average(N,lambda_,N_measurements,N_thermalization,HMC=False, ):
    observable_name = "Correlation"  
    file_name = "Results/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'

    vals = np.load(file_name)
   
    results = np.zeros(N+1)
    errs = np.zeros(N+1)
    avs = np.zeros((N+1,N_measurements))
    
    x = np.arange(0,N_measurements)
    colors = ['red','blue','green','orange','purple','black','brown','pink','yellow','cyan']
    fig,ax = plt.subplots()
    for j in range(N+1):
        color = colors[j]
        for i in range(0,N_measurements):
            
            avs[j,i] = np.mean(vals[j][1:i+2]).copy()
    for j in range(N+1):
        plt.plot(x,avs[j].copy(),label = "|x-y| = "+str(j),color = colors[j])
    for j in range(N+1):
        color = colors[j]
        for i in range(0,N_measurements):
            avs[j][i] = np.mean(2*vals[0][1:i+2]-2*vals[j][1:i+2])
        
        plt.plot(avs[j],label = "Substracted |x-y| = "+str(j),color = color,alpha = 0.5)
    
    plt.xlabel("N measurements")
    ax.spines[['right', 'top']].set_visible(False)

    plt.ylabel("$\langle\phi(x)\phi(y)\\rangle$")
    plt.legend()
    plt.show()



    return results,errs