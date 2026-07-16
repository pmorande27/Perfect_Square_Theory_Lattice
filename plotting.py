import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import kn
import os
def plot_difference(N,lambda_,N_measurements,N_thermalization,HMC=False):
    observable_name = "Difference"  
    file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
    if HMC:
        file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"
    vals = np.load(file_name)
    results = vals[0]#/np.abs(lambda_)
    errs = vals[1]#/np.abs(lambda_)
    axis = plt.subplot(111)

    print(results)
    x = np.arange(1,int(N/2)+1,1)
    print(x)
    axis.errorbar(x,results,yerr=errs,fmt='.',capsize=2,label = "Data")
    axis.plot(x,results,linestyle = 'dashed',color = 'blue',linewidth = 0.75)
    plt.xlabel("$|x-y|$")
    plt.ylabel("$\langle( \phi(x) - \phi (y))^2\\rangle$")
    axis.spines[['right', 'top']].set_visible(False)
    
    plt.title("$\lambda = $ "+str(lambda_)+ " N = "+str(N) + " HMC = "+str(HMC))

    plt.show()
def plot_exponential(N,lambda_,N_measurements,N_thermalization,c,HMC=False,msq = 0):
    observable_name = "Exponential"  

    file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" c = "+str(c)+'.npy'
    if HMC:
        if msq == 0:
            file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" c = "+str(c)+" HMC.npy"
        else:
            file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" c = "+str(c)+f" msq = {msq} HMC.npy"
    vals = np.load(file_name)
    results = np.array(list(vals[0]) + [vals[0][0]])
    errs = np.array(list(vals[1]) + [vals[1][0]])
    axis = plt.subplot(111)


    x = np.arange(0,N+1)
    def model(x,a):
        return a*1/(4*np.pi**2*x**2)
    popt,pcov = curve_fit(model,np.arange(1,2),results[1:2],sigma = errs[1:2])
    
    chisq = np.sum((results[1:2]-model(np.arange(1,2),*popt))**2/errs[1:2]**2)/(1)
    print("Chisq = ", chisq)
    asq = 1/popt[0]
    x = np.arange(0,N+1)
    #axis.errorbar(x,results,yerr=errs,fmt='.',capsize=2,label = "Data")
    m_eff = np.abs(np.log(results/np.roll(results,-1)))
    m_errs = np.sqrt((errs/results)**2 + (np.roll(errs,-1)/np.roll(results,-1))**2)
    #print(m_eff)

    axis.errorbar(x,m_eff,yerr = m_errs,fmt = ".",capsize = 2,label = "Data")
    model = lambda x,a: a
    popt,pcov = curve_fit(model,np.arange(2,N//2),m_eff[2:N//2],sigma = m_errs[2:N//2])
    print("Correlation Length: ",1/popt[0])
    
    plt.plot(x, [popt[0] for i in x],label = f"$\chi^2 = {round(chisq,2)}$",color = 'red')
    #plt.fill_between(x,model(x,*popt)-np.sqrt(np.diag(pcov)),model(x,*popt)+np.sqrt(np.diag(pcov)),alpha = 0.5,color = 'red')
    #print(popt)

    
    #axis.plot(x,results,linestyle = 'dashed',color = 'blue',linewidth = 0.75)
    plt.xlabel("$|x-y|$")
    plt.ylabel("$\langle( \exp{ic(\phi(x)-\phi(y))}\\rangle$")

    #axis.errorbar(x,[results[i]*x[i]**(3/2) for i in range(len(x))],yerr=errs,fmt='.',capsize=2,label = "Data")
    l = 0
    for i in range(len(x)//2):
        l+= results[i]*x[i]**(3/2)*asq*2*np.pi**2

    #axis.plot(x,results,linestyle = 'dashed',color = 'blue',linewidth = 0.75)
    plt.xlabel("$|x-y|$")
    plt.ylabel("$\langle( \exp{ic(\phi(x)-\phi(y))}\\rangle$")
    
    

    axis.spines[['right', 'top']].set_visible(False)
    
    plt.title("$\lambda = $ "+str(lambda_)+ " N = "+str(N) + " HMC = "+str(HMC))
    observable_name = "Exponential Average"
    file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" c = "+str(c)+' HMC.npy'
    #print(file_name)
    if os.path.exists(file_name):
        vals = np.load(file_name)
        result = vals[0]
        err = vals[1]
        x = np.linspace(0,N,100)
        #axis.plot(x,[result for i in x],linestyle = 'dashed',color = 'red',linewidth = 0.75,label = "$\langle hh^\dag \\rangle$")

        #axis.fill_between(x,[result-err for i in x],[result+err for i in x],alpha = 0.5,color = 'red')
    else:
        pass
        #print("No average file")
    observable_name = "Difference"  
    file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
    if HMC:
        file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"
    """def model(x,b,c):
        return  b*np.exp(-c*x)
    t_min = 1
    popt,pcov = curve_fit(model,np.arange(t_min,N//2),results[t_min:N//2],sigma = errs[t_min:N//2])
    hsq,c,m =result,popt[0],popt[1]
    
    err_hsq,(err_c,err_m) = err,np.sqrt(np.diag(pcov))
    print(f"L>= {8.5/m}?")
    print(f"$hsq = {hsq}\pm {err_hsq}, c = {c} \pm {err_c} ,m = {m}\pm {err_m}$")
    chisq = np.sum((results[t_min:N//2]-model(np.arange(t_min,N//2),*popt))**2/errs[t_min:N//2]**2)/(len(results[t_min:N//2])-3)
    x = np.linspace(t_min,N//2,100)
    plt.plot(x,model(x,*popt),label = f"$\chi^2 = {round(chisq,2)}$",color = 'red')"""
    plt.legend()
    if os.path.exists(file_name):
        vals = np.load(file_name)
        results,errs = list(vals[0]),list(vals[1])
        #print(results)
   
        results = [0] + results + np.flip(results[0:N//2-1]).tolist() + [0]
        #print(x,vals)
        results = np.array(results)

        axis = plt.subplot(111)


        x = np.arange(0,N+1)
        c = 4*np.pi
        #print(c**2/2*np.array(results)/np.abs(lambda_))
        #axis.plot(x,np.exp(-c**2/2*np.array(results)/np.abs(lambda_)),label = "Data")
        #axis.plot(x,results,linestyle = 'dashed',color = 'blue',linewidth = 0.75)
    plt.xlabel("$|x-y|$")
    plt.ylabel("$\langle( \exp{ic(\phi(x)-\phi(y))}\\rangle$")
    #plt.yscale('log')
    

    axis.spines[['right', 'top']].set_visible(False)
    
    plt.show()

def plot_exponential_lambdas_volumes(Ns,lambdas,N_measurements,N_thermalization,c,HMC=False):
    colors = ['red','blue','green','orange','purple'] + ['cyan','magenta','yellow','brown','pink','grey','olive','lime','teal','coral','navy','maroon','gold','indigo','tan','wheat','azure','ivory','lavender','plum','crimson','salmon','tomato','khaki','snow','orchid']
    for j in range(len(Ns)):
        N = Ns[j]
        for i,lambda_ in enumerate(lambdas):
            color = colors[i]
            observable_name = "Exponential"  
            file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" c = "+str(c)+'.npy'
            if HMC:
                file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" c = "+str(c)+" HMC.npy"
            vals = np.load(file_name)
            results = np.array(list(vals[0]) + [vals[0][0]])
            errs = np.array(list(vals[1]) + [vals[1][0]])
            axis = plt.subplot(111)
            axis.spines[['right', 'top']].set_visible(False)
        
            plt.title(f"$c = {c}$ "+" HMC = "+str(HMC))


            print(results)
            x = np.arange(0,N+1)
            print(x)
            axis.errorbar(x,results,yerr=errs,fmt='.',capsize=2,label = "$\lambda = $"+ str(lambda_),color = color)
            plt.plot(x,results,linestyle = 'dashed',color = color,linewidth = 0.75)
    plt.ylabel("$\langle( \exp{ic(\phi(x)-\phi(y))}\\rangle$")
    plt.legend()
    plt.show()
def plot_exponential_lambdas_cs(N,lambdas,N_measurements,N_thermalization,cs,HMC=False):
    colors = ['red','blue','green','orange','purple'] + ['cyan','magenta','yellow','brown','pink','grey','olive','lime','teal','coral','navy','maroon','gold','indigo','tan','wheat','azure','ivory','lavender','plum','crimson','salmon','tomato','khaki','snow','orchid']
    for j in range(len(cs)):
        c = cs[j]
        for i,lambda_ in enumerate(lambdas):
            color = colors[j]
            observable_name = "Exponential"  
            file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" c = "+str(c)+'.npy'
            if HMC:
                file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" c = "+str(c)+" HMC.npy"
            vals = np.load(file_name)
            results = np.array(list(vals[0]) + [vals[0][0]])
            errs = np.array(list(vals[1]) + [vals[1][0]])
            axis = plt.subplot(111)
            axis.spines[['right', 'top']].set_visible(False)
            model = lambda x,a,b: a/(x**2) +b
            popt,pcov = curve_fit(model,np.arange(1,N//2),results[1:N//2])
            chisq = np.sum((results[1:N//2]-model(np.arange(1,N//2),*popt))**2/errs[1:N//2]**2)/(N//2-2)
            print(c,chisq)
            xs = np.linspace(1,N//2,100)
            plt.plot(xs,model(xs,*popt),label = f"$c = {c}$ $\chi^2 = {round(chisq,2)}$",color = color)


        
            plt.title(f"$c = {cs}$")


            print(results)
            x = np.arange(0,N+1)
            print(x)
            axis.errorbar(x,results,yerr=errs,fmt='.',capsize=2,color = color)
            plt.plot(x,results,linestyle = 'dashed',color = color,linewidth = 0.75)
    plt.ylabel("$\langle( \exp{ic(\phi(x)-\phi(y))}\\rangle$")
    plt.legend()
    plt.show()


def plot_exponential_lambdas(N,lambdas,N_measurements,N_thermalization,c,HMC=False):
    colors = ['red','blue','green','orange','purple'] + ['cyan','magenta','yellow','brown','pink','grey','olive','lime','teal','coral','navy','maroon','gold','indigo','tan','wheat','azure','ivory','lavender','plum','crimson','salmon','tomato','khaki','snow','orchid']

    for i,lambda_ in enumerate(lambdas):
        color = colors[i]
        observable_name = "Exponential"  
        file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" c = "+str(c)+'.npy'
        if HMC:
            file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" c = "+str(c)+" HMC.npy"
        vals = np.load(file_name)
        results = np.array(list(vals[0]) + [vals[0][0]])
        errs = np.array(list(vals[1]) + [vals[1][0]])
        axis = plt.subplot(111)


        print(results)
        x = np.arange(0,N+1)
        print(x)
        axis.errorbar(x,results,yerr=errs,fmt='.',capsize=2,label = "$\lambda = $"+ str(lambda_),color = color)
        axis.plot(x,results,linestyle = 'dashed',color = color,linewidth = 0.75)
        plt.xlabel("$|x-y|$")
        plt.ylabel("$\langle( \exp\{ic(f(x)-f(y))\}\\rangle$")
        axis.spines[['right', 'top']].set_visible(False)
        plt.legend()
        
        plt.title(f"$c = {c}$ "+ " N = "+str(N) + " HMC = "+str(HMC))


    plt.show()
    return results
def plot_exponential_lambdas_f(N,lambdas,N_measurements,N_thermalization,c,HMC=False):
    colors = ['red','blue','green','orange','purple'] + ['cyan','magenta','yellow','brown','pink','grey','olive','lime','teal','coral','navy','maroon','gold','indigo','tan','wheat','azure','ivory','lavender','plum','crimson','salmon','tomato','khaki','snow','orchid']

    for i,lambda_ in enumerate(lambdas):
        color = colors[i]
        observable_name = "Exponential"  
        file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" c = "+str(c)+'.npy'
        if HMC:
            file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" c = "+str(c)+" HMC.npy"
        vals = np.load(file_name)
        results = np.array(list(vals[0]) + [vals[0][0]])
        errs = np.array(list(vals[1]) + [vals[1][0]])
        axis = plt.subplot(111)


        print(results)
        x = np.arange(0,N+1)
        print(x)
        axis.errorbar(x,results,yerr=errs,fmt='.',capsize=2,label = "$\lambda = $"+ str(lambda_),color = color)
        axis.plot(x,results,linestyle = 'dashed',color = color,linewidth = 0.75)
        plt.xlabel("$|x-y|$")
        plt.ylabel("$\langle( \exp\{ic(f(x)-f(y))\}\\rangle$")
        axis.spines[['right', 'top']].set_visible(False)
        
        plt.title(f"$c = {c}$ "+ " N = "+str(N) + " HMC = "+str(HMC))
        plt.legend( )

    plt.show()
def plot_exponential_lambdas_f_Vs(Ns,lambdas,N_measurements,N_thermalization,c,HMC=False):
    colors = ['red','blue','green','orange','purple'] + ['cyan','magenta','yellow','brown','pink','grey','olive','lime','teal','coral','navy','maroon','gold','indigo','tan','wheat','azure','ivory','lavender','plum','crimson','salmon','tomato','khaki','snow','orchid']
    labels = []
    for N in Ns:
        for i,lambda_ in enumerate(lambdas):
            color = colors[i]
            observable_name = "Exponential"  
            file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" c = "+str(c)+'.npy'
            if HMC:
                file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" c = "+str(c)+" HMC.npy"
            vals = np.load(file_name)
            results = np.array(list(vals[0]) + [vals[0][0]])
            errs = np.array(list(vals[1]) + [vals[1][0]])
            axis = plt.subplot(111)


            print(results)
            x = np.arange(0,N+1)
            print(x)
            if lambda_ not in labels:
                axis.errorbar(x,results,yerr=errs,fmt='.',capsize=2,label = "$\lambda = $"+ str(lambda_),color = color)
                labels.append(lambda_)
            else:
                axis.errorbar(x,results,yerr=errs,fmt='.',capsize=2,color = color)
            axis.plot(x,results,linestyle = 'dashed',color = color,linewidth = 0.75)
            plt.xlabel("$|x-y|$")
            plt.ylabel("$\langle( \exp\{ic(f(x)-f(y))\}\\rangle$")
            axis.spines[['right', 'top']].set_visible(False)
            
        plt.title(f"$c = {c}$ "+ " N = "+str(Ns) + " HMC = "+str(HMC))
        plt.legend( )

    plt.show()
                
def plot_difference_lambdas_volumes_phis(Ns,lambdas,N_measurements,N_thermalization,HMC=False):
    axis = plt.subplot(111)
    colors = ['red','blue','green','orange','purple'] + ['cyan','magenta','yellow','brown','pink','grey','olive','lime','teal','coral','navy','maroon','gold','indigo','tan','wheat','azure','ivory','lavender','plum','crimson','salmon','tomato','khaki','snow','orchid']
    labels = []
    for j in range(len(Ns)):
        N = Ns[j]
        for i,lambda_ in enumerate(lambdas):
            observable_name = "Difference"  
            file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
            if HMC:
                file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"       
            color = colors[i]
            vals = np.load(file_name)
            results = list(vals[0]/np.abs(lambda_))
            errs = list(vals[1]/np.abs(lambda_))
            results = [0] + results + np.flip(results[0:N//2-1]).tolist() + [0]
            errs = [0] + errs + np.flip(errs[0:N//2-1]).tolist() + [0]


            x = np.arange(0,N+1,1)
            print(x.shape,np.array(results).shape)
            print(lambda_,results)
            #results = np.roll(results,-1)-results 
            label = "$\lambda = $"+ str(lambda_)
            if label not in labels:

                axis.errorbar(x,results,yerr=errs,fmt='.',capsize=2,label = "$\lambda = $"+ str(lambda_),color = color)
                labels.append(label)
            else:
                axis.errorbar(x,results,yerr=errs,fmt='.',capsize=2,color = color)

            axis.plot(x,results,linestyle = 'dashed',color = color,linewidth = 0.75)
            plt.xlabel("$|x-y|$")
            plt.ylabel("$\langle( \phi(x) - \phi (y))^2\\rangle$")
            
            axis.spines[['right', 'top']].set_visible(False)
            
            #plt.title(" N = "+str(N) + " HMC = "+str(HMC))
            #analytic =np.array([ 0.0381944,0.0575829,0.0661302, 0.0686041])
        #print(lambda_,results/analytic)

    #analytic =np.array([0, 0.0381944,0.0575829,0.0661302, 0.0686041,0.0661302,0.0575829,0.0381944,0])

    #x = np.arange(0,N+1,1)
    #axis.plot(x,analytic,linestyle = 'dashed',color = 'black',linewidth = 0.75)
    #axis.plot(x,analytic,'o',color = 'black',markersize = 3,label = "Analytic")
    observable_name = "Difference multiple Lambda N"
    if HMC:
        observable_name = "Difference multiple Lambda HMC Multiple Volumes"
    
    file_name = "Results/Processed/Plots/"+observable_name+".pdf"

    plt.legend(loc = "upper right")

    #plt.xlim(0,N+2)
    #plt.ylim(0,0.1)
    plt.savefig(file_name)
    plt.show()
def plot_difference_lambdas(N,lambdas,N_measurements,N_thermalization,HMC=False):
    axis = plt.subplot(111)
    colors = ['red','blue','green','orange','purple'] + ['cyan','magenta','yellow','brown','pink','grey','olive','lime','teal','coral','navy','maroon','gold','indigo','tan','wheat','azure','ivory','lavender','plum','crimson','salmon','tomato','khaki','snow','orchid']

    for i,lambda_ in enumerate(lambdas):
        observable_name = "Difference"  
        file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
        if HMC:
            file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"       
        color = colors[i]
        vals = np.load(file_name)
        results = list(vals[0])#/np.abs(lambda_))
        errs = list(vals[1])#/np.abs(lambda_))
        results = [0] + results + np.flip(results[0:N//2-1]).tolist() + [0]
        errs = [0] + errs + np.flip(errs[0:N//2-1]).tolist() + [0]


        x = np.arange(0,N+1,1)
        print(x.shape,np.array(results).shape)
        print(lambda_,results)
        #results = np.roll(results,-1)-results 
        axis.errorbar(x,results,yerr=errs,fmt='.',capsize=2,label = "$\lambda = $"+ str(lambda_),color = color)
        axis.plot(x,results,linestyle = 'dashed',color = color,linewidth = 0.75)
        plt.xlabel("$|x-y|$")
        plt.ylabel("$\langle( f(x) - f (y))^2\\rangle$")
        
        axis.spines[['right', 'top']].set_visible(False)
        
        plt.title(" N = "+str(N) + " HMC = "+str(HMC))
        #analytic =np.array([ 0.0381944,0.0575829,0.0661302, 0.0686041])
        #print(lambda_,results/analytic)

    #analytic =np.array([0, 0.0381944,0.0575829,0.0661302, 0.0686041,0.0661302,0.0575829,0.0381944,0])

    #x = np.arange(0,N+1,1)
    #axis.plot(x,analytic,linestyle = 'dashed',color = 'black',linewidth = 0.75)
    #axis.plot(x,analytic,'o',color = 'black',markersize = 3,label = "Analytic")
    observable_name = "Difference multiple Lambda"
    if HMC:
        observable_name = "Difference multiple Lambda HMC"
    
    file_name = "Results/Processed/Plots/"+observable_name+".svg"

    plt.legend(loc = "upper right")

    #plt.xlim(0,N+2)
    #plt.ylim(0,0.1)
    plt.savefig(file_name)
    plt.show()
def plot_corr_lambdas(N,lambdas,N_measurements,N_thermalization,HMC=False):
    axis = plt.subplot(111)
    colors = ['red','blue','green','orange','purple'] + ['cyan','magenta','yellow','brown','pink','grey','olive','lime','teal','coral','navy','maroon','gold','indigo','tan','wheat','azure','ivory','lavender','plum','crimson','salmon','tomato','khaki','snow','orchid']

    for i,lambda_ in enumerate(lambdas):
        observable_name = "Correlation"  
        file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
        if HMC:
            file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"       
        color = colors[i]
        vals = np.load(file_name)
        results = list(vals[0])
        errs = list(vals[1])
        print(results)
        #results = results + np.flip(results).tolist()
        print(results)
        #errs = errs + np.flip(errs).tolist() 
        x = np.arange(0,N+1,1)
        #results = np.log(results[1:])
        #x =np.log( np.arange(1,8,1))
        #errs = errs[1:]
        print(x.shape,np.array(results).shape)
        print(lambda_,results)
        axis.errorbar(x,results,yerr=errs,fmt='.',capsize=2,label = "$\lambda = $"+ str(lambda_),color = color)
        axis.plot(x,results,linestyle = 'dashed',color = color,linewidth = 0.75)
        plt.xlabel("$|x-y|$")
        plt.ylabel("$\langle( \phi(x) * \phi (y))\\rangle$")
        
        axis.spines[['right', 'top']].set_visible(False)
        
        plt.title(" N = "+str(N) + " HMC = "+str(HMC))
       
    observable_name = "Correlation multiple Lambda"
    if HMC:
        observable_name = "Correlation multiple Lambda HMC"
    
    file_name = "Results/Processed/Plots/"+observable_name+".pdf"

    plt.legend(loc = "upper right")

    #plt.xlim(0,N+2)
    #plt.ylim(0,0.1)
    plt.savefig(file_name)
    plt.show()
def plot_corr_lambdas_fs(N,lambdas,N_measurements,N_thermalization,HMC=False):
    axis = plt.subplot(111)
    colors = ['red','blue','green','orange','purple'] + ['cyan','magenta','yellow','brown','pink','grey','olive','lime','teal','coral','navy','maroon','gold','indigo','tan','wheat','azure','ivory','lavender','plum','crimson','salmon','tomato','khaki','snow','orchid']

    for i,lambda_ in enumerate(lambdas):
        observable_name = "Correlation"  
        file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
        if HMC:
            file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"       
        color = colors[i]
        vals = np.load(file_name)
        results = list(vals[0])/np.abs(lambda_)
        errs = list(vals[1])/np.abs(lambda_)
        #results = results + np.flip(results).tolist()
        #errs = errs + np.flip(errs).tolist() 
        x = np.arange(0,N+1,1)
        #results = np.log(results[1:])
        #x =np.log( np.arange(1,8,1))
        #errs = errs[1:]
        print(x.shape,np.array(results).shape)
        print(lambda_,results)
        axis.errorbar(x,results,yerr=errs,fmt='.',capsize=2,label = "$\lambda = $"+ str(lambda_),color = color)
        axis.plot(x,results,linestyle = 'dashed',color = color,linewidth = 0.75)
        plt.xlabel("$|x-y|$")
        plt.ylabel("$\langle( f(x) * f (y))\\rangle$")
        
        axis.spines[['right', 'top']].set_visible(False)
        
        plt.title(" N = "+str(N) + " HMC = "+str(HMC))
       
    observable_name = "Correlation multiple Lambda"
    if HMC:
        observable_name = "Correlation multiple Lambda HMC"
    
    file_name = "Results/Processed/Plots/"+observable_name+".svg"

    plt.legend(loc = "upper right")

    #plt.xlim(0,N+2)
    #plt.ylim(0,0.1)
    plt.savefig(file_name)
    plt.show()
def plot_corr_lambdas_log(N,lambdas,N_measurements,N_thermalization,HMC=False):
    axis = plt.subplot(111)
    colors = ['red','blue','green','orange','purple'] + ['cyan','magenta','yellow','brown','pink','grey','olive','lime','teal','coral','navy','maroon','gold','indigo','tan','wheat','azure','ivory','lavender','plum','crimson','salmon','tomato','khaki','snow','orchid']

    for i,lambda_ in enumerate(lambdas):
        observable_name = "Correlation"  
        file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
        if HMC:
            file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"       
        color = colors[i]
        vals = np.load(file_name)
        results = list(vals[0])
        errs = list(vals[1])
        print(results)
        #results = results + np.flip(results).tolist()
        print(results)
        #errs = errs + np.flip(errs).tolist() 
        x = np.arange(0,N+1,1)
        x = np.log(x)
        errs = errs/np.array(results)
        results = np.log(results)
        
        #results = np.log(results[1:])
        #x =np.log( np.arange(1,8,1))
        #errs = errs[1:]
        print(x.shape,np.array(results).shape)
        print(lambda_,results)
        axis.errorbar(x,results,yerr=errs,fmt='.',capsize=2,label = "$\lambda = $"+ str(lambda_),color = color)
        axis.plot(x,results,linestyle = 'dashed',color = color,linewidth = 0.75)
        plt.xlabel("$\log{|x-y|}$")
        plt.ylabel("$\log{\langle( f(x) * f (y))}\\rangle$")
        
        axis.spines[['right', 'top']].set_visible(False)
        
        plt.title(" N = "+str(N) + " HMC = "+str(HMC))
       
    observable_name = "Correlation multiple Lambda log"
    if HMC:
        observable_name = "Correlation multiple Lambda HMC log" + " N = "+str(N) + " N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)
    
    file_name = "Results/Processed/Plots/"+observable_name+".pdf"

    plt.legend(loc = "upper right")

    #plt.xlim(0,N+2)
    #plt.ylim(0,0.1)
    plt.savefig(file_name)
    plt.show()

def plot_difference_corr_sq_lambdas(N,lambdas,N_measurements,N_thermalization,HMC=False):
    axis = plt.subplot(111)
    colors = ['red','blue','green','orange','purple'] + ['cyan','magenta','yellow','brown','pink','grey','olive','lime','teal','coral','navy','maroon','gold','indigo','tan','wheat','azure','ivory','lavender','plum','crimson','salmon','tomato','khaki','snow','orchid']

    for i,lambda_ in enumerate(lambdas):
        observable_name = "Difference"  
        file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
        if HMC:
            file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"       
        color = colors[i]
        vals = np.load(file_name)
        results = list(vals[0]/np.abs(lambda_)).copy()
        errs = list(vals[1]/np.abs(lambda_)).copy()
        results = [0] + results + np.flip(results[0:N//2-1]).tolist() + [0]
        errs = [0] + errs + np.flip(errs[0:N//2-1]).tolist() + [0]
        print(results)


        x = np.arange(0,N+1,1)
        print(x.shape,np.array(results).shape)
        print(lambda_,results)
        axis.errorbar(x,results,yerr=errs,fmt='.',capsize=2,label = "$\lambda = $"+ str(lambda_),color = color)
        axis.plot(x,results,linestyle = 'dashed',color = color,linewidth = 0.75)
        plt.xlabel("$|x-y|$")
        plt.ylabel("$\langle( f(x) - f (y))^2\\rangle$")
        observable_name = "Square"  
        file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
        if HMC:
            file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"       
        color = colors[i]
        values = np.load(file_name)
        results_sq = 2*np.array(list(values[0]/np.abs(lambda_))).copy()
        errs_sq =2* np.array(list(values[1]/np.abs(lambda_))).copy()
        observable_name = "Correlation"  
        file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
        if HMC:
            file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"       
        
        vals_corr = np.load(file_name)
        results_corr = 2*vals_corr[0].copy()/np.abs(lambda_)
        errs_corr = 2*vals_corr[1].copy()/np.abs(lambda_)
        result_total = results_sq-results_corr
        err_total = np.sqrt(errs_sq**2+errs_corr**2)
        errs = [0 for i in err_total]
        plt.errorbar(x,results_corr,yerr=errs,fmt='.',capsize=2,label = "$corr $"+ str(lambda_),color = "green")
        plt.errorbar(x,results_sq,yerr=errs,fmt='.',capsize=2,label = "$sq $"+ str(lambda_),color = "black")

        plt.errorbar(x,result_total,yerr=errs_sq,fmt='.',capsize=2,label = "$total $"+ str(lambda_),color = "blue")
        #print(np.abs(result_total-results))
        


        
        axis.spines[['right', 'top']].set_visible(False)
        
        plt.title(" N = "+str(N) + " HMC = "+str(HMC))
        #analytic =np.array([ 0.0381944,0.0575829,0.0661302, 0.0686041])
        #print(lambda_,results/analytic)

    #analytic =np.array([0, 0.0381944,0.0575829,0.0661302, 0.0686041,0.0661302,0.0575829,0.0381944,0])

    #x = np.arange(0,N+1,1)
    #axis.plot(x,analytic,linestyle = 'dashed',color = 'black',linewidth = 0.75)
    #axis.plot(x,analytic,'o',color = 'black',markersize = 3,label = "Analytic")
    observable_name = "Difference multiple Lambda"
    if HMC:
        observable_name = "Difference multiple Lambda HMC"
    
    file_name = "Results/Processed/Plots/"+observable_name+".svg"

    plt.legend(loc = "upper right")

    #plt.xlim(0,N+2)
    plt.ylim(0,0.1)
    plt.savefig(file_name)
    plt.show()
def plot_difference_corr_sq_lambdas_Ns(N,lambdas,N_measurementss,N_thermalizations,HMC=False):
    alphas = np.linspace(0.1,1,len(N_measurementss))
    axis = plt.subplot(111)
    colors = ['red','blue','green','orange','purple'] + ['cyan','magenta','yellow','brown','pink','grey','olive','lime','teal','coral','navy','maroon','gold','indigo','tan','wheat','azure','ivory','lavender','plum','crimson','salmon','tomato','khaki','snow','orchid']
    for j in range(len(N_measurementss)):
        
        N_measurements = N_measurementss[j]
        N_thermalization = N_thermalizations[j]
        for i,lambda_ in enumerate(lambdas):
            print("Hola")
            observable_name = "Difference"  
            file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
            if HMC:
                file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"       
            color = colors[i]
            vals = np.load(file_name)
            results = list(vals[0]/np.abs(lambda_)).copy()
            errs = list(vals[1]/np.abs(lambda_)).copy()
            results = [0] + results + np.flip(results[0:N//2-1]).tolist() + [0]
            errs = [0] + errs + np.flip(errs[0:N//2-1]).tolist() + [0]
            print(results)


            x = np.arange(0,N+1,1)
           
            #axis.errorbar(x,results,yerr=errs,fmt='.',capsize=2,label = "$\lambda = $"+ str(lambda_),color = color)
            #axis.plot(x,results,linestyle = 'dashed',color = color,linewidth = 0.75)
            plt.xlabel("$|x-y|$")
            plt.ylabel("$\langle( f(x) - f (y))^2\\rangle$")
            observable_name = "Square"  
            file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
            if HMC:
                file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"       
            color = colors[j]
            values = np.load(file_name)
            results_sq = 2*np.array(list(values[0]/np.abs(lambda_))).copy()
            errs_sq =2* np.array(list(values[1]/np.abs(lambda_))).copy()
            observable_name = "Correlation"  
            file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
            if HMC:
                file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"       
            
            vals_corr = np.load(file_name)
            results_corr = 2*vals_corr[0].copy()/np.abs(lambda_)
            errs_corr = 2*vals_corr[1].copy()/np.abs(lambda_)
            result_total = results_sq-results_corr
            print(np.abs(result_total-results))
            err_total = np.sqrt(errs_sq**2+errs_corr**2)
            errs = [0 for i in err_total]
            plt.errorbar(x,results_corr,yerr=errs_corr,fmt='.',capsize=2,label = f"$Corr {N_measurements} $"+ str(lambda_),color = color)
            plt.errorbar(x,results_sq,yerr=errs_sq,fmt='.',capsize=2,label = f"$sq {N_measurements}$"+ str(lambda_),color = color,alpha = 0.3)
            #plt.errorbar(x,result_total,yerr=errs_sq,fmt='.',capsize=2,label = "$total $"+ str(lambda_),color = "blue",alpha = alpha)
            #print(np.abs(result_total-results))
            


            
            axis.spines[['right', 'top']].set_visible(False)
            
            plt.title(" N = "+str(N) + " HMC = "+str(HMC))
        #analytic =np.array([ 0.0381944,0.0575829,0.0661302, 0.0686041])
        #print(lambda_,results/analytic)

    #analytic =np.array([0, 0.0381944,0.0575829,0.0661302, 0.0686041,0.0661302,0.0575829,0.0381944,0])

    #x = np.arange(0,N+1,1)
    #axis.plot(x,analytic,linestyle = 'dashed',color = 'black',linewidth = 0.75)
    #axis.plot(x,analytic,'o',color = 'black',markersize = 3,label = "Analytic")
    observable_name = "Difference multiple Lambda"
    if HMC:
        observable_name = "Difference multiple Lambda HMC"
    
    file_name = "Results/Processed/Plots/"+observable_name+".svg"

    plt.legend(loc = "upper right")

    #plt.xlim(0,N+2)
    #plt.ylim(0,0.1)
    plt.savefig(file_name)
    plt.show()
def plot_difference_lambdas_volumes_phis_corr(Ns,lambdas,N_measurements,N_thermalization,HMC=False):
    axis = plt.subplot(111)
    colors = ['red','blue','green','orange','purple'] + ['cyan','magenta','yellow','brown','pink','grey','olive','lime','teal','coral','navy','maroon','gold','indigo','tan','wheat','azure','ivory','lavender','plum','crimson','salmon','tomato','khaki','snow','orchid']
    labels = []
    for j in range(len(Ns)):
        N = Ns[j]
        for i,lambda_ in enumerate(lambdas):
            observable_name = "Difference"  
            file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
            if HMC:
                file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"       
            color = colors[i]
            vals = np.load(file_name)
            results = list(vals[0])#/np.abs(lambda_))
            errs = list(vals[1])#/np.abs(lambda_))
            results = [0] + results + np.flip(results[0:N//2-1]).tolist() + [0]
            errs = [0] + errs + np.flip(errs[0:N//2-1]).tolist() + [0]


            x = np.arange(0,N+1,1)
            print(x.shape,np.array(results).shape)
            print(lambda_,results)
            #errs = errs/np.array(results)
            results =-0.5* (np.array(results) -np.max(results))
            errs = 0.5*np.array(errs)
            

            #results = np.log(results)
            #x = np.log(x)
            label = "$\lambda = $"+ str(lambda_)
            if label not in labels:

                axis.errorbar(x,results,yerr=errs,fmt='.',capsize=2,label = "$\lambda = $"+ str(lambda_),color = color)
                labels.append(label)
            else:
                axis.errorbar(x,results,yerr=errs,fmt='.',capsize=2,color = color)

            axis.plot(x,results,linestyle = 'dashed',color = color,linewidth = 0.75)
            
            
            plt.xlabel("$|x-y|$")
            plt.ylabel("$\langle \phi(x)  \phi (y)\\rangle$")
            
            axis.spines[['right', 'top']].set_visible(False)
            
            #plt.title(" N = "+str(N) + " HMC = "+str(HMC))
            #analytic =np.array([ 0.0381944,0.0575829,0.0661302, 0.0686041])
        #print(lambda_,results/analytic)

    #analytic =np.array([0, 0.0381944,0.0575829,0.0661302, 0.0686041,0.0661302,0.0575829,0.0381944,0])

    #x = np.arange(0,N+1,1)
    #axis.plot(x,analytic,linestyle = 'dashed',color = 'black',linewidth = 0.75)
    #axis.plot(x,analytic,'o',color = 'black',markersize = 3,label = "Analytic")
    observable_name = "corr multiple Lambda N"
    if HMC:
        observable_name = "corr Lambda HMC Multiple Volumes"
    
    file_name = "Results/Processed/Plots/"+observable_name+".pdf"

    plt.legend(loc = "upper right")

    #plt.xlim(0,N+2)
    #plt.ylim(0,0.1)
    plt.savefig(file_name)
    plt.show()
def plot_difference_lambdas_volumes_phis_corr_fit(Ns,lambdas,N_measurements,N_thermalization,HMC=False):
    axis = plt.subplot(111)
    colors = ['red','blue','green','orange','purple'] + ['cyan','magenta','yellow','brown','pink','grey','olive','lime','teal','coral','navy','maroon','gold','indigo','tan','wheat','azure','ivory','lavender','plum','crimson','salmon','tomato','khaki','snow','orchid']
    labels = []
    for j in range(len(Ns)):
        N = Ns[j]
        for i,lambda_ in enumerate(lambdas):
            observable_name = "Difference"  
            file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+'.npy'
            if HMC:
                file_name = "Results/Processed/"+observable_name+"/"+observable_name+" lambda = " + str(lambda_) + " N = " + str(N) +" N measurements = "  + str(N_measurements)+" N Thermal = "  + str(N_thermalization)+" HMC.npy"       
            color = colors[i]
            vals = np.load(file_name)
            results = list(vals[0])#/np.abs(lambda_))
            errs = list(vals[1])#/np.abs(lambda_))
            results = [0] + results + np.flip(results[0:N//2-1]).tolist() + [0]
            errs = [0] + errs + np.flip(errs[0:N//2-1]).tolist() + [0]


            x = np.arange(0,N+1,1)
            print(x.shape,np.array(results).shape)
            print(lambda_,results)
            #errs = errs/np.array(results)
            results =-0.5* (np.array(results) -np.max(results))*x**2
            errs = 0.5*np.array(errs)*x**2
            

            #results = np.log(results)
            #x = np.log(x)
            label = "$\lambda = $"+ str(lambda_)
            if label not in labels:

                axis.errorbar(x[:N//2+1],results[:N//2+1],yerr=errs[:N//2+1],fmt='.',capsize=2,label = "$\lambda = $"+ str(lambda_),color = color)
                labels.append(label)
            else:
                axis.errorbar(x[:N//2+1],results[:N//2+1],yerr=errs[:N//2+1],fmt='.',capsize=2,color = color)

            axis.plot(x[:N//2+1],results[:N//2+1],linestyle = 'dashed',color = color,linewidth = 0.75)
            def model(x,a,m):
                return a*x*kn(1,m*x)
            min_x = 2
            popt,pcov = curve_fit(model,x[min_x:N//2+1],results[min_x:N//2+1],sigma=errs[min_x:N//2+1])
            chisq = np.sum((results[min_x:N//2+1]-model(x[min_x:N//2+1],*popt))**2/errs[min_x:N//2+1]**2)/(len(results[min_x:N//2+1])-2)
            print(popt,chisq)
            xs = np.linspace(min_x,N//2,100)
            err_a = np.sqrt(pcov[0,0])
            err_m = np.sqrt(pcov[1,1])
            axis.plot(xs,model(xs,*popt),color = color)
            axis.fill_between(xs,model(xs,popt[0]+err_a,popt[1]+err_m),model(xs,popt[0]-err_a,popt[1]-err_m),alpha = 0.2,color = color)
            #plt.plot(x[1:N//2+1],model(x[1:N//2+1],*popt),color = color)

            print(popt[1])
            plt.xlabel("$|x-y|$")
            plt.ylabel("$|x-y|^2\langle \phi(x)  \phi (y)\\rangle$")
            
            axis.spines[['right', 'top']].set_visible(False)
            
            #plt.title(" N = "+str(N) + " HMC = "+str(HMC))
            #analytic =np.array([ 0.0381944,0.0575829,0.0661302, 0.0686041])
        #print(lambda_,results/analytic)

    #analytic =np.array([0, 0.0381944,0.0575829,0.0661302, 0.0686041,0.0661302,0.0575829,0.0381944,0])

    #x = np.arange(0,N+1,1)
    #axis.plot(x,analytic,linestyle = 'dashed',color = 'black',linewidth = 0.75)
    #axis.plot(x,analytic,'o',color = 'black',markersize = 3,label = "Analytic")
    observable_name = "Difference multiple Lambda N"
    if HMC:
        observable_name = "Difference multiple Lambda HMC Multiple Volumes"
    
    file_name = "Results/Processed/Plots/"+observable_name+".svg"

    plt.legend(loc = "upper right")

    #plt.xlim(0,N+2)
    #plt.ylim(0,0.1)
    plt.savefig(file_name)
    plt.show()