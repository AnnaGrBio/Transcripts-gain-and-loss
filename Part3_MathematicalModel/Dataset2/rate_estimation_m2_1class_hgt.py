import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from ete3 import Tree
from scipy.optimize import minimize

#################
################# number of populations that are investigated (hard coded!)
#################
no_pops = 7
add_info = 1

#################
################# read data ---------- dataset 1
#################
df = pd.read_csv("data_m2.txt", sep=",", header=None, engine='python', names=['col' + str(x) for x in range(no_pops+add_info+1) ])

################# population names hard coded!
################# find indices of transcripts present in the different populations
################# produce an array where row = population and column = de-novo transcript number 
dnt_inds = np.zeros((no_pops,len(df)))
for j in range(len(df)):
    if ((df.iloc[j]=='AK5').any() == True):
        dnt_inds[0,j] = 1
    if ((df.iloc[j]=='DK5').any() == True):
        dnt_inds[1,j] = 1
    if ((df.iloc[j]=='GI5').any() == True):
        dnt_inds[2,j] = 1
    if ((df.iloc[j]=='SW5').any() == True):
        dnt_inds[3,j] = 1
    if ((df.iloc[j]=='UM').any() == True):
        dnt_inds[4,j] = 1 
    if ((df.iloc[j]=='YE').any() == True):
        dnt_inds[5,j] = 1
    if ((df.iloc[j]=='Zamb').any() == True):
        dnt_inds[6,j] = 1


#################
################# get empirical frequency spectrum
dnt_presence = np.zeros(len(df))
for j in range(len(df)):
    dnt_presence[j] = int(np.sum(dnt_inds[:,j]))

emp_spec = np.histogram(dnt_presence,bins = no_pops)[0]

#################
################# compute mean number of transcripts per population
mean_transcripts = 0
for j in range(no_pops):
    mean_transcripts += np.sum(dnt_inds[j,:])/no_pops

#################
################# compute number of transcripts per population for a parameter set
def mean(params):
    res = 0
    for i in range(int((len(params)-1)/2)):
        if (params[2*i + 1] != 0):
            res += params[2*i]/params[2*i +1]
    res += params[-1]
    return(res)

#################
################# determining initial estimates for gain and loss rate for the optimization procedure
################# according to Baumdicker, Hess & Pfaffelhuber (2012), Eq. (1) -> uses empirical gene frequency spectrum! (note that in Eqs. G12 is replaced by G12/2) 
G1 = 0
G12 = 0
for i in range(no_pops):
    G1 += (i+1)*emp_spec[i]/no_pops
    G12 += (i+1)*(no_pops-i-1)*emp_spec[i]/(no_pops*(no_pops-1))

loss_init = G12/(2*G1-G12) * 2
gain_init = G1*loss_init

####### setting some boundary values by hand
if (gain_init < 50):
    gain_init = 50
elif (gain_init > 10000):
    gain_init = 10000

if (loss_init < 0.1):
    loss_init = 0.1
elif (loss_init > 20):
    loss_init = 20

# setting essential transcript class + upper boundary
ess_init = emp_spec[-1]/2

if (ess_init>0):
    ess_upp_bound = emp_spec[-1]
else:
    ess_upp_bound = 0

# summarize initial parameters (2 classes of transcripts + essential transcripts)
init_params = [gain_init,loss_init,0.1,ess_init]
bounds = [(0,20000),(0,1000),(0,10),(0,ess_upp_bound)]

#################
#################
################# coalescent tree analysis (Baumdicker et al. (2010); Collins & Higgs (2012))
#################
################# theoretical gene frequency spectrum, e.g. Eq. (1) in Collins & Higgs (2012)
def theo_freq_spec_coal_hgt(params):
    n = no_pops
    gamma = params[-2]
    theo_spec_coal = np.zeros(no_pops)
    #
    for i in range(int((len(params)-2)/2)):
        theta = params[2*i]         # gain rate
        rho = params[2*i + 1]       # loss rate
        
        if (theta > 0 and rho > 0):
            prod_enum = 1
            prod_denom = 1
            #
            for k in range(n):
                prod_enum = prod_enum * (n-k)         # note that running index is k+1 in formula! -> n-(k+1)+1 = n-k 
                prod_denom = prod_denom * (n-(k+1)+rho)
                kfac = 1
                nfac = 1
                mfac = 1
                summ = 0
                for m in range(10000):
                    summ += gamma**m * kfac*(k+m) / (nfac * (n+rho+m) * mfac * (m+1))
                #
                theo_spec_coal[k] += ( theta/(k+1) ) * prod_enum /(prod_denom) * (1+summ)
    
    # add essential transcripts
    theo_spec_coal[-1] += params[-1]
    #
    return(theo_spec_coal)


#################
################# weighted least-squares estimation according to Eq. (24) in Collins & Higgs (2012)
def chisquare_hgt(params):
    gain1 = max(0,min(20000,(mean_transcripts - params[-1])*params[0]))
    #params_aux = np.concatenate(([gain1],params))
    theo_spec = theo_freq_spec_coal_hgt(params)
    
    result = np.sum((emp_spec-theo_spec)**2/theo_spec)
    
    return(result)


#################
################# optimization algorithm (Sequential Least Squares Programming (SLSQP) minimization)
estimate = minimize(chisquare_hgt,init_params,method='SLSQP',bounds=bounds).x

#################
#################
#################   without hgt
#################
#################
################# theoretical gene frequency spectrum, e.g. Eq. (1) in Collins & Higgs (2012)
def theo_freq_spec_coal(params):
    n = no_pops
    theo_spec_coal = np.zeros(no_pops)
    #
    for i in range(int((len(params)-1)/2)):
        theta = params[2*i]         # gain rate
        rho = params[2*i + 1]       # loss rate
        
        if (theta > 0 and rho > 0):
            prod_enum = 1
            prod_denom = 1
            #
            for k in range(n):
                prod_enum = prod_enum * (n-k)         # note that running index is k+1 in formula! -> n-(k+1)+1 = n-k 
                prod_denom = prod_denom * (n-(k+1)+rho)
                #
                theo_spec_coal[k] += ( theta/(k+1) ) * prod_enum / prod_denom
    
    # add essential transcripts
    theo_spec_coal[-1] += params[-1]
    #
    return(theo_spec_coal)


#################
################# weighted least-squares estimation according to Eq. (24) in Collins & Higgs (2012)
def chisquare(params):
    gain1 = max(0,min(20000,(mean_transcripts - params[-1])*params[0]))
    params_aux = np.concatenate(([gain1],params))
    theo_spec = theo_freq_spec_coal(params_aux)
    
    result = np.sum((emp_spec-theo_spec)**2/theo_spec)
    
    return(result)


#################
################# optimization algorithm (Sequential Least Squares Programming (SLSQP) minimization)
# summarize initial parameters (2 classes of transcripts + essential transcripts)
init_params = [loss_init,ess_init]
bounds = [(0,1000),(0,ess_upp_bound)]

estimate_wo = minimize(chisquare,init_params,method='SLSQP',bounds=bounds).x
# add gain rate to estimate
estimateFinal_wo = np.concatenate(([(mean_transcripts - estimate_wo[-1])*estimate_wo[0]],estimate_wo))

#################
################# plot
plt.plot(np.arange(1,no_pops+1,1),theo_freq_spec_coal(estimateFinal_wo),'o',markersize=10)
plt.plot(np.arange(1,no_pops+1,1),theo_freq_spec_coal_hgt(estimate),'^',markersize=10, color="black")
plt.hist(dnt_presence,bins = np.arange(0.5,no_pops+0.6,1),color='C0',alpha=0.5)
plt.ylim((-100,8000))
plt.tick_params(axis='both', which='major', labelsize=12, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.show()

