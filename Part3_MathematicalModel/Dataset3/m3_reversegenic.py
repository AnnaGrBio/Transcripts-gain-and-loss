import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from ete3 import Tree
from ete3 import TreeStyle
from scipy.optimize import minimize

#################
################# number of populations that are investigated (hard coded!)
#################
no_pops = 7
add_info = 1

#################
################# read data ---------- dataset 1
#################
df = pd.read_csv("data_m3.txt", sep=",", header=None, engine='python', names=['col' + str(x) for x in range(no_pops+add_info+1) ])

#################
################# transform data to de-novo counts per pop
#################
names = ['col' + str(x+1) for x in range(no_pops+add_info)]

################# remove certain de-novo type from data
ind = [] 
for j in range(no_pops):
    ind = np.concatenate((ind,np.where(df[names[j+1]].eq('Intergenic').values==True)[0],np.where(df[names[j+1]].eq('Intronic').values==True)[0],np.where(df[names[j+1]].eq('ExonLonger').values==True)[0],np.where(df[names[j+1]].eq('NcRNA').values==True)[0],np.where(df[names[j+1]].eq('Pseudogene').values==True)[0]))

df = df.drop(ind)

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
    if len(params)%2 == 0:
        gain1 = (mean_transcripts - params[1]/params[2])*params[0]
        params = np.concatenate(([gain1],params))
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

loss_init = G12/(2*G1-G12)
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
init_params = [loss_init,gain_init/100,loss_init/100,ess_init]
bounds = [(0,1000),(0,20000),(0,1000),(0,ess_upp_bound)]

#################
#################
################# coalescent tree analysis (Baumdicker et al. (2010); Collins & Higgs (2012))
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
    result = 0
    # gain rate 1 is set, so that the mean number of transcripts per population is correctly estimated
    gain1 = max(0,(mean_transcripts - params[1]/params[2]-params[3])*params[0])
    if gain1 < 0:
        result = 10**9
    else:
        params_aux = np.concatenate(([gain1],params))
        theo_spec = theo_freq_spec_coal(params_aux)
        #
        result += np.sum((emp_spec-theo_spec)**2/theo_spec)
    #
    return(result)


#################
################# optimization algorithm (Sequential Least Squares Programming (SLSQP) minimization)
estimate = minimize(chisquare,init_params,method='SLSQP',bounds=bounds).x
# add gain rate to estimate
estimateFinal = np.concatenate(([(mean_transcripts - estimate[1]/estimate[2]-estimate[3])*estimate[0]],estimate))


#################
#################
################# fixed tree analysis (with least-squares: Collins & Higgs (2012), with likelihood: Baumdicker et al. (2012))
#################
################# load phylogenetic tree of D. melanogaster populations
t = Tree("newick_tree.txt",format=1)
nleaves = len(t.get_leaves())
no_pops = nleaves
nnodes = nleaves + (nleaves-1)

#################
################# create matrix for nodes with entries: time to ancestor, no of descendent leaves/tips, index of daughter 1, index of daughter 2
treeMatrix = np.zeros((nnodes,4))

index = 0
for node in t.traverse("postorder"):            # postorder goes backwards from leafs to root
    # if root -> no time to ancestor; otherwise record time to ancestor
    if (node.is_root() == False):
        treeMatrix[index,0] = node.dist
    #
    # if leaf -> no offspring
    if (node.is_leaf()==False):
        treeMatrix[index,1] = len(node.get_leaves())
        
        # getting names of direct descendants    
        j = 2
        for child in node.get_children():
            treeMatrix[index,j] = child.name    
            j = 3
    # if leaf -> set offspring values to -1 (simplifies identification of leaves)
    else:
        treeMatrix[index,2] = -1
        treeMatrix[index,3] = -1
    #
    node.name = index           # this simplifies cross referring to distances later
    #
    index += 1

#################
################# rescale branch length of tree to fit coalescent time scale
leaves = np.where(treeMatrix[:,1]==0)[0]
inter_times = np.zeros(nleaves-1)

inter_counter = 0                # counts number of recorded intercoalescent times
tM = np.copy(treeMatrix)         # copy of the tree matrix that can be adjusted!
while (inter_counter <= nleaves-2):
    # find smallest coalescent time among the leaves
    current_min = np.max(treeMatrix[:,0])
    for i in range(len(leaves)):
        current_min = min(current_min,tM[leaves[i],0])
    # record inter-coalescence time
    inter_times[inter_counter] = current_min
    # remove leaves from vector
    index_min = np.where(tM[:,0]==current_min)[0]
    ind_parent = np.concatenate((np.where(tM[:,2]==index_min[0])[0],np.where(tM[:,3]==index_min[0])[0]))
    # if minimal coalescent time involves the leaves
    if (len(index_min)>1): indices = index_min
    # otherwise we need to take the daughter nodes of the parent
    else: indices = [tM[ind_parent[0],2],tM[ind_parent[0],3]]
    # remove the leaves
    for i in range(2):
        leaves = np.delete(leaves,np.where(leaves==indices[i])[0][0])
    # adjust remaining branch lengths of the nodes left in the leaves array
    for i in range(len(leaves)):
        tM[leaves[i],0] -= current_min
    # add parent of last coalescent event to the leaves array
    leaves = np.concatenate((leaves,ind_parent))
    #
    inter_counter += 1

# reorder coalescent times from largest to smallest (in theory!) -- here: largest coalescence time is actually the one with 7 populations!
inter_times = inter_times[::-1]

################# scaling factor for branch lengths to transform to coalescent scale; WARNING: our tree is clearly not a coalescent! 
choose_array = np.zeros(len(inter_times))
for i in range(len(choose_array)):
    choose_array[i] = (i+2)*(i+1)/2

scaling_coal = (no_pops-1)/np.sum(inter_times*choose_array)

################# rescale branch lengths to coalescent scaling
treeMatrix_scaled_coal = np.copy(treeMatrix)
treeMatrix_scaled_coal[:,0] = treeMatrix[:,0]*scaling_coal

################ Alternative scaling because the tree is not a coalescent! (we simply divide by the longest branch length and later re-scale the rates)
### scaling the tree to simplify optimization algorithm! otherwise exp(-years) is too small and numerically unfeasible
scaling_alt = np.max(treeMatrix[:,0])
treeMatrix_scaled_alt = np.copy(treeMatrix)
treeMatrix_scaled_alt[:,0] = treeMatrix[:,0]/scaling_alt

#################
################# create frequency spectrum per node (following Collins & Higgs (2012) Eqs. (14-17)) 

################# loss function
def loss(t,params):
    return(1-np.exp(-params[1]*t))    

################# survival function
def surv(t,params):
    return(np.exp(-params[1]*t))      

################# create frequency spectrum per node   
def theo_freq_fixed(params,coal=False):
    theo_freq_spec_fixed = np.zeros(no_pops)
    #
    if (coal==True):
        TM = treeMatrix_scaled_coal
    else:
        TM = treeMatrix_scaled_alt
    #
    # loop over different classes of transcripts
    for l in range(int((len(params)-1)/2)):
        u = params[2*l]             # gain rate
        v = params[2*l + 1]         # loss rate
        #
        nNewDNTs = np.zeros(nnodes)
        # relative contribution to frequency distribution at each node, e.g. leaf node has contribution to only the single pop class, node with 2 offspring has
        # contribution to 0,1,2 classes
        aux_freq = np.zeros((nnodes, no_pops+1))
        #
        # important that treeMatrix has been constructed going backwards through the tree!
        for i in range(nnodes):
            # number of new DNTs at that node
            if (TM[i,0] > 0):
                nNewDNTs[i] = (u/v) * (1- np.exp(-v*TM[i,0]))
            else:   # root exception
                nNewDNTs[i] = u/v
            #
            # compute the auxiliary frequency distribution
            if (TM[i,1] == 0):      # leaf nodes only have contributions to k=1
                aux_freq[i,1] = 1
            else:       # internal nodes have contributions to k=0,...,k=no of offspring leaves
                child1 = int(TM[i,2])
                child2 = int(TM[i,3])
                #
                for k in range(int(TM[i,1]+1)):
                    if (k==0):
                        # loss of dnt along both subtrees (Eq. (17) in Collins & Higgs (2012))
                        aux_freq[i,0] = (loss(TM[child1,0],[u,v]) + surv(TM[child1,0],[u,v])*aux_freq[child1,0]) * (loss(TM[child2,0],[u,v]) + surv(TM[child2,0],[u,v])*aux_freq[child2,0])
                    #
                    else:       # go through terms in Eq. (16) in Collins & Higgs (2012), e.g. survival along child 1 times contribution to k and loss along child 2
                        term1 = surv(TM[child1,0],[u,v])*loss(TM[child2,0],[u,v])*aux_freq[child1,k]
                        term2 = surv(TM[child2,0],[u,v])*loss(TM[child1,0],[u,v])*aux_freq[child2,k]
                        term3 = 0
                        j = 0
                        while (j<=k):
                            term3 += surv(TM[child1,0],[u,v])*surv(TM[child2,0],[u,v])*aux_freq[child1,j]*aux_freq[child2,k-j]
                            j += 1
                        #
                        aux_freq[i,k] = term1 + term2 + term3
        # compute frequency distribution (Eq. (14) in Collins & Higgs (2012), Eq. (6) and adaptation for larger number comparisons in Baumdicker et al. (2012))
        for k in range(no_pops):
            for i in range(nnodes):
                theo_freq_spec_fixed[k] += nNewDNTs[i] * aux_freq[i,k+1]
            #print(theo_freq_spec_fixed[k])
    #
    # add essential transcripts
    theo_freq_spec_fixed[-1] += params[-1]
    #
    return(theo_freq_spec_fixed)


#################
################# chi-square like estimation -- Collins & Higgs (2012) 
def chisquare_fixed(params):
    coal=False                  # if tree does not resemble coalescent, use alternative scaling! -> only to simplify numerical estimation!
    #
    # fix gain rate so that mean number of transcripts is exactly met 
    gain1 = max(0,(mean_transcripts - params[1]/params[2] - params[-1])*params[0])
    params = np.concatenate(([gain1],params))
    #
    theo_spec = theo_freq_fixed(params,coal)
    #
    result = 0
    for k in range(len(theo_spec)):
        if (theo_spec[k]==0):
            result += 10**9
            break
        #
        else:
            result += (emp_spec[k]-theo_spec[k])**2/theo_spec[k]
    return(result)

#################
################# optimization algorithm (Sequential Least Squares Programming (SLSQP) minimization)
estimate_fixed = minimize(chisquare_fixed,init_params,method='SLSQP',bounds=bounds).x
# add gain rate to estimate
estimate_fixedFinal = np.concatenate(([(mean_transcripts - estimate_fixed[1]/estimate_fixed[2]-estimate_fixed[3])*estimate_fixed[0]],estimate_fixed))

#################
################# plot
plt.plot(np.arange(1,no_pops+1,1),theo_freq_fixed(estimate_fixedFinal,False),'o',markersize=10)
plt.plot(np.arange(1,no_pops+1,1),theo_freq_spec_coal(estimateFinal),'^',markersize=10, color="black")
plt.hist(dnt_presence,bins = np.arange(0.5,no_pops+0.6,1),color='C0',alpha=0.5)
plt.ylim((-10,8000))
plt.tick_params(axis='both', which='major', labelsize=12, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.show()