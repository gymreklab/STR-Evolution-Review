import numpy as np
from scipy.stats import geom
import matplotlib.pyplot as plt

########## Function definitions and imports ##############
def GetStepSizeProb(a1, a2, beta, p):
    """
    Get step size probability

    Parameters
    ----------
    a1 : int
        Allele mutating from
    a2: int
        Allele mutating to
    beta: float 
        Strength of mutation size directional bias
    p: float
        Parameterizes the mutation geometric step size distribution
    
    Returns
    -------
    step_size_prob: float
        Given that a1 mutates, the porbability that a1 mutates to a2
    """
    step_size = (a2-a1)
    if abs(step_size)>10: return 0
    up_prob = max([0.001,0.5*(1-beta*p*a1)])
    up_prob = min(up_prob, 0.999) # Maximum value is 0.99 (allow for minimal probability of contraction at small alleles)
    down_prob = 1-up_prob
    if step_size>0: dir_prob = up_prob
    else: dir_prob = down_prob
    step_prob = geom.pmf(abs(step_size), p)
    return dir_prob*step_prob

def GetTransitionMatrix(num_alleles, mu, beta, p, L):
    """
    Build transition matrix 

    Parameters
    ----------
    num_alleles : int
        Size of transition matrix to build. Centered at "0" (most optimal allele in the center)
    mu: float
        Per-generation mutation rate of the central allele
    beta: float 
        Strength of mutation size directional bias
    p: float
        Parameterizes the mutation geometric step size distribution
    L: float
        Length-dependent mutation rate parameter (slope of the increase of mutation rate with allele size)
    
    Returns
    -------
    transition_matrix: numpy.ndarray object
        Matrix of transition probabilities
    """

    # Initialize matrix (optimal=0)
    transition_matrix = np.zeros((num_alleles, num_alleles))

    # Fill in probability to transition from a1 to a2
    for i in range(num_alleles):
        for j in range(num_alleles):
            a1 = -1*int(num_alleles/2)+i
            a2 = -1*int(num_alleles/2)+j
            log_mu_prime = np.log10(mu)+L*a1 # Length-dependent mutation rate
            mu_prime = 10**log_mu_prime
            if mu_prime < 10**-8: mu_prime = 10**-8 
            if mu_prime > 10**-3: mu_prime = 10**-3
            if a1==a2: transition_matrix[i,j] = 1-mu_prime
            else:
                prob = GetStepSizeProb(a1, a2, beta, p)
                transition_matrix[i,j] = mu_prime*prob
        
    # Rescale each row to sum to 1 
    for i in range(num_alleles):
        rowsum = np.sum(transition_matrix[i,:])
        transition_matrix[i,:] = transition_matrix[i,:]/rowsum

    return transition_matrix

def Simulate(num_alleles, N_e, mu, beta, p, L, max_iter, end_samp_n, use_drift=True, set_start_equal=False):
    """Simulate allele frequencies

    Parameters
    ----------
    num_alleles : int
        Number of alleles included in simulation (must be odd)   
    N_e: int
        Effective population size
    mu: float
        Per-generation mutation rate of the central allele
    beta: float
        Strength of mutation size directional bias
    p: float
        Parameterizes the mutation geometric step size distribution
    L: float
        Length-dependent mutation rate parameter (slope of the increase of mutation rate with allele size)
    max_iter: int
        Number of generations to perform forward simulations
    end_samp_n: int
        Sample size for observed allele frequencies for last sampling step
        If set to 0, do not apply last sampling step
        Whether to perform multinomial sampling step at each generation
    set_start_equal: bool
        Whether to set starting vector of allele frequencies as equal
    
    Return
    ------
    allele_freqs: Simulated allele frequencies
    """
    print("Simulating mu=%s beta=%s p=%s L=%s"%(mu, beta, p, L))

    # Set the starting vector of allele frequencies 
    allele_freqs = np.zeros(num_alleles)
    
    if set_start_equal == True:
        # Set starting vector of allele frequencies: All alleles equal frequencies
        for i in range(num_alleles):
            allele_freqs[i] = 1/num_alleles
    else:
        # Set starting vector of allele frequencies: Optimal/ancestral allele frequency 1
        allele_freqs[int(num_alleles/2)] = 1

    # Get transition matrix (constant)
    transition_matrix = GetTransitionMatrix(num_alleles, mu, beta, p, L)
    
    # Transpose transition matrix
    transition_matrix_transpose = transition_matrix.transpose()

    t = 0
    while t < max_iter:
        # Apply mutation 
        allele_freqs = np.matmul(transition_matrix_transpose, allele_freqs)        
        if use_drift == True:
            # Use multinomial sampling
            allele_counts = np.random.multinomial(2*N_e, allele_freqs)
            # Rescale allele_freqs to sum to 1
            rowsum = np.sum(allele_counts)
            allele_freqs = allele_counts/rowsum
        t = t + 1
    
    # End sampling step
    # Use multinomial sampling on smaller sample size
    if end_samp_n > 0:
        allele_counts = np.random.multinomial(end_samp_n, allele_freqs)
        # Rescale allele_freqs to sum to 1
        rowsum = np.sum(allele_counts)
        allele_freqs = allele_counts/rowsum
        
    return allele_freqs

def PlotAfreqs(afreqs, fname=None):
    num_alleles = len(afreqs)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(list(range(-1*int(num_alleles/2), int(num_alleles/2)+1)), afreqs, color="orange", edgecolor="black")
    ax.set_xlabel("Allele (relative to opt.)")
    ax.set_ylabel("Frequency");
    if fname is not None: fig.savefig(fname)