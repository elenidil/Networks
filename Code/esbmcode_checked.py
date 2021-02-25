import time
import random
import ctypes
import numpy as np
import pandas as pd
import numba
import math
from functools import partial
from tqdm import tqdm
from scipy import special
from scipy.special import betaln
#from spicy.special import gamma
import matplotlib.pyplot as plt


# define the functions
# optimization of scipy special betaln into numba
addr = numba.extending.get_cython_function_address('scipy.special.cython_special', 'betaln')
functype = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.c_double)
betaln_simple = functype(addr)


@numba.vectorize('float64(float64, float64)')
def nb_betaln(x, y):
    return betaln_simple(x, y)


addr = numba.extending.get_cython_function_address("scipy.special.cython_special", "gammaln")
functype = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
gammaln_simple = functype(addr)


@numba.vectorize('float64(float64)')
def nb_gammaln(x):
    return gammaln_simple(x)


def vec2mat(clust_lab):
    """ Put cluster labels in binary matrix form
    parameters: clust_lab, ndarray shape (V): clust_lab[v] = h iff node v in cluster h.
    returns: ndarray, shape (V, H): binary array M[v, h] = 1 iff node v in cluster h."""
    num_clusters = max(clust_lab)
    M = np.zeros((len(clust_lab), num_clusters), dtype=float)
    for v, cluster in enumerate(clust_lab):
        M[v, cluster - 1] = 1
    return M


@numba.jit(nopython=True)
def _compute_logp3(m_full, r_v, v_minus, a, b, log_beta_ab):
    # create log_p3
    # shape of term1 and term2 is HxH
    m_bar = np.outer(v_minus, v_minus) - np.diag(0.5 * v_minus * (v_minus + 1)) - m_full
    log_lhds_old = nb_betaln(m_full + a + r_v, m_bar + b + v_minus - r_v) - nb_betaln(m_full + a, m_bar + b)
    log_lhd_new = nb_betaln(r_v + a, v_minus - r_v + b) - log_beta_ab
    log_p3 = np.concatenate((log_lhds_old, log_lhd_new.reshape(1, -1)), axis=0).sum(axis=1)
    return log_p3


@numba.jit(nopython=True, cache=False)
def _update_cluster_sizes(Y, Z, m_full, v):
    if (Z.shape[1] > 1):
        nonempty_v = (Z.sum(axis=0) - Z[v, :]).nonzero()[0]
        # Reduce the dimensions of m_full and remove the empty columns from Z
        Z = Z[:, nonempty_v]
        m_full = m_full[nonempty_v][:, nonempty_v]

    # v_minus is number of nodes in each cluster excluding node v
    # H is the number of active clusters
    # r_v is number of edges from node v to each cluster (no self loops)
    # h_v are the clusters v has been assigned to (possibly empty)
    H = Z.shape[1]
    v_minus = Z.sum(axis=0) - Z[v, :]
    r_v = Z.T @ Y[:, v] - Z[v, :] * Y[v, v]
    h_v = Z[v].nonzero()[0]

    # Compute the m matrix by difference - (if h_v is empty this does nothing)
    if len(h_v) == 1:
        resid1 = np.zeros((H, H))
        resid1[:, h_v] = r_v.reshape(r_v.shape[0], 1)
        resid1[h_v] = r_v
        m_full -= resid1
    return Z, m_full, r_v, h_v, v_minus


@numba.jit(nopython=True, cache=False)
def _update_zeta_rv_mfull(Y, Z, r_v, m_full, new_sample, v, h_v):
    num_vertices = Z.shape[0]
    num_clusters = Z.shape[1]

    if new_sample < num_clusters:
        Z[v, new_sample] = 1
    else:  # todo check

        # Attach an extra row with zeros to Z and an extra
        # column *and* row with zeros to m_mat, this makes
        # initially shapes are Z: VxH, M: HxH
        # after updating shapes are Z: Vx(H+1), M:(H+1)x(H+1)
        # so there is an extra cluster
        Z = np.concatenate((Z, np.zeros((num_vertices, 1))), axis=1)
        Z[v, new_sample] = 1
        m_tmp = np.concatenate((m_full, np.zeros((num_clusters, 1))), axis=1)
        m_full = np.concatenate((m_tmp, np.zeros((1, num_clusters + 1))), axis=0)
        num_clusters += 1
        r_v = Z.T @ Y[:, v] - Z[v, :] * Y[v, v]
    # Updating m_full
    resid2 = np.zeros((num_clusters, num_clusters))
    resid2[:, new_sample] = r_v
    resid2[new_sample] = r_v
    m_full += resid2
    return Z, r_v, m_full


@numba.jit(nopython=True, cache=False)
def _get_probs(log_p2, log_p3, log_addit):
    # Continue computing probabilities
    # log_p2 and log_p3 have the same shape always
    # log_addit is a single integer
    log_p = log_addit + log_p2 + log_p3

    # length of probs is always (1 + H)
    probs = np.exp(log_p - log_p.max())
    probs = probs / probs.sum()
    return probs


@numba.jit(nopython=True, cache=False)
def logp2_DM(x, beta_DM, H_DM):
    # uses DM urn
    size = len(x)
    out = np.append(x + beta_DM, beta_DM * (H_DM - size) * (H_DM > size))
    out = np.log(out)
    return out
@numba.jit(nopython=True, cache=False)
def logp2_GN(x, gamma_GN):
    # uses GN urn
    size = len(x)
    out = np.append((x + 1) * (x.sum() - size + gamma_GN), (size**2) - (size * gamma_GN))
    out = np.log(out)
    return out
@numba.jit(nopython=True, cache=False)
def logp2_DP(x, alpha_PY):
    # uses DP urn
    out = np.append(x, alpha_PY)
    out = np.log(out)
    return out
@numba.jit(nopython=True, cache=False)
def logp2_PY(x, alpha_PY, sigma_PY):
    # uses PY urn
    size = len(x)
    out = np.append(x - sigma_PY, alpha_PY + size * sigma_PY)
    out = np.log(out)
    return out



def esbm(Y, seed, num_iter, prior, z_init=None, a=1, b=1, alpha_PY = None, sigma_PY = None,
         beta_DM = None, H_DM=None, gamma_GN=None, x=None, alpha_xi=None):
    """ Gibbs sampler for the extended stochastic block model

    Arguments:
        Y , ndarray: shape (V, V) symmetric adjacency matrix (boolean)
        z_init, ndarray: shape (V) int vector that holds cluster numbers
        a, int: a parameter of beta prior
        b, int: b parameter of beta prior
        prior, str: One of "DP", "PY", "DM", "GN"
        alpha_PY, sigma_PY: Parameters for PY prior
        beta_DM, H_DM: Parameters for DM prior
        gamma_GN: Parameter for GN prior
        x, ndarray: shape (V) vector of categorical covariates
        alpha_xi, ndarray: (C-vector of parameters for Dirichlet) must be
            specified along with x
    Returns:
        z_post, ndarray: posterior samples of community labels for each node
    """

    assert prior in ["DP", "PY", "DM", "GN"], "Invalid Input to esbm"

    # Selection of the prior distribution to be used
    # (different Urn schemes in article)
    if prior == "DP":
        assert alpha_PY is not None
        _compute_logp2 = partial(logp2_DP, alpha_PY = alpha_PY)

    elif prior == "PY":
        assert (alpha_PY is not None) and (sigma_PY is not None)
        _compute_logp2 = partial(logp2_PY, alpha_PY=alpha_PY, sigma_PY = sigma_PY)

    elif prior == "DM":
        assert (beta_DM is not None) and (H_DM is not None)
        _compute_logp2 = partial(logp2_DM, H_DM=H_DM, beta_DM=beta_DM)

    elif prior == "GN":
        assert gamma_GN is not None
        _compute_logp2 = partial(logp2_GN, gamma_GN=gamma_GN)

    # Pre-processing of the node attributes
    if x is not None:
        assert alpha_xi is not None, "alpha_xi must be set along with x"
        print("Covariates x have been provided")
        X = vec2mat(x)
        alpha0 = sum(alpha_xi)

    # Cluster assignments are encoded in 2 ways
    # 1) Z is a ndarray of shape (V, H) where H is the total number of clusters
    #    Z[v, h] = 1 iff node v is in cluster h (faster)
    #
    # 2) z_init is a vector of length V, with the cluster index of each v (compact)
    # vectors like that are created for all iterations and are packed in a
    # Vxnum_iter matrix z_post, s.t. z_post[v,t]=h if node v is in cluster h at
    # iteration t
    # note that Z is boolean, and Y is also boolean
    #use it for now, might extend it later
    assert (Y == Y.T).all(), "Y is not symmetric"
    Y = Y.astype(float)
    V = Y.shape[1]
    # initialize number of clusters and Z
    #if z_init is None:
    #    Z = np.eye(V, dtype=float)
    #else:
    #    Z = vec2mat(z_init)
    assert z_init is not None
    Z = vec2mat(z_init)
    z_post = []

    # Initialization for np.random.multinomial
    random.seed(seed)

    # m_full is initialized with shape H x H, matrix with block connections
    # note that the number of clusters, H <= V always
    m_full = Z.T @ Y @ Z - np.diag(0.5 * ((Y @ Z) * Z).sum(axis=0))

    log_addit = 0
    log_beta_ab = nb_betaln(a, b)
    for t in tqdm(range(0, num_iter)):
        for v in range(V):
            Z, m_full, r_v, h_v, v_minus = _update_cluster_sizes(Y, Z, m_full, v)
            # Compute the probabilities
            log_p3 = _compute_logp3(m_full, r_v, v_minus, a, b, log_beta_ab)
            #log_p2 = logp2(urn(v_minus)) #this makes the code slow
            log_p2 = _compute_logp2(v_minus)


            # Update log_addit, if x is None log_addit is always 0
            # todo check if I can make the code faster if putting this outside the esbm function

            if x is not None:
                Vx = np.dot(Z.T, X) - np.outer(Z[v, :], X[v, :])
                addit_old = (Vx[:, x[v] - 1] + alpha_xi[x[v] - 1]) / (v_minus + alpha0)
                addit_new = alpha_xi[x[v] - 1] / alpha0
                log_addit = np.log(np.append(addit_old, addit_new))

            probs = _get_probs(log_p2, log_p3, log_addit)
            # Sampling the indicator
            new_sample = random.choices(range(len(probs)), weights=probs)[0]
            if len(h_v) == 1:
                Z[v, h_v] = 0 # todo make sure this is correct
                # todo understand why this is needed only when length(h_v)=1, since h_v can be either empty or 1 or 2 or greater
            Z, r_v, m_full = _update_zeta_rv_mfull(Y, Z, r_v, m_full, new_sample, v, h_v)

        z_post.append(Z @ np.arange(1, Z.shape[1] + 1))
    return np.asarray(z_post).T


def expected_cl_py(n, sigma, theta, H):
    # COMPUTE THE DISTRIBUTION OF H AND THE EXPECTED VALUE UNDER VARIOUS GIBBS TYPE ####
    # takes n, sigma, theta and H and returns the expected value under various gibbs type
    n = int(n)
    assert (n > 0), "n > 0 is not True"
    assert (sigma >= 0), "sigma >= 0 is not True"
    assert (sigma < 1), "sigma < 1 is not True"
    assert (theta > -sigma), "theta > -sigma is not True"
    assert (H > 1), "H > 1 is not True"

    if (H == float("inf")):
        if (sigma == 0):
            out = theta * sum([(1 / (theta + k)) for k in range(n)])
        else:
            out = 1 / sigma * math.exp( nb_gammaln(theta + sigma + n) - nb_gammaln(theta + sigma) - nb_gammaln(theta + n)
                + nb_gammaln(theta + 1)) - theta / sigma
    elif (H < float("inf")):
        if (sigma == 0):
            # index = range(n)
            out = H - H * math.exp(sum([(math.log(k + theta * (1 - 1 / H)) - math.log(theta + k)) for k in range(n)]))

    return out


@numba.jit(nopython=True, cache=False)
def lower(matrix, diagonal): #diagonal takes value False when we don't want it to be included
    # the lower function takes the same time like R lowerTriangle for (100,100) matrix but takes
    # 10times the time of R for (10000,10000) matrix.
    low = []
    for j in range(matrix.shape[1]):
        for i in range(matrix.shape[0]):
            if (diagonal is False):
                if (i > j):
                    low.append(matrix[i,j])
            else:
                if (i >= j):
                    low.append(matrix[i,j])
    return low

#@numba.jit(nopython=True, cache=False)
def melt(matrix, rows, cols):
    #for size (10000,10000) this takes 3.5sec while R takes 2.5 sec
    assert isinstance(rows, np.ndarray)
    assert isinstance(cols, np.ndarray)
    assert isinstance(matrix, np.ndarray)
    X2 = np.repeat(rows, len(cols))
    X1 = np.tile(cols, len(rows))
    return np.stack((X1, X2, matrix.T.ravel())).T

#@numba.jit(nopython=True, cache=False)
def log_py_z(Y, z, a, b): # todo check that is the same with R where you set the colnames
    # takes an adjacency matrix Y, one vector of node labels z, hyperparameters (a,b) of Beta priors for block probabilities
    # returns the logarithm of the marginal likelihood in eq. [3] evaluated at z.
    H = len(set(z.tolist()))
    Y = np.asarray(Y)
    z = np.asarray(z)

    edge_counts = pd.DataFrame(melt(Y, z, z), columns=["X2", "X1", "value"])
    non_edge_counts = pd.DataFrame(melt(1 - Y, z, z), columns=["X2", "X1", "value"])

    edge_sum = edge_counts.groupby(['X2', 'X1'])['value'].sum()
    Edge = np.array(edge_sum).reshape((H,H), order = "F")
    np.fill_diagonal(Edge, np.diagonal(Edge) / 2)
    non_edge_sum = non_edge_counts.groupby(['X2', 'X1'])['value'].sum()
    No_Edge = np.array(non_edge_sum).reshape((H,H), order='F')
    np.fill_diagonal(No_Edge, np.diagonal(No_Edge) / 2)
    lower_Edge = lower(Edge, True)
    lower_No_Edge = lower(No_Edge, True)
    a_n = np.asarray(lower_Edge) + a       #a_n = [k + a for k in lower_Edge]
    b_bar_n = np.asarray(lower_No_Edge) + b    #b_bar_n = [l + b for l in lower_No_Edge]

    return sum(betaln(a_n, b_bar_n)) - (H * (H+1)/2) *betaln(a,b)


# GNEDIN
import math
def nCr(n, r):
    f = math.factorial
    return f(n) // (f(r) * f(n - r))


def HGnedin(V, h, gamma=0.5):
    return math.exp(math.log(nCr(V, h)) + nb_gammaln(h - gamma) - nb_gammaln(1 - gamma)
                    + math.log(gamma) + nb_gammaln(V + gamma - h) - nb_gammaln(V + gamma))

#------------------------------

#-------
#load the data (different Y matrices)

#c_true = pd.read_csv("C:\\Users\\elena\OneDrive\Έγγραφα\c_true.csv")
#Ymat = pd.read_csv(r"C:\\Users\\elena\OneDrive\Έγγραφα\Ymatrix.csv")
#Y = Ymat.drop("Unnamed: 0", axis = 1)
#c_true = c_true.drop("Unnamed: 0", axis = 1)
#Y = Y.to_numpy()


#new Y in github 80x80
#Ymat = pd.read_csv(r"C:\\Users\\elena\Desktop\ESBM\network_1_80.csv")
#Y = Ymat.drop("Unnamed: 0", axis = 1)
#c_true = c_true.drop("Unnamed: 0", axis = 1)
#Y = Y.to_numpy()

#new Y generated by R package
#Y_new = pd.read_csv(r"C:\\Users\\elena\Desktop\ESBM\Y200.csv")
Y_new = pd.read_csv(r"C:\\Users\\elena\Desktop\ESBM\Y100_85.csv")

Y1 = Y_new.drop("Unnamed: 0", axis = 1)
Y1 = Y1.to_numpy()
Y = Y1
# another Y generated from R with size 1000x1000
#Y10 = pd.read_csv(r"C:\\Users\\elena\Desktop\ESBM\Y10.xlsx")
#Y10 = Y10.drop('Unnamed: 0', axis=1)
#Y10 = Y10.to_numpy()

# Y is shape VxV, symmetric
num_iter = 20000
seed = 1
V = Y.shape[0]
print("Y is 100x100")
#Setting the hyperparameters
print("The expected value of H (non empty number of groups) under various Gibbs type priors")
#Dirichlet Multinomial
sigma_DM = 0
H_DM = 50
beta_DM = 3.5 / 50
print("DM", expected_cl_py(V, sigma = sigma_DM, theta = beta_DM*H_DM, H = H_DM))

#Dirichlet process
sigma_DP = 0
H_DP = float("inf")
alpha_DP = 3
print("DP", expected_cl_py(V, sigma = sigma_DP, theta = alpha_DP, H = H_DP))

#Pitman-Yor process
sigma_PY = 0.6
H_PY = float("inf")
alpha_PY = -0.3
print("PY", expected_cl_py(V, sigma = sigma_PY, theta = alpha_PY, H = H_PY))

#Gnedin process
gamma = 0.45
print("GN", sum([k * HGnedin(V, k, gamma = gamma) for k in range(1,V+1)] ) )

#Posterior computation via collapsed Gibbs sampler
z_init = list(range(1,(Y.shape[1]+1)))


#Dirichlet Multinomial
prior = "DM"
start = time.time()
Z_DM = esbm(Y=Y, seed=seed, num_iter=num_iter, z_init= z_init, a= 1, b= 1, prior=prior, beta_DM=3.5/50, H_DM=50)
end = time.time()
print("Time for esbm DM", end - start)

#Dirichlet Process
prior = "DP"
start = time.time()
Z_DP = esbm(Y=Y, seed=seed, num_iter=num_iter, z_init= z_init, a= 1, b= 1, prior=prior, alpha_PY=3, sigma_PY=0)
end = time.time()
print("Time for esbm DP", end - start)

#Pitman Yor Process
prior = "PY"
start = time.time()
Z_PY = esbm(Y=Y, seed=seed, num_iter=num_iter, z_init= z_init, a= 1, b= 1, prior=prior, alpha_PY=-0.3, sigma_PY=0.6)
end = time.time()
print("Time for esbm PY", end - start)

# Gnedin Process
prior = "GN"
start = time.time()
Z_GN = esbm(Y=Y, seed=seed, num_iter=num_iter, z_init= z_init, a= 1, b= 1, prior=prior, gamma_GN=0.45)
end = time.time()
print("Time for esbm GN", end - start)




#compute the logarithm of the marginal likelihoods under different priors

start = time.time()
l_y_DM = np.zeros(num_iter, dtype = object)
for t in tqdm(range(num_iter)):
    l_y_DM[t] = log_py_z(Y, Z_DM[:,t],1,1)
end = time.time()
print(end-start, 'Dirichlet Multinomial evaluate marginal likelihoods')
plt.plot(range(len(l_y_DM)),np.vstack(l_y_DM))
plt.title("Dirichlet Multinomial without node attributes")
plt.ylabel("l_y_DM")
plt.xlabel("x-axis")
plt.show()

start = time.time()
l_y_DP = np.zeros(num_iter, dtype = object)
for t in range(num_iter):
    l_y_DP[t] = log_py_z(Y, Z_DP[:,t],1,1)
end = time.time()
print(end-start, 'Dirichlet Process evaluate marginal likelihoods')
plt.plot(range(len(l_y_DP)),np.vstack(l_y_DP))
plt.title("Dirichlet Process without node attributes")
plt.ylabel("l_y_DP")
plt.xlabel("x-axis")
plt.show()

start = time.time()
l_y_PY = np.zeros(num_iter, dtype = object)
for t in range(num_iter):
    l_y_PY[t] = log_py_z(Y, Z_PY[:,t],1,1)
end = time.time()
print(end-start, 'Pitman Yor Process evaluate marginal likelihoods')
plt.plot(range(len(l_y_PY)),np.vstack(l_y_PY))
plt.title("Pitman Yor Process without node attributes")
plt.ylabel("l_y_PY")
plt.xlabel("x-axis")
plt.show()

start = time.time()
l_y_GN = np.zeros(num_iter, dtype = object)
for t in tqdm(range(num_iter)):
    l_y_GN[t] = log_py_z(Y, Z_GN[:,t],1,1)
end = time.time()
print(end-start, 'Gnedin Process evaluate marginal likelihoods')
plt.plot(range(len(l_y_GN)),np.vstack(l_y_GN))
plt.title("Gnedin Process without node attributes")
plt.ylabel("l_y_GN")
plt.xlabel("x-axis")
plt.show()

#posterior inference under ESBM
#perform estimation, uncertainty quantification and model selection for ESBM
#burn_in = 4000 #about 1/5 of N_iter
#z_0 = [1] * 20 + [2] * 20 + [3]*15 + [4]*15 + [5] *10
#before performing posterior inference, visualize the traceplots for the logarithms of the likelihood


#Now work with node attributes
N_iter = 20000
V = Y.shape[0]
seed = 1
my_z = z_init = list(range(1,(Y.shape[1]+1)))
#define the vector with node attributes
#my_x =  [1] *40 + [2] *40 + [3]*40 + [4]*40 + [5] *40
my_x =  [1] *20 + [2] *20 + [3]*20 + [4]*20 + [5] *20
my_alpha_xi = [1]*5

#Gnedin Process
my_prior = "GN"
start = time.time()
Z_GN_x = esbm(Y, seed = 1, num_iter= N_iter, prior=my_prior, z_init=my_z, a=1, b=1, gamma_GN=0.45, x = my_x, alpha_xi= my_alpha_xi)
end = time.time()
print(end-start, "Gnedin with node attributes - esbm") # this takes 576.95 sec - 10 mins for 50k num_iter

start = time.time()
l_y_GN_x = np.zeros(N_iter, dtype = object)
for t in range(N_iter):
    l_y_GN_x[t] = log_py_z(Y, Z_GN_x[:,t],1,1)
end = time.time()
print(end-start, "Gnedin with node attributes - evaluate marginal likelihoods")  # takes 294.26 sec for 50k = n_iter

plt.plot(range(len(l_y_GN_x)),np.vstack(l_y_GN_x))
plt.title("Gnedin Process with node attributes")
plt.ylabel("l_y_GN_x")
plt.xlabel("x-axis")
plt.show()

#the above traceplots confirm that the Gibbs sampler has satisfactory mixing and rapid convergence.
#Due to the stability of the chains for the quantity in Eq. [1] we can reliably compute the logarithm
#of the marginal likelihoods for the different priors and models via the harmonic mean in Eq.[14]
burn_in = 5000
#Dirichlet Multinomial
l_y_DM_burn = l_y_DM[burn_in:num_iter]
neg_l_y_DM = -l_y_DM_burn
extra = np.array(neg_l_y_DM - max(neg_l_y_DM), dtype = np.float32)
l_y_post_DM = math.log(len(l_y_DM_burn)) - max(neg_l_y_DM) - math.log(sum(np.exp(extra)))
print("DM", l_y_post_DM)

#Dirichlet Process
l_y_DP_burn = l_y_DP[burn_in:num_iter]
neg_l_y_DP = -l_y_DP_burn
extra = np.array(neg_l_y_DP - max(neg_l_y_DP), dtype = np.float32)
l_y_post_DP = math.log(len(l_y_DP_burn)) - max(neg_l_y_DP) - math.log(sum(np.exp(extra)))
print("DP", l_y_post_DP)

#Pitman Yor Process
l_y_PY_burn = l_y_PY[burn_in:num_iter]
neg_l_y_PY = -l_y_PY_burn
extra = np.array(neg_l_y_PY - max(neg_l_y_PY), dtype = np.float32)
l_y_post_PY = math.log(len(l_y_PY_burn)) - max(neg_l_y_PY) - math.log(sum(np.exp(extra)))
print("PY", l_y_post_PY)

#Gnedin Process
l_y_GN_burn = l_y_GN[burn_in:num_iter]
neg_l_y_GN = -l_y_GN_burn
extra = np.array(neg_l_y_GN - max(neg_l_y_GN), dtype = np.float32)
l_y_post_GN = math.log(len(l_y_GN_burn)) - max(neg_l_y_GN) - math.log(sum(np.exp(extra)))
print("GN", l_y_post_GN)


import csv
#Z_DM_df = pd.DataFrame(Z_DM)
#numpy.savetxt("Z_DM.csv", Z_DM, delimiter = ",")
#pd.DataFrame(np_array).to_csv("path/to/file.csv")
pd.DataFrame(Z_DM).to_csv("C:\\Users\\elena\Desktop\ESBM\python exported data\Z_DM_py.csv")
pd.DataFrame(Z_DP).to_csv("C:\\Users\\elena\Desktop\ESBM\python exported data\Z_DP_py.csv")
pd.DataFrame(Z_PY).to_csv("C:\\Users\\elena\Desktop\ESBM\python exported data\Z_PY_py.csv")
pd.DataFrame(Z_GN).to_csv("C:\\Users\\elena\Desktop\ESBM\python exported data\Z_GN_py.csv")
pd.DataFrame(Z_GN_x).to_csv("C:\\Users\\elena\Desktop\ESBM\python exported data\Z_GN_x_py.csv")

