from __future__ import division
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from scipy.linalg import cholesky
from matplotlib import colors
import matplotlib.pyplot as plt

def assembleKE(nu):
    k = np.array([1/2 - nu/6, 1/8 + nu/8, -1/4 - nu/12, -1/8 + 3*nu/8, -1/4 + nu/12, -1/8 - nu/8, nu/6, 1/8 - 3*nu/8])
    KE = 1/(1 - pow(nu, 2)) * np.array([
        [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
        [k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
        [k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
        [k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
        [k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
        [k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
        [k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
        [k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]]]);
    return KE

def performFEA(KE, rho, p, E):
    edofMat = np.zeros((nElements, 8), dtype = int)
    for elx in range(nelx):
        for ely in range(nely):
            el = ely + elx * nely
            n1 = (nely + 1) * elx + ely
            n2 = (nely + 1) * (elx + 1) + ely
            edofMat[el , :] = np.array([2*n1+2, 2*n1+3, 2*n2+2, 2*n2+3, 2*n2, 2*n2+1, 2*n1, 2*n1+1])

    # Construct the index pointers for the coo format
    iK = np.kron(edofMat, np.ones((8, 1))).flatten()
    jK = np.kron(edofMat, np.ones((1, 8))).flatten()
    simpE = rho ** p @ E

    sK = ((KE.flatten()[np.newaxis]).T * simpE.flatten()).flatten(order = 'F')
    K = coo_matrix((sK, (iK, jK)), shape = (nDOFS, nDOFS)).tocsc()

    totalDOFS = np.arange(nDOFS)
    fixedDOFS = np.union1d(totalDOFS[0:2 * (nely + 1):2], np.array([nDOFS - 1]))
    freeDOFS = np.setdiff1d(totalDOFS,fixedDOFS)

    # Solution and RHS vectors
    F = np.zeros((nDOFS, 1))
    U = np.zeros((nDOFS, 1))
    # Set load
    F[1,0] = -1

    # Remove constrained dofs from matrix
    K = K[freeDOFS , :][: , freeDOFS]

    # Solve system 
    U[freeDOFS, 0] = spsolve(K, F[freeDOFS, 0])
    return U, edofMat, simpE

def evaluateObjectiveGradient(rho):
    KE = assembleKE(nu)
    U, edofMat, simpE = performFEA(KE, rho, p, E)
    ce = np.sum((U[edofMat].reshape(nElements, 8) @ KE) * U[edofMat].reshape(nElements, 8), axis = 1)
    objC = (simpE.flatten() * ce).sum()
    gradC = (-1 * p * rho ** (p - 1)) * ce.reshape((nElements, 1)) * E.reshape((1, -1))

    v = np.zeros((nely + 2, nelx + 2, nMaterials - 1))
    v[1:nely + 1, 1:nelx + 1, :] = rho[ : , : nMaterials - 1].reshape(nely, nelx, nMaterials - 1)
    v[0, :, :] = v[1, :, :]
    v[nely + 1, :, :] = v[nely, :, :]
    v[:, 0, :] = v[:, 1, :]
    v[:, nelx + 1, :] = v[:, nelx, :]

    objP1 = (rho ** 2) * (rho - 1) ** 2
    objP2 = (v[1:nely + 1, 1:nelx + 1: ] - v[0:nely, 1:nelx + 1: ]) ** 2 + (v[1:nely + 1, 1:nelx+1, :] - v[1:nely + 1, 0:nelx, :]) ** 2

    objP = (1 / epsilon) * np.sum(objP1) + (epsilon / 8) * np.sum(objP2)
    gradP = (4 * rho ** 3 - 6 * rho ** 2 + 2 * rho) / epsilon

    G = gradC + alpha * gradP
    G = G - G[:, [nMaterials - 1]]
    return objC, objP, G

# Project gradients on appropriate spaces
def projectGradient(b, x):
    b = b.reshape((1, nMaterials))
    m = x.shape[1] - 1
    n = x.shape[0]
    res = np.inf
    mu = np.zeros((m, 1))
    zeros = np.zeros((n, 1))
    ones = np.ones((n, 1))

    while res > 1e-10:
        x_old = x.copy()
        diff = 1
        while diff > 1e-2:
            lambda_u = (1 / m) * np.maximum(np.sum(x[:, :m], axis=1).reshape(n, 1) - 1 - mu.sum(), zeros)
            lambda_l = (1 / m) * np.minimum(np.sum(x[:, :m], axis=1).reshape(n, 1) - mu.sum(), zeros)
            lambda_ = lambda_u + lambda_l
            mu_new = ((np.sum(x[:, :m] - np.tile(lambda_,m)) - b.flatten()[:m]) / n).reshape((m, 1))
            diff = np.linalg.norm(mu - mu_new, np.inf)
            mu = mu_new
        x[:, :m] = np.maximum(zeros, np.minimum(ones, x[:, :m] - np.tile(lambda_,m) - (np.tile(mu, n)).reshape(n, m)))
        res = np.linalg.norm(x - x_old, np.inf)
    return x

def filterDensity(rho, filterRadius):
    nfilter = int(nElements * ((2 *(np.ceil(filterRadius)-1)+1)**2))
    iH = np.zeros(nfilter)
    jH = np.zeros(nfilter)
    sH = np.zeros(nfilter)
    cc = 0
    for i in range(nelx):
        for j in range(nely):
            row = i * nely + j
            kk1 = int(np.maximum(i - (np.ceil(filterRadius)-1), 0))
            kk2 = int(np.minimum(i + np.ceil(filterRadius),nelx))
            ll1 = int(np.maximum(j - (np.ceil(filterRadius)-1), 0))
            ll2 = int(np.minimum(j + np.ceil(filterRadius),nely))
            for k in range(kk1, kk2):
                for l in range(ll1,ll2):
                    col = k * nely + l
                    fac = filterRadius - np.sqrt(((i-k)*(i-k)+(j-l)*(j-l)))
                    iH[cc] = row
                    jH[cc] = col
                    sH[cc] = np.maximum(0.0, fac)
                    cc = cc + 1
    # Finalize assembly and convert to csc format
    H = coo_matrix((sH, (iH, jH)), shape = (nElements, nElements)).tocsc()	
    Hs = H.sum(1)
    return np.asarray(H * rho[np.newaxis].T/Hs)[:,0]

# Main function to put everything together
def topoptPhaseField(nelx, nely, rho, maxIt):
    b = nElements * volFractions
    # Optimization loop
    for iter in range(maxIt):
        objC, objP, G = evaluateObjectiveGradient(rho)
        F = objC + alpha * objP
        rho = rho + stepSize * (projectGradient(b, rho - G) - rho)
        for i in range(nMaterials - 1):
            rho[:, i] = filterDensity(rho[:, i], filterRadius)
        rho[:, nMaterials - 1] = 1 - np.sum(rho[:, :nMaterials - 1], axis=1)
        print(f'Iteration: {iter + 1:5d} Objective: {objC:8.6f} Perimeter: {objP:8.6f}')

    # Initialize plot and plot the initial design
    plt.ion() # Ensure that redrawing is possible
    fig, ax = plt.subplots()
    im = ax.imshow(-rho[:, 0].reshape((nelx, nely)).T, cmap = 'gray', interpolation = 'none', norm = colors.Normalize(vmin = -1, vmax = 0))
    fig.show()
    fig.canvas.draw()
    im.set_array(-rho[:, 0].reshape((nelx, nely)).T)
    input("Press any key...")


# The main function  
if __name__ == "__main__":
    # Default input parameters
    nelx = 120
    nely = 40
    E = np.array([1.0, 1.0e-9]) #e
    volFractions = np.array([0.4, 0.6])
    nMaterials = len(volFractions) #p
    epsilon = 1.0
    stepSize = 0.3
    p = 2.0
    alpha = 0.125
    maxIt = 200
    nu = 0.3
    filterRadius = 2.0


    nElements = nelx * nely
    nDOFS = 2 * (nelx + 1) * (nely + 1)
    # Initialize design densities rho with shape of (nElements, nMaterials)
    rho = np.ones((nElements, nMaterials)) @ np.tile(volFractions, (len(volFractions), 1)) / len(volFractions)
    
    import sys
    if len(sys.argv)>1: nelx    =int(sys.argv[1])
    if len(sys.argv)>2: nely    =int(sys.argv[2])
    if len(sys.argv)>3: epsilon =int(sys.argv[3])
    if len(sys.argv)>4: maxIt   =int(sys.argv[4])

    topoptPhaseField(nelx, nely, rho, maxIt)
