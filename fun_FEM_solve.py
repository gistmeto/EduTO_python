import numpy as np
from scipy.sparse import csr_matrix, coo_matrix
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

import warnings
warnings.simplefilter("ignore")

def fun_FEM_solve(fem, opt):
    # Initializing FEM Variable
    K = coo_matrix((fem["ndof"], fem["ndof"]))
    F = np.zeros(fem["ndof"])

    # Build K Matrix
    
    kelist = np.kron(fem["Emin"] + (fem["E0"] - fem["Emin"]) * opt["erho"] ** opt["penal"], np.ones(144)).T * fem["K_S"]

    K = coo_matrix((kelist, (fem["is"]-1, fem["js"]-1)), shape=(fem["ndof"], fem["ndof"])).tocsc()

    
    # Apply boundary conditions
    K[fem["bcdof"]-1, :] = 0

    K[:, fem["bcdof"]-1] = 0

    
    for i in range(len(fem["bcdof"])):
        K[fem["bcdof"][i]-1, fem["bcdof"][i]-1] = 1

    
    F[fem["bcdof"]-1] = fem["bcval"]

    # Apply traction force
    for i in range(fem["fdof"].shape[0]):
        F[fem["fdof"][i]-1] += fem["fval"][i]
    fem["F"]=F
    # Solve
    U = spsolve(K, F)[:,np.newaxis]

    return (U, K)