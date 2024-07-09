import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

def fun_OPT_fgdfdg(fem,opt):

    opt["f"].append(np.dot(sp.csc_matrix.dot(fem["U"].transpose(),fem["K"]),fem["U"]).squeeze())
    opt["g"].append((np.dot(fem["Ve"], opt["erho"]) / opt["VT"]) - opt["volfrac"])
    # Derivatives
    DerhoDnrho = opt["Ten"].T
    DnrhoDfdv = (1/np.cosh(opt["bt"] * opt["fdv"]))**2 * opt["bt"] / (2 * np.tanh(opt["bt"]))  

    # Evaluation of dfdx
    dkelist = (np.kron((fem["E0"] - fem["Emin"]) * opt["penal"] * np.power(opt["erho"], opt["penal"] - 1), np.ones(144)).T * fem["K_S"]).reshape(1,-1)
    ulist1 = np.reshape(np.kron(np.transpose(fem["U"][fem["edof"]-1])[0], np.ones((1,12))).T, (fem["ne"] * 144, 1)).T
    ulist2 = np.reshape(np.kron(np.transpose(fem["U"][fem["edof"]-1])[0], np.ones((12,1))).T, (fem["ne"] * 144, 1)).T
    DfDerho = -np.sum(np.reshape((ulist1 * dkelist * ulist2), (fem["ne"],144)).T, axis=0).T
    DfDdv = DerhoDnrho.dot(DfDerho) * DnrhoDfdv
    dfdx = DfDdv[opt["dof_dd"]-1][:, np.newaxis]

    # Evaluation of dgdx
    DgDerho = fem["Ve"].T / opt["VT"]
    DgDdv = DerhoDnrho.dot(DgDerho) * DnrhoDfdv
    dgdx = DgDdv[opt["dof_dd"]-1]
    dgdx = dgdx.T[np.newaxis, :]

    # Storing results in dictionary
    
    return (opt["f"],opt["g"],dfdx, dgdx)

