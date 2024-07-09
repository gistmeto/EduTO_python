import numpy as np
from scipy.sparse.linalg import splu
import scipy.sparse as sp
from scipy.sparse import coo_matrix, csr_matrix, csc_matrix, lil_matrix

def fun_Pre_optinit(input, fem):
    opt = {}
    opt["f"] = []
    opt["g"] = []
    # Separation of design and non-design domains
    nd_ele = np.where(fem["IX"][:, 4] == 2)[0]
    opt["dof_nd"] = np.unique(fem["IX"][nd_ele, :4].flatten())
    opt["dof_dd"] = np.setdiff1d(np.arange(1, fem["nn"] + 1), opt["dof_nd"]).flatten()

    # Setting problem parameters
    opt["VT"] = input["VT"]
    opt["volfrac"] = input["volfrac"]
    opt["penal"] = input["penal"]

    # Setting of optimization variables
    opt["dv"] = np.ones((len(opt["dof_dd"]),1)) * input["initdv"]
    opt["nv"] = np.zeros((np.max(opt["dof_dd"]),1))
    opt["nv"][opt["dof_nd"]-1] = 1
    opt["nv"][opt["dof_dd"]-1] = opt["dv"]

    opt["dvold"] = opt["dv"].copy()
    opt["dvolder"] = opt["dv"].copy()
    opt["dvmin"] = opt["dv"] * 0 - 1
    opt["dvmax"] = opt["dv"] * 0 + 1

    opt["iter"] = 1
    opt["deltaf"] = 1.0
    opt["cont_sw"] = 0.
    opt["bt"] = input["bt"]

    # Setting of MMA variables
    MMA = {}
    MMA["a0"] = 1
    MMA["a"] = np.array([[0.]])
    MMA["c"] = input["MMA_c"]
    MMA["d"] = np.array([1.])
    MMA["low"] = opt["dvmin"]
    MMA["upp"] = opt["dvmax"]
    opt["MMA"]=MMA

    # Build filter element stiffness (Kft) and Transformation matrix (Tft)
    nx = np.transpose(np.vstack([fem["X"][fem["IX"][:,0]-1,[0]],
                                 fem["X"][fem["IX"][:,1]-1,[0]],
                                 fem["X"][fem["IX"][:,2]-1,[0]],
                                 fem["X"][fem["IX"][:,3]-1,[0]]]))
    ny = np.transpose(np.vstack([fem["X"][fem["IX"][:,0]-1,[1]],
                                 fem["X"][fem["IX"][:,1]-1,[1]],
                                 fem["X"][fem["IX"][:,2]-1,[1]],
                                 fem["X"][fem["IX"][:,3]-1,[1]]]))
    nz = np.transpose(np.vstack([fem["X"][fem["IX"][:,0]-1,[2]],
                                 fem["X"][fem["IX"][:,1]-1,[2]],
                                 fem["X"][fem["IX"][:,2]-1,[2]],
                                 fem["X"][fem["IX"][:,3]-1,[2]]]))
    
    Kd = (input["rmin"] / (2 * np.sqrt(3)))**2 * np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    NN = np.array([[2, 1, 1, 1], [1, 2, 1, 1], [1, 1, 2, 1], [1, 1, 1, 2]]) / 12
    isf = np.reshape(np.kron(fem["IX"][:, 0:4], np.ones((1, 4), dtype=int)), (16 * fem["ne"], ))
    jsf = np.reshape(np.kron(fem["IX"][:, 0:4], np.ones((4, 1), dtype=int)), (16 * fem["ne"], ))
    
    # Initialization
    Kft = np.zeros((16 * fem["ne"],))
    Tft = np.zeros((16 * fem["ne"],))


    for e in range(fem["ne"]):
        px = np.array([nx[e, 0], nx[e, 1], nx[e, 2], nx[e, 3]]) 
        py = np.array([ny[e, 0], ny[e, 1], ny[e, 2], ny[e, 3]])
        pz = np.array([nz[e, 0], nz[e, 1], nz[e, 2], nz[e, 3]])
    
        Ve = fem["Ve"][e]
        B = np.zeros((3, 4))
    
        for i in range(4):
            ind = np.arange(4)
            ind = np.delete(ind, i)  # Remove i-th index
            pm = np.array([1, -1, 1, -1])
        
            beta  = -pm[i] * np.linalg.det(np.vstack((np.ones(3), py[ind], pz[ind])))
            gamma = pm[i] * np.linalg.det(np.vstack((np.ones(3), px[ind], pz[ind])))
            delta = -pm[i] * np.linalg.det(np.vstack((np.ones(3), px[ind], py[ind])))
        
            B[:, i] = 1 / (6 * Ve) * np.array([beta, gamma, delta])

        Ke = (np.matmul(np.matmul(B.T, Kd), B) + NN) * Ve 
        Kft[e * 16:(e + 1) * 16] = np.reshape(Ke,(16,)) 
        Tft[e * 16:(e + 1) * 16] = np.reshape(NN * Ve,(16,))


    # LU decomposition (sparse)
    Kft_sparse = coo_matrix((Kft, (isf-1, jsf-1))).tocsc()
    opt["Kft_sparse"] = Kft_sparse
    LU = splu(Kft_sparse, permc_spec="NATURAL")
    opt["lu_L_Kft"] = LU.L
    opt["lu_U_Kft"] = LU.U
    opt["Tft"] = coo_matrix((Tft, (isf-1, jsf-1))).tocsc()
    # Build matrix for transformation from nodal to element density (Ten)
    opt["Ten"] = lil_matrix((fem["ne"], fem["nn"]))

    for e in range(fem["ne"]):
        opt["Ten"][e, fem["IX"][e, 0:4]-1] = 1 / 4
    
    return (opt,MMA)