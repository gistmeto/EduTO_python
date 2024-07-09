import numpy as np

def fun_Pre_feminit(input,msh):
    # Separation of FEM domain and exterior domain
    ele_ind_fem = np.where(msh["IX"][:, 4] <= 2)[0]
    node_ind_fem = np.unique(np.hstack((msh["IX"][ele_ind_fem, :4])))
    node_ind_inv = np.zeros((np.max(node_ind_fem),))
    for i in range(len(node_ind_fem)):
        node_ind_inv[node_ind_fem[i]-1] = i+1 
    fem = {}
    fem["extX"] = msh["X"][np.setdiff1d(np.arange(len(msh["X"]))+1, node_ind_fem)-1, :]

    # Setting of FEM variables
    
    fem["IX"] = msh["IX"][ele_ind_fem, :]  # Element connectivity
    fem["IX"] = np.column_stack((node_ind_inv[fem["IX"][:, :4]-1], fem["IX"][:, 4])).astype(int)
    fem["X"] = msh["X"][node_ind_fem - 1, :]# Nodal coordinates
    fem["nn"] = fem["X"].shape[0]
    fem["ndof"] = fem["nn"]*3
    fem["ne"] = fem["IX"].shape[0]
    fem["edof"] = np.array([fem["IX"][:,0]*3-2, fem["IX"][:,0]*3-1, fem["IX"][:,0]*3,
                            fem["IX"][:,1]*3-2, fem["IX"][:,1]*3-1, fem["IX"][:,1]*3,
                            fem["IX"][:,2]*3-2, fem["IX"][:,2]*3-1, fem["IX"][:,2]*3,
                            fem["IX"][:,3]*3-2, fem["IX"][:,3]*3-1, fem["IX"][:,3]*3], dtype=int).T
    
    fem["is"] = np.reshape(np.kron(fem["edof"],np.ones((1,12),dtype=int)),144*fem["ne"])
    fem["js"] = np.reshape(np.kron(fem["edof"],np.ones((12,1),dtype=int)),144*fem["ne"])
    fem["E0"] = input["mprop"][0]
    fem["Emin"] = fem["E0"]*(1E-9)
    fem["D"] = 1 / ((1+input["mprop"][1])*(1-2*input["mprop"][1])) * (
                np.array([[1-input["mprop"][1], input["mprop"][1],   input["mprop"][1],   0,0,0],
                          [input["mprop"][1],   1-input["mprop"][1], input["mprop"][1],   0,0,0],
                          [input["mprop"][1],   input["mprop"][1],   1-input["mprop"][1], 0,0,0],
                          [0,        0,        0,             (1-2*input["mprop"][1])/2    ,0,0],
                          [0,        0,        0,             0,  (1-2*input["mprop"][1])/2  ,0],
                          [0,        0,        0,             0,   0, (1-2*input["mprop"][1])/2]]))

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
    fem["nx"]=nx
    fem["ny"]=ny
    fem["nz"]=nz
    # Build solid stiffnes matrix (K_S) and volume vector ()
    fem["K_S"] = np.zeros((144*fem["ne"]))

    fem["Ve"] = np.zeros(fem["ne"])

    for e in range(fem["ne"]):
        px = np.array([nx[e, 0], nx[e, 1], nx[e, 2], nx[e, 3]]).reshape(-1,1)  
        py = np.array([ny[e, 0], ny[e, 1], ny[e, 2], ny[e, 3]]).reshape(-1,1)  
        pz = np.array([nz[e, 0], nz[e, 1], nz[e, 2], nz[e, 3]]).reshape(-1,1)  
      
        V = np.abs(1./6 * np.linalg.det(np.vstack(([np.ones((1,4)),px.T,py.T,pz.T]))))
        
        B = np.zeros((6,12))
        for i in range(4):
            ind = [0,1,2,3]
            ind.remove(i)

            pm = [1,-1,1,-1]

            beta  = -pm[i]*np.linalg.det(np.vstack(([np.ones((1,3)),py[ind].T,pz[ind].T])))
            gamma =  pm[i]*np.linalg.det(np.vstack(([np.ones((1,3)),px[ind].T,pz[ind].T])))
            delta = -pm[i]*np.linalg.det(np.vstack(([np.ones((1,3)),px[ind].T,py[ind].T])))

            B[:,(i)*3:(i)*3+3] = 1./(6*V) * np.array([[beta, 0, 0],
                                                          [0, gamma, 0],
                                                          [0, 0, delta],
                                                          [gamma, beta, 0],
                                                          [0, delta, gamma],
                                                          [delta, 0, beta]] )

        Ke = np.matmul(np.matmul(np.transpose(B),fem["D"]),B)*V
        Ve = V

        fem["K_S"][e*144:(e+1)*144] = np.reshape(Ke,(144))
        
        fem["Ve"][e] = Ve
    # Setting of FEM boundary condition variables
    
    node_ind_inv = node_ind_inv.astype(int)
    for i in range(len(msh["bc"]["r"]["u"])):
        msh["bc"]["r"]["u"][i] = node_ind_inv[msh["bc"]["r"]["u"][i+1]-1]

    for i in range(len(msh["bc"]["r"]["v"])):
        msh["bc"]["r"]["v"][i] = node_ind_inv[msh["bc"]["r"]["v"][i+1]-1]

    for i in range(len(msh["bc"]["r"]["w"])):
        msh["bc"]["r"]["w"][i] = node_ind_inv[msh["bc"]["r"]["w"][i+1]-1]

    for i in range(len(msh["bc"]["f"])):
        msh["bc"]["f"][i] = node_ind_inv[msh["bc"]["f"][i+1]-1]

    for i in range(len(msh["lc"]["t"])):
        msh["lc"]["t"][i] = node_ind_inv[msh["lc"]["t"][i+1]-1]

    for i in range(len(msh["lc"]["b"])):
        msh["lc"]["b"][i] = node_ind_inv[msh["lc"]["b"][i+1]-1]


    
    fem["bcdof"] = []
    fem["bcval"] = []

    # roller boundary condition (u = 0)
    for i in range(len(msh["bc"]["r"]["u"])-1):  # Iterate over the boundary conditions
        reshaped = (msh["bc"]["r"]["u"][i].T.flatten()) * 3 - 2  # Adjust for 0-based indexing
        fem["bcdof"] = np.concatenate((fem["bcdof"], reshaped))
        fem["bcval"] = np.concatenate((fem["bcval"], np.zeros((msh["bc"]["r"]["u"][i].T.flatten()).shape)))

    for i in range(len(msh["bc"]["r"]["v"])-1):  # Iterate over the boundary conditions
        reshaped = msh["bc"]["r"]["v"][i].T.flatten()* 3 - 1  # Adjust for 0-based indexing
        fem["bcdof"] = np.concatenate((fem["bcdof"], reshaped.flatten()))
        fem["bcval"] = np.concatenate((fem["bcval"], np.zeros((msh["bc"]["r"]["v"][i].T.flatten()).shape)))

    for i in range(len(msh["bc"]["r"]["w"])-1):  # Iterate over the boundary conditions
        reshaped = msh["bc"]["r"]["w"][i].T.flatten() * 3  # Adjust for 0-based indexing
        fem["bcdof"] = np.concatenate((fem["bcdof"], reshaped.flatten()))
        fem["bcval"] = np.concatenate((fem["bcval"], np.zeros((msh["bc"]["r"]["w"][i].T.flatten()).shape)))
    
    for i in range(len(msh["bc"]["f"])-1):  # Iterate over the boundary conditions
        reshaped = (msh["bc"]["f"][i].T.flatten()) * 3 - 2
        reshaped1 = (msh["bc"]["f"][i].T.flatten()) * 3 - 1
        reshaped2 = (msh["bc"]["f"][i].T.flatten()) * 3  # Adjust for 0-based indexing
        fem["bcdof"] = np.concatenate((fem["bcdof"], reshaped,reshaped1,reshaped2))
        fem["bcval"] = np.concatenate((fem["bcval"], np.zeros(3*((msh["bc"]["f"][i].T.flatten())).shape[0])))
    fem["bcdof"] = fem["bcdof"].astype(int)
    fem["bcval"] = fem["bcval"].astype(int)

    
    fem["fdof"] = []
    fem["fval"] = []

    # Setting of FEM traction force variables
    for i in range(len(msh["lc"]["t"])-1):
        for j in range(np.size(msh["lc"]["t"][i], 0)):
            vec_a = np.zeros(3)
            vec_b = np.zeros(3)
            for k in range(3):
                vec_a[k] = fem["X"][msh["lc"]["t"][i][j, 2]-1, k] - fem["X"][msh["lc"]["t"][i][j, 0]-1, k]
                vec_b[k] = fem["X"][msh["lc"]["t"][i][j, 1]-1, k] - fem["X"][msh["lc"]["t"][i][j, 0]-1, k]
            
   

            S = np.linalg.norm(np.cross(vec_a, vec_b)) / 2
            
            for k in range(3):
                fem["fdof"] = np.concatenate((fem["fdof"], [msh["lc"]["t"][i][j, k] * 3 - 2, msh["lc"]["t"][i][j, k] * 3 - 1, msh["lc"]["t"][i][j, k] * 3]))
                fem["fval"] = np.concatenate((fem["fval"], [S / 3 * input["f"][0], S / 3 * input["f"][1], S / 3 * input["f"][2]]))
    fem["fdof"]= fem["fdof"].astype(int)

    return (fem)