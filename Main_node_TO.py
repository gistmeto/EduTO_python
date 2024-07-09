import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from mma import mmasub
from fun_Pre_inputsload import fun_Pre_inputsload
from fun_Pre_mshload import fun_Pre_mshload
from fun_Pre_feminit import fun_Pre_feminit
from fun_Pre_optinit import fun_Pre_optinit

from fun_FEM_solve import fun_FEM_solve

from fun_OPT_fgdfdg import fun_OPT_fgdfdg
from fun_OPT_plot import fun_OPT_plot
from fun_post import fun_post
import os

#Pre Processing
modelname= 'model/Ex1_Beam_Coarse'
(input) = fun_Pre_inputsload(modelname)
(msh) = fun_Pre_mshload(input)
(fem) = fun_Pre_feminit(input, msh) 
(opt,MMA) = fun_Pre_optinit(input, fem)

#main processing - Topology Optimization
while  opt["bt"]<input["bt_fn"]:
    opt["fdv"] = spsolve(opt["Kft_sparse"],(sp.csc_matrix.dot(opt["Tft"],opt["nv"])))
    
    opt["nrho"] = np.maximum(np.minimum(np.tanh(np.dot(opt["bt"],opt["fdv"])) / (2 * np.tanh(opt["bt"])) + 0.5, 1), -1)   
    opt["erho"] = opt["Ten"].dot(opt["nrho"])
    
    (fem["U"],fem["K"]) = fun_FEM_solve(fem,opt)
   
    (opt["f"], opt["g"], opt["dfdx"],opt["dgdx"]) = fun_OPT_fgdfdg(fem,opt)

    if opt["iter"] > 1:
        opt["deltaf"] = np.abs( (opt["f"][-1]-opt["f"][-2])/opt["f"][-2] )

    print("iter : %3d\tf : %.4f\tVolume : %.4f\tdeltaf : %.5f\tbeta : %.2f"%(opt["iter"],opt["f"][-1],opt["g"][-1]+input["volfrac"],opt["deltaf"],opt["bt"]))
    fun_OPT_plot(opt,fem)

    (opt["dvnew"], ymma, zmma, lam, xsi, eta, mu, zet, s, MMA["low"], MMA["upp"]) = mmasub(1,len(opt["dv"]),opt["iter"],
                                                                                           opt["dv"],opt["dvmin"],opt["dvmax"],opt["dvold"],opt["dvolder"],
                                                                                           opt["f"][-1],opt["dfdx"],opt["g"][-1],opt["dgdx"],
                                                                                           MMA["low"],MMA["upp"],MMA["a0"],MMA["a"],MMA["c"],MMA["d"],1) 

    opt["iter"] += 1
    opt["dvolder"] = opt["dvold"]
    opt["dvold"] = opt["dv"]
    opt["dv"] = opt["dvnew"]
    opt["nv"][opt["dof_dd"]-1] = opt["dv"]

    if  (opt["cont_sw"]==0) and (opt["deltaf"] < input["conv"]):
        opt["cont_sw"]=1
        opt["cont_iter"] = 0
        print("Continuation start")
    elif (opt["cont_sw"]==1):
        opt["cont_iter"] += 1
        if (np.mod(opt["cont_iter"],input["bt_ns"])==1):
            opt["bt"] *= input["bt_ic"]

    
    
#post processing

fun_post(fem, opt, input)
print('finish')


