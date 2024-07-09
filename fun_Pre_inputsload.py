def fun_Pre_inputsload(modelname):
    input={}
    input["modelname"] = modelname
    input["penal"] = 3
    input["rmin"] = 0.03
    input["initdv"] = -0.5
    input["VT"] = 0.0011780970000000001
    input["volfrac"] = 0.15
    input["mprop"] = [1.463E+9,    0.37]
    input["f"] = [0, 0, -1000000]
    input["conv"] = 0.001
    input["bt"] = 0.1
    input["bt_ic"] = 1.5
    input["bt_ns"] = 4
    input["bt_fn"] = 10
    #resol, xmir, ymir, zmir
    input["MMA_c"] = 10
    input["resol"] = 0.005
    input["xmir"] = 1
    input["ymir"] = 1
    input["zmir"] = 0
    return input