import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import  LinearNDInterpolator
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage.measure import marching_cubes
from stl import mesh
import os
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def fun_post(fem, opt, input):
    file_name = f"Result_{input['modelname']}.npz"
    dir_name = os.path.dirname(file_name)
    if dir_name and not os.path.exists(dir_name):
        os.makedirs(dir_name)
    np.savez(f'Result_{input["modelname"]}.npz', fem=fem,opt=opt,input=input)

    data = np.load(f'Result_{input["modelname"]}.npz', allow_pickle=True)


    fem = data['fem'].item()  
    opt = data['opt'].item()
    input = data['input'].item()


    
    box = np.array([np.min(fem['extX'], axis=0), np.max(fem['extX'], axis=0)])
    Xgd, Ygd, Zgd = np.meshgrid(np.arange(box[0, 0], box[1, 0], input['resol']),
                                np.arange(box[0, 1], box[1, 1]+input['resol'], input['resol']),
                                np.arange(box[0, 2], box[1, 2]+input['resol'], input['resol']),
                                indexing='ij')

    fdv_values = np.concatenate((opt['fdv'], -1 * np.ones(fem['extX'].shape[0])))
    points = np.concatenate((fem['X'], fem['extX']), axis=0)

    interpolator = LinearNDInterpolator(points, fdv_values)
    Rlfdy = interpolator((Xgd, Ygd, Zgd))

    if input['xmir'] == 1:
        Xgd = np.concatenate((np.flip(-Xgd, axis=0), Xgd), axis=0)
        Ygd = np.concatenate((Ygd, Ygd), axis=0)
        Zgd = np.concatenate((Zgd, Zgd), axis=0)
        Rlfdy = np.concatenate((np.flip(Rlfdy, axis=0), Rlfdy), axis=0)

    if input['ymir'] == 1:
        Xgd = np.concatenate((Xgd, Xgd), axis=1)
        Ygd = np.concatenate((np.flip(-Ygd, axis=1), Ygd), axis=1)
        Zgd = np.concatenate((Zgd, Zgd), axis=1)
        Rlfdy = np.concatenate((np.flip(Rlfdy, axis=1), Rlfdy), axis=1)

    if input['zmir'] == 1:
        Xgd = np.concatenate((Xgd, Xgd), axis=2)
        Ygd = np.concatenate((Ygd, Ygd), axis=2)
        Zgd = np.concatenate((np.flip(-Zgd, axis=2), Zgd), axis=2)
        Rlfdy = np.concatenate((np.flip(Rlfdy, axis=2), Rlfdy), axis=2)

    verts, faces, _, _ = marching_cubes(Rlfdy, level=0)
    fig = plt.figure(4)
    plt.clf()
    ax = fig.add_subplot(111, projection='3d')
    mesh3D = Poly3DCollection(verts[faces], alpha=0.1)
    mesh3D.set_facecolor((0.1, 0.1, 0.1))
    ax.set_xlim([0, box[1, 0] * 400])
    ax.set_ylim([0, box[1, 1] * 400])
    ax.set_zlim([0, box[1, 2] * 400])
    ax.grid(True)
    ax.view_init(elev=10, azim=100)
    ax.add_collection3d(mesh3D)
    plt.show()

    stl_data = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            stl_data.vectors[i][j] = verts[f[j], :]

    stl_data.save(f'Result_{input["modelname"]}.stl')





