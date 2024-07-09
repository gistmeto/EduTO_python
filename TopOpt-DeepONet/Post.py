import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import tensorflow as tf

######################### Function definition #####################
def post_losses(loss_list):
    #loss_list = np.load("loss_list.npy")
    fig, ax= plt.subplots(1,1)
    ax.semilogy(loss_list)
    ax.set_xlabel('epochs')
    ax.set_ylabel('$log \mathcal{L}$')
    ax.set_title("history", fontsize=14)
    plt.savefig('Loss_history')
    plt.close()
    return True

def post_fields(DeepONet,case_indices):
    Y = DeepONet.Y
    x=DeepONet.Y[:, 0]
    y=DeepONet.Y[:, 1]
    triang=tri.Triangulation(x, y)

    for i in case_indices:
        rho_data_i = tf.gather(DeepONet.rho_data, i, axis=1)

        Y_i = DeepONet.Y
        U_i = tf.expand_dims(tf.gather(DeepONet.U_data, i, axis=1),axis=0)   
        rho_i_NN = DeepONet.cal_rho_N(U_i,Y_i)

        fig, ax = plt.subplots(2,1)
        ax[0].set_aspect('equal')
        ax[0].tripcolor(triang, tf.squeeze(rho_data_i), shading='flat', vmin = 0, vmax = 1, cmap = 'gray_r')
        ax[0].set_title('Rho_Data %3d'%(i), fontsize=10)

        ax[1].set_aspect('equal')
        ax[1].tripcolor(triang, tf.squeeze(rho_i_NN), shading='flat', vmin = 0, vmax = 1, cmap = 'gray_r')
        ax[1].set_title('Rho_NN %3d'%(i), fontsize=10)
        plt.savefig('Rho_%03d'%(i))
        plt.close()
    return True

def post_R_line(DeepONet,case_indices):
    Y = DeepONet.Y
    for i in case_indices:
        rho_data_i = tf.gather(DeepONet.rho_data, i, axis=1)

        Y_i = DeepONet.Y
        U_i = tf.expand_dims(tf.gather(DeepONet.U_data, i, axis=1),axis=0)   
        rho_i_NN = DeepONet.cal_rho_N(U_i,Y_i)

        fig1, ax1 = plt.subplots()
        ax1.set_aspect('equal')
        plt.scatter(rho_data_i,rho_i_NN)
        plt.plot([-0.2,1.2],[-0.2,1.2],color="black")
        plt.xlabel("Ground truth")
        plt.ylabel("Predicted")
        plt.ylim([-0.2,1.2])
        plt.ylim([-0.2,1.2])
        ax1.set_title('Rho correlration %3d'%(i))
        plt.savefig('Rho_Correlation_%03d'%(i))
        plt.close()
    return True
