import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import time

def fun_OPT_plot(opt,fem):

    plt.figure(0)
    plt.clf()
    ax = plt.gca()
    ax.set_title("Objective")
    plt.plot(np.array(opt["f"]),'k')
    plt.grid(True)
    plt.draw()
    plt.pause(0.001)

    plt.figure(1)
    plt.clf()
    ax = plt.gca()
    ax.set_title("Constraint")
    plt.plot(np.array(opt["g"]+np.ones((len(opt["g"],)))*opt["volfrac"]),'b')
    plt.grid(True)
    plt.draw()
    plt.pause(0.001)

    # Figure 생성 및 크기 설정
    fig = plt.figure(2)
    fig.clf()
    fig.set_size_inches(6, 3.5)
    # 3D 산점도 그리기
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title("Density distribution")
    scatter = ax.scatter(fem["X"][:, 0], fem["X"][:, 1], fem["X"][:, 2], c=-opt["nrho"], s=3, cmap='gray_r')
    # 시점 설정
    ax.view_init(elev=20, azim=-30)
    # 축 설정
    ax.grid(True)
    # 각 축의 범위 설정하여 비율 유지
    x_range = np.ptp(fem["X"][:, 0])
    y_range = np.ptp(fem["X"][:, 1])
    z_range = np.ptp(fem["X"][:, 2])
    max_range = max(x_range, y_range, z_range)
    ax.set_xlim([np.mean(fem["X"][:, 0]) - max_range / 2, np.mean(fem["X"][:, 0]) + max_range / 2])
    ax.set_ylim([np.mean(fem["X"][:, 1]) - max_range / 2, np.mean(fem["X"][:, 1]) + max_range / 2])
    ax.set_zlim([np.mean(fem["X"][:, 2]) - max_range / 2, np.mean(fem["X"][:, 2]) + max_range / 2])
    # 컬러맵 설정
    fig.colorbar(scatter, pad=0.1)
    # 그림 업데이트
    plt.draw()
    plt.pause(0.001)

    if opt["iter"] % 10 == 0:
        print(f"계산 시간: {time.time()} 초")
    return