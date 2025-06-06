import matplotlib.pyplot as plt
from matplotlib import animation
from pylab import *
import numpy as np

def plot_distribution_x(w_dir, number, file_name):

    F = loadtxt(w_dir + '/F' + str(number) + '.dat')
    xgrid = loadtxt(w_dir + '/xgrid.dat')
    ygrid = loadtxt(w_dir + '/ygrid.dat')
    zgrid = loadtxt(w_dir + '/zgrid.dat')
    pgrid = loadtxt(w_dir + '/pgrid.dat')

    Nx = len(xgrid)
    Ny = 1
    Nz = 1
    Np = pgrid.shape[0]

    F1 = np.zeros([Nz, Ny, Nx, Np])
    F2 = np.zeros([Nz, Ny, Nx, Np])

    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                for l in range(Np):
                    F1[k][j][i][l] = F[Np*Nx*Ny*k + Np*Nx*j + Np*i + l, 1]
                    F2[k][j][i][l] = F[Np*Nx*Ny*k + Np*Nx*j + Np*i + l, 2]

    Fx1 = np.zeros([Nx])
    Fx2 = np.zeros([Nx])
    Fxa = np.zeros([Nx])

    for i in range(Nx):
        Fx1[i] = F1[int(Nz/2), int(Ny/2), i, 2]
        Fx2[i] = F2[int(Nz/2), int(Ny/2), i, 2]

    for i in range(int(Nx/2), Nx):
        Fxa[i] = Fx2[int(Nx/2)]

    for i in range(int(Nx/2)):
        Fxa[i] = Fx2[int(Nx/2)]*exp((xgrid[i] - xgrid[int(Nx/2)])/50)

    plt.rcParams.update({'font.size': 40})
    #plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 0.5
    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    ax.set_xlabel(r'$x$', fontsize=40, fontweight='bold')
    ax.set_ylabel(r'$F(x)$', fontsize=40, fontweight='bold')
    ax.set_yscale("log")
    #ax.set_ylim([6E-15, 5E-10])
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    ax.minorticks_on()

    plt.plot(xgrid, Fx1, label = 'explicit')
    plt.plot(xgrid, Fx2, label = 'implicit')
    plt.plot(xgrid, Fxa, label = 'analytic')

    ax.legend(fontsize = "20")
    plt.savefig(file_name + '.png', bbox_inches='tight')
    plt.close()