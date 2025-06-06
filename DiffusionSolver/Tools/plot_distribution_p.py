import matplotlib.pyplot as plt
from matplotlib import animation
from pylab import *
import numpy as np

def plot_distribution_p(w_dir, number, file_name):

    F = loadtxt(w_dir + '/F' + str(number) + '.dat')
    xgrid = loadtxt(w_dir + '/xgrid.dat')
    ygrid = loadtxt(w_dir + '/ygrid.dat')
    zgrid = loadtxt(w_dir + '/zgrid.dat')
    pgrid = loadtxt(w_dir + '/pgrid.dat')

    Nx = xgrid.shape[0]
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

    Fp1 = np.zeros([Np])
    Fp2 = np.zeros([Np])
    Fpa = np.zeros([Np])

    for i in range(Np):
        Fp1[i] = F1[int(Nz/2), int(Ny/2), int(Nx/2), i]
        Fp2[i] = F2[int(Nz/2), int(Ny/2), int(Nx/2), i]

    for i in range(Np):
        Fpa[i] = Fp2[2]*pgrid[2]/pgrid[i]


    plt.rcParams.update({'font.size': 40})
    #plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 0.5
    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    ax.set_xlabel(r'$p$', fontsize=40, fontweight='bold')
    ax.set_ylabel(r'$F(p)$', fontsize=40, fontweight='bold')
    ax.set_yscale("log")
    ax.set_xscale("log")
    #ax.set_ylim([6E-15, 5E-10])
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    ax.minorticks_on()

    plt.plot(pgrid, Fp1, label = 'explicit')
    plt.plot(pgrid, Fp2, label = 'implicit')
    plt.plot(pgrid, Fpa, label = 'analytic')

    ax.legend(fontsize = "20")
    plt.savefig(file_name + '.png', bbox_inches='tight')