import sys
import numpy as np
import os
from scipy.io import savemat

def compute_Cvv_func(sim_name):
    print('begin compute_Cvv_func')
    sim_name = str(sim_name)
    pathname = os.getcwd() + '/'
    #pathname = '/Users/yitongzheng/Documents/Corey/tumor3D/'
    fstr = os.path.join(pathname, 'pos', sim_name + '.pos')

    Nt = 10**5
    Np = 1500

    finfo = os.stat(fstr)
    print('-- Reading in {}'.format(os.path.basename(fstr)))
    print('-- File size = {:.2f} MB'.format(finfo.st_size / 1e6))

    with open(fstr, 'r') as file:
        VX = np.zeros((Nt, Np))
        VY = np.zeros((Nt, Np))
        VZ = np.zeros((Nt, Np))
        Cvv = np.ones(Nt)

        line = file.readline()
        line = file.readline()
        NCELL, tN = map(int, line.split()[1:3])
        aN = NCELL - tN
        nV = 42

        for _ in range(4):
            file.readline()
        Lx, Ly, wpos = map(float, file.readline().split()[1:4])
        L = [Lx-wpos, Ly, Ly]
        particle_index = -1
        vid = 1
        nFrame = 0
        xyz = np.zeros((tN, 3))

        for line in file:
            if line.startswith('VINFO'):
                if vid > tN:
                    continue
                parts = line.split()
                x, y, z = map(float, parts[7:10])
                xyz[vid-1, :] = [x, y, z]
                vid += 1
            elif line.startswith('NEWFR'):
                for _ in range(5):
                    file.readline()
                vid = 1
            elif line.startswith('ENDFR'):
                nFrame += 1
                print(nFrame)
                VX[nFrame-1, :] = xyz[:tN, 0]
                VY[nFrame-1, :] = xyz[:tN, 1]
                VZ[nFrame-1, :] = xyz[:tN, 2]

        VX = VX[3:nFrame-1, :]
        VY = VY[3:nFrame-1, :]
        VZ = VZ[3:nFrame-1, :]
        Cvv = Cvv[3:nFrame-1]
        
        N = len(Cvv)
        for dt in range(1,N):
            Cvv[dt] = np.mean(VX[:-dt,:]*VX[dt:,:] + VY[:-dt,:]*VY[dt:,:] + VZ[:-dt,:]*VZ[dt:,:]) / \
                      np.mean(VX[:-dt,:]**2 + VY[:-dt,:]**2 + VZ[:-dt,:]**2)

        your_folder = os.path.join(pathname, 'Cvv')
        if not os.path.exists(your_folder):
            os.mkdir(your_folder)
        data_name = os.path.join(your_folder, sim_name + '_Cvv.mat')
        savemat(data_name, {'Cvv': Cvv})


if __name__ == "__main__":
    # Check if at least one argument is provided
    if len(sys.argv) > 1:
        argument = sys.argv[1]
        compute_Cvv_func(argument)
    else:
        print("No argument provided")
