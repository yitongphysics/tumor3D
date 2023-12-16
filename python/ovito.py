import numpy as np

# Read the input file
pth1 = "/Users/yitongzheng/Documents/Corey/tumor3D/0425_inv4/"
pth2 = "/Users/yitongzheng/Documents/Corey/tumor3D/"
filename = '11.pos'
filename = pth1 + filename

with open(filename, 'r') as fileID:
    # Skip the first two lines
    next(fileID)

    # Read NUMCL line
    line = next(fileID)
    NCELL, tN = map(int, line.split()[1:3])
    nV = 42

    # Skip three lines
    for _ in range(3):
        next(fileID)

    # Read BOXSZ line
    line = next(fileID)
    Lx, Ly, wpos = map(float, line.split()[1:4])
    IMG = 1

# Open the output file for writing
output_filename = '0_data.txt'
output_filename = pth2 + output_filename

with open(output_filename, 'w') as output_fileID:

    # Initialize variables
    pid = 1
    particle_index = 1
    frame = 0

    with open(filename, 'r') as fileID:
        for line in fileID:
            if line.startswith('NEWFR'):
                pid = 1
                line = next(fileID)
                NCELL, tN = map(int, line.split()[1:3])
                nV = 42

                for _ in range(3):
                    next(fileID)

                line = next(fileID)
                Lx, Ly, wpos = map(float, line.split()[1:4])
                output_fileID.write(f'{((NCELL-tN)*nV+tN)*IMG}\n')
                output_fileID.write(f'Lattice="{Lx-wpos} 0 0 0 {Ly} 0 0 0 {Ly}" Properties=species:S:1:pos:R:3:radius:R:1\n')

            elif line.startswith('VINFO'):
                parts = line.split()
                particle_index_new = int(parts[1])
                x, y, z, radius = map(float, parts[3:7])
                x = Lx - x

                if particle_index_new < tN:
                    particle_index = 1
                    y = y % Ly
                    z = z % Ly

                    for dLy in range(1):
                        for dLz in range(1):
                            particle_index_prt = particle_index
                            output_fileID.write(f'{particle_index_prt} {x} {y + dLy * Ly} {z + dLz * Ly} {radius}\n')
                            pid += 1
                elif particle_index_new == particle_index:
                    dy = y - y_old
                    dz = z - z_old
                    while dy>Ly/2:
                        dy -= Ly
                    while dy<-Ly/2:
                        dy += Ly
                    while dz>Ly/2:
                        dz -= Ly
                    while dz<-Ly/2:
                        dz += Ly
                    y = y_old + dy
                    z = z_old + dz
                    X.append(x)
                    Y.append(y)
                    Z.append(z)
                else:
                    if particle_index_new > tN:
                        while sum(Y)/len(Y)<0:
                            for kk in range(len(Y)):
                                Y[kk] += Ly
                        while sum(Y)/len(Y)>Ly:
                            for kk in range(len(Y)):
                                Y[kk] -= Ly
                        while sum(Z)/len(Z)<0:
                            for kk in range(len(Z)):
                                Z[kk] += Ly
                        while sum(Z)/len(Z)>Ly:
                            for kk in range(len(Z)):
                                Z[kk] -= Ly

                        for kk in range(len(X)):
                            for dLy in range(1):
                                for dLz in range(1):
                                    particle_index_prt = particle_index
                                    if particle_index_prt > 1:
                                        while particle_index_prt%9 == 1:
                                            particle_index_prt += int(particle_index_prt/9)
                                    output_fileID.write(f'{particle_index_prt} {X[kk]} {Y[kk] + dLy * Ly} {Z[kk] + dLz * Ly} {radius}\n')
                                    pid += 1

                    particle_index = particle_index_new
                    X=[]
                    Y=[]
                    Z=[]

                    if y<0:
                        y += Ly
                    elif y>Ly:
                        y -= Ly
                    if z<0:
                        z += Ly
                    elif z>Ly:
                        z -= Ly
                    X.append(x)
                    Y.append(y)
                    Z.append(z)


                y_old, z_old = y, z

            elif line.startswith('ENDFR'):
                while sum(Y)/len(Y)<0:
                    for kk in range(len(Y)):
                        Y[kk] += Ly
                while sum(Y)/len(Y)>Ly:
                    for kk in range(len(Y)):
                        Y[kk] -= Ly
                while sum(Z)/len(Z)<0:
                    for kk in range(len(Z)):
                        Z[kk] += Ly
                while sum(Z)/len(Z)>Ly:
                    for kk in range(len(Z)):
                        Z[kk] -= Ly

                for kk in range(len(X)):
                    for dLy in range(1):
                        for dLz in range(1):
                            particle_index_prt = particle_index
                            if particle_index_prt > 1 and particle_index_prt%9 == 1:
                                particle_index_prt += int(particle_index_prt/9)
                            output_fileID.write(f'{particle_index_prt} {X[kk]} {Y[kk] + dLy * Ly} {Z[kk] + dLz * Ly} {radius}\n')
                            pid += 1
                frame += 1
                print(frame)
                #break
