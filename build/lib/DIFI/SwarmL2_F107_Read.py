def SwarmL2_F107_Read(filename):
    # [t_f107,f107] = SwarmL2_F107_Read(filename)

    # Reads f107 values from DBL data file (Swarm L2 format).

    # Input:    filename        data filename

    # Output:   t_f107            time vector
    #           f107              data

    # A. Chulliat, 2016-09-19
    # (from earlier version: read_f107_new.m)

    t_f107 = []
    f107 = []
    fid = open(filename, 'r')

    for tline in fid:
        # skip comment lines
        if tline[0] == '#':
            print("comment line")
            continue
        # read data lines
        f = [float(i) for i in tline.split()]
        t_f107.append(f[0])
        f107.append(f[1])

    fid.close()

    return t_f107, f107
