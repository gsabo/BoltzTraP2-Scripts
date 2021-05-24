#!/usr/bin/python
# Tested on Ubuntu 20.04.2 LTS with Python 3.8.5

eps3 = 1.0e-3
inf6 = 1.0e6
rydberg = 13.60569253

def main(argv = None):
    if argv is None:
        argv = sys.argv
    if len(argv) < 5 or len(argv) > 7:
        self = '/' + argv[0]
        self = self[len(self)-self[::-1].find('/'):]
        print("")
        print("    Version 1.0.3")
        print("    Converts the output of Quantum Espresso 6.7")
        print("    to the input of BoltzTraP2. Modified version of qe2boltz.py for BoltzTraP2")
        print("    from that orginal written by Georgy Samsonidze,")
        print("    An Li, Daehyun Wee, Bosch Research (October 2011).")
        print("")
        print("    Usage: %s prefix format efermi nbnd_exclude [fn_pw [fn_energy]]" % self)
        print("")
        print("  * prefix = name of the system, same as in espresso input file")
        print("  * format = pw | bands | inteqp (read energies from the output of")
        print("    espresso/PW/pw.x, espresso/PP/bands.x or BerkeleyGW/BSE/inteqp.flavor.x)")
        print("  * efermi = Fermi energy in eV (if abs(efermi) > 1e6 and format = pw | bands")
        print("    efermi is read from fn_inp)")
        print("  * nbnd_exclude = number of the lowest energy bands to exclude in the output")
        print("  * fn_pw = name of a file containing the output of pw.x (the default is")
        print("    prefix.nscf.out, must be run with verbosity = high and ibrav = 0)")
        print("  * fn_energy = name of a file containing the output of bands.x or")
        print("    inteqp.flavor.x (the default is bands.out or bandstructure.dat,")
        print("    not used if format = pw)")
        print("")
        print("    Creates files prefix.energy and prefix.structure.")
        print("")
        return 1

    prefix = argv[1]
    ftype_inp = argv[2].lower()
    if ftype_inp != 'pw' and ftype_inp != 'bands' and ftype_inp != 'inteqp':
        print("\n    Error: unknown format.\n")
        return 2
    efermi = float(argv[3]) / rydberg
    nbnd_exclude = int(argv[4])
    if nbnd_exclude < 0:
        print("\n    Error: invalid nbnd_exclude.\n")
        return 2
    if len(argv) > 5:
        fname_pw = argv[5]
    else:
        fname_pw = prefix + '.nscf.out'
    if len(argv) > 6:
        fname_inp = argv[6]
    else:
        if ftype_inp == 'bands':
            fname_inp = 'bands.out'
        elif ftype_inp == 'inteqp':
            fname_inp = 'bandstructure.dat'

    fname_def = 'BoltzTraP.def'
    fname_intrans = prefix + '.intrans'
    fname_energy = prefix + '.energy'
    fname_struct = prefix + '.structure'

    deltae = 0.0005
    ecut = 0.4
    lpfac = 5
    efcut = 0.15
    tmax = 800.0
    deltat = 50.0
    ecut2 = -1.0
    dosmethod = 'TETRA'

    f = open(fname_pw, 'r')
    f_pw = f.readlines()
    f.close()

    if ftype_inp != 'pw':
        f = open(fname_inp, 'r')
        f_inp = f.readlines()
        f.close()

    i = 0
    natom = 0
    efermi_scf = (inf6 + eps3) / rydberg
    avec = []
    idxsym = []
    idxbnd = []
    spin = False
    symbolstring = ''
    for line in f_pw:
        if 'lattice parameter (alat)  =' in line:
            alat = float(line.split()[4])
        elif ' a(' in line:
            atext = line[23:57].split()
            avec.append([float(atext[0]) * alat, float(atext[1]) * alat, float(atext[2]) * alat])
        elif 'cryst.   s' in line:
            idxsym.append(i)
        elif 'the Fermi energy is' in line:
            efermi_scf = float(line.split()[4]) / rydberg
        elif 'highest occupied, lowest unoccupied level' in line:
            efermi_scf = (float(line.split()[6]) + float(line.split()[7])) / (2.0 * rydberg)
        elif 'number of atoms' in line:
            natom = int(line.split()[4])
        elif 'number of electrons' in line:
            nelec = float(line.split()[4])
        elif 'Sym.Ops.' in line or 'Sym. Ops.' in line:
            nsym = int(line.split()[0])
        elif 'positions (cryst. coord.)' in line:
            idatom = i + 1
        elif 'No symmetry found' in line:
            nsym = 1
        elif 'number of k points=' in line:
            nkpt = int(line.split()[4])
        elif 'number of Kohn-Sham states=' in line:
            nbnd = int(line.split()[4])
        elif ' cryst. coord.' in line:
            idxkpt = i + 1
        elif 'band energies (ev)' in line or 'bands (ev)' in line:
            idxbnd.append(i + 2)
        elif 'SPIN' in line:
            spin = True
        i += 1

    if abs(efermi) > inf6 / rydberg and (ftype_inp == 'pw' or ftype_inp == 'bands'):
        if abs(efermi_scf) > inf6 / rydberg:
            print("\n Error: Fermi energy not found.\n")
            return 2
        else:
            efermi = efermi_scf

    if spin:
        nelec -= nbnd_exclude
        nspin = 2
    else:
        nelec -= 2 * nbnd_exclude
        nspin = 1

    symbolstring = []
    PosCartCoord = []
    f_pw1 = []
    f_pw2 = []
    for ia in range(natom):
            f_pw1 = f_pw[idatom + ia][15:24]
            f_pw2 = f_pw[idatom + ia][38:73] 
            atmtext1 = f_pw1.split()
            atmtext2 = f_pw2.split()
            symbolstring.append(atmtext1[0])
            PosCartCoord.append([float(atmtext2[0]),float(atmtext2[1]),float(atmtext2[2])])

    kpoint = []
    for ik in range(nkpt):
        ktext = f_pw[idxkpt + ik][20:56].split()
        kpoint.append([float(ktext[0]), float(ktext[1]), float(ktext[2])])

    if ftype_inp == 'pw':
        energy = []
        ncol = 8
        nrow = nbnd // ncol
        if nbnd % ncol != 0:
            nrow += 1
        for ispin in range(nspin):
            for ik in range(nkpt):
                energy.append([])
                nelem = ncol
                for ir in range(nrow):   
                
                    # Store each energy line from case.nscf.out
                    eline = f_pw[idxbnd[ik] + ir]
                
                    # Set loop counter for etext to zero
                    ei = 0
                    # "range(2," -> Skip the first 2 blank characters
                    # eline[ei:ei+9] -> Split energy line every 9 characters to the energy values
                    etext = [eline[ei:ei+9] for ei in range(2, len(eline), 9)]
                                                
                    if ir == nrow - 1:
                        nelem = nbnd - ncol * (nrow - 1)
                    for ie in range(nelem):  
                        energy[ik].append(float(etext[ie]) / rydberg)
    elif ftype_inp == 'bands':
        energy = []
        ncol = 10
        nrow = nbnd / ncol
        if nbnd % ncol != 0:
            nrow += 1
        for ik in range(nkpt):
            energy.append([])
            nelem = ncol
            for ir in range(nrow):
                etext = f_inp[ik * (nrow + 1) + ir + 2].split()
                if ir == nrow - 1:
                    nelem = nbnd - ncol * (nrow - 1)
                for ie in range(nelem):
                    energy[ik].append(float(etext[ie]) / rydberg)
    elif ftype_inp == 'inteqp':
        nhead = 2
        ntot = len(f_inp) - nhead
        bndmin = int(f_inp[nhead].split()[1]) - 1
        bndmax = int(f_inp[nhead + ntot - 1].split()[1]) - 1
        nbnd = bndmax + 1
        energy = []
        for ik in range(nkpt):
            energy.append([])
            for ib in range(bndmin):
                energy[ik].append(0.0)
            for ib in range(bndmax - bndmin + 1):
                energy[ik].append(float(f_inp[nhead + ik + ib * nkpt].split()[6]) / rydberg)

    f_energy = prefix + '\n'
    f_energy += str(nkpt) + ' ' + str(nspin) + ' ' + str(efermi) + '\n'
    for ispin in range(nspin):
        for ik in range(nkpt):
            f_energy += str(kpoint[ik][0]) + ' ' + str(kpoint[ik][1]) + ' ' + str(kpoint[ik][2]) + ' ' + str(nbnd - nbnd_exclude) + '\n'
            for ib in range(nbnd_exclude, nbnd):
                f_energy += str(energy[ik][ib]) + '\n'

    f = open(fname_energy, 'w')
    f.write(f_energy)
    f.close()

    f_struct = prefix + '\n'
    for i in range(3):
        f_struct += str(avec[i][0]) + ' ' + str(avec[i][1]) + ' ' + str(avec[i][2]) + '\n'
    f_struct += str(natom) + '\n'
    for ia in range(natom):
        f_struct += str(symbolstring[ia]) + ' '
        f_struct += str(PosCartCoord[ia][0]) + ' ' + str(PosCartCoord[ia][1]) + ' ' + str(PosCartCoord[ia][2]) + '\n'

    f = open(fname_struct, 'w')
    f.write(f_struct)
    f.close()

    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
