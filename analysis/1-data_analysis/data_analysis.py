import numpy as np
import json
import os
import time
import pandas as pd
from copy import deepcopy
from typing import List

from pymatgen.core import Molecule
from pymatgen.analysis.molecule_matcher import MoleculeMatcher


def analyze_molinfo(data_folder: List):
    start_time = time.time()
    # Input data
    elemental_energy = {'H': -16.04216696, 'C': -1038.35838704, 'O': -2044.20177304, 'F': -2714.17577144}
    ID = []  # Molecular EP ID
    BE = []  # Binding energy, eV
    FE_sol = []  # Formation energy in solvent, eV/atom
    FE_vac = []  # Formation energy in vaccum, eV/atom
    LUMO_sol = []  # Unit in eV
    HOMO_sol = []  # Unit in eV
    LUMO_cls = []  # Unit in eV
    HOMO_cls = []  # Unit in eV
    DE_sol = []  # Dielectric constant of solvents
    DE_ely = []  # Dielectric constant of electrolytes
    Vis_sol = []  # Viscosity of solvents
    Vis_ely = []  # Viscosity of electrolytes
    IC = []  # Ionic conductivity
    D_Li = []  # Diffusivity of Li
    D_FSI = []  # Diffusivity of FSI
    D_sol = []  # Diffusivity of solvents
    Similarity = []  # RMS displacement between two structures;
    # Molecules are placed in a box with a=b=c=20 A; rms displacement normalized by (Vol / nsites) ** (1/3)

    num = 0
    error_list = []
    sm = MoleculeMatcher(tolerance=10)  # Molecule Matcher to get the rmsd between structures

    for x, data_folder_x in enumerate(data_folder):
        for root, dis, files in os.walk(data_folder_x):
            for file in files:
                try:
                    num += 1
                    if num % 100 == 0:
                        print(num)
                    with open(os.path.join(root, file), 'r') as fin:
                        data = json.loads(fin.read())
                    # print('Loading ', data['name'])

                    ID.append(data['name'])
                    try:
                        DE_sol.append(data['dielectric_constant_RESP2'])
                    except:
                        DE_sol.append(data['dielectric_constant_REPS2'])
                    Vis_sol.append(float(data['viscosity_RESP2(mPas)'][0].split('=')[-1]))
                    try:
                        DE_ely.append(data['dielectric_constant_ely'])
                        Vis_ely.append(float(data['viscosity_ely(mPas)'][0].split('=')[-1]))
                    except:
                        DE_ely.append(data['dielectric_constant_ely_RESP2'])
                        Vis_ely.append(float(data['viscosity_ely_RESP2(mPas)'][0].split('=')[-1]))

                    # Calculating binding energy
                    try:
                        BE.append(27.21138 * (min(list(map(float, data['HFenergy_cls'])))
                                              - float(data['HFenergy_RESP2'])
                                              - float(data['HFenergy_Li+_RESP2'])))
                    except:
                        BE.append(27.21138 * (-float(data['HFenergy_cls'].split('-')[-1])
                                              - float(data['HFenergy_RESP2'])
                                              - float(data['HFenergy_Li+_RESP2'])))

                    # Calculating formation energy
                    comp = data['HOMOLUMO_RESP2'][0][0]
                    if len(data['HOMOLUMO_RESP2']) > 1:
                        for i in range(1, len(data['HOMOLUMO_RESP2'])):
                            comp = comp + data['HOMOLUMO_RESP2'][i][0]
                    elemetal = []
                    e_num = []
                    for i in range(len(comp)):
                        if comp[i].isdigit() == False:
                            elemetal.append(comp[i])
                            num_i = []
                            if i + 1 < len(comp) and comp[i + 1].isdigit() == True:
                                num_i = comp[i + 1]
                                ii = deepcopy(i)
                                if i + 2 < len(comp):
                                    while comp[ii + 2].isdigit() == True:
                                        num_i = num_i + comp[ii + 2]
                                        if ii + 3 < len(comp):
                                            ii = ii + 1
                                        else:
                                            break
                                e_num.append(float(num_i))
                            else:
                                e_num.append(float(1))
                    fe = 0
                    for i, ele_i in enumerate(elemetal):
                        fe = fe + elemental_energy[ele_i] * e_num[i]

                    FE_sol.append((27.21138 * float(data['HFenergy_RESP2']) - fe) / sum(e_num))
                    FE_vac.append((27.21138 * float(data['HFenergy_vum']) - fe) / sum(e_num))

                    # Loading diffusivity; calculating conductivity and transference number (TBA)
                    try:
                        for i, di_i in enumerate(data['diffusivity(m^2/s)']):
                            if 'cation' in di_i or 'Li' in di_i:
                                D_Li.append(float(di_i.split('=')[-1]))
                            elif 'anion' in di_i or 'FSI' in di_i:
                                D_FSI.append(float(di_i.split('=')[-1]))
                            elif 'ep' in di_i:
                                D_sol.append(float(di_i.split('=')[-1]))
                    except:
                        for i, di_i in enumerate(data['diffusivity_RESP2(m^2/s)']):
                            if 'cation' in di_i or 'Li' in di_i:
                                D_Li.append(float(di_i.split('=')[-1]))
                            elif 'anion' in di_i or 'FSI' in di_i:
                                D_FSI.append(float(di_i.split('=')[-1]))
                            elif 'ep' in di_i:
                                D_sol.append(float(di_i.split('=')[-1]))

                    # Loading LUMO & HOMO of solvents
                    HL1s = data['HOMOLUMO_RESP2']
                    energy_sol_unoccupied = []
                    energy_sol_occupied = []
                    for i, HL_i in enumerate(HL1s):
                        for j, HL_j in enumerate(HL_i):
                            if j > 0:
                                if float(HL_j.split(' ')[
                                             -4]) < 0.5:  # Differing occupied and unoccupied MOs by occupancy number
                                    energy_sol_unoccupied.append(float(HL_j.split(' ')[-1]))
                                else:
                                    energy_sol_occupied.append(float(HL_j.split(' ')[-1]))
                    LUMO_sol.append(27.21138 * min(energy_sol_unoccupied))
                    HOMO_sol.append(27.21138 * max(energy_sol_occupied))

                    # Loading LUMO & HOMO of solvents in clusters
                    HL2s = data['HOMOLUMO_cls']
                    energy_cls_unoccupied = []
                    energy_cls_occupied = []

                    # Selecting the HOMO&LUMO of the most stable cluster
                    try:
                        tol_energy_cls = list(map(float, data['HFenergy_cls']))
                        min_index = tol_energy_cls.index(min(tol_energy_cls))
                        HL2s_selected = HL2s[min_index]
                        cls = data['coordinate_cls'][min_index]  # The structure of the most stable cluster
                    except:
                        HL2s_selected = HL2s
                        cls = data['coordinate_cls']

                    for i, HL_i in enumerate(HL2s_selected):
                        if 'Li' not in HL_i[0]:  # Exluding the parts of Li
                            for j, HL_j in enumerate(HL_i):
                                if j > 0:
                                    if float(HL_j.split(' ')[
                                                 -4]) < 0.5:  # Differing occupied and unoccupied MOs by occupancy number
                                        energy_cls_unoccupied.append(float(HL_j.split(' ')[-1]))
                                    else:
                                        energy_cls_occupied.append(float(HL_j.split(' ')[-1]))

                    LUMO_cls.append(27.21138 * min(energy_cls_unoccupied))
                    HOMO_cls.append(27.21138 * max(energy_cls_occupied))

                    # Calculating molecule similarity before and after interacting with a Li+
                    mol = data['coordinate_RESP2']  # The pristine molecule structure
                    # cls is read in the HOMO&LUMO part; The Li-Solvent complex structure

                    species = []
                    coords = []
                    for i in range(len(mol)):
                        mol_i = mol[i].split(' ')
                        species.append(str(mol_i[0]))
                        coords.append([float(mol_i[3]), float(mol_i[6]), float(mol_i[9])])
                    solvent = Molecule(species=species, coords=coords)

                    species = []
                    coords = []
                    for i in range(len(cls)):
                        mol_i = cls[i].split(' ')
                        species.append(str(mol_i[0]))
                        coords.append([float(mol_i[3]), float(mol_i[6]), float(mol_i[9])])
                    cluster = Molecule(species=species, coords=coords)
                    cluster.remove_species(['Li'])
                    # Get RMSD between two molecule with arbitrary atom order
                    Similarity.append(sm.get_rmsd(solvent, cluster))
                    # print('Finished')
                    num = num + 1
                except:
                    error_list.append(file)
    print(num, ' molecules have been loaded')
    print(f'--- {time.time() - start_time:.4f} seconds ---')

    """
    plt.figure()
    ax=plt.gca()
    dLUMO = []
    dHOMO = []
    for i in range(len(LUMO_cls)):
        dLUMO.append(LUMO_cls[i]-LUMO_sol[i])
        dHOMO.append(HOMO_cls[i]-HOMO_sol[i])
    plt.scatter(BE,dLUMO,s=8.0,color='r')
    front1={'family':'Arial','weight': 'normal','size':16}
    front2={'family':'Arial','weight': 'normal','size':14}
    plt.xlabel('Binding energy (eV)',front1)
    plt.ylabel('LUMO energy change (eV)',front1)
    plt.tick_params(axis='both',which='major',length=6,width=2,direction='in',labelsize=14)#设置主坐标轴刻度大小
    plt.tick_params(axis='both',which='minor',length=3,width=1,direction='in',labelsize=14)#设置次坐标轴刻度大小
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2)
    plt.xlim(-6,2)
    plt.ylim(-12,6)
    x_major_locator=MultipleLocator(2)
    y_major_locator=MultipleLocator(3)
    ax.xaxis.set_major_locator(x_major_locator)
    ax.yaxis.set_major_locator(y_major_locator)
    plt.show()

    plt.figure()
    ax=plt.gca()
    plt.scatter(BE,dHOMO,s=3.0,color='r')
    front1={'family':'Arial','weight': 'normal','size':16}
    front2={'family':'Arial','weight': 'normal','size':14}
    plt.xlabel('Binding energy (eV)',front1)
    plt.ylabel('HOMO energy change (eV)',front1)
    plt.tick_params(axis='both',which='major',length=6,width=2,direction='in',labelsize=14)#设置主坐标轴刻度大小
    plt.tick_params(axis='both',which='minor',length=3,width=1,direction='in',labelsize=14)#设置次坐标轴刻度大小
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2)
    plt.xlim(-6,2)
    plt.ylim(-12,6)
    x_major_locator=MultipleLocator(2)
    y_major_locator=MultipleLocator(3)
    ax.xaxis.set_major_locator(x_major_locator)
    ax.yaxis.set_major_locator(y_major_locator)
    plt.show()

    plt.figure()
    ax=plt.gca()
    plt.scatter(LUMO_sol,HOMO_sol,s=10.0,color='r')
    front1={'family':'Arial','weight': 'normal','size':16}
    front2={'family':'Arial','weight': 'normal','size':14}
    plt.xlabel('LUMO energy (eV)',front1)
    plt.ylabel('HOMO energy (eV)',front1)
    plt.tick_params(axis='both',which='major',length=6,width=2,direction='in',labelsize=14)#设置主坐标轴刻度大小
    plt.tick_params(axis='both',which='minor',length=3,width=1,direction='in',labelsize=14)#设置次坐标轴刻度大小
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2)
    plt.xlim(-6,12)
    plt.ylim(-12,3)
    x_major_locator=MultipleLocator(2)
    y_major_locator=MultipleLocator(3)
    ax.xaxis.set_major_locator(x_major_locator)
    ax.yaxis.set_major_locator(y_major_locator)
    plt.show()

    plt.figure()
    ax=plt.gca()
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    plt.axes(xscale = 'log', yscale = 'log')
    plt.scatter(D_Li,D_FSI,s=10.0,color='r')
    front1={'family':'Arial','weight': 'normal','size':16}
    front2={'family':'Arial','weight': 'normal','size':14}
    plt.xlabel('Diffusivity of cations (m^2/s)',front1)
    plt.ylabel('Diffusivity of cations (m^2/s)',front1)
    plt.tick_params(axis='both',which='major',length=6,width=2,direction='in',labelsize=14)#设置主坐标轴刻度大小
    plt.tick_params(axis='both',which='minor',length=3,width=1,direction='in',labelsize=14)#设置次坐标轴刻度大小
    plt.xlim(1e-15,1e-5)
    plt.ylim(1e-15,1e-5)
    #x_major_locator=MultipleLocator(1e)
    #y_major_locator=MultipleLocator(3)
    #ax.xaxis.set_major_locator(x_major_locator)
    #ax.yaxis.set_major_locator(y_major_locator)
    plt.show()
    """

    # Save as csv
    x = {
        'EP ID': ID,
        'Binding energy (eV)': BE,
        'Formation energy in solvent (eV/atom)': FE_sol,
        'Formation energy in vaccum (eV/atom)': FE_vac,
        'LUMO_sol (eV)': LUMO_sol,
        'HOMO_sol (eV)': HOMO_sol,
        'LUMO_cls (eV)': LUMO_cls,
        'HOMO_cls (eV)': HOMO_cls,
        'LUMO change (eV)': np.array(LUMO_cls) - np.array(LUMO_sol),
        'HOMO change (eV)': np.array(HOMO_cls) - np.array(HOMO_sol),
        'Dielectric constant of solvents': DE_sol,
        'Dielectric constant of electrolytes': DE_ely,
        'Viscosity of solvents (mPa s)': Vis_sol,
        'Viscosity of electrolytes (mPa s)': Vis_ely,
        'Diffusivity of Li (m^2/s)': D_Li,
        'Diffusivity of FSI (m^2/s)': D_FSI,
        'Diffusivity of solvens (m^2/s)': D_sol,
        'Transference number of Li': np.array(D_Li) / (np.array(D_Li) + np.array(D_FSI)),
        'Similarity of solvents': Similarity
    }
    print(error_list)
    data = pd.DataFrame(x)
    data = data.sort_values(axis=0, by='EP ID', ascending=True)
    data.to_csv(f'ep_data-{time.strftime("%Y%m%d", time.localtime())}.csv', index=False)


def main():
    analyze_molinfo(['20230909'])


if __name__ == '__main__':
    main()

