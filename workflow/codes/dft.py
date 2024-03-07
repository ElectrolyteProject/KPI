#coding:utf-8
import os
import math
import numpy as np
from copy import deepcopy
from pymatgen.core import Structure, Molecule, Lattice
from pymatgen.analysis.structure_matcher import StructureMatcher,ElementComparator

import global_config as gc
import assist
import local_fs, custom_shell

from typing import Any, Dict, List, Tuple


def _extract_coor(log_path: str) -> Tuple[str, List[str], int]:
    '''Extract radius, coordinates and imaginary freq of a molecule from a log file.

    Args:
    - log_path: Path to the log file.

    Returns:
    - radius: Molecule radius.
    - coordinates: Atom coordinates.
    - imaginary: Whether(0/1) has imaginary frequency.
    '''
    with open(log_path, 'r') as fin:
        info = fin.read()
    infos = info.splitlines()
    file_name = os.path.basename(log_path)
    if not 'Normal termination of Gaussian' in infos[-1]:
        raise RuntimeError(file_name + ' is NOT terminated normally!')
    if 'imaginary' in info:
        print('*******************************')
        print(file_name, 'exists imaginary frequency!')
        print('*******************************')
        img = 1
    else:
        img = 0
    so_position = []
    sb_position = []
    a0_position = []
    for i in range(len(infos)):
        if 'Standard orientation' in infos[i]:
            so_position.append(i)
        elif 'Standard basis' in infos[i]:
            sb_position.append(i)
        elif 'Recommended a0' in infos[i]:
            a0_position.append(i)
    so = max(so_position)
    sb = max(sb_position)
    if len(a0_position) > 0:
        a0 = max(a0_position)
        r = infos[a0].split()[6]
    else:
        r = None
    atom = assist.number2element()
    mol_coor = []
    for j in range(so + 5, sb - 2):
        infoss = infos[j].split()
        element_num = int(infoss[1])
        infoss[1] = atom[element_num]
        mol_coor.append(infoss[1] + '   ' + infoss[3] + '   ' + infoss[4] + '   ' + infoss[5])
    return r, mol_coor, img


def _extract_elapsedtime(log_path: str) -> float:
    '''Extract Gaussian running time.

    Args:
    - log_path: Path to the log file.

    Returns:
    - totaltime: Gaussian running time in unit hour.
    '''
    with open(log_path, 'r+') as f:
        f.seek(0, 0)
        infos = f.readlines()
    hour = 0
    minute = 0
    second = 0
    for info in infos:
        if 'Elapsed time' in info:
            hour += float(info.split()[-6])
            minute += float(info.split()[-4])
            second += float(info.split()[-2])
    totaltime = round(hour + minute / 60 + second / 3600, 4)
    return totaltime


def _judge_spin(log_path: str) -> bool:
    '''Return the spin state of the system.'''
    with open(log_path, 'r') as fin:
        fin.seek(0, 0)
        info = fin.read()
    if 'alpha spin orbitals' in info:
        spin = True
    else:
        spin = False
    return spin


def _nbo_charge(log_path: str, spin: str) -> List[str]:
    '''Extract Natural Bond Orbital charges of atoms from a log file.

    Args:
    - log_path: Path to the log file.
    - spin: Spin state of the system, `True`(open shell) or `False`(closed shell).

    Returns:
    - chagres: NBO charge per atom.
    '''
    with open(log_path, 'r+') as fin:
        infos = fin.readlines()
    a = []
    b = []
    NBOcharge = []
    for i in range(len(infos)):
        if 'Summary of Natural Population Analysis' in infos[i]:
            a.append(i)
        elif 'Natural Electron Configuration' in infos[i]:
            b.append(i)
    if spin:
        A = a[-3]
        B = b[-3]
    else:
        A = a[-1]
        B = b[-1]
    for j in range(A + 6, B - 11):
        infoss = infos[j].split()
        infossj = '   '.join(infoss[0:3])
        NBOcharge.append(infossj)
    return NBOcharge


def _mk_charge_spinpop(log_path: str) -> List[str]:
    '''Extract Mulliken charge and spin population of each atom from a log file.
    This function is only available for open shell systems.

    Args:
    - log_path: Path to the log file.

    Returns:
    - charges: MK charge and spin population per atom.
    '''
    with open(log_path, 'r') as fin:
        infos = fin.readlines()
    a = []
    charge_spinpop = []
    for i in range(len(infos)):
        if 'spin densities' in infos[i]:
            a.append(i)
    A = a[-2]
    B = a[-1]
    for j in range(A + 2, B - 1):
        infoss = infos[j].split()
        infossj = '   '.join(infoss[0:4])
        charge_spinpop.append(infossj)
    return charge_spinpop


def _energies(dir_log: str) -> Tuple[float, List[str]]:
    '''Extract energy and Gibbs energy information from a log file.

    Args:
    - log_path: Path to the log file.

    Returns:
    - hf_energy: Hatree-Fock energy.
    - thermo_energy: Kinds of free energies.
    '''
    with open(dir_log, 'r+') as f1:
        f1.seek(0, 0)
        infos = f1.readlines()
    for l, line in enumerate(infos):
        if r'1\1\GINC' in line:
            b = l
        elif r'@' in line:
            e = l
        elif 'Zero-point correction' in line:
            thermo = infos[l:l+8]
    tojoin = [i.strip(' \n') for i in infos[b:e+1]]
    infopart = ''.join(tojoin)
    _energy = infopart.split('\HF=')[-1]
    energy = eval(_energy.split('\\')[0])
    for j in range(len(thermo)):
        thermo[j] = thermo[j].strip('\n')
    return energy, thermo


def _homo_lumo(dir_log, spin):
    '''Extract HOMO & LUMO information from a log file depending on the spin state of the system.

    Args:
    - log_path: Path to the log file.
    - spin: Spin state of the system, `True`(open shell) or `False`(closed shell).

    Returns:
    - HOMOLUMO2
    '''
    with open(dir_log, 'r') as f1:
        infos = f1.readlines()
    HOMOLUMO = []
    NBO_position = []
    if spin:
        for i in range(len(infos)):
            if 'Natural Bond Orbitals (Summary)' in infos[i]:
                NBO_position.append(i)
        N_a = NBO_position[-2]  # alpha spin
        N_b = NBO_position[-1]  # beta spin
        for i in range(N_a, len(infos)):
            if ' Sorting of NBOs:' in infos[i]:
                N_a_end = i
                break
        for i in range(N_b, len(infos)):
            if ' Sorting of NBOs:' in infos[i]:
                N_b_end = i
                break
        HOMOLUMO1 = _homo_lumo_detail(infos, N_a, N_a_end, HOMOLUMO)
        HOMOLUMO2 = _homo_lumo_detail(infos, N_b, N_b_end, HOMOLUMO1)
    else:
        Lewis_position = []
        for i in range(len(infos)):
            if 'Natural Bond Orbitals (Summary)' in infos[i]:
                NBO_position.append(i)
            elif 'Total Lewis' in infos[i]:
                Lewis_position.append(i)
        N = max(NBO_position)
        M = max(Lewis_position)
        HOMOLUMO2 = _homo_lumo_detail(infos, N, M+1, HOMOLUMO)
    return HOMOLUMO2


def _homo_lumo_detail(infos, ranges, rangee, homo_lumo) -> None:
    '''Extract HOMO & LUMO information from the given log text.

    Args:
    - infos: Lines of log.
    - ranges: Starting line number of NBO orbital information.
    - rangee: Ending line number of NBO orbital information.
    - homo_lumo: A list to append extracted HOMO & LUMO information.

    Returns:
    - homo_lumo: The input paramter `homo_lumo`.
    '''
    content = []  # Storing the energy level of the molecular obirtals
    MU_position = []  # Molecular unit
    Lewis_position2 = []
    for i in range(ranges + 5, rangee):
        if 'Molecular unit' in infos[i]:
            MU_position.append(i)
        elif 'Total Lewis' in infos[i]:
            Lewis_position2.append(i)
    for mu in range(len(MU_position)):
        content3_tot = []
        for i in range(MU_position[mu] + 1, Lewis_position2[mu] - 1):
            if '. ' in infos[i]:
                if 'RY*' not in infos[i]:  # Excluding the Rydberg orbitals
                    content1 = infos[i]
                    content1 = content1.split(' ')
                    content2 = []
                    content3 = infos[i]
                    content3_tot.append(content3)
                    for cc in content1:
                        if '.' in cc:
                            content2.append(cc)
                    if len(content2) != 3:
                        raise 'Error in NBO analysis!'
                    content.append(content2)
        occupied = []
        unoccupied = []
        for j in range(len(content)):
            if float(content[j][1]) > 0.5:  # modified to 0.5, spin
                occupied.append(float(content[j][2]))
            else:
                unoccupied.append(float(content[j][2]))
        occupied.sort()
        unoccupied.sort()
        if len(occupied) > 0:
            HOMO = occupied[-1]
        else:
            HOMO = 0
        if len(occupied) > 1:
            HOMOd1 = occupied[-2]
        else:
            HOMOd1 = 0
        if len(unoccupied) > 0:
            LUMO = unoccupied[0]
        else:
            LUMO = 0
        if len(unoccupied) > 1:
            LUMOu1 = unoccupied[1]
        else:
            LUMOu1 = 0
        unit = []
        unit.append(infos[MU_position[mu]].split(' ')[-1].split('\n')[0].split('(')[-1].split(')')[0])
        content3_tot_clean = []
        for k in range(len(content)):
            content3_tots = content3_tot[k].split()
            if '(' in content3_tots[-1]:
                del content3_tots[-1]
            content3_totsj = '   '.join(content3_tots)
            content3_tot_clean.append(content3_totsj)
        for k in range(len(content)):
            if float(content[k][2]) == HOMOd1:
                unit.append(content3_tot_clean[k])
        for k in range(len(content)):
            if float(content[k][2]) == HOMO:
                unit.append(content3_tot_clean[k])
        for k in range(len(content)):
            if float(content[k][2]) == LUMO:
                unit.append(content3_tot_clean[k])
        for k in range(len(content)):
            if float(content[k][2]) == LUMOu1:
                unit.append(content3_tot_clean[k])
        homo_lumo.append(unit)  # Information of all molecular units
        content.clear()
    return homo_lumo


def _dipole_moment(dir_log: str) -> float:
    with open(dir_log, 'r+') as f1:
        f1.seek(0, 0)
        infos = f1.readlines()
    for l, line in enumerate(infos):
        if 'Dipole moment (field-independent basis, Debye)' in line:
            dm = eval(infos[l + 1].split()[-1])   # unit: Debye
    return dm


def analyze_log(log_path: str) -> Dict[str, Any]:
    '''Analyze a Gaussian log file and extract required properties.

    Args:
    - log_path: Path to the log file.

    Returns:
    - extracted: A mapping from property to value.
    '''
    r, mol_coor, img = _extract_coor(log_path)
    spin = _judge_spin(log_path)
    NBOcharge = _nbo_charge(log_path, spin)
    MKchargespinpop = _mk_charge_spinpop(log_path) if spin else ['NONE']
    HFenergy, thermoenergy = _energies(log_path)
    HOMOLUMO = _homo_lumo(log_path, spin)
    dipole_moment = _dipole_moment(log_path)
    if log_path.endswith('_vum.log') and '+Li+' not in os.path.basename(log_path):
        extracted = {
            'radius': r,
            'coordinate_vum': mol_coor,
            'HFenergy_vum': str(HFenergy),
            'thermochemistry_vum': thermoenergy,
            'HOMOLUMO_vum': HOMOLUMO,
            'NBOcharge_vum': NBOcharge,
            'MKcharge_spinpop_vum': MKchargespinpop,
            'dipole_moment_vum(Debye)': dipole_moment
        }
        if img == 1:
            extracted['img_vum'] = str(img)
    elif log_path.endswith('_sol.log'):
        extracted = {
            'coordinate_sol': mol_coor,
            'HFenergy_sol': str(HFenergy),
            'thermochemistry_sol': thermoenergy,
            'HOMOLUMO_sol': HOMOLUMO,
            'NBOcharge_sol': NBOcharge,
            'MKcharge_spinpop_sol': MKchargespinpop,
            'dipole_moment_sol(Debye)': dipole_moment,
        }
        if img == 1:
            extracted['img_sol'] = str(img)
    elif log_path.endswith('+Li+_resp2.log'):
        extracted = {
            'coordinate_cls': mol_coor,
            'HFenergy_cls': str(HFenergy),
            'thermochemistry_cls': thermoenergy,
            'HOMOLUMO_cls': HOMOLUMO,
            'NBOcharge_cls': NBOcharge,
            'MKcharge_spinpop_cls': MKchargespinpop,
            'img_cls': str(img),
            'dipole_moment_cls(Debye)': dipole_moment,
        }
    elif log_path.endswith('+Li+_vum.log'):
        extracted = {
            'coordinate_cls_vum': mol_coor,
            'HFenergy_cls_vum': str(HFenergy),
            'thermochemistry_cls_vum': thermoenergy,
            'HOMOLUMO_cls_vum': HOMOLUMO,
            'NBOcharge_cls_vum': NBOcharge,
            'MKcharge_spinpop_cls_vum': MKchargespinpop,
            'img_cls_vum': str(img),
            'dipole_moment_cls_vum(Debye)': dipole_moment,
        }
    else:
        extracted = {
            'coordinate_RESP2': mol_coor,
            'HFenergy_RESP2': str(HFenergy),
            'thermochemistry_RESP2': thermoenergy,
            'HOMOLUMO_RESP2': HOMOLUMO,
            'NBOcharge_RESP2': NBOcharge,
            'MKcharge_spinpop_RESP2': MKchargespinpop,
            'dipole_moment_RESP2(Debye)': dipole_moment,
        }
        if img == 1:
            extracted['img_RESP2'] = str(img)
    return extracted


def analyze_log_li(log_path: str) -> Dict[str, Any]:
    '''Analyze a Gaussian log file of Li+ and extract required properties.

    Args:
    - log_path: Path to the log file.

    Returns:
    - extracted: A mapping from property to value.
    '''
    HFenergy, thermoenergy = _energies(log_path)
    if log_path.endswith('_sol.log'):
        extracted = {
            'HFenergy_Li+_sol': str(HFenergy),
            'thermochemistry_Li+_sol': thermoenergy,
        }
    else:
        extracted = {
            'HFenergy_Li+_RESP2': str(HFenergy),
            'thermochemistry_Li+_RESP2': thermoenergy,
        }
    return extracted


def analyze_log_cls(log_path: str) -> Dict[str, Any]:
    '''Analyze a Gaussian log file and extract required properties.

    Args:
    - log_path: Path to the log file.

    Returns:
    - extracted: A mapping from property to value.
    '''
    r, mol_coor, img = _extract_coor(log_path)
    spin = _judge_spin(log_path)
    NBOcharge = _nbo_charge(log_path, spin)
    MKchargespinpop = _mk_charge_spinpop(log_path) if spin else ['NONE']
    HFenergy, thermoenergy = _energies(log_path)
    HOMOLUMO = _homo_lumo(log_path, spin)
    if log_path.endswith('+Li+_resp2.log'):
        extracted = {
            'coordinate_cls': mol_coor,
            'HFenergy_cls': str(HFenergy),
            'thermochemistry_cls': thermoenergy,
            'HOMOLUMO_cls': HOMOLUMO,
            'NBOcharge_cls': NBOcharge,
            'MKcharge_spinpop_cls': MKchargespinpop
        }
        if img == 1:
            extracted['img_cls'] = str(img)
    return extracted


def keywords(job_type: int) -> str:
    '''Decide the Gaussian keywords for the DFT job according to the job_type.

    Args:
    - job_type: The job type number.

    Returns:
    - keywords: Gaussian keywords.
    '''
    path = os.path.join(gc.kPublicFilesDir, 'keywords')
    keywords = local_fs.load_json_file(path)
    return keywords[str(job_type)]


def job(dir_run: str, mol_info: Dict, chg_spmul: str, keywords: str,
        solva_info: Dict = {}, sc: str = 'zghpc', ncores: int = 8) -> List[str]:
    '''Make a Gaussian input file and execute it in a blocking way.

    Args:
    - dir_run: The path where to run the Gaussian job.
    - mol_info: Molecule info containing name and coordinate.
    - chg_spmul: Charge and spin multiplicity of molecule.
    - keywords: Keywords to define the Gaussian job.
    - solva_info: Solvation parameters.
    - sc: Supercomputer name: zghpc, zghpc_gpu, ts1000 (up to 20230816)
    - ncores: Available cores.

    Returns:
    - files: Output files of Gaussian. By sequence, they are the log file and the checkpoint file.
    '''
    try:
        dc = solva_info['dielectric_constant_RESP2']
        soltype = 'resp2'
    except:
        try:
            dc = solva_info['dielectric_constant']
            soltype = 'cm1a'
        except:
            soltype = None
    if 'scrf' in keywords:
        if soltype != 'resp2':
            suffix = '_sol'
        else:
            suffix = '_resp2'
    else:
        suffix = '_vum'
    gjf_name = mol_info['name'] + suffix
    gjf_path = os.path.join(dir_run, gjf_name + '.gjf')
    ckpt_path = os.path.join(dir_run, gjf_name + '.chk')

    # Generate the input file
    mol_coord = mol_info['coordinate']
    with open(gjf_path, 'w') as fout:
        fout.write('%chk=' + ckpt_path + '\n')
        fout.write('%nproc=' + '%s\n' % ncores)
        fout.write(keywords + '\n\n')
        fout.write(gjf_name + '\n\n')
        fout.write(str(chg_spmul) + '\n')
        for line in mol_coord:
            fout.write(line + '\n')
        fout.write('\n')
        for i in range(len(mol_coord)):
            fout.write('%s' % (i + 1) + '\n')
        fout.write('\n')
        if 'scrf' in keywords:
            fout.write('solventname=%s\n' % solva_info['name'])
            fout.write('RSOLV=%s\n' % solva_info['radius'])
            fout.write('EPS=%s\n' % dc)
            fout.write('EpsInf=%s\n' % str(1.758))
            fout.write('\n')
        if 'nbo' in keywords:
            fout.write('$nbo bndidx archive NRT BOAO $end\n\n')

    dir_subscript = os.path.join(dir_run, 'gaussian.sh')
    with open(dir_subscript, 'w+') as f2:
        if sc == 'zghpc':
            f2.write('#!/bin/bash\n'
                     '#SBATCH -J g16\n'
                     '#SBATCH -N 1\n'
                     '#SBATCH -n %s\n' % ncores +
                     '#SBATCH -o stdout.%j\n'
                     '#SBATCH -e stderr.%j\n'
                     '\n' +
                     'module load Gaussian\n' +
                     'g16 %s.gjf\n' % gjf_name)
        elif sc == 'zghpc_gpu':
            f2.write('#!/bin/bash\n'
                     '#SBATCH -J g16\n'
                     '#SBATCH -p gpu\n'
                     '#SBATCH -N 1\n'
                     '#SBATCH -n %s\n' % ncores +
                     '#SBATCH -o stdout.%j\n'
                     '#SBATCH -e stderr.%j\n'
                     '\n' +
                     'module load Gaussian\n' +
                     'g16 %s.gjf\n' % gjf_name)
        elif sc == 'ts1000':
            f2.write('#!/bin/bash\n'
                     '#SBATCH -J g16\n'
                     '#SBATCH -p cnall\n'
                     '#SBATCH -N 1\n'
                     '#SBATCH -n %s\n' % ncores +
                     '#SBATCH -o stdout.%j\n'
                     '#SBATCH -e stderr.%j\n'
                     '\n' +
                     'g16 %s.gjf\n' % gjf_name)
        else:
            raise ValueError('Supercomputer not supported!')
    return [os.path.join(dir_run, gjf_name + '.log'), os.path.join(dir_run, gjf_name + '.chk')]


def formchk(dir_chk):
    """
    Perform formchk for .chk file and remove .chk file.
    :param dir_moldft: the path of the DFT job
    :return:
    """
    sh_command = 'formchk %s && rm %s' % (dir_chk, dir_chk)
    custom_shell.execute(sh_command, os.path.dirname(dir_chk))


def cal_resp2(chk_vum_path: str, chk_sol_path: str):
    '''Calculate the RESP2 charge by Multiwfn program and calcRESP2.sh (modified) from Tian Lu.

    Args:
    - chk_vum_path: Path to the formatted checkpoint of G16 on the molecule in vacum.
    - chk_sol_path: Path to the formatted checkpoint of G16 on the molecule in solvent.

    Returns:
    - charge_path: Path to the file containing charge assignment.
    '''
    spath = os.path.join(gc.kPublicFilesDir, 'calcRESP2.sh')
    dchk_vum, fchk_vum = os.path.split(os.path.abspath(chk_vum_path))
    dchk_sol, fchk_sol = os.path.split(os.path.abspath(chk_sol_path))
    sh_command = 'sh %s %s %s %s %s' % (spath, fchk_vum, dchk_vum, fchk_sol, dchk_sol)
    custom_shell.execute(sh_command, dchk_sol)
    charge_path = os.path.join(dchk_sol, 'RESP2.chg')
    assert os.path.exists(charge_path), 'No charge file found at: ' + charge_path
    return charge_path


def _mol_to_str(mol, L=20):
    """
    Convert the molecule structure into the crystal structure
    L: The lattice of the cell
    """
    lattice_str = Lattice.from_parameters(a=L,b=L,c=L,alpha=90,beta=90,gamma=90)
    mol_coords = []
    for i in range(len(mol)):
        mol_coords.append(list( mol.sites[i].coords))
    s = Structure(lattice=lattice_str, species=mol.species, coords=mol_coords,coords_are_cartesian=True)
    return s


def _str_to_mol(struct):
    """
    Convert the crystal structure into the molecule structure
    L: The lattice of the cell
    """
    mol_coords = []
    for i in range(len(struct)):
        mol_coords.append(list(struct.sites[i].coords))
    mol = Molecule(species=struct.species, coords=mol_coords)
    return mol


def _dedup(slist,stol=0.05):
    # Deleting the same structures
    ulist = []
    sm = StructureMatcher(stol=stol, ltol=0.05, angle_tol=1, comparator=ElementComparator())
    for i, d in enumerate(slist):
        unique = True
        for j, ud in enumerate(ulist):
            if sm.fit(d, ud):
                unique = False
                break
        if unique:
            ulist.append(d)
    return ulist


def _coords_cal_2(c0,c1,d=2.0):
    """
    c0: The coordinates of the atom to be bound
    c1: The coordinates of the neighbor atom of c0
    d: The distant between the new atom and the atom to be bonded. Default value=2.0 Å
    """
    d01 = np.linalg.norm(np.array(c0)-np.array(c1)) # ord=2, distance
    c_d = np.array(c0)+(np.array(c0)-np.array(c1))*d/d01
    return c_d  #The coordinates of the new added site


def _coords_cal_3(c0,c1,c2,d=2.0):
    """
    c0: The coordinates of the atom to be bound
    c1 & c2: The coordinates of the neighbor atom of c0
    d: The distant between the new atom and the atom to be bonded. Default value=2.0 Å
    """
    c3 = np.array(c1) + np.array(c2) - np.array(c0)
    c_d = _coords_cal_2(c0,c3,d)
    return c_d  #The coordinates of the new added site


def _coords_cal_dual_1(c0,c1,c2,d=2.35):
    """
    c0: The coordination of the neighbor atom of c1 & c2
    c1 & c2: The coordination of the atom to be bound
    d: The distant between the new atom and the atom to be bond. Defalt value=2.23 Å
    """
    c3= 3*np.array(c0) - np.array(c1) - np.array(c2)
    c_d = _coords_cal_2(c0,c3,d)
    return c_d  #The coordination of the new added site


def _coords_cal_dual_2(c0,c1,c2,c3,d=2.6):
    """
    c0: The coordination of the neighbor atom of c1 & c2 & c3
    c1 & c2 & c3: The coordination of the atom to be bound
    d: The distant between the new atom and the atom to be bond. Defalt value=2.6 Å
    """
    c3= np.array(c1) + np.array(c2) + np.array(c3) - 2*np.array(c0)
    c_d = c0 + (c3 - np.array(c0)) * d/np.linalg.norm(c3- np.array(c0))
    return c_d  #The coordination of the new added site


def _coords_cal_dual_3(ci,cj,cki,ckj,d=2.0):
    """
    ci & cj: The coordination of the birdging atoms
    cki & ckj: The coordination of the atom to be bound
    d: The distant between the new atom and the atom to be bond. Defalt value=2 Å
    """
    ckij = (np.array(cki) + np.array(ckj))/2  #The center of cki and ckj
    d_ki_kj = np.linalg.norm(np.array(cki) - np.array(ckj))
    dd = math.sqrt(d*d - d_ki_kj*d_ki_kj/4)
    c_d = _coords_cal_3(ckij,ci,cj, dd)
    return c_d  #The coordination of the new added site


def _o_sites(mol,d=2.0):
    """
    mol: the pristine molecule
    d: the distance between Li and O
    return: the coordinates of Li sites
    """
    O_index=[] #The oxygen sites in the molecule
    for i,site_i in enumerate(mol):
        if 'O' in site_i:
            O_index.append(i)
    Li_coords=[]
    if len(O_index)>0:
        for i in range(len(O_index)):
            index=O_index[i]
            c_oxy=mol[index].coords #The coordinates of the oxygen atom
            #Getting the neighbors of the atom
            neighbors=mol.get_neighbors(mol.sites[index],r=1.8) # from 1.6 to 1.8, such as Si-O bond length exceeds 1.6
            if len(neighbors) > 2:
                print('Warning! The oxygen is bound with more than two atoms')
            elif len(neighbors) == 2:
                c0=neighbors[0].coords
                c1=neighbors[1].coords
                Li_coord = _coords_cal_3(c_oxy,c0,c1,d)
                Li_coords.append(Li_coord)
            elif len(neighbors) == 1:
                c0=neighbors[0].coords
                Li_coord = _coords_cal_2(c_oxy,c0,d)
                Li_coords.append(Li_coord)
            else:
                print('Warning! There is not neighbors found with 1.6 Å of the oxygen')
    return Li_coords


def _f_sites(mol,d=2.0):
    """
    mol: the pristine molecule
    d: the distance between Li and F
    return: the coordinates of Li sites
    """
    F_index=[] #The fluorine sites in the molecule
    for i,site_i in enumerate(mol):
        if 'F' in site_i:
            F_index.append(i)
    Li_coords=[]
    if len(F_index)>0:
        for i in range(len(F_index)):
            index=F_index[i]
            c_flu=mol[index].coords #The coordinates of the fluorine atom
            #Getting the neighbors of the atom
            neighbors=mol.get_neighbors(mol.sites[index],r=1.6)
            if len(neighbors) > 1:
                print('Warning! The fluorine is bound with more than one atoms')
            elif len(neighbors) == 1:
                c0=neighbors[0].coords
                Li_coord = _coords_cal_2(c_flu,c0,d)
                Li_coords.append(Li_coord)
            else:
                print('Warning! There is not neighbors found with 1.6 Å of the fluorine')
    return Li_coords


def _n_sites(mol,d=2.0):
    """
    mol: the pristine molecule
    d: the distance between Li and N
    return: the coordinates of Li sites
    """
    N_index=[] #The fluorine sites in the molecule
    for i,site_i in enumerate(mol):
        if 'N' in site_i:
            N_index.append(i)
    Li_coords=[]
    if len(N_index)>0:
        for i in range(len(N_index)):
            index=N_index[i]
            c_flu=mol[index].coords #The coordinates of the fluorine atom
            #Getting the neighbors of the atom
            neighbors=mol.get_neighbors(mol.sites[index],r=1.6)
            if len(neighbors) > 1:
                print('Warning! The fluorine is bound with more than one atoms')
            elif len(neighbors) == 1:
                c0=neighbors[0].coords
                Li_coord = _coords_cal_2(c_flu,c0,d)
                Li_coords.append(Li_coord)
            else:
                print('Warning! There is not neighbors found with 1.6 Å of the fluorine')
    return Li_coords


def _multi_interaction_sites(mol_with_H):
    """
    mol_with_H: the pristine molecule
    return: the coordinations of Li sites
    Li will be bound to at least two O/F/N atoms.
    """
    Li_coords = []
    mol = deepcopy(mol_with_H)
    mol.remove_species('H')  # Deleting the H atoms

    # The interaction sites are bound to a same atom, e.x. -CF2-
    for i in range(len(mol)):
        # Getting the neighbors of the atom, except for H atoms
        neighbors = mol.get_neighbors(mol.sites[i], r=1.6)
        neighbors_index = []  # The index of the neighbors
        for j, n_j in enumerate(neighbors):
            if 'O' in n_j or 'F' in n_j or 'N' in n_j:
                for k, site_k in enumerate(mol):
                    if n_j.distance(site_k) < 0.0001:
                        neighbors_index.append(k)  # The neighboring O, F, N atom index

        # Atom_i is bound with two O/F/N atoms
        if len(neighbors_index) == 2:
            c0 = mol[i].coords
            c1 = mol[neighbors_index[0]].coords
            c2 = mol[neighbors_index[1]].coords
            Li_coord = _coords_cal_dual_1(c0, c1, c2, d=2.35)  # The d value should be checked
            Li_coords.append(Li_coord)

        # Atom_i is bound with three O/F/N atoms
        if len(neighbors_index) == 3:
            c0 = mol[i].coords
            c1 = mol[neighbors_index[0]].coords
            c2 = mol[neighbors_index[1]].coords
            c3 = mol[neighbors_index[2]].coords
            Li_coord = _coords_cal_dual_2(c0, c1, c2, c3, d=2.6)  # The d value should be checked
            Li_coords.append(Li_coord)

            # Interact with only two atoms of the three neighbors
            Li_coords.append(_coords_cal_dual_1(c0, c1, c2, d=2.4))
            Li_coords.append(_coords_cal_dual_1(c0, c1, c3, d=2.4))
            Li_coords.append(_coords_cal_dual_1(c0, c3, c2, d=2.4))

    # The interaction sites are bound to two neighboring atoms, e.x. -CFH-C(O)-
    for i in range(len(mol)):
        for j in range(i + 1, len(mol)):
            if mol.get_distance(i, j) < 1.6:  # Getting the two neighboring sites
                neighbors_i = mol.get_neighbors(mol.sites[i], r=1.6)
                neighbors_j = mol.get_neighbors(mol.sites[j], r=1.6)
                neighbors_index_i = []
                neighbors_index_j = []
                for ii, n_i in enumerate(neighbors_i):
                    if 'O' in n_i or 'F' in n_i or 'N' in n_i:
                        for k, site_k in enumerate(mol):
                            if n_i.distance(site_k) < 0.0001 and k != j:
                                # The neighboring O, F, N atom index, but excluding the site_j
                                neighbors_index_i.append(k)
                for jj, n_j in enumerate(neighbors_j):
                    if 'O' in n_j or 'F' in n_j or 'N' in n_j:
                        for k, site_k in enumerate(mol):
                            if n_j.distance(site_k) < 0.0001 and k != i:
                                # The neighboring O, F, N atom index, but excluding the site_i
                                neighbors_index_j.append(k)

                                # Getting the ki-i-j-kj case and d(ki,kj) < 3.6 A
                # If d(ki,kj) > 4 A, the Coords_cal_dual_3 function will poduce an error.
                # Besides (i,kj) and (ki,j) should not be bound
                for ni, ki in enumerate(neighbors_index_i):
                    for nj, kj in enumerate(neighbors_index_j):
                        if mol.get_distance(ki, kj) < 3.6 \
                                and mol.get_distance(i, kj) > 1.6 \
                                and mol.get_distance(j, ki) > 1.6:
                            ci = mol[i].coords
                            cj = mol[j].coords
                            cki = mol[ki].coords
                            ckj = mol[kj].coords
                            Li_coord = _coords_cal_dual_3(ci, cj, cki, ckj, d=2.4)  # The d value should be checked
                            Li_coords.append(Li_coord)
    return Li_coords


def _distance_checking(mol,d_Li_H = 1.5, d_Li_X = 1.8):
    """
    This funtion aims to check the distance between Li and the rest sites
    mol: the pristine molecule with a Li site
    d_Li_H: the distance between Li and H, default = 1.5 A
    d_Li_X: the distance between Li and X, X =C/O/F/N, default = 1.8 A
    return: True or false
    """
    site_Li = []
    site_H = []
    site_rest = []
    check = True
    for i,site_i in enumerate(mol):
        if 'H' in site_i:
            site_H.append(i)
        elif 'Li' in site_i:
            site_Li.append(i)
        else:
            site_rest.append(i)
    for i,site_i in enumerate(site_Li):
        for j,site_j in enumerate(site_H):
            if mol.get_distance(site_i,site_j) < d_Li_H:
                check =False
        for k,site_k in enumerate(site_rest):
            if mol.get_distance(site_i,site_k) < d_Li_X:
                check =False
    return check


def _append_li_sites(mol):
    """
    Add a Li atom to pristine molecules
    Currently, only N, O and F sites are considered as the binding sites #20230817 N added
    Return: The lists of new molecules with a Li
    """
    Li_O = _o_sites(mol, d=2.0)
    Li_F = _f_sites(mol, d=2.0)
    Li_N = _n_sites(mol, d=2.0)
    Li_multi = _multi_interaction_sites(mol)
    Li_sites = Li_O + Li_F + Li_N + Li_multi
    new_mols = []
    structs = []
    # The molecule is first converted into structure and dedup
    # and then converted back to molecule
    for i, site_i in enumerate(Li_sites):
        new_mol = deepcopy(mol)
        new_mol.append(species='Li', coords=site_i)
        if _distance_checking(new_mol) == True:
            # Ensure the structure is reasonable:
            structs.append(_mol_to_str(new_mol))
    structs_dep = _dedup(structs, stol=0.05)
    for i, struct_i in enumerate(structs_dep):
        new_mols.append(_str_to_mol(struct_i))
    return new_mols


def li_coor_sampling(mol_info):
    mol_coor = mol_info['coordinate']
    name = mol_info['name']
    dir_tmpxyz = 'mol_coor-%s.xyz' % name
    with open(dir_tmpxyz, 'w') as f:
        f.write(str(len(mol_coor)) + '\n')
        f.write('temporary xyz file from code' + '\n')
        for j in mol_coor:
            f.write(j + '\n')
    mol = Molecule.from_file(dir_tmpxyz)
    os.system('rm %s' % dir_tmpxyz)
    new_mols = _append_li_sites(mol)
    mol_coor_collection = []
    for i, mol_i in enumerate(new_mols):
        # mol_i.to(filename='%s+Li-%s.xyz' % (name, i))
        ss = _str_to_mol(mol_i)
        species = ss.species
        for s in range(len(species)):
            species[s] = str(species[s]).split()[-1]
        # print(species)
        coords = ss.cart_coords
        mol_coords = []
        for count, coord in enumerate(coords):
            coord_str = '     '.join(str(round(c, 4)) for c in coord)
            specie_coord = '   '.join([species[count], coord_str])
            # print(specie_coord)
            mol_coords.append(specie_coord)
        mol_coor_collection.append(mol_coords)
    return mol_coor_collection


__all__ = [
    'analyze_log',
    'analyze_log_li',
    'cal_resp2',
    'job',
    'keywords',
    'li_coor_sampling'
]
