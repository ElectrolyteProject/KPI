#coding:utf-8
import functools as ft
import glob
import logging
import re
import multiprocessing as mp
import shutil
import numpy as np
import os
import pandas as pd
import random
import time
from scipy.integrate import trapz
import json
import string
from typing import Dict, List, Tuple

from pymatgen.io.lammps.data import LammpsData, CombinedData
import MDAnalysis as mda
import MDAnalysis.analysis.msd as msd
import mdtraj as mdt

import global_config as gc
import assist
import custom_shell, local_fs


def boxsize(mol_info: Dict, mol_num: int) -> float:
    '''Estimate the box size of MD model according to the radius calculated by DFT_vum job.

    Args:
    - mol_info: A dict containing molecule metadata.
    - mol_num: Number of molecules in box.

    Returns:
    - box_size: Length of the cubic simulation box.
    '''
    r = float(mol_info['radius'])
    mol_volume = 4/3*np.pi*r*r*r
    box_size = round(pow(mol_volume*mol_num, 1/3), 2)
    return float(box_size)


def def_model(mol_name: str, mol_num: str) -> Tuple[str, List[str], List[str]]:
    '''Define the model (i.e., detailed compositions) for the MD job.

    Args:
    - mol_name: the name/id of the molecule/cluster (for cluster, this is usually '%s Li')
    - mol_num: the number of the molecule/each molecule in cluster (for cluster, this is usually '1 1')

    Returns:
    - model: Combined name of the MD model.
    - mol_names: List of molecule names.
    - mol_nums: Count of molecules.
    '''
    mol_names = str(mol_name).split()
    mol_nums = str(mol_num).split()
    model = mol_nums[0] + mol_names[0]
    for i in range(1, len(mol_names)):
        model = model + '+' + mol_nums[i] + mol_names[i]
    return model, mol_names, mol_nums


def mkfile_xyz(dir_run, mol_info: Dict) -> str:
    '''Make the .xyz file of a molecule/cluster to generate the .xyz file of a model for the MD job.

    Args:
    - mol_info: A dict containing molecule metadata.

    Returns:
    - path: Path to the `.xyz` file.
    '''
    xyz_path = os.path.join(dir_run, mol_info['name'] + '.xyz')
    mol_coor = mol_info['coordinate']
    with open(xyz_path, 'w') as fout:
        fout.write(str(len(mol_coor)) + '\n')
        fout.write('PDB File\n')
        for line in mol_coor:
            fout.write(line + '\n')
    return xyz_path


def mkfile_lmp(dir_run, mol_info: Dict) -> str:
    '''Make the molecular force field file (i.e., .lmp file) for the MD job by BOSS program and LigParGen code.

    Args:
    - mol_info: A dict containing molecule metadata.

    Returns:
    - lmp_path: Path to the `.lmp` file.
    '''
    mol_name = mol_info['name']
    mol_coor = mol_info['coordinate']

    # Generate .pdb file
    pdb_path = os.path.join(dir_run, mol_name + '.pdb')
    with open(pdb_path, 'w') as fout:
        fout.write('REMARK   1 File created by mkfile_lmp\n')
        for i in range(len(mol_coor)):
            line = mol_coor[i].split()
            line[1] = '%.3f' % float(line[1])  # %.3f means the number is rounded to 3 decimal places
            line[2] = '%.3f' % float(line[2])
            line[3] = '%.3f' % float(line[3])
            if i + 1 < 10:
                line.insert(0, 'HETATM    %s' % (i + 1))
            else:
                line.insert(0, 'HETATM   %s' % (i + 1))
            # Character length of element symbol = 1 or 2 (e.g., H, Li)
            if len(line[1]) == 1:
                line.insert(2, '        0    ')
            else:
                line.insert(2, '       0    ')
            linej0 = '  '.join(line[0:2])
            list0 = [linej0, line[2]]
            linej = '   '.join(list0)
            for j in range(3, 6):
                list = [linej, line[j]]
                if '-' in line[j]:
                    linej = '  '.join(list)
                else:
                    linej = '   '.join(list)
            fout.write(str(linej) + '\n')
        fout.write('END\n')

    # Generate .lmp file
    #sh_command = 'LigParGen -r %s -p %s.pdb -c 0 -pa .' % (mol_name, mol_name)
    # Generate .lmp file by latest ligpargen code from Israel
    sh_command = 'ligpargen -r %s -i %s -c 0 -o 0 -cgen CM1A -p %s' % (mol_name,
                                                                       pdb_path,
                                                                       os.path.join(dir_run, 'lpgtmp'))
    try:
        custom_shell.execute(sh_command, dir_run)
    except custom_shell.ShellError:
        custom_shell.execute(sh_command, dir_run)
    ff_fname = os.path.join(dir_run, mol_name + '.lmp')
    shutil.copyfile(os.path.join(os.path.join(dir_run, 'lpgtmp'), mol_name + '.lammps.lmp'), ff_fname)
    count = 10  # Wait 10 seconds as LigParGen works in an ASYNC-like style
    while count > 0 and not os.path.exists(ff_fname) and count:
        time.sleep(1)
        count -= 1
    shutil.rmtree(os.path.join(dir_run, 'lpgtmp'), ignore_errors=True)
    # check_lmp(ff_fname)
    return ff_fname


def check_lmp(lmp_path):
    with open(lmp_path, 'r') as f:
        infos = f.readlines()
    blank_types = []
    for info in infos[8:13]:
        if info.split()[0] == '0':
            blank_types.append(info.split()[-2])
    if len(blank_types) > 0:
        del_indexes = []
        for i, info in enumerate(infos):
            for blank_type in blank_types:
                if info.startswith(blank_type.capitalize()) \
                        or info.startswith(blank_type.capitalize() + 's'):
                    del_indexes.append(i)
        for del_index in reversed(del_indexes):
            infos.pop(del_index)
        with open(lmp_path, 'w') as f:
            f.writelines(infos)


def modfile_lmp(dir_run, mol_name: str, charge_path: str, lmpog_path: str) -> str:
    '''Modify the CM1A charge in .lmp file to RESP2 charge.

    Args:
    - mol_name: the name/id of the molecule
    - charge_path: Path to the file containing charge assignment.
    - lmpog_path: Path to the original `.lmp` file.

    Returns:
    - lmp_path: Path to the modified `.lmp` file.
    '''
    try:
        with open(charge_path, 'r') as f:
            infos = f.readlines()
    except:
        print('*******************************')
        print('RESP2.chg file NOT found!')
        print('*******************************')
        raise FileExistsError
    chglist1 = []
    atomlist = []
    for i in infos:
        chg1 = round(float(i.split()[-1]), 4)
        chglist1.append(chg1)
        atomlist.append(i.split()[0])
    sumchg = round(sum(chglist1), 4)
    exatoms = ['F', 'O', 'N', 'S', 'H']
    if sumchg == 0:
        chglist = chglist1
    else:
        for index, atom in enumerate(atomlist):
            if atom not in exatoms:
                chglist1[index] = chglist1[index] - sumchg
                break
        chglist = chglist1
    with open(lmpog_path, 'r') as f:
        infos = f.readlines()
    for l in range(len(infos)):
        if 'Atoms' in infos[l]:
            b = l + 2  # beginning line of "Atoms" part
        elif 'Bonds' in infos[l]:
            e = l - 2  # ending line of "Atoms" part
    count = 0
    for a in range(b, e + 1):
        atominfo = infos[a].split()
        atominfo[3] = str(chglist[count])  # replace charge
        atominfo.append('\n')
        atominfo = '\t\t'.join(atominfo)
        infos[a] = atominfo
        count += 1
    lmp_path = os.path.join(dir_run, mol_name + '.lmp')
    with open(lmp_path, 'w') as f:
        f.writelines(infos)
    return lmp_path


def cpfile(file: str) -> str:
    '''Copy file from PublicFilesDir.

    Args:
    - file: File name of the file to be copied.

    Returns:
    - file: Path to the copied file.
    '''
    shutil.copyfile(os.path.join(gc.kPublicFilesDir, os.path.basename(file)), file)
    return file


def mkfile_modelxyz(dir_run, model: str, mol_paths: List[str], mol_nums: List[str], box_size: float, seed: int = None) -> str:
    '''Make the .xyz file of a model containing multiple molecules for the MD job by Packmol program.

    Args:
    - model: Model name.
    - mol_paths: Paths of `.xyz` to used molecules.
    - mol_nums: Numbers of each used molecule.
    - box_size: Length of the cubic simulation box.
    - seed: Seed keyword in packmol.

    Returns:
    - xyz_path: Path to the `.xyz` file containing lots of molecules.
    '''
    pkmin, pkmout = 'packmol.inp', 'packmol.oup'
    xyz_path = os.path.join(dir_run, model + '.xyz')
    if not seed:
        seed = 1001
    with open(os.path.join(dir_run, pkmin), 'w') as fout:
        fout.write('seed %s\n' % str(seed))
        fout.write('tolerance 2.0\n')
        fout.write('filetype xyz\n')
        fout.write('output ' + xyz_path + '\n')
        for i in range(len(mol_paths)):
            fout.write('structure ' + mol_paths[i] + '\n')
            fout.write('  number ' + mol_nums[i] + '\n')
            fout.write('  inside box 0.0 0.0 0.0 {l} {l} {l}\n'.format(l=box_size))
            fout.write('end structure\n')
    sh_command = 'packmol <%s >%s' % (pkmin, pkmout)
    custom_shell.execute(sh_command, dir_run)
    with open(os.path.join(dir_run, pkmout), 'r') as f1:
        info = f1.read()
    assert 'Success!' in info, 'Packmol reported failure!'
    return xyz_path


def mkfile_data(dir_run, model: str, model_path: str,
                mol_names: List[str], mol_paths: List[str], mol_nums: List[str]) -> str:
    '''Make a .data file for the MD job.

    Args:
    - model: Model name.
    - model_path: Path to `.xyz` containing packed molecules.
    - mol_names: List of molecule names.
    - mol_paths: Paths of `.lmp` to used molecules.
    - mol_nums: Numbers of each used molecule.

    Returns:
    - data_path: Path to the LAMMPS data file.
    '''
    mols = []
    list_of_numbers = []
    for path, num in zip(mol_paths, mol_nums):
        mols.append(LammpsData.from_file(path))
        list_of_numbers.append(int(num))
    coordinates = CombinedData.parse_xyz(model_path)
    combined = CombinedData.from_lammpsdata(mols, mol_names, list_of_numbers, coordinates)
    data_path = os.path.join(dir_run, model + '.data')
    combined.write_file(data_path)
    return data_path


def modfile_data(dir_run, model: str, lmp_path: str, opdata_path: str) -> str:
    '''Modify the CM1A charge in op.data file (from molmd1) to RESP2 charge.

    Args:
    - model: Model name.
    - lmp_path: Path to the RESP2 `.lmp` file.
    - opdata_path: Path to the `op.data` file.

    Returns:
    - data_path: Path to the modified `.data` file.
    '''
    # extract charge from .lmp file
    with open(lmp_path, 'r') as flmp:
        infos_lmp = flmp.readlines()
    for l, line in enumerate(infos_lmp):
        if line.startswith('Atoms'):
            ab = l + 2  # beginning line of "Atoms" part
        elif line.startswith('Bonds'):
            ae = l - 2  # ending line of "Atoms" part
    charge = []
    for a in range(ab, ae + 1):
        charge.append(infos_lmp[a].split()[3])
    # replace charge in op.data file according to atom types
    with open(opdata_path, 'r') as fdata:
        infos_data = fdata.readlines()
    for l, line in enumerate(infos_data):
        if line.startswith('Atoms'):
            ab = l + 2  # beginning line of "Atoms" part
        elif line.startswith('Velocities'):
            ae = l - 2  # ending line of "Atoms" part
            vb = l  # beginning line of "Velocities" part
        elif line.startswith('Bonds'):
            ve = l  # ending line of "Velocities" part
    for a in range(ab, ae + 1):
        atominfo = infos_data[a].split()
        del atominfo[-3:]
        atomtype = int(atominfo[2])
        atominfo[3] = charge[atomtype - 1]
        atominfo.append('\n')
        infos_data[a] = ' '.join(atominfo)
    del infos_data[vb:ve]
    data_path = os.path.join(dir_run, model + '.data')
    with open(data_path, 'w') as f:
        f.writelines(infos_data)
    # os.remove(opdata_path)
    return data_path


def simulation(simulation_fname: str) -> Tuple[float, float, List[int]]:
    '''Tailor the simulation temperature and runtime of the MD job.

    Args:
    - simulation_fname: Name of the json file containing temperature and runtime information.

    Returns:
    - temperature1: Temperature low bound.
    - temperature2: Temperature high bound.
    - runtime: Steps to run in each stage.
    '''
    path = os.path.join(gc.kPublicFilesDir, 'staged', simulation_fname)
    simul_info = local_fs.load_yaml_file(path)
    temperature1 = simul_info['TEMPERATURE1']
    temperature2 = simul_info['TEMPERATURE2']
    runtime = simul_info['RUNTIME']
    return temperature1, temperature2, runtime


def lammps_fix_commands(compute_type: str):
    '''
    Get related `compute` & `fix` commands of the given type.

    Args:
    - compute_type: The compute keywords of the MD job, na:None/dm:dipole moment/press:pressure/msd:msd/com:com.

    Returns:
    - commands: Lammps commands.
    '''
    compute_command_info = {
        "na": [
            "\n"
        ],
        "dm": [
            "compute dpall1 all chunk/atom molecule\n",
            "compute dpall2 all dipole/chunk dpall1\n",
            "fix dpall all ave/time 1000 1 1000 c_dpall2[*] file $$ModelName_dipole.out mode vector\n"
        ],
        "press": [
            "variable pxy equal pxy\n",
            "variable pxz equal pxz\n",
            "variable pyz equal pyz\n",
            "fix pressure all ave/time 1 1 1 v_pxy v_pxz v_pyz file $$ModelName_pressure.out\n"
        ],
        "msd": [
            "group cation type $$CationType\n",
            "compute cation1 cation chunk/atom molecule\n",
            "compute cation2 cation msd/chunk cation1\n",
            "group anion type $$AnionType\n",
            "compute anion1 anion chunk/atom molecule\n",
            "compute anion2 anion msd/chunk anion1\n",
            "fix msdcation cation ave/time 1 1 1000 c_cation2[*] file Li_msd.out mode vector\n",
            "fix msdanion anion ave/time 1 1 1000 c_anion2[*] file $$AnionName_msd.out mode vector\n"
        ],
        "com": [
            "TOBEDETERMINED"
        ]
    }
    compute_command_list = []
    for c in compute_type.split():
        if c in compute_command_info:
            ccjoin = ''.join(compute_command_info[c])
            compute_command_list.append(ccjoin)
        else:
            logging.warning('Compute type %s NOT defined in compute_commands file!' % c)
            compute_command_list.append('\n')
    return '\n'.join(compute_command_list)


class LammpsTemplate(string.Template):
    """A string class for supporting $$-PascalStyleToken."""
    delimiter = '$$'
    idpattern = '[a-zA-Z][a-zA-Z0-9]*'


def md_judge(dir_run):
    """
    Judge whether the MD job terminates normally.
    :param dir_molmd: the path of the MD job
    :return: None
    """
    with open(os.path.join(dir_run, 'log.lammps'), 'r') as f:
        infos = f.readlines()
    if 'Total wall time' in infos[len(infos) - 1]:
        files = []
        lammpss = glob.glob(os.path.join(dir_run, '*.lammps'))
        files.extend(lammpss)
        frames = glob.glob(os.path.join(dir_run, '*.data'))
        files.extend(frames)
        fixes = glob.glob(os.path.join(dir_run, '*.out'))
        files.extend(fixes)
        tractories = glob.glob(os.path.join(dir_run, '*.lammpstrj'))
        files.extend(tractories)
        return files
    else:
        print('*******************************')
        print(dir_run.split(os.sep)[-1], 'is NOT terminated normally!')
        print('*******************************')
        raise Exception


def job(dir_run: str, in_template: str, model: str, compute_type: str, ncores: int = 4) -> List[str]:
    '''
    Modify `in.lammps` files for the MD job.

    Args:
    - dir_run: Path of the folder to run LAMMPS.
    - in_template: Name of the in.lammps template file.
    - model: Model name.
    - compute_type: The compute keywords of the MD job, 0:None/1:dipole/2:pressure/3:msd.
    '''
    init = local_fs.load_json_file(os.path.join(gc.kPublicFilesDir, 'lammps.init'))
    compute_commands = lammps_fix_commands(compute_type)
    temperature1, temperature2, runtime = init['temperature'][0], init['temperature'][1], init['runtime']
    to_replace = {
        'ModelName': model,
        'Temperature1': str(temperature1),
        'Temperature2': str(temperature2)
    }
    to_replace['ComputeType'] = LammpsTemplate(compute_commands).substitute(**to_replace)
    for idx, item in enumerate(runtime):
        to_replace['Runtime' + str(idx)] = str(item)
    intemp = os.path.join(gc.kPublicFilesDir, in_template)
    with open(intemp, 'r') as fin:
        raw_content = fin.read()
    filled_content = LammpsTemplate(raw_content).substitute(**to_replace)
    in_ = os.path.join(dir_run, 'in.lammps')
    with open(in_, 'w') as fout:
        fout.write(filled_content)
    dir_sub = os.path.join(dir_run, 'lammps.sh')
    with open(dir_sub, 'w+') as f2:
        f2.write('#!/bin/bash\n'
                 '#SBATCH -J lmp\n'
                 '#SBATCH -N 1\n'
                 '#SBATCH -n %s\n' % ncores +
                 '#SBATCH -o lammps.log\n'
                 '#SBATCH -e stderr.%j\n'
                 '\n'
                 'module load LAMMPS/23Jun2022\n'
                 'mpirun lmp -sf omp -pk omp 0 neigh yes -in in.lammps\n')


def caldc(in_path: str, log_path: str, dipole_path: str,
          molnumtot: int, molnumsol: int, freq=2, cutoff=5):
    '''Calculate the dielectric constant from file `dipole.out`.

    Args:
    - in_path: Path of `in.lammps`.
    - log_path: Path of `log.lammps`.
    - dipole_path: Path of `dipole.out`.
    - molnumtot: Total number of molecules including salts in box.
    - molnumsol: Total number of molecules excluding salts in box.
    - freq: Statistical frequency for calculating the dielectric constant.
    - cutoff: Cutoff distance between the cation and anion in salt.

    Returns:
    - dc: Dielectric constant.
    - T: Target temperature.
    - V: Volume of the box.
    '''
    data = []
    with open(dipole_path, 'r') as f1:
        for (l, line) in enumerate(f1.readlines()):
            if l > 2:
                ls = line.split()
                number = map(float, ls)
                data.append(number)
    data_df = pd.DataFrame(data)

    # dataframe to dielectric constant
    dpxyz = []
    co = []
    n = 0
    cycle = molnumtot + 1   # +1 represents the step information line
    for j in data_df.index:
        if j > 0:
            if (j + 1) % cycle == 0:
                dpxyz_sol = data_df.iloc[n * cycle + 1:n * cycle + 1 + molnumsol, 1:4]    # sol
                dpxyz_sol_pd = pd.DataFrame(dpxyz_sol)
                dpxyz_salt = data_df.iloc[n * cycle + 1 + molnumsol:j + 1, 1:5]  # salt
                dpxyz_salt_pd = pd.DataFrame(dpxyz_salt)
                # only calculate dipole moment of CIP salt
                for jco in dpxyz_salt_pd.index:
                    if dpxyz_salt_pd.loc[jco, 4] <= cutoff:  # cutoff = distance between cation and anion
                        co.append(dpxyz_salt_pd.loc[jco, 0:3])
                co_df = pd.DataFrame(co)
                dpxyz_all = pd.concat([dpxyz_sol_pd, co_df], axis=0)
                sum_dpxyz = np.sum(dpxyz_all, axis=0)  # 0=sum of column
                dpxyz.append(sum_dpxyz)
                n = n + 1
                co.clear()
    dpxyz_df = pd.DataFrame(dpxyz)  # dipole moment of every single step (sum of dp of all molecules in a step)
    Mxyz = []
    f = freq
    for k in dpxyz_df.index:
        if k >= (f - 1):
            dpxyz_avg = dpxyz_df.iloc[k - (f - 1):k + 1, 0:3]
            dpxyz_avg_pd = pd.DataFrame(dpxyz_avg)
            avg_dpxyz_avg = np.mean(dpxyz_avg_pd, axis=0)
            Mxyz.append(avg_dpxyz_avg)
    Mxyz_df = pd.DataFrame(Mxyz)  # average dipole moment of every freq steps
    Mxyz2 = Mxyz_df * Mxyz_df
    Mavg2 = Mxyz2.sum(axis=1)
    Mavg2_df = pd.DataFrame(Mavg2)
    dpxyz2 = dpxyz_df * dpxyz_df
    M2 = dpxyz2.sum(axis=1)
    M2avg = []
    for m in M2.index:
        if m >= (f - 1):
            M2avgm = M2.iloc[m - (f - 1):m + 1]
            M2avgm_pd = pd.DataFrame(M2avgm)
            avg1_M2avgm = np.mean(M2avgm_pd, axis=0)
            M2avg.append(avg1_M2avgm)
    M2avg_df = pd.DataFrame(M2avg)

    # dielectric constant
    M_M = M2avg_df - Mavg2_df
    T, V = logals(in_path, log_path)
    dielectric_const = 1 + (2.56e-58 * M_M) / (3 * 1.38e-23 * 8.85e-12 * T * 1e-30 * V)  # temperature, volume
    dc_mean = dielectric_const.mean()
    dc = str(round(float(dc_mean.iloc[0]), 4))
    return dc, T, V


def logals(in_path: str, log_path: str) -> Tuple[float, float]:
    '''Read `in.lammps` and `log.lammps` of the MD job to get the temperature and volume.

    Args:
    - in_path: Path of `in.lammps`.
    - log_path: Path of `log.lammps`.

    Returns:
    - temperature: Target temperature.
    - volume: Box volume.
    '''
    # Read the temperature of nvt run
    with open(in_path, 'r') as f1:
        for line in f1.readlines():
            if 'nvt' in line:
                T = float(line.split()[-2])
    # Read the volume of the cell
    with open(log_path, 'r') as f2:
        loop_num = []
        for (l_n, line) in enumerate(f2.readlines()):
            if 'Loop time' in line:
                loop_num.append(l_n)
        f2.seek(0, 0)
        for (l_m, line) in enumerate(f2.readlines()):
            if l_m == max(loop_num) - 1:
                V = float(line.split()[-1]) # PLEASE MIND THE THERMO_DUMP INFO OF LOG.LAMMPS
    return T, V


def calvis(dir_run, mol_name: str, pressure_path: str, T: float, V: float,
           freqs: List[int], interval: int, split_parts: int, ncores: int = 8) -> Tuple[list, str]:
    '''Calculate the viscosity from the pressure file.

    Args:
    - mol_name: Molecule name.
    - pressure_path: Path to the pressure file.
    - T: Target temperature.
    - V: Volume of the box.
    - freqs: Statistical frequencies (unit=fs) for averaging the pressure.
    - interval: Interval (unit=fs) to split the pressure.out (usually the NVT run).
    - split_parts: Total split parts, interval*split_parts = NVT run time.
    - ncores: Available cores.

    Returns:
    - viscosities: Viscosity of each frequency.
    - excelpath: Path to stats.
    '''
    logging.info("start cal_viscosity")
    excelpath = os.path.join(dir_run,
                             'viscosity-nvt%sps-%s-stdavg.xlsx' % (int(interval / 1000), mol_name))
    writer = pd.ExcelWriter(excelpath)
    data = pd.read_csv(pressure_path, sep=' ', skiprows=[0, 1], header=None)
    allvis = []
    avgall = []
    stdall = []
    for freq in freqs:
        for split_id in range(split_parts):
            datas = data.iloc[split_id * interval:(split_id + 1) * interval, -3:]
            avpres = datas.rolling(window=freq, axis=0).mean()
            avpres = avpres.drop(index=avpres.index[:freq-1], axis=0).reset_index().iloc[:, -3:]
            deltasq = (avpres * freq) * (avpres * freq)
            avdeltasq = deltasq.expanding(min_periods=1).mean()
            kB = 1.3806504e-23  # [J/K] Boltzmann
            atm2Pa = 101325.0
            A2m = 1.0e-10
            fs2s = 1.0e-15
            convert = atm2Pa * atm2Pa * fs2s * A2m * A2m * A2m
            viscosity_xyz = V * avdeltasq * convert / (2 * kB * T * freq)
            viscosityPas = viscosity_xyz.T.mean()
            viscosity = viscosityPas * 1000
            avgall.append(round(viscosity.mean(), 8))
            stdall.append(round(viscosity.std(), 8))
        pd.DataFrame(avgall).to_excel(writer, sheet_name='%s-avg' % freq, index=False, header=None)
        pd.DataFrame(stdall).to_excel(writer, sheet_name='%s-std' % freq, index=False, header=None)
        sum = 0
        for i in avgall:
            sum += i
        avg_avg = round((sum / len(avgall)), 8)
        stdcutofflist = [round(0.1 * avg_avg, 8)]
        for stdcutoff in stdcutofflist:
            dataall = []
            for i in range(len(stdall)):
                std = stdall[i]
                if std <= stdcutoff:
                    dataall.append(avgall[i])
            data_df = pd.DataFrame(dataall)
            data_df_sort = data_df.sort_values(by=0, axis=0)
            data_df_sort.to_excel(writer, sheet_name='%s-std%s-sort' % (freq, stdcutoff), index=False, header=None)
            midposition = int(len(dataall) / 2)
            mid = float(data_df_sort.iloc[midposition])
            allvis.append('freq%s-std%s = %s' % (freq, stdcutoff, str(mid)))
    writer.save()
    logging.info("end cal_viscosity")
    return allvis, excelpath


def _cal_viscosity(split_id: int,
                  data: pd.DataFrame, itv: int, freq: int, V: float, T: float) -> Tuple[float, float]:
    '''Calculate the viscosity.

    Args:
    - split_id: Which part to handle.
    - data: Tracked pressures.
    - itv: Interval (unit=fs) to split the pressure.out (usually the NVT run).
    - freq: Statistical frequency (unit=fs) for averaging the pressure.
    - V: Volume of the box.
    - T: Target temperature.

    Returns:
    - avg: Viscosity averaged value.
    - std: Viscosity standard error.
    '''
    print("start _cal_viscosity")
    datas = data.iloc[split_id * itv:(split_id + 1) * itv]
    datasplit = datas.reset_index()
    avpress = []
    for j in datasplit.index:
        if j >= (freq - 1):
            press_freq = datasplit.iloc[j - (freq - 1):j + 1, -3:]
            press_freq_pd = pd.DataFrame(press_freq)
            avg_press_freq = np.mean(press_freq_pd, axis=0)  # 0=sum of column
            avpress.append(avg_press_freq)
    avpress_df = pd.DataFrame(avpress)
    deltasq = (avpress_df * freq) * (avpress_df * freq)
    avdeltasq = []
    for i in deltasq.index:
        deltasq_freq = deltasq.iloc[0:i + 1]
        deltasq_freq_pd = pd.DataFrame(deltasq_freq)
        avg_deltasq_freq = np.mean(deltasq_freq_pd, axis=0)  # 0=sum of column
        avdeltasq.append(avg_deltasq_freq)
    avdeltasq_df = pd.DataFrame(avdeltasq)
    kB = 1.3806504e-23  # [J/K] Boltzmann
    atm2Pa = 101325.0
    A2m = 1.0e-10
    fs2s = 1.0e-15
    convert = atm2Pa * atm2Pa * fs2s * A2m * A2m * A2m
    viscosity_xyz = V * avdeltasq_df * convert / (2 * kB * T * freq)
    viscosityPas = viscosity_xyz.T.mean()
    viscosity = viscosityPas * 1000
    avgi = round(viscosity.mean(), 8)
    stdi = round(viscosity.std(), 8)
    print("end _cal_viscosity")
    return avgi, stdi


def mkfile_psf(dir_run, model: str) -> str:
    '''Make the .psf file (topology file) by VMD program and writePSF.tcl code.

    Args:
    - model: Model name

    Returns:
    - psf_path: Path to the `.psf` file
    '''

    tcl_path = os.path.join(gc.kPublicFilesDir, 'writePSF.tcl')
    sh_command = 'vmd -dispdev none -e %s -args %s > /dev/null 2>&1' % (tcl_path, model)
    custom_shell.execute(sh_command, dir_run)
    psf_path = os.path.join(dir_run, model + '.psf')
    return psf_path


def _atomnum(dir_run, mol_names: List[str]) -> List[int]:
    '''Obtain the list of atom numbers of compositions in MD model.

    Args:
    - mol_names: List of molecule names.

    Returns:
    - atom_nums: List of atom numbers
    '''
    atom_nums = []
    for mol in mol_names:
        xyz_path = os.path.join(dir_run, mol + '.xyz')
        with open (xyz_path, 'r') as fxyz:
            infos = fxyz.readlines()
        atom_nums.append(int(len(infos) - 2))
    return atom_nums


def _least_squares(x: np.array, y: np.array) -> Tuple[np.array, np.array]:
    '''Linearly fit the curve based on the Least Squares method.
    (e.g. fitting the mean square displacement (MSD) to calculate the diffusivity)
    Args:
    - x: x data.
    - y: y data

    Returns:
    - a: slope
    - b: intercept
    '''
    x_ = x.mean()
    y_ = y.mean()
    m = np.zeros(1)
    n = np.zeros(1)
    for i in np.arange(len(y)):
        k = (x[i]-x_)*(y[i]-y_)
        m += k
        p = np.square(x[i]-x_)
        n = n + p
    a = m/n
    b = y_ - a*x_
    return a, b


def calmsd(dir_run, model: str, mol_names: List[str]) -> Tuple[list, str]:
    '''Calculate the msd and diffusivity of each composition in the electrolyte.

    Args:
    - model: Model name.
    - mol_names: List of the name of each composition

    Returns:
    - diffs: List of the diffusivity of each composition
    - excelpath: Excel file of the msd results
    '''
    shutil.copyfile(os.path.join(dir_run, '%s_unwrapped.lammpstrj' % model),
                    os.path.join(dir_run, '%s_unwrapped.lammpsdump' % model))
    u = mda.Universe(os.path.join(dir_run, '%s_unwrapped.lammpsdump' % model))
    # a = u.select_atoms('type x', periodic=True)
    an = _atomnum(dir_run, mol_names)
    mol_type = 'type %s' % an[0]
    anion_type = 'type %s' % (sum(an) - 1)
    cation_type = 'type %s' % sum(an)
    types = [mol_type, anion_type, cation_type]
    names = [mol_names[0], 'anion', 'cation']
    excelpath = os.path.join(dir_run,
                             'mda_msd-%s.xlsx' % mol_names[0])
    writer = pd.ExcelWriter(excelpath)
    diffs = []
    for index, type in enumerate(types):
        MSD = msd.EinsteinMSD(u, select=type, msd_type='xyz', fft=True)
        MSD.run()
        msd_array = MSD.results.timeseries
        msd_df = pd.DataFrame(msd_array) / 100  # convert Ang^2 to nm^2
        msd_df.to_excel(writer, '%s' % names[index], index=False, header=None)
        # print(msd_df)
        for colnum in range(msd_df.shape[1]):
            value = msd_df.iloc[:, colnum]
            # print(value)
            rownum = value.shape[0]
            x = np.linspace(0, rownum - 1, rownum)
            y = np.array(value)
            a, b = _least_squares(x, y)
            slope = float(a[0])
            D = slope * pow(10, -6) / 6  # unit = m^2/s
            diffs.append('%s = %s' % (names[index], D))
    writer.save()
    os.remove(os.path.join(dir_run, '%s_unwrapped.lammpsdump' % model))
    return diffs, excelpath


def calrdfCN(dir_run, model, mol_names, mol_nums, V, f1=1, f2=5000, r=10, dr=0.02):
    '''Calculate the rdf and CN from the .lammpstrj file.

    Args:
    - model:
    - mol_names: List of the name of each composition
    - mol_nums: List of the number of each composition
    - V: Model volume
    - f1: First frame to calculate rdf and CN
    - f2: Last frame to calculate rdf and CN
    - r: Maximum r
    - dr: Delta r

    Returns:
    - rdfCN_list: List of rdf and CN results
    '''
    data_path = os.path.join(dir_run, model + '.data')
    an = _atomnum(dir_run, mol_names)
    mass2ele = assist.mass2element()
    with open(data_path, 'r') as fdata:
        infos = fdata.readlines()
    for index, info in enumerate(infos):
        if 'Masses' in info:
            b = index + 2
        elif 'Pair Coeffs' in info:
            e = index - 2
            break
    atmass = [infos[j].split()[1] for j in range(b, e + 1)]
    sel1 = []
    for a in range(sum(an)):
        if mass2ele[atmass[a]] == 'Li':
            sel1.append(str(a + 1))
    sel1 = [' '.join(sel1)]
    sel2 = []
    sel2_i = []
    Atom = []
    AtomNum = []
    count_element = ['N', 'O', 'F', 'P', 'S']
    for o in range(len(an)):
        start = sum(an[0:o])
        stop = sum(an[0:o + 1])
        count_list = [0, 0, 0, 0, 0]
        for k in range(len(count_element)):
            for a in range(start, stop):
                if mass2ele[atmass[a]] == count_element[k]:
                    sel2_i.append(str(a + 1))   # atom type = a + 1
                    count_list[k] += 1
            if sel2_i:
                sel2.append(' '.join(sel2_i))
                sel2_i = []
        # Atoms with the same element are classified into one group.
        for l in range(len(count_list)):
            if count_list[l] > 0:
                Atom.append('%s_%s' % (mol_names[o], count_element[l]))
                AtomNum.append('%s' % (int(mol_nums[o]) * count_list[l]))
    tj = os.path.join(dir_run, model + '.lammpstrj')
    tp = mkfile_psf(dir_run, model)
    traj = mdt.load(tj, top=tp)
    traj_extract = traj[f1:f2 + 1]
    topo = mdt.load_psf(tp)
    rdfCN_list = []
    for i in range(len(sel1)):
        sl1 = sel1[i]
        for j in range(len(sel2)):
            sl2 = sel2[j]
            select1 = topo.select('resSeq %s' % sl1)  # 'resSeq' represents the type of atom
            select2 = topo.select('resSeq %s' % sl2)
            pr1 = topo.select_pairs(selection1=select1, selection2=select2)
            rdf1 = mdt.compute_rdf(traj_extract, pr1, r_range=(0.0, r / 10),
                                  bin_width=dr / 10)  # default r=10 A deltar=0.05 A
            rdf1_df = pd.DataFrame(rdf1)
            rdf1_df_stack = rdf1_df.stack().unstack(0)
            rdf = pd.concat([rdf1_df_stack[0] * 10, rdf1_df_stack[1]], axis=1,
                            ignore_index=True)  # change unit from nm to A
            varr = np.linspace(0.01, r - 0.01, int(r / dr))
            intgr = np.array(4 * np.pi * varr ** 2 * rdf[1] * int(AtomNum[j]) / V)
            CN = []
            for i in range(len(varr)):
                CN.append(float(round(trapz(intgr[:i + 1], varr[:i + 1]), 4)))
            rdfCN = rdf[0].round(4).map(str) + ' ' + rdf[1].round(4).map(str) + ' ' + pd.DataFrame(CN)[0].map(str)
            rdfCN.loc[0] = '%s' % str(Atom[j])
            rdfCN_json = rdfCN.to_json(orient='values', indent=1)
            rdfCN_list.append(json.loads(rdfCN_json))
    return rdfCN_list


def extract_walltime(log_path: str) -> float:
    '''Extract the simulation time from log.lammps.

    Args:
    - log_path: Path to the log.lammps file.

    Returns:
    - time: Simulation time in unit hour
    '''
    with open(log_path, 'r') as flog:
        infos = flog.readlines()
    hms = [0] * 3
    for info in infos:
        if info.startswith('Total wall time'):
            for id in range(3):
                hms[id] += int(re.split('[:\n ]', info)[id - 4])
    time = round(hms[0] + hms[1] / 60 + hms[2] / 3600, 4)
    return time


__all__ = [

]