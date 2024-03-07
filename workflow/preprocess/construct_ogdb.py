import os
import shutil
import time
import pandas as pd
import json
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
import random
import re


def convert(src, dst):
    """
    Construct the original molecular database.

    Args:
    - src: the file containing id and smile of molecules
    """
    info = pd.read_csv(src, header=None)
    error_file = os.path.join(os.path.abspath(os.path.dirname(dst)), "convert_error.txt")
    with open(error_file, "a") as f:
        f.write(
            "==============================\n" +
            time.strftime("%Y.%m.%d %H:%M:%S", time.localtime(time.time())) + "\n" +
            f"FULL LIST of {os.path.basename(dst)}" + "\n" +
            str(info) + "\n" +
            f"ERROR LIST of {os.path.basename(dst)}" + "\n" +
            "==============================\n"
        )
    #print(info)
    for i in info.index:
        infoi = info.iloc[i, 0].split()
        # print(infoi)
        id = infoi[0]
        smi = infoi[1]
        mol_woH = Chem.MolFromSmiles(smi)
        mol = Chem.AddHs(mol_woH)  # add H atoms
        try:
            AllChem.EmbedMolecule(mol, randomSeed=1)
            AllChem.MMFFOptimizeMolecule(mol)
        except Exception:
            try:
                AllChem.EmbedMolecule(mol, useRandomCoords=True)
                AllChem.MMFFOptimizeMolecule(mol)
                # with open('random.txt', 'a') as fr:
                #    fr.write(id + '\n')
            except Exception:
                try:
                    AllChem.EmbedMolecule(mol, maxAttempts=5000, useRandomCoords=True)
                    AllChem.MMFFOptimizeMolecule(mol)
                    # with open('random+max.txt', 'a') as fr:
                    #    fr.write(id + '\n')
                except Exception:
                    with open(error_file, 'a') as f:
                        f.write(id + ' ' + smi + '\n')
                    continue
        a = Chem.MolToXYZBlock(mol)  # string type information
        if list(a)[1] != '\n':
            b = a[4:]  # exclude the first 4 characters in a if the total atom number > 9
            c = b.split('\n')
        else:
            b = a[3:]  # exclude the first 3 characters in a if the total atom number < 10
            c = b.split('\n')
        del c[-1]
        # del c[0]
        with open(os.path.join(dst, f'{str(id)}.mol'), 'w') as f:
            f.write(json.dumps({'SMILE': smi, 'name': str(id), 'coordinate': c}, indent=1))


def convert_new(src, dst, num):
    """
    Construct the original molecular database.

    Args:
    - src: The path of input .dat files.
    - dst: The path where to write the .mol files.
    - num: The number of molecules in a group.
    """
    global error_file
    error_file = os.path.join(src, "convert_error.txt")
    epids = []
    for dat in os.listdir(src):
        if dat.endswith('.dat'):
            info = pd.read_csv(os.path.join(src, dat), header=None)
            with open(error_file, "a") as f:
                f.write(
                    "==============================\n" +
                    time.strftime("%Y.%m.%d %H:%M:%S", time.localtime(time.time())) + "\n" +
                    f"FULL LIST of {os.path.basename(dst)}" + "\n" +
                    str(info) + "\n" +
                    f"ERROR LIST of {os.path.basename(dst)}" + "\n" +
                    "==============================\n"
                )
            for i in info.index:
                infoi = info.iloc[i, 0].split()
                success = 0
                failure = 0
                while success == 0 and failure < 10:  # try 10 different seeds
                    try:
                        epid = smi2coords(infoi, random.randint(1, 1000), dst)
                        epids.append(epid)
                        success += 1
                    except:
                        failure += 1
    epids.sort()
    write_config(src, os.path.basename(dst), epids, num)


def smi2coords(infoi, seed, dst):
    id = infoi[0]
    smi_ = infoi[1]
    # Standardize SMILES
    mol_ = Chem.MolFromSmiles(smi_)
    smi = Chem.MolToSmiles(mol_)
    mol_woH = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol_woH)  # add H atoms
    try:
        AllChem.EmbedMolecule(mol, randomSeed=seed)
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        try:
            AllChem.EmbedMolecule(mol, useRandomCoords=True)
            AllChem.MMFFOptimizeMolecule(mol)
            # with open('random.txt', 'a') as fr:
            #    fr.write(id + '\n')
        except Exception:
            try:
                AllChem.EmbedMolecule(mol, maxAttempts=5000, useRandomCoords=True)
                AllChem.MMFFOptimizeMolecule(mol)
                # with open('random+max.txt', 'a') as fr:
                #    fr.write(id + '\n')
            except Exception:
                with open(error_file, 'a') as f:
                    f.write(id + ' ' + smi + '\n')
                raise Exception
    a = Chem.MolToXYZBlock(mol)  # string type information
    if list(a)[1] != '\n':
        b = a[4:]  # exclude the first 4 characters in a if the total atom number > 9
        c = b.split('\n')
    else:
        b = a[3:]  # exclude the first 3 characters in a if the total atom number < 10
        c = b.split('\n')
    del c[-1]
    # del c[0]
    with open(os.path.join(dst, f'{str(id)}.mol'), 'w') as f:
        f.write(json.dumps({'SMILE': smi, 'name': str(id), 'coordinate': c}, indent=1))
    return id


def write_subscript(config_name, startid):
    #name = f'{prefix}_{suffix}'
    dir_sub = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'SUB', f'sub.{re.split("-", config_name, 1)[1]}')
    with open(dir_sub, 'w+') as f:
        if args.sc == 'zghpc':
            f.write('#!/bin/bash\n'
                    '#SBATCH -J py\n'
                    f'#SBATCH -w {args.w}\n'
                    '#SBATCH -N 1\n'
                    '#SBATCH -n 2\n'
                    f'#SBATCH -o {re.split("-", config_name, 1)[1]}.%j.stdout\n'
                    #f'#SBATCH -e {re.split("-", config_name, 1)[1]}.%j.stderr\n'
                    '\n'
                    'module load Gaussian\n'
                    'module load BOSS\n'
                    'module load Multiwfn\n'
                    'module load Packmol\n'
                    'module load VMD\n'
                    f'python {os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "codes", "RUN_start.py")} '
                    f'-config {config_name}.json '
                    f'-start {startid} '
                    '2>&1\n')
        elif args.sc == 'zghpc_gpu':
            f.write('#!/bin/bash\n'
                    '#SBATCH -J py\n'
                    '#SBATCH -p gpu\n'
                    '#SBATCH -N 1\n'
                    '#SBATCH -n 2\n'
                    f'#SBATCH -o {re.split("-", config_name, 1)[1]}.%j.stdout\n'
                    #f'#SBATCH -e {re.split("-", config_name, 1)[1]}.%j.stderr\n'
                    '\n'
                    'module load Gaussian\n'
                    'module load BOSS\n'
                    'module load Multiwfn\n'
                    'module load Packmol\n'
                    'module load VMD\n'
                    f'python {os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "codes", "RUN_start.py")} '
                    f'-config {config_name}.json '
                    f'-start {startid} '
                    '2>&1\n')
        elif args.sc == 'ts1000':
            f.write('#!/bin/bash\n'
                    '#SBATCH -J py\n'
                    '#SBATCH -p cnall\n'
                    '#SBATCH -N 1\n'
                    '#SBATCH -n 2\n'
                    f'#SBATCH -o {re.split("-", config_name, 1)[1]}.%j.stdout\n'
                    #f'#SBATCH -e {re.split("-", config_name, 1)[1]}.%j.stderr\n'
                    '\n'
                    'module load soft/anaconda3/config\n'
                    'source /apps/soft/anaconda3/etc/profile.d/conda.sh\n'
                    'conda activate yaonan\n'
                    'module load soft/vmd/vmd-1.9.3\n'
                    'module load soft/packmol\n'
                    f'python {os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "codes", "RUN_start.py")} '
                    f'-config {config_name}.json '
                    f'-start {startid} '
                    '2>&1\n')


def write_batchsub(config_name, count):
    dir_batchsub = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'SUB', f'batchsub.{config_name}')
    with open(dir_batchsub, 'w+') as f:
        f.write(
            '#!/bin/bash\n'
            f'prefix={config_name}\n'
            f'for id in $(seq 0 {count}); do\n'
            ' echo $id\n'
            ' sbatch sub.${prefix}_${id}\n'
            'done\n'
        )


def write_config(src, config_name, molecules, num):
    """
    Write the config file.

    Args:
    - src: The path where to write the config file.
    - config_name: The config name.
    - molecules: The list of molecules.
    - num: The number of molecules in a group.
    """
    ogdb = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'ogdb', config_name)
    moldb = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'moldb', config_name)
    index = 0
    remaining = len(molecules)
    count = 0
    if args.salt is None:
        raise ValueError('Salt type cannot be None!')
    while remaining > 0:
        with open(os.path.join(src, f"epoch_config-{config_name}_{count}.json"), 'w') as f:
            f.write(json.dumps(
                {
                    "ogdb": ogdb,
                    "moldb": moldb,
                    "molecules": molecules[:num],
                    "index": index,
                    "salt": args.salt,
                    "supercomputer": args.sc
                },
                indent=1))
        write_subscript(f'epoch_config-{config_name}_{count}', args.taskid)
        shutil.copyfile(os.path.join(src, f"epoch_config-{config_name}_{count}.json"),
                        os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'SUB', f"epoch_config-{config_name}_{count}.json"))
        del molecules[:num]
        count += 1
        remaining = len(molecules)
    write_batchsub(config_name, count - 1)


def main(args):
    if not os.path.exists(args.o):
        os.mkdir(args.o)
    convert_new(args.i, args.o, args.num)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help="The path of input .dat files.", type=str, default=None)
    parser.add_argument('-o', help="The path to output converted files.", type=str, default=None)
    parser.add_argument('--num', help="The number of molecules in a group.", type=int, default=None)
    parser.add_argument('-n', help="The config name.", type=str, default=None)
    parser.add_argument('--salt', help="Salt type", type=str, default=None)
    parser.add_argument('--sc', help="Supercomputer name", type=str, default=None)
    parser.add_argument('-w', help="The node list.", type=str, default=None)
    parser.add_argument('--taskid', help="The starting task id.", type=int, default=1)
    args = parser.parse_args()
    main(args)

