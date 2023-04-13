import os
import time
import pandas as pd
import json
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
import random


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
    smi = infoi[1]
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


def write_config(src, config_name, molecules, num):
    """
    Write the config file.

    Args:
    - src: The path where to write the config file.
    - config_name: The config name.
    - molecules: The list of molecules.
    - num: The number of molecules in a group.
    """
    ogdb = f"/home/yaonan/Electrolyte_Project/ogdb/{config_name}"
    moldb = f"/home/yaonan/Electrolyte_Project/moldb/{config_name}"
    index = 0
    salt = "LiFSI"
    remaining = len(molecules)
    count = 0
    while remaining > 0:
        with open(os.path.join(src, f"epoch_config-{config_name}_{count}.json"), 'w') as f:
            f.write(json.dumps(
                {
                    "ogdb": ogdb,
                    "moldb": moldb,
                    "molecules": molecules[:num],
                    "index": index,
                    "salt": salt,
                },
                indent=1))
        del molecules[:num]
        count += 1
        remaining = len(molecules)


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
    args = parser.parse_args()
    main(args)

