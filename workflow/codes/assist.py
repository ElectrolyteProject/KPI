import string
import os
import time
import pandas as pd

from typing import Dict

import custom_shell

class LammpsTemplate(string.Template):
    """A string class for supporting $$-PascalStyleToken."""
    delimiter = '$$'
    idpattern = '[a-zA-Z][a-zA-Z0-9]*'


def number2element() -> Dict[int, str]:
    '''Return a mapping from atomic number to element symbol.'''
    return {
        1:'H',
        2: 'He',
        3: 'Li',
        4: 'Be',
        5: 'B',
        6: 'C',
        7: 'N',
        8: 'O',
        9: 'F',
        10: 'Ne',
        11: 'Na',
        12: 'Mg',
        13: 'Al',
        14: 'Si',
        15: 'P',
        16: 'S',
        17: 'Cl',
        18: 'Ar',
        19: 'K',
        20: 'Ca',
        21: 'Sc',
        22: 'Ti',
        23: 'V',
        24: 'Cr',
        25: 'Mn',
        26: 'Fe',
        27: 'Co',
        28: 'Ni',
        29: 'Cu',
        30: 'Zn',
        31: 'Ga',
        32: 'Ge',
        33: 'As',
        34: 'Se',
        35: 'Br',
        36: 'Kr',
        37: 'Rb',
        38: 'Sr',
        39: 'Y',
        40: 'Zr',
        41: 'Nb',
        42: 'Mo',
        43: 'Tc',
        44: 'Ru',
        45: 'Rh',
        46: 'Pd',
        47: 'Ag',
        48: 'Cd',
        49: 'In',
        50: 'Sn',
        51: 'Sb',
        52: 'Te',
        53: 'I',
        54: 'Xe',
        55: 'Cs',
        56: 'Ba',
        57: 'La',
        58: 'Ce',
        59: 'Pr',
        60: 'Nd',
        61: 'Pm',
        62: 'Sm',
        63: 'Eu',
        64: 'Gd',
        65: 'Tb',
        66: 'Dy',
        67: 'Ho',
        68: 'Er',
        69: 'Tm',
        70: 'Yb',
        71: 'Lu',
        72: 'Hf',
        73: 'Ta',
        74: 'W',
        75: 'Re',
        76: 'Os',
        77: 'Ir',
        78: 'Pt',
        79: 'Au',
        80: 'Hg',
        81: 'Tl',
        82: 'Pb',
        83: 'Bi',
        84: 'Po',
        85: 'At',
        86: 'Rn',
        87: 'Fr',
        88: 'Ra',
        89: 'Ac',
        90: 'Th',
        91: 'Pa',
        92: 'U',
        93: 'Np',
        94: 'Pu',
        95: 'Am',
        96: 'Cm',
        97: 'Bk',
        98: 'Cf',
        99: 'Es',
        100: 'Fm',
        101: 'Md',
        102: 'No',
        103: 'Lr',
        104: 'Rf',
        105: 'Db',
        106: 'Sg',
        107: 'Bh',
        108: 'Hs',
        109: 'Mt',
        110: 'Ds',
        111: 'Rg',
        112: 'Cn',
        113: 'Nh',
        114: 'Fl',
        115: 'Mc',
        116: 'Lv',
        117: 'Ts',
        118: 'Og'
    }


def mass2element() -> Dict[str, str]:
    '''Return a mapping from atomic mass to element symbol.'''
    return {"1.008": "H",
     "4.002": "He",
     "6.941": "Li",
     "9.012": "Be",
     "10.811": "B",
     "10.812": "B",
     "12.011": "C",
     "14.007": "N",
     "15.999": "O",
     "18.998": "F",
     "20.179": "Ne",
     "22.989": "Na",
     "24.305": "Mg",
     "26.982": "Al",
     "28.086": "Si",
     "30.974": "P",
     "30.977": "P",
     "32.065": "S",
     "32.067": "S",
     "35.450": "Cl",
     "39.948": "Ar",
     "39.098": "K",
     "40.078": "Ca",
     "44.955": "Sc",
     "47.867": "Ti",
     "50.941": "V",
     "51.996": "Cr",
     "54.938": "Mn",
     "55.845": "Fe",
     "58.933": "Co",
     "58.693": "Ni",
     "63.546": "Cu",
     "65.390": "Zn",
     "69.723": "Ga",
     "72.64": "Ge",
     "74.921": "As",
     "78.96": "Se",
     "79.904": "Br",
     "83.798": "Kr",
     "85.467": "Rb",
     "87.62": "Sr",
     "88.905": "Y",
     "91.224": "Zr",
     "92.906": "Nb",
     "95.94": "Mo",
     "97.9072": "Tc",
     "101.07": "Ru",
     "102.905": "Rh",
     "106.42": "Pd",
     "107.868": "Ag",
     "112.411": "Cd",
     "114.818": "In",
     "118.71": "Sn",
     "121.76": "Sb",
     "127.6": "Te",
     "126.904": "I",
     "131.293": "Xe",
     "132.905": "Cs",
     "137.327": "Ba",
     "138.905": "La",
     "140.116": "Ce",
     "140.907": "Pr",
     "144.242": "Nd",
     "145": "Pm",
     "150.36": "Sm",
     "151.964": "Eu",
     "157.25": "Gd",
     "158.925": "Tb",
     "162.5": "Dy",
     "164.93": "Ho",
     "167.259": "Er",
     "168.934": "Tm",
     "173.04": "Yb",
     "174.967": "Lu",
     "178.49": "Hf",
     "180.947": "Ta",
     "183.84": "W",
     "186.207": "Re",
     "190.23": "Os",
     "192.217": "Ir",
     "195.084": "Pt",
     "196.966": "Au",
     "200.59": "Hg",
     "204.383": "Tl",
     "207.2": "Pb",
     "208.98": "Bi",
     "208.982": "Po",
     "209.987": "At",
     "222.017": "Rn",
     "223": "Fr",
     "226": "Ra",
     "227": "Ac",
     "232.038": "Th",
     "231.035": "Pa",
     "238.028": "U",
     "238.8486": "Np",
     "242.8798": "Pu",
     "244.8594": "Am",
     "246.911": "Cm",
     "248.9266": "Bk",
     "252.9578": "Cf",
     "253.9656": "Es",
     "259.0046": "Fm",
     "260.0124": "Md",
     "261.0202": "No",
     "264.0436": "Lr",
     "269.0826": "Rf",
     "270.0904": "Db",
     "273.1138": "Sg",
     "274.1216": "Bh",
     "272.106": "Hs",
     "278.1528": "Mt",
     "283.1918": "Ds",
     "282.184": "Rg",
     "287.223": "Cn",
     "286.2152": "Nh",
     "291.1964": "Fl",
     "290.1888": "Mc",
     "295.2268": "Lv",
     "293.2116": "Ts",
     "299.2572": "Og"}


def monitor_job(dir_mol, sleeptime=60):
    """
    Monitor the job status and judge whether the job is finished.
    :param dir_mol: the path of the job
    :param sleeptime: the sleeping time (unit: second) interval between two judgements
    :return: None
    """
    job = dir_mol.split(os.sep)[-2]
    job2 = dir_mol.split(os.sep)[-3]
    if job == 'DFT':
        subscript = 'gaussian.sh'
    elif job == 'MD':
        subscript = 'lammps.sh'
    else:
        if job2 == 'DFT':
            subscript = 'gaussian.sh'
        elif job2 == 'MD':
            subscript = 'lammps.sh'
        else:
            print('*******************************')
            print('NEITHER DFT NOR MD job!')
            print('*******************************')
            raise Exception
    sh_command = 'sbatch %s > %s' % (subscript, 'jobID.txt')
    custom_shell.execute(sh_command, dir_mol)
    #os.chdir(dir_mol)
    # dir_subscript = os.path.join(dir_mol, subscript)
    # dir_jobID = os.path.join(dir_mol, 'jobID.txt')
    #os.system('chmod +x %s' % subscript)
    #os.system('sbatch %s > %s' % (subscript, 'jobID.txt'))
    dir_jobID = os.path.join(dir_mol, 'jobID.txt')
    jobIDfile = 0
    while jobIDfile == 0:
        if os.path.isfile(dir_jobID):
            try:
                with open(dir_jobID, 'r') as f:
                    jobID = int(f.read().split()[-1])
                # print('jobID = %s' % jobID)
                jobIDfile = 1
            except Exception as e:
                print('*******************************')
                print('JobID NOT found in jobID file!')
                print('*******************************')
                raise e
        else:
            # print('The jobID file is not generated!')
            time.sleep(1)
    # print('Job submitted!\n')
    # print('In progress...\n')
    job = 0
    while job == 0:
        mon = os.system('squeue | grep %s > /dev/null' % jobID)
        if mon != 0:
            # print('Next job!')
            job = 1
            os.remove(dir_jobID)
            os.remove(os.path.join(dir_mol, subscript))
        else:
            # print('The job is not finished!')
            time.sleep(sleeptime)  # second
            job = 0
    # print('Job finished!\n')


def job_status(dir_jobstatus, groupnum):
    """
    Monitor the job status.
    :param groupnum: the group number
    :param dirjob_status: the path where to write the job_status file
    :return:
    """
    dir_js = os.path.join(dir_jobstatus, 'job_status-%s.csv' % groupnum)
    column = ['name/id',
              'moldft1',
              'molmd1',
              'moldft2',
              'molmd_RESP2',
              'moldft_RESP2',
              'moldft_Li',
              'moldft3_0',
              'moldft3_1',
              'moldft3_2',
              'molmd2',
              'moldft4']
    job_status_df = pd.DataFrame(columns=column)
    job_status_df.to_csv(dir_js, index=False, header=True)
    return dir_js


__all__ = [
    'LammpsTemplate',
    'number2element'
]
