import json
import os
import re
import pandas as pd


def smile_analysis(dir_csv):
    data = pd.read_csv(dir_csv)
    smis = data['SMILES']
    #print(smis)
    for smi in smis:
        C = smi.count('C')
        O = smi.count('O')
        F = smi.count('F')
        #print(C+O+F)
        CO = C + O
        if CO > 9 and F == 0:
            print(smi)


def screening_mols(dir_csv):
    # Molecule screening
    data = pd.read_csv(dir_csv)
    dict = {}
    num = []
    print(f'Initial number of molecules: {len(data)}\n')
    num.append(len(data))
    selected = list(data['EP ID'] + ' ' + data['SMILES'])
    selected.sort()
    dict['0-total'] = selected

    #data = data.loc[(data['CAS'] != 'None')]
    num.append(len(data))
    print(f'0-CAS:')
    print(f'Rest molecules: {len(data):5}; Excluding: {1 - len(data) / num[0]:.2%}\n')
    selected = list(data['EP ID'] + ' ' + data['SMILES'])
    selected.sort()
    dict[f'0-CAS'] = selected

    criteria = {
        'Formation energy in vaccum (eV/atom)': '<-0.2',
        'Binding energy (eV)': '>-2.5;<0',
        'Viscosity of solvents (mPa s)': '<2',
        'Viscosity of electrolytes (mPa s)': '<5',
        'Transference number of Li': '',
        'Diffusivity of Li (m^2/s)': '>10e-11',
        'LUMO_sol (eV)': '>4',
        'LUMO_cls (eV)': '>0',
        'HOMO_sol (eV)': '<-5',
        'HOMO_cls (eV)': '<-8'
    }

    count = 0
    for key, value in criteria.items():
        count += 1
        if '<' in value and '>' not in value:
            upbound = eval(re.split('[<>]', value)[-1])
            data = data.loc[(data[key] < upbound)]
        elif '>' in value and '<' not in value:
            lowbound = eval(re.split('[<>]', value)[-1])
            data = data.loc[(data[key] > lowbound)]
        elif '<' in value and '>' in value:
            uplow = value.split(';')
            if '<' in uplow[0]:
                upbound = eval(re.split('[<>]', uplow[0])[-1])
                lowbound = eval(re.split('[<>]', uplow[1])[-1])
            else:
                upbound = eval(re.split('[<>]', uplow[1])[-1])
                lowbound = eval(re.split('[<>]', uplow[0])[-1])
            data = data.loc[(data[key] < upbound) & (data[key] > lowbound)]
        num.append(len(data))
        #print(f'{count}-{key} {value}:')
        #print(f'Rest molecules: {len(data):5}; Excluding: {1 - len(data) / num[0]:.2%}\n')
        print(f'{1 - len(data) / num[0]:.2%} ({len(data)})')
        selected = list(data['EP ID'] + ' ' + data['SMILES'])
        selected.sort()
        dict[f'{count}-{key} {value}'] = selected

    epid_smiles = {}
    for s in selected:
        epid_smiles[s.split()[0]] = s.split()[1]
    with open(os.path.join(os.path.dirname(dir_csv),
                           f'{os.path.basename(dir_csv).split(".")[0]}-selected.json'), 'w') as f:
        f.write(json.dumps(epid_smiles, indent=1))
    data.sort_values(by=['EP ID']).to_csv(
        os.path.join(os.path.dirname(dir_csv), f'{os.path.basename(dir_csv).split(".")[0]}-selected.csv'), index=False)


def main():
    screening_mols(r'E:\Electrolyte_Project\analysis\8-20230228-data\ep_data-20230909-wSMILES.csv')


if __name__ == '__main__':
    main()

