import os
import glob
import assist
import pandas as pd


def tasks():
    return {
        1: 'moldft1',
        2: 'molmd1',
        3: 'moldft2',
        4: 'molmd_RESP2',
        5: 'moldft_RESP2',
        6: 'moldft_Li',
        7: 'moldft3',
        8: 'molmd2'
    }


def check_stage(mol_name):
    df_js = pd.read_csv(dir_js)
    df_js.at[len(df_js), 'name/id'] = mol_name
    df_js.to_csv(dir_js, index=False, header=True)

    dir_mol = os.path.join(dir_db, mol_name)
    dir_moldft = os.path.join(dir_mol, 'DFT')
    dir_molmd = os.path.join(dir_mol, 'MD')
    taskid_name = tasks()
    for taskid in range(1, 9):
        taskname = taskid_name[taskid]
        if taskid != 7:
            if 'dft' in taskname:
                dir_run = os.path.join(dir_moldft, taskname)
            else:
                dir_run = os.path.join(dir_molmd, taskname)
            if not os.path.exists(dir_run) or len(os.listdir(dir_run)) == 0:
                break
            molinfo = glob.glob(os.path.join(dir_run, '*.molinfo'))
            sub_script = glob.glob(os.path.join(dir_run, '*.sh'))
            if len(molinfo) > 0:
                df_js.at[len(df_js) - 1, taskname] = 0
                df_js.to_csv(dir_js, index=False, header=True)
            else:
                if len(sub_script) > 0:
                    break
                else:
                    df_js.at[len(df_js) - 1, taskname] = 1
                    df_js.to_csv(dir_js, index=False, header=True)
                    break
        else:
            success = 0
            for sampid in range(3):
                dir_run = os.path.join(dir_moldft, taskname + f'_{sampid}')
                if not os.path.exists(dir_run) or len(os.listdir(dir_run)) == 0:
                    break
                molinfo = glob.glob(os.path.join(dir_run, '*.molinfo'))
                fchk = glob.glob(os.path.join(dir_run, '*.fchk*'))
                sub_script = glob.glob(os.path.join(dir_run, '*.sh'))
                if len(molinfo) > 0 and len(fchk) > 0:
                    df_js.at[len(df_js) - 1, taskname + f'_{str(sampid)}'] = 0
                    df_js.to_csv(dir_js, index=False, header=True)
                    success += 1
                else:
                    if len(sub_script) > 0:
                        break
                    else:
                        df_js.at[len(df_js) - 1, taskname + f'_{str(sampid)}'] = 1
                        df_js.to_csv(dir_js, index=False, header=True)
            if success == 0:
                break


def main():
    global dir_db
    global dir_js
    dir_db = '/home/yaonan/Electrolyte_Project/moldb'
    dir_js = assist.job_status(os.getcwd(), 'all')
    count = 0
    mols = [mol for mol in os.listdir(dir_db) if mol.startswith('ep-')]
    mols.sort()
    for mol in mols:
        mol_name = mol.split('.')[0]
        check_stage(mol_name)
        count += 1
        if count % 100 == 0:
            print(f'{count}')


if __name__ == '__main__':
    main()

