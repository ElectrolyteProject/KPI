import os
import dft
import local_fs
import assist
import time
import argparse
import logging
import traceback
import bzip2_files
import re
import pandas as pd
import multiprocessing as mp


def gjf_modifier(dir_run, chg_spmul, keywords, solva_info, ncores: int = 8):
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
    gjf_path = False
    for file in os.listdir(dir_run):
        if file.endswith('.gjf'):
            gjf_name = file.split('.')[0]
            gjf_path = os.path.join(dir_run, gjf_name + '.gjf')
            ckpt_path = os.path.join(dir_run, gjf_name + '.chk')
    if gjf_path:
        with open(gjf_path, 'r') as f1:
            f1.seek(0, 0)
            infos = f1.readlines()
        infos[1] = keywords + '\n'
        infos[3] = gjf_name + '\n'
        infos[5] = str(chg_spmul) + '\n'  # charge and spin multiplicity
        infos.insert(1, f'%nproc={str(ncores)}\n')
        for l, line in enumerate(infos):
            if '%chk=' in line:
                infos[l] = '%chk=' + ckpt_path + '\n'
        if 'scrf' in keywords and 'water' not in keywords:  # Please mind the letter case
            infos.insert(-1, '\nsolventname=%s\n' % solva_info['name'])
            infos.insert(-1, 'RSOLV=%s\n' % solva_info['radius'])
            infos.insert(-1, 'EPS=%s\n' % dc)
            infos.insert(-1, 'EpsInf=%s\n' % str(1.758))
        if 'nbo' in keywords:
            infos.insert(-1, '\n$nbo bndidx archive NRT BOAO $end\n')
        with open(gjf_path, 'w') as f2:
            f2.writelines(infos)
        dir_subscript = os.path.join(dir_run, 'gaussian.sh')
        with open(dir_subscript, 'w+') as f2:
            f2.write('#!/bin/bash\n'
                     '#SBATCH -J g16\n'
                     '#SBATCH -N 1\n'
                     '#SBATCH -n %s\n' % ncores +
                     '#SBATCH -o stdout.%j\n'
                     '#SBATCH -e stderr.%j\n'
                     '\n'
                     'module load Gaussian\n' +
                     'g16 %s.gjf\n' % gjf_name)
    else:
        raise FileExistsError
    return [os.path.join(dir_run, gjf_name + '.log'), os.path.join(dir_run, gjf_name + '.chk')]


def only_sub_dft(mol_name, sampid=0):
    #job_status_df = pd.read_csv(dir_js)
    #job_status_df.at[len(job_status_df), 'name/id'] = mol_name
    #index = job_status_df[job_status_df['name/id'] == mol_name].index.tolist()[0]
    #job_status_df.to_csv(dir_js, index=False, header=True)
    dir_mol = os.path.join(dir_db, mol_name)
    dir_moldft = os.path.join(dir_mol, 'DFT')
    if args.start == 1:
        try:
            start_time = time.time()
            logging.info('--- Task1 on %s: DFT optimization in vacuum ---', mol_name)
            dir_moldft1 = os.path.join(dir_moldft, 'moldft1')
            if not os.path.exists(dir_moldft1):
                os.mkdir(dir_moldft1)
            metadata = local_fs.load_json_file(os.path.join(dir_ogdb, mol_name + '.mol'))
            mol_name = metadata['name']
            # Call Gaussian to optimize molecule geometry
            keywords = dft.keywords(1)
            #keywords = keywords.replace('pop=(nbo,savenbo)', 'pop=nbo')
            #files = dft.job(dir_moldft1, mol_info, '0 1', keywords, ncores=ncores_dft)
            files = gjf_modifier(dir_moldft1, '0 1', keywords, solva_info=metadata, ncores=ncores_dft)
            #gjf_name = mol_info['name'] + '_vum'
            #files = [os.path.join(dir_moldft1, gjf_name + '.log'),
            #         os.path.join(dir_moldft1, gjf_name + '.chk')]
            assist.monitor_job(dir_moldft1, wait_dft)
            extracted = dft.analyze_log(files[0])
            metadata.update(extracted)
            dft.formchk(files[1])
            # Update time
            metadata.update({'time(hour)': round((time.time() - start_time) / 3600, 4)})

            local_fs.write_json_file(os.path.join(dir_moldft1, mol_name + '.molinfo'), metadata)
            compressed_files = bzip2_files.bzip2_shell([files[0]])
            logging.info('--- Task1 on %s: %.4f s ---', mol_name, time.time() - start_time)

            #job_status_df.at[len(job_status_df) - 1, 'moldft1'] = 0
            #job_status_df.to_csv(dir_js, index=False, header=True)
        except Exception as e:
            #job_status_df.at[len(job_status_df) - 1, 'moldft1'] = 1
            #job_status_df.to_csv(dir_js, index=False, header=True)
            raise e

    if args.start == 3:
        try:
            start_time = time.time()
            logging.info('--- Task3 on %s: DFT optimization in solvent environment ---', mol_name)
            dir_moldft2 = os.path.join(dir_moldft, 'moldft2')
            if not os.path.exists(dir_moldft2):
                os.mkdir(dir_moldft2)
            dir_molmd1 = os.path.join(dir_mol, 'MD', 'molmd1')
            metadata = local_fs.load_json_file(os.path.join(dir_molmd1, mol_name + '.molinfo'))
            mol_name = metadata['name']
            # Call Gaussian to optimize molecule geometry
            keywords = dft.keywords(2)
            files = gjf_modifier(dir_moldft2, '0 1', keywords, solva_info=metadata, ncores=ncores_dft)
            assist.monitor_job(dir_moldft2, wait_dft)
            extracted = dft.analyze_log(files[0])
            metadata.update(extracted)
            dft.formchk(files[1])
            # Update time
            metadata.update({'time(hour)': round((time.time() - start_time) / 3600, 4)})

            local_fs.write_json_file(os.path.join(dir_moldft2, mol_name + '.molinfo'), metadata)
            compressed_files = bzip2_files.bzip2_shell([files[0]])
            logging.info('--- Task3 on %s: %.4f s ---', mol_name, time.time() - start_time)

            #job_status_df.at[len(job_status_df) - 1, 'moldft2'] = 0
            #job_status_df.to_csv(dir_js, index=False, header=True)
        except Exception as e:
            #job_status_df.at[len(job_status_df) - 1, 'moldft2'] = 1
            #job_status_df.to_csv(dir_js, index=False, header=True)
            raise e

    if args.start == 5:
        try:
            start_time = time.time()
            logging.info('--- Task5 on %s: DFT optimization in solvent environment (RESP2) ---', mol_name)
            dir_moldft_RESP2 = os.path.join(dir_moldft, 'moldft_RESP2')
            if not os.path.exists(dir_moldft_RESP2):
                os.mkdir(dir_moldft_RESP2)
            dir_molmd_RESP2 = os.path.join(dir_mol, 'MD', 'molmd_RESP2')
            metadata = local_fs.load_json_file(os.path.join(dir_molmd_RESP2, mol_name + '.molinfo'))
            mol_name = metadata['name']
            # Call Gaussian to optimize molecule geometry
            keywords = dft.keywords(2)
            files = gjf_modifier(dir_moldft_RESP2, '0 1', keywords, solva_info=metadata, ncores=ncores_dft)
            assist.monitor_job(dir_moldft_RESP2, wait_dft)
            extracted = dft.analyze_log(files[0])
            for key, value in extracted.items():
                if 'coordinate' in key:
                    mol_info_samp = {'name': mol_name, 'coordinate': value}
                    mol_coor_collection = dft.li_coor_sampling(mol_info_samp)
                    metadata.update({'coordinate_samp': mol_coor_collection})
            metadata.update(extracted)
            dft.formchk(files[1])
            # Update time
            if metadata.get('time(hour)'):
                metadata['time(hour)'] += (time.time() - start_time) / 3600
                round(metadata['time(hour)'], 4)
            else:
                metadata.update({'time(hour)': round((time.time() - start_time) / 3600, 4)})

            local_fs.write_json_file(os.path.join(dir_moldft_RESP2, mol_name + '.molinfo'), metadata)
            compressed_files = bzip2_files.bzip2_shell([files[0]])
            logging.info('--- Task5 on %s: %.4f s ---', mol_name, time.time() - start_time)

            #job_status_df.at[len(job_status_df) - 1, 'moldft_RESP2'] = 0
            #job_status_df.to_csv(dir_js, index=False, header=True)
        except Exception as e:
            #job_status_df.at[len(job_status_df) - 1, 'moldft_RESP2'] = 1
            #job_status_df.to_csv(dir_js, index=False, header=True)
            raise e

    if args.start == 6:
        try:
            start_time = time.time()
            logging.info('--- Task6 on %s: DFT optimization of Li+ in solvent environment (RESP2) ---', mol_name)
            dir_moldft_Li = os.path.join(dir_moldft, 'moldft_Li')
            if not os.path.exists(dir_moldft_Li):
                os.mkdir(dir_moldft_Li)
            dir_moldft_RESP2 = os.path.join(dir_moldft, 'moldft_RESP2')
            metadata = local_fs.load_json_file(os.path.join(dir_moldft_RESP2, mol_name + '.molinfo'))
            mol_name = metadata['name']
            mol_info = {'name': '%s-Li+' % mol_name, 'coordinate': ["Li   0.0   0.0   0.0"]}
            # Call Gaussian to optimize molecule geometry
            keywords = dft.keywords(2)
            files = dft.job(dir_moldft_Li, mol_info, '1 1', keywords, solva_info=metadata, ncores=ncores_dft)
            assist.monitor_job(dir_moldft_Li, wait_dft)
            extracted = dft.analyze_log_li(files[0])
            metadata.update(extracted)
            dft.formchk(files[1])

            # Update time
            if metadata.get('time(hour)'):
                metadata['time(hour)'] += (time.time() - start_time) / 3600
                round(metadata['time(hour)'], 4)
            else:
                metadata.update({'time(hour)': round((time.time() - start_time) / 3600, 4)})

            # Move generated files to the output folder
            local_fs.write_json_file(os.path.join(dir_moldft_Li, mol_name + '.molinfo'), metadata)
            compressed_files = bzip2_files.bzip2_shell([files[0]])
            logging.info('--- Task6 on %s: %.4f s ---', mol_name, time.time() - start_time)

            #job_status_df.at[len(job_status_df) - 1, 'moldft_Li'] = 0
            #job_status_df.to_csv(dir_js, index=False, header=True)
        except Exception as e:
            #job_status_df.at[len(job_status_df) - 1, 'moldft_Li'] = 1
            #job_status_df.to_csv(dir_js, index=False, header=True)
            raise e

    if args.start == 7:
        try:
            start_time = time.time()
            logging.info('--- Task7_%s on %s: DFT optimization of multiple ion-solvent clusters (RESP2) ---', sampid,
                        mol_name)
            dir_moldft3 = os.path.join(dir_moldft, 'moldft3_%s' % sampid)
            try:
                metadata = local_fs.load_json_file(
                    os.path.join(dir_moldft3, mol_name + '.molinfo'))
            except:
                metadata = local_fs.load_json_file(
                    os.path.join(dir_moldft, 'moldft3_%s-error' % sampid, mol_name + '.molinfo'))
            keywords = dft.keywords(2)
            keywords = keywords.replace('opt=calcfc', '')
            files = gjf_modifier(dir_moldft3, '1 1', keywords, solva_info=metadata, ncores=ncores_dft)
            assist.monitor_job(dir_moldft3, wait_dft)
            extracted = dft.analyze_log(files[0])
            dft.formchk(files[1])
            clskeys = [
                'coordinate_cls',
                'HFenergy_cls',
                'thermochemistry_cls',
                'HOMOLUMO_cls',
                'NBOcharge_cls',
                'MKcharge_spinpop_cls',
                'img_cls',
                'dipole_moment_cls(Debye)'
            ]
            for clskey in clskeys:
                if metadata.get(clskey):
                    metadata[clskey].append(extracted[clskey])
                else:
                    listing = {
                        clskey: [extracted[clskey]]
                    }
                    metadata.update(listing)
            try:
                coordinate_samp = metadata['coordinate_samp']
                coordinate_samp.pop(0)
                if len(coordinate_samp) == 0:
                    metadata.pop('coordinate_samp')
                else:
                    metadata.update({'coordinate_samp': coordinate_samp})
            except:
                pass
            # Update time
            if metadata.get('time(hour)'):
                metadata['time(hour)'] += (time.time() - start_time) / 3600
                round(metadata['time(hour)'], 4)
            else:
                metadata.update({'time(hour)': round((time.time() - start_time) / 3600, 4)})
            local_fs.write_json_file(os.path.join(dir_moldft3, mol_name + '.molinfo'), metadata)
            compressed_files = bzip2_files.bzip2_shell([files[0]])
            logging.info('--- Task7_%s on %s: %.4f s ---', sampid, mol_name, time.time() - start_time)

            #job_status_df.at[index, 'moldft3_%s' % sampid] = 0
            #job_status_df.to_csv(dir_js, index=False, header=True)
        except Exception as e:
            #job_status_df.at[index, 'moldft3_%s' % sampid] = 1
            #job_status_df.to_csv(dir_js, index=False, header=True)
            logging.info(e)
            logging.info(traceback.format_exc())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-config', help="The name of the config file", type=str, default=None)
    parser.add_argument('-start', help="The task id", type=int, default=None)
    args = parser.parse_args()
    zero_time_tot = time.time()

    config = local_fs.load_json_file(args.config)
    dir_ogdb = config['ogdb']
    dir_db = config['moldb']

    # formal params
    ncores_dft = 64#16
    ncores_md = 64#24
    dcfreq = 18000
    visfreq = 50000
    visitv = 100000
    sp = 200
    wait_dft = 60
    wait_md = 600
    """
    # test params
    ncores_dft = 64
    ncores_md = 64
    dcfreq = 2
    visfreq = 50
    visitv = 100
    sp = 200
    wait_dft = 60
    wait_md = 60
    """
    dir_jobstatus = os.path.join(os.path.dirname(dir_db), 'job_status')
    if not os.path.exists(dir_jobstatus):
        os.mkdir(dir_jobstatus)
    dir_js = assist.job_status(dir_jobstatus, re.split('[-.]', args.config)[1])

    allmols = config['molecules']
    allmols.sort()

    pool = mp.Pool(10)
    pool.map(only_sub_dft, allmols)
    pool.close()
    pool.join()
    """
    for mol in allmols:
        try:
            only_sub_dft(mol)
        except Exception as e:
            logging.info(e)
            logging.info(traceback.format_exc())
            continue
    """
    logging.info('==================== END %s s ====================' % float(round((time.time() - zero_time_tot), 4)))

