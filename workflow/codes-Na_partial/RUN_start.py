import os
import re
import pandas as pd
import shutil
import argparse
import time
import dft
import md
import assist
import local_fs
import bzip2_files
import traceback
import logging
import glob
import global_config as gc
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)-15s [%(filename)s:%(lineno)d] %(levelname)s %(message)s'
)


def process_molinfo(data):
    keys_tbr = ['HFenergy_Li+_RESP2',
                'thermochemistry_Li+_RESP2',
                'coordinate_samp',
                'coordinate_cls',
                'HFenergy_cls',
                'thermochemistry_cls',
                'HOMOLUMO_cls',
                'NBOcharge_cls',
                'MKcharge_spinpop_cls',
                'coordinate_cls_vum',
                'HFenergy_cls_vum',
                'thermochemistry_cls_vum',
                'HOMOLUMO_cls_vum',
                'NBOcharge_cls_vum',
                'MKcharge_spinpop_cls_vum',
                'composition_ely',
                'volume_ely(Ang^3)',
                'dielectric_constant_ely',
                'viscosity_ely(mPas)',
                'diffusivity(m^2/s)',
                'conductivity(S/cm)',
                'conductivity-fitting_betas',
                'rdfCN']
    return {key: val for key, val in data.items() if key not in keys_tbr}


def workflow(mol_name):
    zero_time = time.time()
    logging.info(f'=== {mol_name} ===')
    dir_mol = os.path.join(dir_db, mol_name)
    dir_moldft = os.path.join(dir_mol, 'DFT')
    dir_molmd = os.path.join(dir_mol, 'MD')
    dirs = [dir_mol, dir_moldft, dir_molmd]
    for dir in dirs:
        if not os.path.exists(dir):
            os.mkdir(dir)
    df_js = pd.read_csv(dir_js)
    df_js.at[len(df_js), 'name/id'] = mol_name
    df_js.to_csv(dir_js, index=False, header=True)

    ##############################   6   ##############################
    try:
        start_time = time.time()
        logging.info('--- Task6 on %s: DFT optimization of Na+ in solvent environment (RESP2) ---', mol_name)
        dir_moldft_Na = os.path.join(dir_moldft, 'moldft_Na')
        if not os.path.exists(dir_moldft_Na):
            os.mkdir(dir_moldft_Na)
        metadata = local_fs.load_json_file(os.path.join(MOLINFODIR, mol_name + '.molinfo'))
        metadata = process_molinfo(metadata)
        for key in list(metadata.keys()):
            if key == 'coordinate_RESP2':
                mol_info_samp = {'name': mol_name, 'coordinate': metadata[key]}
                mol_coor_collection = dft.na_coor_sampling(mol_info_samp)
                metadata.update({'coordinate_samp': mol_coor_collection})
        mol_name = metadata['name']
        mol_info = {'name': '%s-Na+' % mol_name, 'coordinate': ["Na   0.0   0.0   0.0"]}
        molinfo = glob.glob(os.path.join(dir_moldft_Na, '*.molinfo'))
        if args.start <= 6 and len(molinfo) == 0:
            # Call Gaussian to optimize molecule geometry
            keywords = dft.keywords(2)
            files = dft.job(dir_moldft_Na, mol_info, '1 1', keywords, solva_info=metadata, sc=config['supercomputer'], ncores=ncores_dft)
            assist.monitor_job(dir_moldft_Na, wait_dft)
            extracted = dft.analyze_log_na(files[0])
            metadata.update(extracted)
            dft.formchk(files[1])

            # Update time
            if metadata.get('time(hour)'):
                metadata['time(hour)'] += (time.time() - start_time) / 3600
                round(metadata['time(hour)'], 4)
            else:
                metadata.update({'time(hour)': round((time.time() - start_time) / 3600, 4)})

            # Move generated files to the output folder
            local_fs.write_json_file(os.path.join(dir_moldft_Na, mol_name + '.molinfo'), metadata)
            compressed_files = bzip2_files.bzip2_parallel([files[0]])
            logging.info('--- Task6 on %s: %.4f s ---', mol_name, time.time() - start_time)

        df_js.at[len(df_js) - 1, 'moldft_Na'] = 0
        df_js.to_csv(dir_js, index=False, header=True)
    except Exception as e:
        df_js.at[len(df_js) - 1, 'moldft_Na'] = 1
        df_js.to_csv(dir_js, index=False, header=True)
        raise e

    ##############################   7   ##############################
    success = 0
    for sampid in range(3):
        try:
            dir_moldft3 = os.path.join(dir_moldft, 'moldft3_%s' % sampid)
            if not os.path.exists(dir_moldft3):
                os.mkdir(dir_moldft3)
            if sampid == 0:
                metadata = local_fs.load_json_file(os.path.join(dir_moldft_Na, mol_name + '.molinfo'))
            else:
                metadata = local_fs.load_json_file(os.path.join(dir_moldft, 'moldft3_%s' % (sampid - 1), mol_name + '.molinfo'))
            if metadata.get('coordinate_samp'):
                start_time = time.time()
                logging.info('--- Task7_%s on %s: DFT optimization of multiple ion-solvent clusters (RESP2) ---', sampid, mol_name)
                # metadata = local_fs.load_json_file(mol_name)
                mol_name = metadata['name']
                mol_info = {'name': '%s+Na+' % mol_name, 'coordinate': metadata['coordinate_samp'][0]}
                molinfo = glob.glob(os.path.join(dir_moldft3, '*.molinfo'))
                if args.start <= 7 and len(molinfo) == 0:
                    # Call Gaussian to optimize molecule geometry
                    # Use try-except in case gaussian not normally terminated and still remove the calculated samp coordinate
                    try:
                        keywords = dft.keywords(2)
                        #keywords = keywords.replace('pop=(nbo,savenbo)', 'pop=nbo')
                        files = dft.job(dir_moldft3, mol_info, '1 1', keywords, solva_info=metadata, sc=config['supercomputer'], ncores=ncores_dft)
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
                    except Exception as e:
                        coordinate_samp = metadata['coordinate_samp']
                        coordinate_samp.pop(0)
                        if len(coordinate_samp) == 0:
                            metadata.pop('coordinate_samp')
                        else:
                            metadata.update({'coordinate_samp': coordinate_samp})
                        if metadata.get('time(hour)'):
                            metadata['time(hour)'] += (time.time() - start_time) / 3600
                            round(metadata['time(hour)'], 4)
                        else:
                            metadata.update({'time(hour)': round((time.time() - start_time) / 3600, 4)})
                        local_fs.write_json_file(os.path.join(dir_moldft3, mol_name + '.molinfo'), metadata)
                        raise e

                    # Remove the calculated samp coordinate (the first item in metadata['coordinate_samp'] list)
                    coordinate_samp = metadata['coordinate_samp']
                    coordinate_samp.pop(0)
                    if len(coordinate_samp) == 0:
                        metadata.pop('coordinate_samp')
                    else:
                        metadata.update({'coordinate_samp': coordinate_samp})

                    # Update time
                    if metadata.get('time(hour)'):
                        metadata['time(hour)'] += (time.time() - start_time) / 3600
                        round(metadata['time(hour)'], 4)
                    else:
                        metadata.update({'time(hour)': round((time.time() - start_time) / 3600, 4)})

                    # Move generated files to the output folder
                    local_fs.write_json_file(os.path.join(dir_moldft3, mol_name + '.molinfo'), metadata)
                    compressed_files = bzip2_files.bzip2_parallel([files[0]])
                    logging.info('--- Task7_%s on %s: %.4f s ---', sampid, mol_name, time.time() - start_time)

                df_js.at[len(df_js) - 1, 'moldft3_%s' % sampid] = 0
                df_js.to_csv(dir_js, index=False, header=True)
                success += 1
            else:
                moldft3_fchks = glob.glob(os.path.join(dir_moldft, 'moldft3_*', '*.fchk*'))
                if len(moldft3_fchks) == 0:
                    raise KeyError('No ion-solvent clusters to be sampled!')
        except Exception:
            local_fs.write_json_file(os.path.join(dir_moldft3, mol_name + '.molinfo'), metadata)
            df_js.at[len(df_js) - 1, 'moldft3_%s' % sampid] = 1
            df_js.to_csv(dir_js, index=False, header=True)
            continue
    moldft3_fchks = glob.glob(os.path.join(dir_moldft, 'moldft3_*', '*.fchk*'))
    if success == 0 and len(moldft3_fchks) == 0:
        raise RuntimeError('No ion-solvent cluster sampling succeeds!')

    ##############################   8   ##############################
    try:
        start_time = time.time()
        logging.info('--- Task8 on %s: MD simulation of electrolyte ---', mol_name)
        dir_molmd2 = os.path.join(dir_molmd, 'molmd2')
        if not os.path.exists(dir_molmd2):
            os.mkdir(dir_molmd2)
        for sampid in range(2, -1, -1):
            moldft3_fchks = glob.glob(os.path.join(dir_moldft, f'moldft3_{sampid}', '*.fchk*'))
            if len(moldft3_fchks) > 0:
                dir_moldft3 = os.path.join(dir_moldft, f'moldft3_{sampid}')
                break
        mol_info = local_fs.load_json_file(os.path.join(dir_moldft3, mol_name + '.molinfo'))
        mol_name = mol_info['name']
        mol_names = '%s %s' % (mol_name, config['salt'])
        mol_info_mod = {'name': mol_name, 'coordinate': mol_info['coordinate_sol']}
        kMoleculeCount = 200
        kSaltCount = 20
        kTotalCount = kMoleculeCount + kSaltCount
        mol_nums = '%s %s' % (kMoleculeCount, kSaltCount)
        molinfo = glob.glob(os.path.join(dir_molmd2, '*.molinfo'))
        if args.start <= 8 and len(molinfo) == 0:
            # Prepare LAMMPS input files
            box_size = md.boxsize(mol_info, kTotalCount)
            model, mol_names, mol_nums = md.def_model(mol_names, mol_nums)
            assert len(mol_nums) == 2, 'Task 8 consumes 1 molecule type and 1 salt type!'
            xyz_path = md.mkfile_xyz(dir_molmd2, mol_info_mod)  # To pack with other molecules

            lmp_path = os.path.join(MOLINFODIR, mol_name + '.lmp')

            sxyz_path = md.cpfile(os.path.join(dir_molmd2, config['salt'] + '.xyz'))
            slmp_path = md.cpfile(os.path.join(dir_molmd2, config['salt'] + '.lmp'))
            model_path = md.mkfile_modelxyz(dir_molmd2, model, [xyz_path, sxyz_path], mol_nums,
                                            round(box_size * pow(0.7, 1 / 3), 2))
            data_path = md.mkfile_data(dir_molmd2, model, model_path, mol_names, [lmp_path, slmp_path], mol_nums)
            # assert data_path == model + '.data', 'Name of the LAMMPS data file [%s] is not consistent!' % data_path

            # Call LAMMPS to simulate and analyze result
            md.job(dir_molmd2, 'in.lammps', model, 'dm press', sc=config['supercomputer'], ncores=ncores_md)
            assist.monitor_job(dir_molmd2, wait_md)
            files = md.md_judge(dir_molmd2)
            in_path = [item for item in files if item.endswith('in.lammps')][0]
            log_path = [item for item in files if item.endswith('log.lammps')][0]
            dipole_path = [item for item in files if item.endswith('_dipole.out')][0]
            pressure_path = [item for item in files if item.endswith('_pressure.out')][0]
            dc, T, V = md.caldc(in_path, log_path, dipole_path,
                                molnumtot=kTotalCount, molnumsol=kMoleculeCount, freq=dcfreq, cutoff=5)
            vis, xls = md.calvis(dir_molmd2, mol_name, pressure_path, T, V,
                                 freqs=[visfreq], interval=visitv, split_parts=sp)
            diffs, msdxls, cond, betas = md.calmsd(dir_molmd2, model, mol_names)
            rdfs = md.calrdfCN(dir_molmd2, model, mol_names, mol_nums, V)
            files.extend([xls, msdxls, lmp_path])
            mol_info.update({
                'composition_ely': model,
                'volume_ely(Ang^3)': V,
                'dielectric_constant_ely': dc,
                'viscosity_ely(mPas)': vis,
                'diffusivity(m^2/s)': diffs,
                'conductivity(S/cm)': cond,
                'conductivity-fitting_betas': betas,
                'rdfCN': rdfs
            })
            # Update time and core
            if mol_info.get('time(hour)'):
                mol_info['time(hour)'] += md.extract_walltime(log_path)
            else:
                mol_info.update({'time(hour)': md.extract_walltime(log_path)})
            if mol_info.get('ncores_md'):
                pass
            else:
                mol_info.update({'ncores_md': ncores_md})

            local_fs.write_json_file(os.path.join(dir_molmd2, mol_name + '.molinfo'), mol_info)
            compressed_files = bzip2_files.bzip2_parallel(files)
            logging.info('--- Task8 on %s: %.4f s ---', mol_name, time.time() - start_time)
        df_js.at[len(df_js) - 1, 'molmd2'] = 0
        df_js.to_csv(dir_js, index=False, header=True)
    except Exception as e:
        df_js.at[len(df_js) - 1, 'molmd2'] = 1
        df_js.to_csv(dir_js, index=False, header=True)
        raise e

    ##############################   9   ##############################
    try:
        start_time = time.time()
        logging.info(f'--- Task9 on {mol_name}: DFT optimization of ion-solvent cluster in vacuum ---')
        dir_moldft4 = os.path.join(dir_moldft, 'moldft4')
        if not os.path.exists(dir_moldft4):
            os.mkdir(dir_moldft4)
        metadata = local_fs.load_json_file(os.path.join(dir_molmd2, mol_name + '.molinfo'))
        # Select the most stable structure coordinate
        energy_cls = [eval(energy) for energy in metadata['HFenergy_cls']]
        min_index = energy_cls.index(min(energy_cls))
        mol_name = metadata['name']
        mol_info = {'name': '%s+Na+' % mol_name, 'coordinate': metadata['coordinate_cls'][min_index]}
        molinfo = glob.glob(os.path.join(dir_moldft4, '*.molinfo'))
        if args.start <= 9 and len(molinfo) == 0:
            # Call Gaussian to optimize molecule geometry
            keywords = dft.keywords(1)
            files = dft.job(dir_moldft4, mol_info, '1 1', keywords, sc=config['supercomputer'], ncores=ncores_dft)
            assist.monitor_job(dir_moldft4, wait_dft)
            extracted = dft.analyze_log(files[0])
            metadata.update(extracted)
            dft.formchk(files[1])
            # Update time
            if metadata.get('time(hour)'):
                metadata['time(hour)'] += (time.time() - start_time) / 3600
                round(metadata['time(hour)'], 4)
            else:
                metadata.update({'time(hour)': round((time.time() - start_time) / 3600, 4)})
            # Output
            local_fs.write_json_file(os.path.join(dir_moldft4, mol_name + '.molinfo'), metadata)
            compressed_files = bzip2_files.bzip2_parallel([files[0]])
            logging.info(f'--- Task9 on {mol_name}: {time.time() - start_time}.4f s ---')
        df_js.at[len(df_js) - 1, 'moldft4'] = 0
        df_js.to_csv(dir_js, index=False, header=True)
    except Exception as e:
        df_js.at[len(df_js) - 1, 'moldft4'] = 1
        df_js.to_csv(dir_js, index=False, header=True)
        raise e

    fchks = glob.glob(os.path.join(dir_moldft, '*', '*.fchk'))
    bzip2_files.bzip2_parallel(fchks)
    logging.info('=== %s: %s s ===\n' % (mol_name, float(round((time.time() - zero_time), 4))))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-config', help="The name of the config file", type=str, default=None)
    parser.add_argument('-start', help="Start from task N", type=int, default=1)
    args = parser.parse_args()
    zero_time_tot = time.time()

    config = local_fs.load_json_file(args.config)
    dir_ogdb = config['ogdb']
    dir_db = config['moldb']
    if not os.path.exists(dir_db):
        os.mkdir(dir_db)

    MOLINFODIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'MOLINFO')
    if config['supercomputer'] == 'zghpc':
        # formal params
        ncores_dft = 64
        ncores_md = 64
        dcfreq = 18000
        visfreq = 50000
        visitv = 100000
        sp = 200
        wait_dft = 60
        wait_md = 1800
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
    elif config['supercomputer'] == 'ts1000':
        # ts1000
        # formal params
        ncores_dft = 56
        ncores_md = 56
        dcfreq = 18000
        visfreq = 50000
        visitv = 100000
        sp = 200
        wait_dft = 60
        wait_md = 1800
        """
        # test params
        ncores_dft = 56
        ncores_md = 56
        dcfreq = 2
        visfreq = 50
        visitv = 100
        sp = 200
        wait_dft = 60
        wait_md = 60
        """
    dir_jobstatus = os.path.join(dir_ogdb.split('ogdb')[0], 'job_status')
    if not os.path.exists(dir_jobstatus):
        os.mkdir(dir_jobstatus)
    dir_js = assist.job_status(dir_jobstatus, re.split('[-.]', args.config)[-2])

    allmols = config['molecules']
    allmols.sort()
    for mol in allmols:
        try:
            workflow(mol)
        except Exception as e:
            logging.info(e)
            logging.info(traceback.format_exc())
            continue
    logging.info('==================== END %s s ====================' % float(round((time.time() - zero_time_tot), 4)))
