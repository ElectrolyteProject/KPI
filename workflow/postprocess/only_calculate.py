import os
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
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)-15s [%(filename)s:%(lineno)d] %(levelname)s %(message)s'
)


def calculate_task(mol_name, sampid=0):
    logging.info(f'{mol_name} task{args.start}')
    dir_mol = os.path.join(dir_db, mol_name)
    dir_moldft = os.path.join(dir_mol, 'DFT')
    dir_molmd = os.path.join(dir_mol, 'MD')
    if args.start == 2:
        ##############################   2   ##############################
        dir_moldft1 = os.path.join(dir_moldft, 'moldft1')
        dir_molmd1 = os.path.join(dir_molmd, 'molmd1')
        kMoleculeCount = 200
        mol_info = local_fs.load_json_file(os.path.join(dir_moldft1, mol_name + '.molinfo'))
        files = md.md_judge(dir_molmd1)
        in_path = [item for item in files if item.endswith('in.lammps')][0]
        log_path = [item for item in files if item.endswith('log.lammps')][0]
        dipole_path = [item for item in files if item.endswith('_dipole.out')][0]
        pressure_path = [item for item in files if item.endswith('_pressure.out')][0]
        dc, T, V = md.caldc(in_path, log_path, dipole_path,
                            molnumtot=kMoleculeCount, molnumsol=kMoleculeCount, freq=dcfreq, cutoff=5)
        vis, xls = md.calvis(dir_molmd1, mol_name, pressure_path, T, V,
                             freqs=[visfreq], interval=visitv, split_parts=sp)
        xyz_path = os.path.join(dir_molmd1, mol_name + '.xyz')
        lmp_path = os.path.join(dir_molmd1, mol_name + '.lmp')
        files.extend([xls, xyz_path, lmp_path])
        mol_info.update({
            'volume(Ang^3)': V,
            'dielectric_constant': dc,
            'viscosity(mPas)': vis
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

        local_fs.write_json_file(os.path.join(dir_molmd1, mol_name + '.molinfo'), mol_info)
        for item in files:
            if 'op.data' in item:
                files.remove(item)
        compressed_files = bzip2_files.bzip2_shell(files)

    if args.start == 3:
        ##############################   3   ##############################
        dir_molmd1 = os.path.join(dir_molmd, 'molmd1')
        dir_moldft2 = os.path.join(dir_moldft, 'moldft2')
        metadata = local_fs.load_json_file(os.path.join(dir_molmd1, mol_name + '.molinfo'))
        logs = glob.glob(os.path.join(dir_moldft2, '*.log'))
        chks = glob.glob(os.path.join(dir_moldft2, '*.chk'))
        extracted = dft.analyze_log(logs[0])
        metadata.update(extracted)
        dft.formchk(chks[0])
        # Update time
        #metadata.update({'time(hour)': round((time.time() - start_time) / 3600, 4)})

        local_fs.write_json_file(os.path.join(dir_moldft2, mol_name + '.molinfo'), metadata)
        compressed_files = bzip2_files.bzip2_parallel([logs[0]])

    elif args.start == 4:
        ##############################   4   ##############################
        dir_moldft2 = os.path.join(dir_moldft, 'moldft2')
        dir_molmd_RESP2 = os.path.join(dir_molmd, 'molmd_RESP2')
        kMoleculeCount = 200
        mol_info = local_fs.load_json_file(os.path.join(dir_moldft2, mol_name + '.molinfo'))
        files = md.md_judge(dir_molmd_RESP2)
        in_path = [item for item in files if item.endswith('in.lammps')][0]
        log_path = [item for item in files if item.endswith('log.lammps')][0]
        dipole_path = [item for item in files if item.endswith('_dipole.out')][0]
        pressure_path = [item for item in files if item.endswith('_pressure.out')][0]
        dc, T, V = md.caldc(in_path, log_path, dipole_path,
                            molnumtot=kMoleculeCount, molnumsol=kMoleculeCount, freq=dcfreq, cutoff=5)
        vis, xls = md.calvis(dir_molmd_RESP2, mol_name, pressure_path, T, V,
                             freqs=[visfreq], interval=visitv, split_parts=sp)
        lmp_path = os.path.join(dir_molmd_RESP2, mol_name + '.lmp')
        files.extend([xls, lmp_path])
        mol_info.update({
            'volume_RESP2(Ang^3)': V,
            'dielectric_constant_RESP2': dc,
            'viscosity_RESP2(mPas)': vis
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

        local_fs.write_json_file(os.path.join(dir_molmd_RESP2, mol_name + '.molinfo'), mol_info)
        compressed_files = bzip2_files.bzip2_shell(files)

    elif args.start == 7:
        dir_moldft_Li = os.path.join(dir_moldft, f'moldft_Li')
        dir_moldft3 = os.path.join(dir_moldft, f'moldft3_{sampid}')
        if sampid == 0:
            metadata = local_fs.load_json_file(os.path.join(dir_moldft_Li, mol_name + '.molinfo'))
        else:
            metadata = local_fs.load_json_file(
                os.path.join(dir_moldft, 'moldft3_%s' % (sampid - 1), mol_name + '.molinfo'))
        files = []
        for file in os.listdir(dir_moldft3):
            if file.endswith('.log'):
                files.append(os.path.join(dir_moldft3, file))
        extracted = dft.analyze_log(files[0])
        #dft.formchk(files[1])
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
        coordinate_samp = metadata['coordinate_samp']
        if len(coordinate_samp) == 0:
            metadata.pop('coordinate_samp')
        else:
            coordinate_samp.pop(0)
            metadata.update({'coordinate_samp': coordinate_samp})
        """
        # Update time
        if metadata.get('time(hour)'):
            metadata['time(hour)'] += (time.time() - start_time) / 3600
            round(metadata['time(hour)'], 4)
        else:
            metadata.update({'time(hour)': round((time.time() - start_time) / 3600, 4)})
        """
        # Move generated files to the output folder
        local_fs.write_json_file(os.path.join(dir_moldft3, mol_name + '.molinfo'), metadata)
        compressed_files = bzip2_files.bzip2_parallel([files[0]])
    elif args.start == 8:
        ##############################   8   ##############################
        dir_moldft3 = os.path.join(dir_moldft, 'moldft3_2')
        dir_molmd2 = os.path.join(dir_molmd, 'molmd2')
        kMoleculeCount = 200
        kSaltCount = 20
        kTotalCount = kMoleculeCount + kSaltCount
        mol_names = str('%s %s' % (mol_name, config['salt'])).split()
        mol_nums = '%s %s' % (kMoleculeCount, kSaltCount)
        mol_info = local_fs.load_json_file(os.path.join(dir_moldft3, mol_name + '.molinfo'))
        files = md.md_judge(dir_molmd2)
        in_path = [item for item in files if item.endswith('in.lammps')][0]
        log_path = [item for item in files if item.endswith('log.lammps')][0]
        dipole_path = [item for item in files if item.endswith('_dipole.out')][0]
        pressure_path = [item for item in files if item.endswith('_pressure.out')][0]
        model = os.path.basename(pressure_path).split('_')[0]
        dc, T, V = md.caldc(in_path, log_path, dipole_path,
                            molnumtot=kTotalCount, molnumsol=kMoleculeCount, freq=dcfreq, cutoff=5)
        vis, xls = md.calvis(dir_molmd2, mol_name, pressure_path, T, V,
                             freqs=[visfreq], interval=visitv, split_parts=sp)
        diffs, msdxls, cond, betas = md.calmsd(dir_molmd2, model, mol_names)
        rdfs = md.calrdfCN(dir_molmd2, model, mol_names, mol_nums, V)
        lmp_path = os.path.join(dir_molmd2, mol_name + '.lmp')
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
        compressed_files = bzip2_files.bzip2_shell(files)

        fchks = glob.glob(os.path.join(dir_moldft, '*', '*.fchk'))
        bzip2_files.bzip2_shell(fchks)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-config', help="The name of the config file", type=str, default=None)
    parser.add_argument('-start', help="Start from task N", type=int, default=1)
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
    allmols = config['molecules']
    allmols.sort()
    for mol in allmols:
        try:
            calculate_task(mol)
        except Exception as e:
            logging.info(e)
            logging.info(traceback.format_exc())
            continue
    logging.info('==================== END %s s ====================' % float(round((time.time() - zero_time_tot), 4)))

