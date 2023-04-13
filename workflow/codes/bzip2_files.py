#coding:utf-8
import bz2
import multiprocessing as mp
import os
import custom_shell

from typing import List


def bzip2_shell(source_files):
    for source_file in source_files:
        sh_command = 'rm -rf %s && bzip2 %s' % (os.path.basename(source_file) + '.bz2', os.path.basename(source_file))
        custom_shell.execute(sh_command, os.path.dirname(source_file))
    return None


def bzip2(source_file: str) -> str:
    '''Compress a file using bz2.

    Args:
    - source_file: Path to the source file to be compressed.

    Returns:
    - target_file: Path to the target file already compressed.
    '''
    target_file = source_file + '.bz2'
    with open(target_file, 'wb') as tgt:
        tgt.write(bz2.compress(open(source_file, 'rb').read(), compresslevel=9))
    os.remove(source_file)
    return target_file


def bzip2_parallel(source_files: List[str], ncores: int = 4) -> List[str]:
    '''Compress a list of files using bz2 in parallel

    Args:
    - source_files: List of source files to be compressed.
    - ncores: Available cores.

    Returns:
    - target_files: List of target files already compressed.
    '''
    pool = mp.Pool(ncores)
    target_files = pool.map(bzip2, source_files)
    pool.close()
    pool.join()
    return target_files


__all__ = [
    'bzip2',
    'bzip2_parallel'
]
