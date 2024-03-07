#coding:utf-8
import json
import os
import shutil
import subprocess as sp

from typing import Dict, List, Union


def load_json_file(path: str) -> Union[Dict, List]:
    '''Read a text file and parse as JSON.

    Args:
    - path: Path to the Json file.

    Returns:
    - data: A Python dict of list of parsed data.
    '''
    with open(path, 'r') as fin:
        return json.load(fin)


def load_yaml_file(path: str) -> Union[Dict, List]:
    '''Read a text file and parse as YAML.

    Args:
    - path: Path to the Json file.

    Returns:
    - data: A Python dict of list of parsed data.
    '''
    import yaml
    with open(path, 'r') as fin:
        return yaml.safe_load(fin)


def make_empty_folder(path: str) -> None:
    '''Create an empty folder at the given path. If already existed, delete it first.

    Args:
    - path: Path to the folder to create.
    '''
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=False)


def write_json_file(path: str, data: Union[Dict, List]):
    '''Write Jsoniazable data to the given path.
    '''
    with open(path, 'w') as fout:
        json.dump(data, fout, indent=2, sort_keys=True)


def touch(path: str, timeout: int = 60) -> None:
    '''Create an empty file at the given path.

    Args:
    - path: Path where to be touched.
    - timeout: Time out in seconds.
    '''
    command = 'touch %s' % path
    sp.check_output(command, shell=True, timeout=timeout)


def test(path: str, flag: str, timeout: int = 60) -> bool:
    '''Test the existence of the given path.

    Args:
    - path: Path to test.
    - flag: `-d`(folder), `-f`(file) or `-e`(both).
    - timeout: Time out in seconds.
    '''
    command = 'test %s %s' % (flag, path)
    try:
        sp.check_output(command, shell=True, timeout=timeout)
        return True
    except:
        return False


__all__ = [
    'load_json_file',
    'load_yaml_file',
    'make_empty_folder',
    'write_json_file'
]
