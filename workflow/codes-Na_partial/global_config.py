#coding:utf-8
import copy
import enum
import json
import os

from typing import Dict, List, Optional, Union


# ====== Paths in Docker image =====
kG16InitEnvFile = '/opt/tiger/g16_root/g16/bsd/g16.login'
kPublicFilesDir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'assets')
kCshInterpreter = '/bin/csh'

# ===== Paths in HDFS =====
kSuccessFlag = '_SUCCESS'
kBugFlag = '_BUG'


SERVERPORT = 2345


class StatusEnum(enum.Enum):
    """base class for enum to mark the status."""
    def __str__(self):
        return self.name

    @classmethod
    def names(cls) -> List[str]:
        """
        :return: a list of names of all members of this enum
        """
        return [e.name for e in cls]

    @classmethod
    def values(cls) -> list:
        """
        :return: a list of values of all members of this enum
        """
        return [e.value for e in cls]

    def __eq__(self, other: Union['StatusEnum', str]) -> bool:
        if isinstance(other, str):
            return self.name == other
        elif isinstance(other, StatusEnum):
            return self.name == other.name
        return False


class TaskType(enum.Enum):
    '''Tell the task type.'''
    DFT1 = 'moldft1'
    MD1 = 'molmd1'
    DFT2 = 'moldft2'
    MD_RESP2  = 'molmd_RESP2'
    DFT_RESP2 = 'moldft_RESP2'
    DFT_Li = 'moldft_Li'
    DFT3 = 'moldft3'
    MD2 = 'molmd2'

    @staticmethod
    def prev_task_type(current: 'TaskType') -> Optional['TaskType']:
        '''Get the previous task type.'''
        found = False
        previous = None
        for _, candidate in enumerate(TaskType):
            if candidate == current:
                found = True
                break
            previous = candidate
        assert found, 'Illegal current task type: ' + current.type
        return previous

    @classmethod
    def next(cls, current: 'StatusEnum') -> Optional['StatusEnum']:
        """ return next item of the enum"""
        item_list = list(cls)
        assert current in item_list, f"{current.name} not in {item_list}"
        if item_list.index(current) == len(item_list) - 1:
            return None
        current_index = item_list.index(current)
        return item_list[current_index + 1]

    @classmethod
    def previous(cls, current: 'StatusEnum') -> Optional['StatusEnum']:
        """ return previous item of the enum"""
        item_list = list(cls)
        assert isinstance(current, cls), f"Type error, need {type(cls)}, found {type(current)}"
        assert current in item_list, f"{current.name} not in {item_list}"
        if item_list.index(current) == 0:
            return None
        current_index = item_list.index(current)
        return item_list[current_index - 1]


class TaskAction(enum.Enum):
    '''Tell the next action to handle task.'''
    COMPLETE = 'COMPLETE'  # Finished successfully
    RECOVER = 'RECOVER'  # Failed due to transient issue
    NOTICE = 'NOTICE'  # Failed due to permanent issue
    DEBUG = 'DEBUG'   # Failed due to business logic


class WorkerStatus(StatusEnum):
    """Tell the status of the worker at the server side"""
    IDLE = "IDLE"  # No job is assigned
    BUSY = "BUSY"  # A job is assigned
    CLOSED = "CLOSED"  # Worker is closed, no longer accepting new jobs
    NOT_FOUND = "NOT_FOUND"  # Worker is not found in the server


class JobStatus(StatusEnum):
    """
        The jobs status in the workflow.
        Not related to the hdfs flag.
    """
    QUEUE = "QUEUE"  # Job has no status yet
    RUNNING = "RUNNING"  # Running in a worker
    SUCCESS = "SUCCESS"  # Job return with code 0
    DEBUG = "DEBUG"  # Job return with code non-zero
    NOTICE = "NOTICE"  # Job not running as there is no worker
    NOT_FOUND = "NOT_FOUND"  # Job is not found in the request field
    UNKNOWN = "UNKNOWN"  # Unknown status


class ServerStatus(StatusEnum):
    """the stutus to RPC response to the client"""
    NOT_FOUND = "NOT_FOUND"  # ConnectionRefusedError
    NOT_INITIALIZED = "NOT_INITIALIZED"  # Server has no config input
    SUCCESS = "SUCCESS"  # response nornally
    CLOSING = "CLOSING"  # Server is closing
    DEBUG = "DEBUG"  # Server response not normally, but still accept the request
    RUNNING = "RUNNING"  # Server is running
    NO_SENDING = "NO_SENDING"  # the request is closed by the client
    UNKNOWN = "UNKNOWN"  # Unknown status, normally bad requests or internal error


class EnvFlag(StatusEnum):
    """A simple flag to identify the env"""
    SERVER = "SERVER"
    CLIENT = "CLIENT"


class EpochConfig(object):
    '''Describe workload in this epoch.'''

    def __init__(self, path) -> None:
        '''Create a `EpochConfig` instance.'''
        with open(path, 'r', encoding="utf-8") as fin:
            self._raw_data = json.load(fin)

    def copy(self, molecules: List[str] = None, index: int = None) -> 'EpochConfig':
        '''Return a deep copy.

        Args:
        - molecules: Update property `molecules` with the given value.
        - index: Update property `index` with the given value.

        Returns:
        - copied: A new config with updated properties.
        '''
        copied = copy.deepcopy(self)
        if molecules is not None:
            copied.molecules = molecules
        if index is not None:
            copied.index = index
        return copied

    def to_file(self, path: str) -> None:
        '''Write a JSON-format representation into the given file path.

        Args:
        - path: Where to write or overwrite.
        '''
        with open(path, 'w', encoding="utf-8") as fout:
            json.dump(self._raw_data, fout, indent=2)

    @property
    def remote_moldef_dir(self) -> str:
        '''Return the path of HDFS folder containing defintion files of all molecules.'''
        return self._raw_data['remote_moldef_dir']

    @property
    def molecules(self) -> List[str]:
        '''Return a list of molecules as solvent.'''
        return self._raw_data['molecules']

    @molecules.setter
    def molecules(self, mols: list) -> None:
        self._raw_data['molecules'] = mols

    @property
    def index(self) -> int:
        '''Return the sampling index of task7 moldft3.'''
        return self._raw_data['index']

    @index.setter
    def index(self, index: int) -> None:
        self._raw_data['index'] = index

    @property
    def salt(self) -> str:
        '''Return the salt as solute.'''
        return self._raw_data['salt']

    @property
    def remote_epoch_dir(self) -> str:
        '''Return the path of HDFS base folder containing all task results.'''
        return self._raw_data['remote_epoch_dir']

    @property
    def remote_tmp_dir(self) -> str:
        '''Return the path of HDFS temp folder under the epoch folder.'''
        return os.path.join(self.remote_epoch_dir, '_tmp')

    @property
    def max_server_count(self) -> int:
        """
            Return the max server count from config file. If not set, return 9999.
        """
        return self._raw_data.get("max_server_count", 9999)

    def remote_output_dir(self, molecule: str, task: TaskType) -> str:
        '''Return the HDFS path of task folder containing a task result.

        Args:
        - molecule: Molecule name.
        - task: Which task to run.

        Returns:
        - hdfs_dir: HDFS folder path of the given task.
        '''
        name = task.value
        subdir = 'MD' if name.startswith('molmd') else 'DFT'
        return os.path.join(self.remote_epoch_dir, molecule, subdir, name)

    def remote_tmp_output_dir(self, molecule: str, task: TaskType) -> str:
        """Return the temp HDFS path of task folder containing a task result.
        As the folder is used in mutiple places, we create a common function here.
        """
        return self.remote_output_dir(molecule, task) + '.tmp'

    def remote_success_path(self, molecule: str, task: TaskType, samp: int = None) -> str:
        """Return the HDFS path of success flag of a task.
        """
        if samp is None:
            suffix = ''
        else:
            suffix = '_' + str(samp)
        return os.path.join(self.remote_output_dir(molecule, task) + suffix, kSuccessFlag)

    def remote_bug_path(self, molecule: str, task: TaskType, samp: int = None) -> str:
        """Return the HDFS path of success flag of a task.
        """
        if samp is None:
            suffix = ''
        else:
            suffix = '_' + str(samp)
        return os.path.join(self.remote_output_dir(molecule, task) + suffix, kBugFlag)

    def remote_topology_path(self, molecule: str) -> str:
        '''Return the HDFS path of molecule topology definition.'''
        return os.path.join(self.remote_moldef_dir, molecule + '.mol')

    @property
    def arnold_kwargs(self) -> Dict:
        '''Return Arnold configuration to manage all submitted tasks & trials.

        Required fields are:
        - token: `str`
        - job_id: `int`
        - group_ids: `list` of `int`
        - cluster_id: `int`
        - preemptible: `bool`
        '''
        return self._raw_data['arnold_kwargs']

    @property
    def task_type(self) -> TaskType:
        """
            Return the task type with TaskType wrapper.
        """
        return TaskType(self._raw_data.get("task_type"))

    @task_type.setter
    def task_type(self, task_type: str) -> None:
        self._raw_data["task_type"] = task_type

__all__ = [
    'kPublicFilesDir',
    'kSuccessFlag',
    'EpochConfig',
    'TaskType',
    'TaskAction'
]


if __name__ == "__main__":
    ... 
