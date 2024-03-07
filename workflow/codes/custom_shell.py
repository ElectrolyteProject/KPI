#coding:utf-8
import subprocess as sp
import sys

from typing import Dict


class ShellError(Exception):
    '''Describe the error of a Shell command execution.'''

    def __init__(self, sh_command, exit_code) -> None:
        '''Create a `ShellError` instance.'''
        super(ShellError, self).__init__()
        self._sh_command = sh_command
        self._exit_code = exit_code

    @property
    def sh_command(self) -> int:
        return self._sh_command

    @property
    def exit_code(self) -> int:
        return self._exit_code

    @property
    def message(self) -> str:
        return self.args[0]

    def __str__(self) -> str:
        return '''Failed to execute: %s\nExit code:%d''' % (self.sh_command, self.exit_code)


def check_output(sh_command: str, working_dir: str = None, executable: str = None) -> str:
    '''Execute a Shell command and return its STDOUT. STDERR is shared with the caller.

    Args:
    - sh_command: Shell command to execute.
    - working_dir: Working folder to use.
    - executable: Path to shell interpreter.

    Returns:
    - message: Text in STDOUT.
    '''
    try:
        binary = sp.check_output(args=sh_command,
                                 bufsize=0,
                                 executable=executable,
                                 stderr=sys.stderr,
                                 shell=True,
                                 cwd=working_dir)
        return binary.decode()
    except sp.CalledProcessError as e:
        raise ShellError(sh_command, e.returncode)


def execute(sh_command: str, working_dir: str = None, env: Dict = None, executable: str = None) -> None:
    '''Execute a Shell command in a blocking way. STDOUT and STDERR are shared with the caller.

    Args:
    - sh_command: Shell command to execute.
    - working_dir: Working folder to use.
    - env: Defines the environment variables for the new process.
    - executable: Path to shell interpreter.
    '''
    process = sp.Popen(args=sh_command,
                       bufsize=0,
                       executable=executable,
                       stdout=sys.stdout,
                       stderr=sys.stderr,
                       shell=True,
                       cwd=working_dir,
                       env=env)
    process.communicate(timeout=60)
    returncode = process.wait(timeout=5)
    if returncode != 0:
        raise ShellError(sh_command, returncode)


__all__ = [
    'ShellError',
    'check_output',
    'execute'
]
