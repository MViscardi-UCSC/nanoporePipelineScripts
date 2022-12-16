"""
220928_testingRichForConsoleLogging.py
Marcus Viscardi,    September 28, 2022


"""
import subprocess
from rich.console import Console

def live_cmd_call(command):
    with subprocess.Popen(command, stdout=subprocess.PIPE, shell=True,
                          bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end="")
    if p.returncode != 0:
        raise subprocess.CalledProcessError(p.returncode, p.args)

if __name__ == '__main__':
    live_cmd_call("w")
