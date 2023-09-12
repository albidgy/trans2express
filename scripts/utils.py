import os
import subprocess
import sys


def write_logs(stdout, stderr, common_output_dir):
    with open(f'{common_output_dir}trans2express.log', 'a') as stdout_f:
        stdout_f.write(stdout)

    with open(f'{common_output_dir}trans2express.log', 'a') as stderr_f:
        stderr_f.write(stderr)


def run_external_tool(cmd, common_output_dir):
    process = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, close_fds=True, universal_newlines=True)
    stdout, stderr = process.communicate()
    write_logs(stdout, stderr, common_output_dir)
    if process.returncode != 0:
        print(f'The tool returned non-zero. See {os.path.abspath(common_output_dir)}/trans2express.log file', file=sys.stderr)
        sys.exit(1)


class Logger:
    def __init__(self, common_output_dir):
        self.terminal = sys.stdout
        self.log = open(f'{common_output_dir}trans2express.log', 'a', buffering=1)

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        pass
