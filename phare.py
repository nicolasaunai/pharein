#!/usr/bin/env python3

import sys
import os
import subprocess
import shutil


class MiniPHARE:

    def __init__(self):
        self.simu = None
        self.miniphare_path = shutil.which('miniphare')
        if self.miniphare_path is None:
            raise RuntimeError("miniphare not found")

    def ini_path(self):
        return 'phare.ini'

    def command(self):
        return [self.miniphare_path, self.ini_path()]

    def run(self, simu):
        self.simu = simu
        print(' '.join(self.command()))
        os.chdir(simu.path)
        subprocess.call(self.command())

    def run_from_path(self, simu_path):
        os.chdir(simu_path)
        subprocess.call(self.command())


if __name__ == '__main__':

    simu_path = sys.argv[1]
    MiniPHARE().run_from_path(simu_path)