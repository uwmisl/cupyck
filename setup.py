from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import subprocess, os

class build(build_py):
    def run(self):

        nvcc_missing = False
        try:
            subprocess.check_call(["nvcc", "--version"])
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                nvcc_missing = True

        if not nvcc_missing:
            subprocess.check_call("make")

        return build_py.run(self)


setup(
  name = 'cupyck',
  version = "1.0.0",
  packages=find_packages(),
  package_data = {'cupyck':['cupyck.so', 'parameters/*']},
  include_package_data = True,
  cmdclass = {"build_py": build}
)
