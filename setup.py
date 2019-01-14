from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import subprocess

class build(build_py):
    def run(self):
        subprocess.check_call("make")
        build_py.run(self)

setup(
  name = 'cupyck',
  version = "1.0.0",
  packages=find_packages(),
  package_data = {'cupyck':['cupyck.so', 'parameters/*']},
  include_package_data = True,
  cmdclass = {"build_py": build}
)
