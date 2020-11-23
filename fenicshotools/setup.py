#!/usr/bin/env python
import sys, os
from setuptools import setup, Command
from setuptools.command.test import test as TestCommand

class PyTest(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        import pytest
        os.chdir("tests/unit")
        errno = pytest.main(self.pytest_args)
        sys.exit(errno)

class Clean(Command):
    """
    Cleans *.pyc so you should get the same copy as is in the VCS.
    """

    description = "remove build files"
    # FIXME: What should these be used for?
    user_options = [("all","a","the same")]

    def initialize_options(self):
        self.all = None

    def finalize_options(self):
        pass

    def run(self):
        import os
        os.system("utils/clean-files")

class BuildDocs(Command):
    """
    Build documentation
    """

    description = "build documentation"
    user_options = []

    def initialize_options(self):
        self.all = None

    def finalize_options(self):
        pass

    def run(self):
        import os
        curdir = os.path.abspath(".")
        os.chdir("docs")
        os.system("make build_demos")
        os.system("make html latex")
        os.chdir("_build/latex")
        os.system("make")
        os.chdir(curdir)

setup(name = 'fenicshotools',
      version = '0.3',
      description = """
        Higher order elements tools for FEniCs using gmsh.
      """,

      # Packages that should be installed
      packages = [ "fenicshotools",
                   "fenicshotools.gmsh",
                   "fenicshotools.io",
                   "fenicshotools.utils",
                   "fenicshotools.scripts",
                   "sphinxarg"
                   ],

      # We rely on pytest
      tests_require=['pytest'],

      # Additional build targets
      cmdclass = {"test" : PyTest,
                  "clean" : Clean,
                  "build_docs" : BuildDocs},

      # Make functions in the scripts callable from the command line
      # through named scripts (platform independent!)
      entry_points={
          'console_scripts': [
              'gmsh2dolfin = fenicshotools.scripts.gmsh2dolfin:main_func',
              'geo2dolfin = fenicshotools.scripts.geo2dolfin:main_func',
              ]}
       )
