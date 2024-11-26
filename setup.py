#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Note: To use the 'upload' functionality of this file, you must:
#   $ pipenv install twine --dev

import io
import os
import sys
import subprocess
from shutil import rmtree

from setuptools import find_packages, setup, Command

# Package meta-data.
NAME = 'nemophoto'
DESCRIPTION = 'Absorption, fluorescence, phosphorescence spectrum simulations and ISC rate calculations.'
URL = 'https://github.com/LeonardoESousa/NEMO'
EMAIL = 'ledso@dtu.dk'
AUTHOR = 'Leonardo Evaristo de Sousa'
REQUIRES_PYTHON = '>=3.6.0'

# What packages are required for this module to be executed?
REQUIRED = ['numpy', 'scipy', 'pandas', 'LeoX>=0.9.0','wheel']

# What packages are optional?
EXTRAS = {
    # 'fancy feature': ['django'],
}

# The rest you shouldn't have to touch too much :)
# ------------------------------------------------
# Except, perhaps the License and Trove Classifiers!
# If you do change the License, remember to change the Trove Classifier for that!

here = os.path.abspath(os.path.dirname(__file__))

# Import the README and use it as the long-description.
# Note: this will only work if 'README.md' is present in your MANIFEST.in file!
try:
    with io.open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION

# Load the package's __version__.py module as a dictionary.
about = {}
project_slug = 'nemo'
with open(os.path.join(here, project_slug, '__version__.py'),encoding='utf-8') as f:
    exec(f.read(), about)

class UploadCommand(Command):
    """Support setup.py upload."""

    description = 'Build and publish the package.'
    user_options = []

    @staticmethod
    def status(s):
        """Prints things in bold."""
        print(f'\033[1m{s}\033[0m')

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        try:
            self.status('Removing previous builds…')
            rmtree(os.path.join(here, 'dist'))
        except OSError:
            pass

        self.status('Building Source and Wheel (universal) distribution…')
        subprocess.run([sys.executable, 'setup.py', 'sdist', 'bdist_wheel', '--universal'], check=True)

        self.status('Uploading the package to PyPI via Twine…')
        subprocess.run(['twine', 'upload', 'dist/*'], check=True)

        sys.exit()


# Where the magic happens:
setup(
    name=NAME,
    version=about['__version__'],
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    packages=find_packages(exclude=["tests", "*.tests", "*.tests.*", "tests.*"]),
    # If your package is a single module, use this instead of 'packages':
    #py_modules=['mypackage'],
    entry_points={'console_scripts': ["nemo = nemo.__main__:main", "nemo_batch_run = nemo.batch:run_batch", "nemo_tuning = nemo.omega:main"]},
    install_requires=REQUIRED,
    extras_require=EXTRAS,
    include_package_data=True,
    license='MIT',
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
    ],
    # $ setup.py publish support.
    cmdclass={
        'upload': UploadCommand,
    },
)
