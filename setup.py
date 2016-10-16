from setuptools import setup, Extension

setup(
  name='disembl',
  packages=['disembl'],
  version='2.0rc0',
  description='A refactoring/reimplementation of DisEMBL',
  author='Shyam Saladi',
  author_email='saladi@caltech.edu',
  url='https://github.com/smsaladi/disembl',
  download_url='https://github.com/smsaladi/disembl/tarball/2.0',
  keywords=['protein', 'disorder', 'sequence', 'bioinformatics'],
  license='GNU General Public License v2 (GPLv2)',
  ext_modules=[Extension('disembl.libdisembl',
                         sources=['disembl/libdisembl.c'])],
  scripts=['disembl/scripts/DisEMBL.py'],
  setup_requires=['pytest-runner'],
  tests_require=['pytest'],
)
