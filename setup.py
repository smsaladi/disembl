from distutils.core import setup, Extension
from setuptools.command.test import test as TestCommand

class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_suite = True
    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        import pytest
        pytest.main()

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
  license='Non-commercial Academic Use License',
  intsall_requires=['numpy', 'ctypes', 'scipy', 'biopython', 'pandas'],
  ext_modules=[Extension('disembl.libdisembl',
                         sources=['disembl/libdisembl.c'])],
  scripts=['disembl/scripts/DisEMBL.py'],
  zip_safe=True,
  tests_require=['pytest'],
  cmdclass = {'test': PyTest},
)
