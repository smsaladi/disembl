from distutils.core import setup, Extension
setup(
  name='disembl',
  packages=['disembl'],
  version='2.0',
  description='A refactoring/reimplementation of DisEMBL',
  author='Shyam Saladi',
  author_email='saladi@caltech.edu',
  url='https://github.com/smsaladi/disembl',
  download_url='https://github.com/smsaladi/disembl/tarball/2.0',
  keywords=['protein', 'disorder', 'sequence', 'bioinformatics'],
  license='Non-commercial Academic Use License',
  intsall_requires=['numpy', 'ctypes', 'scipy', 'biopython', 'pandas'],
  ext_modules=[Extension('disembl.libdisembl', sources=['disembl/libdisembl.c'])]
)
