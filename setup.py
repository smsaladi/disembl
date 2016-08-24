from distutils.core import setup
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
  license='Non-commercial Academic Use License'
  intsall_requires=['numpy', 'ctypes','scipy','biopython', 'pandas']
  ext_modules=[Extension('_disembl',
                             ['disembl.c'],
                             extra_link_args=['-Wall', '-Werror', '-O2',
                                                '-fpic', '-shared', '-o',])]
)
