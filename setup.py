from setuptools import setup

setup(name='mpanalysis',
      version='0.5',
      description='Pathway analyses based on mpdiag and tobac',
      url='http://github.com/mheikenfeld/mpanalysis',
      author='Max Heikenfeld',
      author_email='max.heikenfeld@physics.ox.ac.uk',
      license='GNU',
      packages=['mpanalysis'],
      install_requires=['mpdiag','tobac'],
      zip_safe=False)
