from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='oskar',
      version='0.0.1',
      description='process data from oskar experiments',
      url='',
      author='Adam Deller',
      author_email='a.deller@ucl.ac.uk',
      license='BSD',
      packages=['oskar'],
      install_requires=[
          'scipy>=0.14','numpy>=1.10','pandas>=0.17', 'h5py>=2.5', 'tqdm>=3.1.4', 'sspals>=0.0.1'
      ],
	  entry_points = {
        'console_scripts': ['oskar_dset=oskar.dset:main',
                            'oskar_info=oskar.info:main',
                            'oskar_average=oskar.average:main',
                            'oskar_vrange=oskar.vrange:main',
                            'oskar_sspals=oskar.sspals_:main',
                            'oskar_count=oskar.count:main'],
	  },
      include_package_data=True,
      zip_safe=False)
