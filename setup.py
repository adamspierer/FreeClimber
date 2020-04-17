from setuptools import setup

'''Setup for FreeClimber
'''

from os import path
from codecs import open
from setuptools import setup, find_packages

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'DESCRIPTION.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='free_climber',
      version='0.5',
      description='A background-adjusted, particle detector used to quantify climbing performance for groups of Drosophila in a negative geotaxis assay',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='http://github.com/adamspierer/free_climber',
      
      author='Adam N. Spierer',
      author_email='anspierer+free_climber@gmail.com',
      license='MIT',

      classifiers=[
      'Development Status :: 4 - Beta',
      'Intended Audience :: Developers',
      'License :: OSI Approved :: MIT License',
      'Programming Language :: Python :: 3.7'],

      keywords='Drosophila climbing particle detection local linear regression',

      packages=setuptools.find_packages(),
      packages=['free_climber'],
      install_requires=['os','sys','argparse','subprocess','time','pandas','numpy','scipy',
                        'pip','matplotlib==3.1.3','wxPython==4.0.4','trackpy==0.4.2'],
      zip_safe=False)
