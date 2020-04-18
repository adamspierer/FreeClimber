from setuptools import setup

'''Setup for FreeClimber
'''

from os import path
from codecs import open
from setuptools import setup, find_packages

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='FreeClimber',
      version='0.1',
      description='A background-adjusted, particle detector used to quantify climbing performance for groups of Drosophila in a negative geotaxis assay',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='http://github.com/adamspierer/FreeClimber',
      
      author='Adam Spierer',
      author_email='anspierer+Github_setup_py@gmail.com',
      license='MIT',

      classifiers=[
      'Development Status :: 4 - Beta',
      'Intended Audience :: Developers',
      'License :: OSI Approved :: MIT License',
      'Programming Language :: Python :: 3.7'],

      keywords='Drosophila climbing negative geotaxis background subtraction particle detection local linear regression',

      packages=find_packages(),
#       packages=['FreeClimber'],
      install_requires=['ffmpeg==1.4',"argparse==1.1"
                        'pandas','numpy','scipy','pip','matplotlib==3.1.3',
                        'wxPython==4.0.4','trackpy==0.4.2'],
      zip_safe=False)
