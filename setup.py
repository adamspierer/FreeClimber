from setuptools import setup, find_packages

# Get the long description from the relevant file
with open('README.md','r') as f:
    long_description = f.read()

setup(name='FreeClimber-FNG',
      version='0.4.0',
      description='FreeClimber-FNG is a fork of FreeClimber adding fast negative geotaxis (FNG) event detection to the Python-based background subtraction and climbing velocity estimation pipeline.',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/jordanvasu/FreeClimber-FNG',

      author='Jordan Vasu',
      author_email='anspierer+Github_setup_py@gmail.com',
      license='MIT',

      classifiers=[
      'Development Status :: 4 - Beta',
      'Intended Audience :: Developers',
      'License :: OSI Approved :: MIT License',
      'Programming Language :: Python :: 3.6',
      'Programming Language :: Python :: 3.7'],

      keywords='Drosophila melanogaster climbing negative geotaxis assay background subtraction particle detection local linear regression high-throughput high throughput automated behavior quantification velocity',

      packages=find_packages(),
#       packages=['FreeClimber'],
      install_requires=['ffmpeg-python==0.2.0',"argparse==1.1",
                        'pandas','numpy','scipy','pip','matplotlib==3.1.3',
                        'wxPython==4.0.4','trackpy==0.4.2'],
      zip_safe=False)
