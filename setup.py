from setuptools import setup, find_packages

setup(name='smfplot',
      version='0.1.0',
      packages=find_packages(),
      description='A python package to visualize single locus single molecule footprinting data.',
      url='https://github.com/jkmckenna/smfplot',
      author="Joseph McKenna",
      author_email="jkmckenna@berkeley.edu",
      license='MIT',
      packages=['smfplot'],
      install_requires=[
          'os',
          'datetime',
          'numpy',
          'scipy',
          'pandas',
          'anndata',
          'matplotlib',
          'seaborn',
          'scanpy',
          'sklearn'
      ]
)