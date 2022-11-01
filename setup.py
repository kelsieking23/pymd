from setuptools import setup, find_packages

setup(
    name='pymd',
    version='0',
    license='GNU GENERAL PUBLIC LICENSE',
    author='Kelsie King',
    author_email='kelsieking23@vt.edu',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    url='https://github.com/kelsieking23/pymd',
    keywords='pymd',
    install_requires=[
          'pandas',
          'numpy',
          'matplotlib',
          'argparse'
      ],
)