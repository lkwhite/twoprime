#!/usr/bin/env python

'''
twoprime: analysis of twoprime-seq data

'''

from ez_setup import use_setuptools
use_setuptools()

from setuptools import find_packages, setup

__version__ = '0.01a'

entry_points = """
[console_scripts]
twoprime-process-signals = twoprime.process_signals:main
"""

install_requires = ["genomedata>1.3.1",  "numpy"]

if __name__ == '__main__':

    setup(name='twoprime',
          version='0.01a',
          description='Analysis of twoprime sequencing data',
          author='Jay Hesselberth',
          author_email='jay.hesselberth@gmail.com',
          packages=['twoprime'],
          install_requires=install_requires,
          entry_points=entry_points
         )
