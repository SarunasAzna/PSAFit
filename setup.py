#!/usr/bin/python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='PSAFit',
    author='BPTI',
    packages=[
        'DCanalysis'
        'DCanalysis.tests',
    ],
)
