#!/usr/bin/env python
from distutils.core import setup

setup(
    name='shiroin',
    version='1.0.0',
    packages=['shiroin',],
    install_requires=['scipy','sympy'],
    author='Grzegorz Adamski',
    author_email='grzes-a1@o2.pl',
    python_requires='>=3.5',
    url='https://github.com/urojony/shiroin',
    license='BSD',
    description='Software for proving inequalities',
    long_description=open('README.md').read(),
)
