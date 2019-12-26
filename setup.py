#!/usr/bin/env python
from setuptools import setup,find_packages

setup(
    name='shiroin',
    version='1.0.1',
    author='Grzegorz Adamski',
    author_email='grzes-a1@o2.pl',
    python_requires='>=3.5',
    url='https://github.com/urojony/shiroin',
    license='BSD',
    description='Software for proving inequalities',
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    packages=find_packages(),
    install_requires=['scipy>=0.15','sympy'],
    classifiers=[
    "License :: OSI Approved :: BSD License",
    "Operating System :: OS Independent",
    "Intended Audience :: Science/Research",
    "Programming Language :: Python :: 3",
    "Topic :: Scientific/Engineering :: Mathematics"
    ],
)
