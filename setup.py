# -*- coding: utf-8 -*-
from setuptools import find_packages, setup

long_description = """MELODIST is an open-source toolbox written in Python for
disaggregating daily meteorological time series to hourly time steps. The
software framework consists of disaggregation functions for each variable
including temperature, humidity, precipitation, shortwave radiation, and wind
speed. These functions can simply be called from a station object, which
includes all relevant information about site characteristics. The data
management of time series is handled using data frame objects as defined in the
pandas package. In this way, input and output data can be easily prepared and
processed."""

setup(
    name='melodist',

    version='0.1.4',

    description='MELODIST: MEteoroLOgical observation time series DISaggregation Tool',
    long_description=long_description,

    url='https://github.com/kristianfoerster/melodist',

    author='Kristian Förster, Florian Hanzer, Benjamin Winter, Thomas Marke, Siling Chen',
    author_email='foerster@iww.uni-hannover.de',

    license='GPLv3',

    classifiers=[
        'Development Status :: 4 - Beta',

        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Atmospheric Science',

        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],

    packages=find_packages(exclude=['docs']),

    install_requires=[
        'numpy',
        'scipy',
        'pandas',
    ],

    extras_require={
        'test': ['pytest'],
    },
)
