#! /usr/bin/env python

from setuptools import setup


setup(
    name="sippv",
    version='0.1.0',
    description=("Add PV corrections to FITS headers "
                 "with existing SIP corrections"),
    author="Evert Rol / GOTO observatory",
    author_email="evert.rol@monash.edu",
    url="http://goto-observatory.org",
    license="ISC",
    scripts=["scripts/wcs-addpv"],
    packages=["sippv"],
    package_dir={'gototile': "gototile"},
    install_requires=["numpy", "scipy", "astropy"],
    classifiers=['Development Status :: 5 - Production/Stable',
                 'License :: OSI Approved :: ISC License (ISCL)',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python',
                 'Intended Audience :: Science/Research',
                 'Topic :: Scientific/Engineering :: Astronomy',
                 'Topic :: Software Development :: Libraries :: Python Modules',],
    keywords=['astronomy', 'astrophysics', 'wcs', 'world coordinate system',
              'fits'],
)
