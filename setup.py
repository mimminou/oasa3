#!/usr/bin/env python

import setuptools


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
      name             = 'oasa3',
      version          = '0.14.1',
      description      = "OASA is a free cheminformatics library written in Python",
      author           = "Beda Kosata",
      author_email     = "beda@zirael.org",
      maintainer       = "Reinis Danne",
      maintainer_email = "rei4dan@gmail.com",
      project_urls={
      'Website': 'http://bkchem.zirael.org/oasa_en.html',
      'Git repository': 'https://gitlab.com/oasa/oasa',},
      license          = "GNU GPL",
      platforms        = ["Unix", "Windows", "hopefully other OSes able to run Python"],
      install_requires = ['pycairo',],
      long_description_content_type="text/markdown",
      long_description = long_description,
      classifiers=[
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: Unix",
      ],
      package_dir={"": "src"},
      packages=setuptools.find_packages(where="src"),
      python_requires=">=2.6",
     )

