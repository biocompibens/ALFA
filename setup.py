#!/usr/bin/env python

from distutils.core import setup

setup(name = "ALFA",
      py_modules = ["ALFA"],
      version = "0.27.6",
      description = "A simple software to get a quick overview of features composing NGS dataset(s).",
      author = "Mathieu Bahin",
      author_email = "mathieu.bahin@biologie.ens.fr",
      maintainer = "Mathieu Bahin",
      maintainer_email = "mathieu.bahin@biologie.ens.fr",
      url = "https://github.com/biocompibens/ALFA",
      scripts=["ALFA"],
      long_description = open("README").read(),
      install_requires=["numpy>=1.15,<1.16", "pysam>=0.15,<0.16", "pybedtools>=0.8,<0.9", "matplotlib>=3.0,<3.1", "progressbar2>=3.37,<3.40"],
      license = "MIT"
)
