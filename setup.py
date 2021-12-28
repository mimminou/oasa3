
import setuptools


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
  name = 'oasa3',
  version = '0.1.1',
  description = "OASA Chemoinformatics library ported to python 3",
  author = "Abdelaziz Amine",
  author_email = "mimminou@gmail.com",
  url = "https://github.com/mimminou/oasa3",
  license = "GNU GPL",
  platforms = ["Unix", "Windows"],
  long_description = long_description,
  long_description_content_type="text/markdown",
  classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: Unix",
  ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
  )

