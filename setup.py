from setuptools import setup, find_packages
from glob import glob
import os

##############################
def find_scripts():
    scripts = []
    # All python files in scripts_electrolytes/scripts
    pyfiles = glob(os.path.join('scripts_electrolytes', 'scripts', "*.py"))
    scripts.extend(pyfiles)
    return scripts

##############################
requirements = ["numpy>=1.18",
                "abipy>=0.9",
                "ase>=3.21",
                "pandas>=1.5"]

setup(
    name="scripts_electrolytes",
    version="1.0.0",
    author="Véronique Brousseau-Couture",
    author_email="veronique.brousseauc@gmail.com",
    description="Utilities for the scripts_electrolytes repo",
    packages=find_packages(),
    scripts=find_scripts(),
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3.9",
    ],
)
