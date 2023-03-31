from setuptools import setup, find_packages

requirements = ["numpy>=1.18"]

setup(
    name="utils",
    version="1.0.0",
    author="Véronique Brousseau-Couture",
    author_email="veronique.brousseauc@gmail.com",
    description="Utilities for the scripts_electrolytes repo",
    packages=find_packages(),
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3.9",
    ],
)
