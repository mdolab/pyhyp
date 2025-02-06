from setuptools import setup, find_packages
import re
from os import path

__version__ = re.findall(
    r"""__version__ = ["']+([0-9\.]*)["']+""",
    open("pyhyp/__init__.py").read(),
)[0]

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

with open("doc/requirements.txt") as f:
    docs_require = f.read().splitlines()

setup(
    name="pyhyp",
    version=__version__,
    description="pyHyp uses hyperbolic volume mesh marching schemes to extrude structured surface meshes into volume meshes.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    keywords="meshing mesh-generation optimization",
    author="",
    author_email="",
    url="https://github.com/mdolab/pyhyp",
    license="Apache License Version 2.0",
    packages=find_packages(
        where=".",
        include=["pyhyp"],
    ),
    package_data={"pyhyp": ["*.so"]},
    install_requires=[
        "numpy>=1.21,!=1.24,!=1.24.1,!=1.24.2",
        "mpi4py>=3.1.5",
        "mdolab-baseclasses>=1.3",
        "cgnsutilities>=2.5",
        "tabulate",
    ],
    extras_require={
        "docs": docs_require,
        "testing": ["pygeo>=1.5", "testflo"],
    },
    classifiers=["Operating System :: Linux", "Programming Language :: Python, Fortran"],
)
