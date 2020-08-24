from setuptools import setup
import re

__version__ = re.findall(
    r"""__version__ = ["']+([0-9\.]*)["']+""",
    open('pyhyp/__init__.py').read(),
)[0]

setup(name='pyhyp',
      version=__version__,


      description="pyHyp uses hyperbolic volume mesh marching schemes to extrude structured surface meshes into volume meshes. ",
      long_description="""
      ## Documentation

      Please see the [documentation](http://mdolab.engin.umich.edu/docs/packages/pyhyp/doc/index.html) for installation details and API documentation.

      To locally build the documentation, enter the `doc` folder and enter `make html` in terminal.
      You can then view the built documentation in the `_build` folder.


      ## Citation

      pyHyp is based on the theory presented in "Enhancements of a three-dimensional hyperbolic grid generation scheme."
      For more background, theory, and figures, see the original [journal article](https://doi.org/10.1016/0096-3003(92)90073-A).

      """,
      long_description_content_type="text/markdown",
      keywords='meshing mesh-generation optimization',
      author='',
      author_email='',
      url='https://github.com/mdolab/pyhyp',
      license='Apache License Version 2.0',
      packages=[
          'pyhyp',
      ],
      package_data={
          'pyhyp': ['*.so']
      },
      install_requires=[
            'numpy>=1.16',
            'mpi4py>=3.0',
      ],
      classifiers=[
        "Operating System :: Linux",
        "Programming Language :: Python, Fortran"]
      )
