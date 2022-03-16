# pyHyp
[![Build Status](https://dev.azure.com/mdolab/Public/_apis/build/status/mdolab.pyhyp?branchName=main)](https://dev.azure.com/mdolab/Public/_build/latest?definitionId=13&branchName=main)
[![Documentation Status](https://readthedocs.com/projects/mdolab-pyhyp/badge/?version=latest)](https://mdolab-pyhyp.readthedocs-hosted.com/en/latest)
[![codecov](https://codecov.io/gh/mdolab/pyhyp/branch/main/graph/badge.svg?token=QJ9DFTGPYX)](https://codecov.io/gh/mdolab/pyhyp)

pyHyp uses hyperbolic volume mesh marching schemes to extrude structured surface meshes into volume meshes.
pyHyp is used as a preprocessing step in the geometry and mesh-creation process prior to an optimization.

## Documentation

Please see the [documentation](https://mdolab-pyhyp.readthedocs-hosted.com/en/latest) for installation details and API documentation.

To locally build the documentation, enter the `doc` folder and enter `make html` in terminal.
You can then view the built documentation in the `_build` folder.


## Citation

If you use pyHyp in any publication for which you find it useful, please cite this paper.

N. Secco, G. K. W. Kenway, P. He, C. A. Mader, and J. R. R. A. Martins, “Efficient Mesh Generation and Deformation for Aerodynamic Shape Optimization”, AIAA Journal, 2021. [doi:10.2514/1.J059491](https://doi.org/10.2514/1.J059491)

```
@article{Secco2021,
    title = {Efficient Mesh Generation and Deformation for Aerodynamic Shape Optimization},
    author = {Ney Secco and Gaetan K. W. Kenway and Ping He and Charles A. Mader and Joaquim R. R. A. Martins},
    doi = {10.2514/1.J059491},
    journal = {AIAA Journal},
    year = {2021}
}
```

## How pyHyp fits within MACH

pyHyp takes structured surface meshes and extrudes them into structured volume meshes.
This is done as a pre-processing step.
Generally, the surface meshes come from ICEM.
The pyHyp-generated volume meshes are then used in [ADflow](https://github.com/mdolab/adflow) to perform CFD.
An example [XDSM](https://github.com/mdolab/pyXDSM) for an optimization setup that uses pyHyp is shown below.

![pySurf XDSM diagram](doc/images/pysurf_xdsm.png)

## License

Copyright 2019 MDO Lab. See the LICENSE file for details.
