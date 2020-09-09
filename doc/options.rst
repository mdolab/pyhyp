.. _pyhyp_options:

Options
=======

Here are the options currently available in pyHyp.

===================  ==========  =======================================================================================
Parameter              Type      Description
===================  ==========  =======================================================================================
``inputFile``         ``char``   Name of the file that contains the surface mesh.
                                 This is a file that has been generated in an external meshing program, typically ICEMCFD.

``fileType``          ``char``   Type of the input file. Use either ``Plot3d`` or ``CGNS``.

``N``                 ``int``    Number of grid levels to march.
                                 This determines the grid dimension in the off-wall direction.
                                 Typically this should be a "multi-grid" friendly number.

``s0``               ``float``   Initial off-wall (normal) spacing of grid.
                                 This is taken to be constant across the entire geometry.
                                 The units are consistent with the rest of the geometry.

``rMin``             ``float``   Relative distance in the normal direction to march.
                                 It is specified as a multiple of the radius of the sphere enclosing the initial surface geometry.
                                 If symmetry is specified, the full mirrored geometry is used to compute the sphere's radius.
                                 Most wing geometries will have ``rMin`` between 10 and 20, that is the farfield boundary is 20 spans away from the geometry.

``cMax``             ``float``   The maximum permissible ratio of marching direction length to the any other in-plane edge.
                                 This parameter effectively operates as a CFL-type limit.
                                 If a step would require a step which would result in a ratio ``c`` greater than ``cMax``, the step is automatically split internally to respect this user-supplied limit.
                                 Typical values of ``cMax`` are around 6-8.
                                 Increased robustness can be achieved at the expense of computational cost by lowering ``cMax``.

``nonLinear``        ``float``   Use the nonlinear formulation.
                                 This is experimental and not currently recommended and may not work at all.

``slExp``            ``float``   Exponent for the :math:`S_l` computation.
                                 The :math:`S_l` value serves the same purpose as found in Chan et al. but the computation is different.
                                 The :math:`S_l` computation in Chan is given as :math:`\sqrt{\frac{N-1}{l-1}}` for :math:`l > 2`.
			       
``ps0``              ``float``   Initial pseudo offwall spacing.
                                 This spacing **must** be less than or equal to ``s0``.
                                 This is actual spacing the hyperbolic scheme uses.
                                 The solver may take many pseudo steps before the first real grid level at ``s0``.

``pGridRatio``       ``float``   The ratio between successive levels in the pseudo grid.
                                 This will be typically somewhere between ~1.05 for large grids to 1.2 for small grids.
                                 This number is **not** the actual grid spacing of the final grid; that spacing ratio is computed and displayed at the beginning of a calculation.
                                 The ``pGridRatio`` **must** be smaller than that number.

``epsE``             ``float``   The explict smoothing parameter.
                                 See the :ref:`Theory<pyhyp_theory>` section for more information.
                                 Typical values are approximately 1.0. Increasing the explicit smoothing may result in a smoother grid, at the expense of orhtogonality.
                                 If the geometry is very sharp corners, too much explicit smoothing will cause the solver to rapidly "soften" the corner and the grid will fold back on itself.
                                 In concave corners, additional smoothing will prevent lines from crossing (avoiding negative cells).

``epsI``             ``float``   Implicit smoothing parameter.
                                 See the :ref:`Theory<pyhyp_theory>` section for more information.
                                 Typical values are from 2.0 to 6.0.
                                 Generally increasing the implicit coefficient results in a more stable solution procedure.
                                 Usually this value should be twice the explicit smoothing parameter.

``theta``            ``float``   Kinsley-Barth coefficient See the :ref:`Theory<pyhyp_theory>` section for more information.
                                 Only a single theta value is used for both directions.
                                 Typical values are ~2.0 to ~4.0.

``volCoef``          ``float``   Coefficient used in point-Jacobi local volume smoothing algorithm.
                                 Typically this value is 0.16 and need not be modified.
                                 Use more ``volSmoothIter`` for stronger local smoothing.
			    
``volBlend``         ``float``   The global volume blending coefficient.
                                 See the :ref:`Theory<pyhyp_theory>` section for more information.
                                 This value will typically be very small, especially if you widely varying cell sizes.
                                 Typically values are from ~0 to 0.001.
                                 Default is 0.0001.

``volSmoothIter``    ``int``     The number of point-Jacobi local volume smoothing iterations to perform.
                                 Typical values are ~5 to ~25.
                                 Default is 10.

``kspRelTol``        ``float``   Tolerance for the solution of the linear system at each iteration.
                                 Typically :math:`1\times 10^{-8}` is sufficient.
                                 Very difficult cases may benefit from a tighter convergence tolerance.

``kspMaxIts``        ``int``     Maximum number of iterations to perform for each step.
                                 Default is 500 which should be sufficient for most cases.

``preConLag``        ``int``     Lag the update of the preconditioner by this number of iterations.
                                 The default value of 10 will typically not need to be changed.

``kspSubspaceSize``  ``int``     Size of the ksp subspace.
                                 Default is 50.
                                 Very large and difficult problems may befefit from a larger subspace size.

``writeMetrics``     ``bool``    Flag to write the mesh gradients to the solution file.
                                 This option should only be used for debugging purposes.
===================  ==========  =======================================================================================