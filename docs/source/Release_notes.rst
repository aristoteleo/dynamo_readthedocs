==================
Release notes
==================


Dynamo Ver 1.3.0
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
Feature Changes
~~~~~~~~~~~~~~~~~~~~~~~~~~
- The preprocessing module has been refactored:

  - Class *Preprocessor* is recommended for most preprocessing methods and recipes. ``pp.recipe_monocle,``
    ``pp.recipe_velocyto`` has been deprecated (`PR 497 <https://github.com/aristoteleo/dynamo-release/pull/497>`_
    `PR 500 <https://github.com/aristoteleo/dynamo-release/pull/500>`_).
    Check the tutorials here for more instructions.
  - Normalization has been refactored (`PR 474 <https://github.com/aristoteleo/dynamo-release/pull/474>`_
    `PR 475 <https://github.com/aristoteleo/dynamo-release/pull/475>`_): ``pp.normalize_cell_expr_by_size_factors``
    has been deprecated, and new APIs are:

    - ``pp.normalize_cell_expr_by_size_factors`` -> ``pp.calc_sz_factor, pp.normalize``.

  - Gene selection has been refactored (`PR 474 <https://github.com/aristoteleo/dynamo-release/pull/474>`_). Now support
    genes selected by fano factors. APIs are ``pp.select_genes_monocle`` and ``pp.select_genes_by_seurat_recipe``.
  - PCA has been refactored (`PR 469 <https://github.com/aristoteleo/dynamo-release/pull/469>`_). ``dyn.pp.pca_monocle``
    has been deprecated. The new API is:

    - ``pp.pca_monocle`` -> ``pp.pca``.

  - Multiple new features added, includes genes selection by fano factors
    (`PR 474 <https://github.com/aristoteleo/dynamo-release/pull/474>`_), external data integration methods
    (`PR 473 <https://github.com/aristoteleo/dynamo-release/pull/473>`_) and ``pp.regress_out``
    (`PR 470 <https://github.com/aristoteleo/dynamo-release/pull/470>`_
    `PR 483 <https://github.com/aristoteleo/dynamo-release/pull/483>`_
    `PR 484 <https://github.com/aristoteleo/dynamo-release/pull/484>`_).
  - Created more tests for preprocessing module (`PR 485 <https://github.com/aristoteleo/dynamo-release/pull/485>`_).
  - Other deprecated APIs include: ``pp.calc_sz_factor_legacy, pp.filter_cells_legacy``,
    ``pp.filter_genes_by_outliers_legacy, pp.select_genes_monocle_legacy, pp.select_genes_by_dispersion_general``,
    ``pp.cook_dist, pp.normalize_cell_expr_by_size_factors``. More information can be found on our preprocessing
    tutorials.

- Debug:

  - Fixed the bug that save_show_or_return flags not working
    (`PR 414 <https://github.com/aristoteleo/dynamo-release/pull/414>`_).
  - Enabled the leiden algorithm to accept the resolution parameters
    (`PR 441 <https://github.com/aristoteleo/dynamo-release/pull/441>`_).
  - Fixed the wrong attribute name of anndata object in `utils_dimensionReduction.py`
    (`PR 458 <https://github.com/aristoteleo/dynamo-release/pull/458>`_)`
  - Fixed the dimensionality issue in `moments.py`
    (`PR 461 <https://github.com/aristoteleo/dynamo-release/pull/461>`_).
  - Fixed part of the bug that h5ad file cannot be saved correctly
    (`PR 467 <https://github.com/aristoteleo/dynamo-release/pull/467>`_).
  - Fixed the bug that `pca_mean` will be `None` under some circumstances
    (`PR 482 <https://github.com/aristoteleo/dynamo-release/pull/482>`_).
  - Removing warning message for nxviz
    (`PR 489 <https://github.com/aristoteleo/dynamo-release/pull/489>`_).
  - Corrected the norm log-likelihood function
    (`PR 495 <https://github.com/aristoteleo/dynamo-release/pull/495>`_).
  - Removed deprecated parameters in gseapy functions
    (`PR 496 <https://github.com/aristoteleo/dynamo-release/pull/496>`_).
  - Fixed the bugs that functions will raise error when no fixed points are found in vector field by sampling
    (`PR 501 <https://github.com/aristoteleo/dynamo-release/pull/501>`_).
  - Removed unwanted operations in dimension reduction
    (`PR 502 <https://github.com/aristoteleo/dynamo-release/pull/502>`_).


Tutorial Updates on Readthedocs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Documentation, Tutorials, and readthedocs update:

  - Update requirements for readthedocs (`PR 466 <https://github.com/aristoteleo/dynamo-release/pull/466>`_).
  - Update readme (`PR 479 <https://github.com/aristoteleo/dynamo-release/pull/479>`_).
  - Fixed documentation error caused by importing Literal
    (`PR 486 <https://github.com/aristoteleo/dynamo-release/pull/486>`_).
  - Fixed readthedocs error caused by the new version of urllib3
    (`PR 488 <https://github.com/aristoteleo/dynamo-release/pull/488>`_).


Other Changes
~~~~~~~~~~~~~~~~~~~~~~~~~~
- Docstring and type hints update:

  - Updated docstring and type hints for tools module
    (`PR 419 <https://github.com/aristoteleo/dynamo-release/pull/419>`_).
  - Updated docstring and type hints for vector field module
    (`PR 434 <https://github.com/aristoteleo/dynamo-release/pull/434>`_).
  - Updated the docstring and type hints for simulation and predicting module
    (`PR 457 <https://github.com/aristoteleo/dynamo-release/pull/457>`_).
  - Update the docstring and type hints for hzplot
    (`PR 456 <https://github.com/aristoteleo/dynamo-release/pull/456>`_).



Dynamo Ver 1.1.0
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
Feature Changes
~~~~~~~~~~~~~~~~~~~~~~~~~~
- Following new function are added, exported or documented in API / class page: 
  
  - *Preprocessing*: ``pp.convert2symbol, pp.filter_cells, pp.filter_gene,`` 
    ``pp.filter_genes_by_pattern, pp.normalize_cells, pp.scale, pp.log1p, pp.pca``
  - *Kinetic parameters and RNA/protein velocity*: ``tl.recipe_deg_data, tl.recipe_kin_data,``
    ``tl.recipe_mix_kin_deg_data, tl.recipe_one_shot_data, tl.velocity_N``
  - *Labeling Velocity recipes*: ``tl.infomap, tl.leiden, tl.louvain, tl.scc``
  - *Clustering*: ``tl.run_scvelo, tl.run_velocyto, tl.vlm_to_adata``
  - *Converter and helper*: ``vf.graphize_vecfld, vf.vector_field_function``
  - *Vector field reconstruction*: ``vf.FixedPoints, vf.VectorField2D, vf.assign_fixedpoints``
  - *Beyond RNA velocity*: ``vf.jacobian, vf.sensitivity``
  - *Vector field ranking*: ``vf.rank_cells, vf.rank_genes, vf.rank_expression_genes,``
    ``vf.rank_jacobian_genes, vf.rank_s_divergence_genes, vf.rank_sensitivity_genes``
  - *Vector field clustering and graph*: ``vf.cluster_field, vf.streamline_clusters``
  - *Prediction* ``pd.andecestor, pd.get_init_path, pd.least_action, pd.perturbation,``
    ``pd.rank_perturbation_cell_clusters, pd.rank_perturbation_cells, pd.rank_perturbation_genes,``
    ``pd.state_graph, pd.tree_model``
  - *Preprocessing plot*: ``pl.biplot, pl.loading, pl.highest_frac_genes, pl.bubble``
  - *Space plot*: ``pl.space``
  - *Kinetics plot*: ``pl.sensitivity_kinetics``
  - *Vector field plots*: ``pl.cell_wise_vectors_3d, pl.plot_fixed_points_2d``
  - *differential geometry plots*: ``pl.acceleration``
  - *Regulatory network plots* ``pl.arcPlot, pl.circosPlot, pl.circosPlotDeprecated, pl.hivePlot``
  - *fate plots* ``pl.fate``
  - *heatmap plots* ``pl.causality, pl.comb_logic, pl.plot_hill_function, pl.response``
  - *Predictions plots* ``pl.lap_min_time``
  - *External functionality* ``ext.normalize_layers_pearson_residuals,``
    ``ext.select_genes_by_pearson_residuals, ext.sctransform``

- More differential geometry analyses

  - include the `switch` mode in rank_jacobian_genes
  - added calculation of `sensitivity` matrix and relevant ranking 

- most probable path and *in silico* perturbation prediction

  - implemented least action path optimization (can be done in high dimensional space) with analytical Jacobian 
  - include genetic perturbation prediction by either changing the vector field function or simulate genetic perturbation via analytical Jacobian

- preprocessor class implementation

  - extensible modular preprocess steps 
  - support following recipes: monocle (dynamo), seurat (seurat V3 flavor), sctransform (seurat), pearson residuals and pearson residuals for feature selection, combined with monocle recipe (ensure no negative values)
  -  following recipes tested on zebrafish dataset to make implemetation results consistent:
    - monocle, seurat, pearson residuals
- CDlib integration

  - leiden, louvain, infomap community detection for cell clustering 
  - wrappers in ``dyn.tl.*`` for computing clusters
  - wrappers in ``dyn.pl.*`` for plotting


Tutorial Updates on Readthedocs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* human HSC hematopoiesis RNA velocity analysis tutorials
* *in silico* perturbation and least action path (LAP) predictions tutorials on HSC dataset
- differential geometry analysis on HSC dataset

  - Molecular mechanism of megakaryocytes
  - Minimal network for basophil lineage commitment
  - Cell-wise analyses: dominant interactions
* gallery: Pancreatic endocrinogenesis differential geometry


Sample Dataset Updates
~~~~~~~~~~~~~~~~~~~~~~~~~~


CI/CD Updates
~~~~~~~~~~~~~~~~~~~~~~~~~~
- update dynamo testing and pytest structure
- test building workflow on 3.7, 3.8, 3.9 (3.6 no longer tested on github building CI)


Performance Improvements
~~~~~~~~~~~~~~~~~~~~~~~~~~


API Changes
~~~~~~~~~~~~~~~~~~~~~~~~~~
- preprocess

 - ``pp.pca`` -> ``pca.pca_monocle``
* Native implementation of various graphical calculus using Numpy without using igraph. 


Other Changes
~~~~~~~~~~~~~~~~~~~~~~~~~~
* **general code refactor and bug fixing**
* **pl.scatters** refactor

