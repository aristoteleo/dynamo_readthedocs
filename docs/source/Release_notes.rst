==================
Release notes
==================


Dynamo Ver 1.1.0
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
Feature Changes
~~~~~~~~~~~~~~~~~~~~~~~~~~
- Following new function are added, exported or documented in API / class page: 
  - *Preprocessing*: ``pp.convert2symbol, pp.filter_cells, pp.filter_gene, pp.filter_genes_by_pattern, pp.normalize_cells, pp.scale, pp.log1p, pp.pca``
  - *Kinetic parameters and RNA/protein velocity*: ``tl.recipe_deg_data, tl.recipe_kin_data, tl.recipe_mix_kin_deg_data, tl.recipe_one_shot_data, tl.velocity_N``
  - *Labeling Velocity recipes*: ``tl.infomap, tl.leiden, tl.louvain, tl.scc``
  - *Clustering*: ``tl.run_scvelo, tl.run_velocyto, tl.vlm_to_adata``
  - *Converter and helper*: ``vf.graphize_vecfld, vf.vector_field_function``
  - *Vector field reconstruction*: ``vf.FixedPoints, vf.VectorField2D, vf.assign_fixedpoints``
  - *Beyond RNA velocity*: ``vf.jacobian, vf.sensitivity``
  - *Vector field ranking*: ``vf.rank_cells, vf.rank_genes, vf.rank_expression_genes, vf.rank_jacobian_genes, vf.rank_s_divergence_genes, vf.rank_sensitivity_genes``
  - *Vector field clustering and graph*: ``vf.cluster_field, vf.streamline_clusters``
  - *Prediction* ``pd.andecestor, pd.get_init_path, pd.least_action, pd.perturbation, pd.rank_perturbation_cell_clusters, pd.rank_perturbation_cells, pd.rank_perturbation_genes, pd.state_graph, pd.tree_model``
  -  *Preprocessing plot*: ``pl.biplot, pl.loading, pl.highest_frac_genes, pl.bubble``
  -  *Space plot*: ``pl.space``
  - *Kinetics plot*: ``pl.sensitivity_kinetics``
  - *Vector field plots*: ``pl.cell_wise_vectors_3d, pl.plot_fixed_points_2d``
  - *differential geometry plots*: ``pl.acceleration``
  - *Regulatory network plots* ``pl.arcPlot, pl.circosPlot, pl.circosPlotDeprecated, pl.hivePlot``
  - *fate plots* ``pl.fate``
  - *heatmap plots* ``pl.causality, pl.comb_logic, pl.plot_hill_function, pl.response``
  - *Predictions plots* ``pl.lap_min_time``
  - *External functionality* ``ext.normalize_layers_pearson_residuals, ext.select_genes_by_pearson_residuals, ext.sctransform``

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

