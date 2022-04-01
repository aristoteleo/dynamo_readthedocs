==================
Release notes
==================


Dynamo Ver 1.1.0
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
Feature Changes
~~~~~~~~~~~~~~~~~~~~~~~~~~
- preprocessor class implementation

  - extensible modular preprocess steps 
  - support following recipes: monocle (dynamo), seurat (seurat V3 flavor), sctransform (seurat), pearson residuals and monocle recipe combined with pearson residuals to prevent negative values
  -  following recipes tested on zebrafish dataset to make implemetation results consistent:
    - monocle, seurat, pearson residuals
- CDlib integration

  - leiden, louvain, infomap community detection
  - wrappers in ``dyn.tl.*`` for computing clusters
  - wrappers in ``dyn.pl.*`` for plotting


Tutorial Updates on Readthedocs
~~~~~~~~~~~~~~~~~~~~~~~~~~
* human HSC raw data process tutorials
* gene perturbation analysis and least action path (LAP) tutorials on HSC dataset
- analysis application tutorials on HSC dataset

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
* new modular APIs in graph_operators authored by Yan Zhang to replace vector field submodule's outdated implementation


Other Changes
~~~~~~~~~~~~~~~~~~~~~~~~~~
* **general code refactor and bug fixing**
* **pl.scatters** refactor


Contributors in lexicographical order
---------------------------------------
* Ke
* Xiaojie, Qiu
* Yan, Zhang