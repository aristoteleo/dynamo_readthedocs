|upload| |conda| |download| |star| |build| |documentation|
|upload_python_package| |test|

**Dynamo**: Mapping Vector Field of Single Cells
===================================================

Inclusive model of expression dynamics with metabolic labeling based
scRNA-seq / multiomics, vector field reconstruction, potential landscape
mapping, differential geometry analyses, and most probably paths / *in
silico* perturbation predictions.

`Installation <https://dynamo-release.readthedocs.io/en/latest/ten_minutes_to_dynamo.html#how-to-install>`__
- `Ten minutes to
dynamo <https://dynamo-release.readthedocs.io/en/latest/ten_minutes_to_dynamo.html>`__
-
`Tutorials <https://dynamo-release.readthedocs.io/en/latest/notebooks/Differential_geometry.html>`__
- `API <https://dynamo-release.readthedocs.io/en/latest/API.html>`__ -
`Citation <https://www.sciencedirect.com/science/article/pii/S0092867421015774?via%3Dihub>`__
-
`Theory <https://dynamo-release.readthedocs.io/en/latest/notebooks/Primer.html>`__

.. image:: https://user-images.githubusercontent.com/7456281/152110270-7ee1b0ed-1205-495d-9d65-59c7984d2fa2.png
   :align: center

Single-cell (sc)RNA-seq, together with RNA velocity and metabolic
labeling, reveals cellular states and transitions at unprecedented
resolution. Fully exploiting these data, however, requires kinetic
models capable of unveiling governing regulatory functions. Here, we
introduce an analytical framework dynamo, which infers absolute RNA
velocity, reconstructs continuous vector fields that predict cell fates,
employs differential geometry to extract underlying regulations, and
ultimately predicts optimal reprogramming paths and perturbation
outcomes. We highlight dynamo’s power to overcome fundamental
limitations of conventional splicing-based RNA velocity analyses to
enable accurate velocity estimations on a metabolically labeled human
hematopoiesis scRNA-seq dataset. Furthermore, differential geometry
analyses reveal mechanisms driving early megakaryocyte appearance and
elucidate asymmetrical regulation within the PU.1-GATA1 circuit.
Leveraging the least-action-path method, dynamo accurately predicts
drivers of numerous hematopoietic transitions. Finally, in silico
perturbations predict cell-fate diversions induced by gene
perturbations. Dynamo, thus, represents an important step in advancing
quantitative and predictive theories of cell-state transitions.

.. |upload| image:: https://img.shields.io/pypi/v/dynamo-release?logo=PyPI
   :target: https://pypi.org/project/dynamo-release/
.. |conda| image:: https://img.shields.io/conda/vn/conda-forge/dynamo-release.svg
   :target: https://anaconda.org/conda-forge/dynamo-release
.. |download| image:: https://static.pepy.tech/badge/dynamo-release
   :target: https://pepy.tech/project/dynamo-release
.. |star| image:: https://img.shields.io/github/stars/aristoteleo/dynamo-release?logo=GitHub&color=red
   :target: https://github.com/aristoteleo/dynamo-release/stargazers
.. |build| image:: https://github.com/aristoteleo/dynamo-release/actions/workflows/python-package.yml/badge.svg
   :target: https://github.com/aristoteleo/dynamo-release/actions/workflows/python-package.yml
.. |documentation| image:: https://readthedocs.org/projects/dynamo-release/badge/?version=latest
   :target: https://dynamo-release.readthedocs.io/en/latest/
.. |upload_python_package| image:: https://github.com/aristoteleo/dynamo-release/actions/workflows/python-publish.yml/badge.svg
   :target: https://github.com/aristoteleo/dynamo-release/actions/workflows/python-publish.yml
.. |test| image:: https://github.com/aristoteleo/dynamo-release/actions/workflows/python-plain-run-test.yml/badge.svg
   :target: https://github.com/aristoteleo/dynamo-release/actions/workflows/python-plain-run-test.yml


Highlights of dynamo
====================
#. Robust and accurate estimation of RNA velocities for regular scRNA-seq datasets:
    * Three methods for the velocity estimations (including the new negative binomial distribution based approach)
    * Improved kernels for transition matrix calculation and velocity projection
    * Strategies to correct RNA velocity vectors (when your RNA velocity direction is problematic)
#. Inclusive modeling of time-resolved metabolic labeling based scRNA-seq:
    * Explicitly model RNA metabolic labeling, in conjunction with RNA bursting, transcription, splicing and degradation
    * Comprehensive RNA kinetic rate estimation for one-shot, pulse, chase and mixture metabolic labeling experiments
#. Move beyond RNA velocity to continuous vector field function for functional and predictive analyses of cell fate transitions:
    * Dynamical systems approaches to identify stable cell types (fixed points), boundaries of cell states (separatrices), etc
    * Calculate RNA acceleration (reveals early drivers), curvature (reveals master regulators of fate decision points), divergence (stability of cell states) and RNA Jacobian (cell-state dependent regulatory networks)
    * Various downstream differential geometry analyses to rank critical regulators/effectors,  and visualize regulatory networks at key fate decision points
#. Non-trivial vector field prediction of cell fate transitions:
    * Least action path approach to predict the optimal paths and transcriptomic factors of cell fate reprogrammings
    * In silico perturbation to predict the gene-wise perturbation effects and cell fate diversion after genetic perturbations


News
==========
.. _Cell: https://www.sciencedirect.com/science/article/pii/S0092867421015774#tbl1

-  5/30/2023: dynamo 1.3.0 released!
-  3/1/2023: We welcome @Sichao25 to join the dynamo develop team!
-  1/28/2023: We welcome @Ukyeon to join the dynamo develop team!
-  15/12/2022: *Thanks for @elfofmaxwell and @MukundhMurthy’s
   contribution*. dynamo 1.2.0 released
-  11/11/2022: the continuing development of dynamo and the Aristotle
   ecosystem will be supported by CZI. See
   `here <https://chanzuckerberg.com/eoss/proposals/predictive-modeling-of-single-cell-multiomics-over-time-and-space/>`__
-  4/14/2022: dynamo 1.1.0 released!
-  3/14/2022: Since today dynamo has its own logo! Here the arrow
   represents the RNA velocity vector field, while the helix the RNA
   molecule and the colored dots RNA metabolic labels (4sU labeling).
   See
   `readthedocs <https://dynamo-release.readthedocs.io/en/latest/index.html>`__
-  2/15/2022: primers and tutorials on least action paths and in silico
   perturbation are released.
-  2/1/2022: after 3.5+ years of perseverance, our dynamo paper is
   finally online in
   `Cell <https://www.sciencedirect.com/science/article/pii/S0092867421015774#tbl1>`__
   today!

Discussion
==========
Please use github issue tracker to report coding related `issues`_ of dynamo. For community discussion of novel usage cases, analysis tips and biological interpretations of dynamo, please join our public slack workspace: `dynamo-discussion`_ (Only a working email address is required from the slack side).

Contribution
============
If you want to contribute to the development of dynamo, please check out CONTRIBUTION instruction: `Contribution`_




.. toctree::
   :caption: Introduction
   :maxdepth: 1
   :hidden:

   notebooks/Introduction
   notebooks/Primer
   notebooks/lap_box_introduction
   notebooks/perturbation_introduction_theory.rst
   .. auto_examples/index.rst

.. toctree::
   :caption: Contents
   :maxdepth: 3
   :hidden:

   ten_minutes_to_dynamo
   Dynamics
   API
   Class
   FAQ
   Release_notes
   References
   Acknowledgement


.. toctree::
   :hidden:
   :caption: Preprocessing
   :maxdepth: 1


   Preprocessor_tutorial


.. toctree::
   :hidden:
   :caption: Conventional scRNA-seq
   :maxdepth: 1
   

   notebooks/zebrafish


.. toctree::
   :hidden:
   :caption: Labeling scRNA-seq
   :maxdepth: 3

   notebooks/tutorial_hsc_velocity
   notebooks/scNT_seq_readthedocs
   notebooks/scEU_seq_rpe1_analysis_kinetic
   notebooks/scEU_seq_organoid_analysis_kinetic


.. toctree::
   :caption: Differential geometry
   :maxdepth: 1
   :hidden:

   notebooks/hsc_differential_geometry_index
   notebooks/Differential_geometry.rst

.. toctree::
   :caption: Vector field predictions
   :maxdepth: 1
   :hidden:

   notebooks/lap_tutorial/lap_tutorial
   notebooks/perturbation_tutorial/perturbation_tutorial
   Shiny

.. toctree::
   :caption: Gallery
   :hidden:

   gallery/index.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _`dynamo`: https://github.com/aristoteleo/dynamo-release
.. _`issues`: (https://github.com/aristoteleo/dynamo-release/issues)
.. _`dynamo-discussion`: https://join.slack.com/t/dynamo-discussionhq/shared_invite/zt-ghve9pzp-r9oJ9hSQznWrDcx1fCog6g
.. _`Contribution`: https://github.com/aristoteleo/dynamo-release/blob/master/CONTRIBUTING.md

