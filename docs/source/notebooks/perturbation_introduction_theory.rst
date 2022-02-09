.. _perturbation_theory_tutorial:

*In-silico* perturbation
=========================================================================================

:math:`\newcommand{\pdv}[2]{\dfrac{\partial #1}{\partial #2}} \newcommand{\trp}{\mathsf{T}}`

We leverage the analytical Jacobian of the reconstructed vector field
function to make :math:`\textit{in silico}` genetic perturbation
and predict cell-fate outcomes after the perturbation, showcased in :ref:`this figure<dynamo_fig7_a>`. 

.. _dynamo_fig7_a:
.. figure:: dynamo_paper_figures/fig7_a.png
    :align: center

    In silico genetic perturbation of the velocity vector field. i) In silico perturbation can predict the gene-wise response. ii) In silico perturbation can predict the cell fate trajectory after genetic perturbation by integrating the displacement of velocities across cells.

Mathematically, for gene :math:`i` in any cell, the genetic perturbation effects or
change in its velocity (or more accurately, the vector field) w.r.t. to
small perturbations in the expression of all genes in the network
(encoded by the Jacobian matrix :math:`\boldsymbol J`),
:math:`\mathrm dx_1`, :math:`\mathrm dx_2`,…, :math:`\mathrm dx_n`, can
be calculated with the :math:`\textit{exact differential}`:

.. math::
    \begin{align*}  \mathrm d f_i = \pdv{f_i}{x_1}\mathrm dx_1 + \pdv{f_i}{x_2}\mathrm dx_2 + ... + \pdv{f_i}{x_n}\mathrm dx_n. \end{align*}

In vectorized form:

.. math::
    \begin{align*}  \begin{bmatrix}  \mathrm df_1 \\[1.5ex] \mathrm df_2 \\[1.5ex] \dots \\[1.5ex] \mathrm df_n  \end{bmatrix} =  \begin{bmatrix}  \pdv{f_1}{x_1} \ &\pdv{f_1}{x_2} \ &\dots \ &\pdv{f_1}{x_n} \\[2ex]  \pdv{f_2}{x_1} \ &\pdv{f_2}{x_2} \ &\dots \ &\pdv{f_2}{x_n} \\[2ex]  \dots \ &\dots \ &\dots \ &\dots \\[2ex]  \pdv{f_n}{x_1} \ &\pdv{f_n}{x_2} \ &\dots \ &\pdv{f_n}{x_n}  \end{bmatrix}  \begin{bmatrix}  \mathrm dx_1 \\[1.5ex] \mathrm dx_2 \\[1.5ex] \dots \\[1.5ex] \mathrm dx_n  \end{bmatrix}. \end{align*}

The matrix on the right hand side is the Jacobian of the vector field.
Replacing infinitesimal changes with finite perturbations, the above
equation becomes:

.. math::
    \begin{align*}  \Delta \boldsymbol f = \boldsymbol J \Delta \boldsymbol x. \end{align*}


In practice, a proportionality constant :math:`c` is often added to the
perturbation :math:`\Delta \boldsymbol x` to amplify the response
:math:`\Delta \boldsymbol f`. Furthermore, because vector fields are
often learned in the PCA space, the perturbations in the
:math:`d`-dimensional gene space are first transformed to the
:math:`k`-dimensional PCA space by:
:math:`\begin{align*}  \Delta \boldsymbol x = \boldsymbol Q^\trp (\Delta \boldsymbol y - \boldsymbol \mu). \end{align*}`
where :math:`\boldsymbol Q` is the :math:`d`-by-:math:`k` PCA loading
matrix, and :math:`\boldsymbol \mu` is the mean of the PCA-transformed
data. The response :math:`\Delta \boldsymbol f` can be transformed back
to the PCA space:
:math:`\begin{align*}  \Delta \boldsymbol g = \boldsymbol Q \Delta \boldsymbol f + \boldsymbol \mu. \end{align*}`
One can then use :math:`\Delta \boldsymbol f`, a gene by cell matrix, to
identify the strongest positive or negative responders of the genetic
perturbation across cells.

Importantly, because :math:`\Delta \boldsymbol f` implies how each cell
state will be affected after genetic perturbations, we can predict the
cell fate trajectory under genetic perturbations by integrating the
perturbation effects across cells over gene expression space, To
visualize the cell fate trajectory, pairs of :math:`\boldsymbol x` and
:math:`\Delta \boldsymbol g` are used in the same vein as the gene
expression and RNA velocity vector to be further projected onto the UMAP
or other low dimensional embeddings using the transition matrix
:cite:p:`Bergen2020-kx, La_Manno2018-vp` and then plotted with
streamlines.