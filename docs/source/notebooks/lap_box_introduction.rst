.. _lap_theory_tutorial:

Optimal cell fate transitions via most probably path
====================================================


The least action path (LAP) principle, first proposed as early as 1744
by Maupertuis :cite:t:`terrall` and famously advocated by Feynman with
his reformulation of quantum mechanics via the path integral of the
classical Hamilton action :cite:t:`fey65`, has previously been used
in predicting the optimal transition path of cell fate transition for
simplistic and designed systems
:cite:t:`Qiu2012-yt, Wang2014-zc, Wang8257`. We reason that with the
reconstructed continuous, differentiable vector field, we can extend the
LAP approach to real datasets in transcriptomic space to computationally
explore optimal paths for differentiation and reprogramming
(dedifferentiation and transdifferentiation), which then helps identify
key transcription factors whose expression levels vary strongest along
these paths. The grand problem of prediction OPCs is summarized :ref:`here<lap_theory_dynamo_paper_fig6_a>`. 

.. _lap_theory_dynamo_paper_fig6_a:
.. figure:: dynamo_paper_figures/fig6_a.png
   :align: center
   :width: 400

   The grand problem of predicting OPtimal cell-fate Conversions(OPCs).


The transcriptomic vector field encodes dynamical information of
pathways connecting different cell types. :ref:`This figure from dynamo paper<lap_theory_dynamo_paper_fig6_b>` :cite:p:`QIU2022` shows developmental path and transition matrix among cell types in hematopoietic stem cell lineages. For the development, the
developmental paths, connecting progenitors and stable cell types, such
as HSCs and megakaryocytes, are
characterized by vector field streamlines, where cells need to overcome
little to no dynamical barrier. On the contrary, the reversed process,
dedifferentiation, requires cells to migrate against the streamline,
overcome the developmental barrier to eventually become progenitors and
reacquire multipotency. Stable cell types are attractors that are
separated by attractor barriers, and during transdifferentiation
processes, cells from one stable attractor overcome these barriers and
transverse into another stable attractor. Overcoming the
dedifferentiation/transdifferentiation barriers driven by stochasticity
are :math:`\textit{rare transitions}`, as evidenced by extremely
low experimental reprogramming efficiency, as low as 0.001–0.01%
:cite:t:`Merkl2013`. This is why reprogramming factors are
generally needed, which reshape the developmental landscape in favor of
the reprogramming transitions. In the context of cell state transitions,
there are two seemingly similar but fundamentally different concepts
worth additional clarification:

-  transition time: The expected waiting time for a cell to initiate and
   finish the transition between two states, regardless of the path it
   takes. This corresponds to the experimentally measured time for one
   cell type to commit into another.

-  traversal time: The time the cell spends traveling along a specific
   path. Theoretically, this is the time for a single cell to complete
   the cell type conversion once the cell has decided on the commitment.

.. _lap_theory_dynamo_paper_fig6_b:
.. figure:: dynamo_paper_figures/fig6_b.png

  Predicting OPCs for hematopoietic cell types. i) The developmental tree, known dedifferentiation and transdifferentiation events previously reported for the six cell types observed in our data. ii) Matrix representation of subpanel i). The optimal paths for hematopoietic transitions can be found by identifying the LAPs between the fixed points
  that correspond to each stable cell type.


:ref:`An example of LAP transition from state0 (HSC) to state1 (Meg)<lap_theory_dynamo_paper_fig6_c>` More explanation regarding math terms shown in this figure can be found in this introduction below.

.. _lap_theory_dynamo_paper_fig6_c:
.. figure:: dynamo_paper_figures/fig6_c.png
   :align: center
   :width: 400

   The optimal paths for hematopoietic transitions can be found by identifying the LAPs between the fixed points
that correspond to each stable cell type.



Given the vector field function, :math:`\boldsymbol f`, optimal pathways
of cell fate conversion can be mathematically analyzed by least action
paths (LAPs)
:cite:t:`freidlin2012random, onsager1953, Maier1997`. The
action is defined as:

.. math::
  \begin{align*}  \newcommand{\trp}{\mathsf{T}}  S_T(\boldsymbol x) = \frac{1}{2} \int_{0}^{T}\mathrm dt {\Big(\boldsymbol v(t) - \boldsymbol f\big(\boldsymbol x(t)\big)\Big)}^\trp \boldsymbol{D}^{-1}  \Big(\boldsymbol v(t) - \boldsymbol f\big(\boldsymbol x(t)\big)\Big), \end{align*}`

:math:`\boldsymbol x` is a path and :math:`\boldsymbol v` its tangential
velocity (the path is parametrized by time :math:`t`, so
:math:`\boldsymbol v(t) = \dot{\boldsymbol x}(t)`).
:math:`\boldsymbol{D}` is the diffusion coefficient matrix accounting
for the stochasticity of gene expression, and for simplicity here we
assume it to be a constant. :math:`T` is the time needed for a cell to
traverse the path. By this definition, a path that strictly follows a
streamline of a vector field whose tangential velocity also equals the
evaluated velocity of the vector field has zero action, whereas any
deviation increases action. In other words, developmental processes are
(mostly) a spontaneous process and driven by intrinsic cell states,
whereas dedifferentiation requires external forces such as ectopic
expression of exogenous TFs or specific chemical inductions.

Computationally, given the starting and end cell states
:math:`\boldsymbol x_0` and :math:`\boldsymbol x_{n}`, such as HSCs and
megakaryocytes, and a specific traversal
time :math:`T`, the LAP can be found by discretizing the path as a
sequence of points
:math:`P=\{\boldsymbol x_0, \boldsymbol x_1, \dots, \boldsymbol x_n\}`,
which forms :math:`n` line segments. For each line segment, the discrete
tangential velocity can be calculated as
:math:`\boldsymbol v_k = (\boldsymbol x_k-\boldsymbol x_{k-1})/\Delta t`,
where :math:`\Delta t = T/n`. The action along the discrete path is
defined as :cite:t:`Perez-Carrasco2016, Tang2017`:

.. math:
  \begin{align*}  S_T(P) = \frac{1}{2D}\sum_{k=1}^{n} \Big(\boldsymbol v_k - \boldsymbol f(\boldsymbol y_k)\Big)^2\Delta t , \end{align*}

| where :math:`y_k` are the middle points of the line segments, i.e.,
  :math:`\boldsymbol y_k = (\boldsymbol x_{k-1} + \boldsymbol x_k)/2`.
  Given a traversal time :math:`T`, the LAP is a path such that:

.. math::
    \begin{align*}  P^* = \underset{P}{\operatorname{argmin}}\ S_T(P) = \underset{P}{\operatorname{argmin}}\ \frac{1}{2D}\sum_{k=1}^{n} \Big(\boldsymbol v_k - \boldsymbol f(\boldsymbol y_k)\Big)^2\Delta t . \end{align*} 
| To obtain the global LAP, the optimal traversal time :math:`T^*` is
  determined as:

.. math::
  \begin{align*}  T^* = \underset{T}{\operatorname{argmin}}\ S_T(P) \end{align*}

The algorithm discretizes the path as a sequence of points,
:math:`P=\{\boldsymbol x_0, \boldsymbol x_1, \dots, \boldsymbol x_n\}`,
which forms :math:`n` line segments. For each line segment, the discrete
tangential velocity can be calculated as
:math:`\boldsymbol v_k=(\boldsymbol x_k - \boldsymbol x_{k-1})/\Delta t`,
where :math:`\Delta t` is the time step for the cell to move from
:math:`\boldsymbol x_{k-1}`. In addition to the deterministic vector
field, we also assume a certain degree of stochasticity in the system:

.. math::
    \begin{align*}  \dot{\boldsymbol x} = \boldsymbol f(\boldsymbol x) + \sigma \boldsymbol\eta(t), \end{align*}

| where :math:`\boldsymbol\eta(t)` is a stochastic white noise and
  :math:`\boldsymbol\sigma` the size of it. The action :math:`S` along
  the discrete path is defined as (Perez-Carrasco et al., 2016):

.. math::
  \begin{align*}  S(P, \Delta t) = \frac{1}{2D}\sum_{k=1}^{n}\Big(\boldsymbol v_k - \boldsymbol f(\boldsymbol y_k)\Big)^2\Delta t, \end{align*}

| where :math:`\boldsymbol y_k` are the middle points of the line
  segments, i.e.,
  :math:`\boldsymbol y_k = (\boldsymbol x_{k-1} + \boldsymbol x_k)/2`.
  We have also assumed the diffusion matrix to be a constant :math:`D`,
  such that :math:`D=\sigma^2/2`. It is intuitive that a path whose
  tangential velocities :math:`\boldsymbol v` align with the vector
  field has smaller action than paths that do not. The LAP is a path
  such that:

.. math::
  \begin{align*}  P^* = \underset{P, \Delta t}{\operatorname{argmin}} S(P, \Delta t) = \underset{P, \Delta t}{\operatorname{argmin}}\frac{1}{2D}\sum_{k=1}^{n}\Big(\boldsymbol v_k - \boldsymbol f(\boldsymbol y_k)\Big)^2\Delta t, \end{align*}

| The algorithm for finding the LAP therefore consists of two steps:

-  Minimization of the action by varying the time step. The optimal time
   step given a fixed path is a simple univariate least square
   minimization, i.e.:

.. math::
  \begin{align*}  \Delta t^* = \underset{\Delta t}{\operatorname{argmin}}\frac{1}{2D}\sum_{k=1}^{n}\Big(\frac{\boldsymbol x_k - \boldsymbol x_{k-1}}{\Delta t} - \boldsymbol f(\boldsymbol y_k)\Big)^2\Delta t,  \end{align*}

-  Minimization of the action by varying the path without moving the
   starting and end points. The optimal path given a fixed time step is
   found by:

.. math::
  \begin{align*}  P^* = \underset{\{\boldsymbol x_1, \boldsymbol x_2, \dots, \boldsymbol x_{n-1}\}}{\operatorname{argmin}}\frac{1}{2D}\sum_{k=1}^{n}\Big(\frac{\boldsymbol x_k - \boldsymbol x_{k-1}}{\Delta t} - \boldsymbol f\big(\frac{\boldsymbol x_{k-1} + \boldsymbol x_k}{2}\big)\Big)^2\Delta t, \end{align*}

For a :math:`d`-dimensional vector field, the number of variables in
the above optimization problem is :math:`d\times n`. To mitigate the
computational cost, the Jacobian of the action w.r.t. the path (more
specifically, the a-th component of the :math:`k`-th point) is
analytically computed:

.. math::
  \begin{align*} \frac{\partial{S}}{\partial{x_k^a}} =& \frac{1}{D}\Big(v_k^a - v_{k+1}^a + f^a(\boldsymbol y_{k+1}) - f^a(\boldsymbol y_k)\Big)\\  &-\frac{1}{2D}\Big(\big(\boldsymbol v_{k+1} - \boldsymbol f(\boldsymbol x_{k+1})\big) \cdot \frac{\partial{f}}{\partial{x^a}}\Big|_{\boldsymbol x_{k+1}} + \big(\boldsymbol v_k - \boldsymbol f(\boldsymbol x_k)\big)\cdot\frac{\partial f}{\partial{x^a}}\Big|_{\boldsymbol x_k}\Big)  \end{align*}

| Note that the partial derivative of the vector field is the
  :math:`a`-th row of the Jacobian of the vector field. With the
  analytical Jacobian, the computation efficiency of the LAP
  optimization improves tremendously, making the LAP calculation
  feasible to operate in high-dimensional space, such as the top 30 PCs.

The LAP is found by iterating between the two steps, and empirically we
found that the path converges in two or three iterations. By default,
the LAP optimization is initialized with the interpolated shortest path
on the kNN graph of cells.

Notably, when LAPs are calculated in the PCA space, we can transform
them back to the original gene expression space to predict the full
transcriptomic kinetics along the optimal path, inspect waves of those
kinetics along the path, and do so in absolute time units when the
vector field used is based on tscRNA-seq.

For rare transitions with :math:`S_{T^*} \gg 0` (e.g., dedifferentiation
and transdifferentiation), the transition rate (number of transitions
per unit time) is proportional to the exponential of actions of all
paths. The Freidlin–Wentzell theorem dictates that the LAP with the
minimal traversal time (which will be referred to as the optimal path
below) contributes the most to this transition rate
:cite:t:`freidlin2012random, onsager1953, Maier1997, Aurell2002`:

.. math::
  \begin{align*}  R(A\rightarrow B) \approx C\exp(-S_{T^*}), \end{align*}

| where :math:`A` and :math:`B` are two cell types, :math:`S_{T^*}` the
  action of the optimal path, and :math:`C` a proportional factor.
  Furthermore, the transition time, or more specifically the mean first
  passage time (MFPT), is related to the transition rate:

.. math::
  \begin{align*}  \mathrm{MFPT} = \frac{1}{R(A\rightarrow B)} \end{align*}

| Therefore, the action of the optimal path predicts both the likelihood
  and transition time for such rare transitions. Again, most
  reprogramming experiments take a few weeks or months, depending on the
  exact initial and terminal cell states
  :cite:t:`takahashi2006induction`.

For natural transitions between points that are connected by the vector
field streamlines (e.g., from a repulsor to an adjacent attractor), the
actions of LAPs, within a certain range of :math:`T`, are all zero,
because a path following the streamline downstream is a LAP with zero
action. The above approximation that the LAP contributes the most to the
transition rate no longer applies. Differentiation processes are often
close to such natural transitions, and the action of a differentiation
LAP cannot tell us any information on the transition rate. However, LAPs
are still the most probable paths for cells to take, as they are
optimized to follow the streamline of the vector field. The waiting time
for the cell to initiate the transition is negligible in this case, so
the transition time can be approximated by the traversal time of the
LAP.

In addition to the computation of transition time and traversal time,
analyzing gene expression variations along LAPs provides essential
information on regulatory genes, and their dynamics, during cell fate
transitions. We calculate the mean squared displacement (MSD) for every
gene :math:`i` along the optimal path:

.. math::
  \begin{align*}  \mathrm{MSD}_i = \sum_{t=0}^{T} \big(y_i(t) - y_i(0)\big)^2 \end{align*}

| Genes with large MSD are potentially genes that regulate the
  corresponding transitions.