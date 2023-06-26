==================
Dynamics
==================

Dynamics combines parameters estimation and velocity calculation, establishing the foundation for vector field
reconstruction and subsequent analyses. This paragraph provides a concise overview of the diverse approaches we have
implemented to cater to different experiment types. More details can be found in the "Method" section in our
`paper <https://www.sciencedirect.com/science/article/pii/S0092867421015774#sec5>`_.

Spliced data:
~~~~~~~~~~~~~~~~~~~~~~~~~~
When the dataset only contains spliced and unspliced RNA data, the following methods can be applied.

-   **Conventional**: The kinetics of RNA is transcription, splicing, and degradation is represented as follows:

    .. math::
       \dot{u} = \alpha - \beta u \\
       \dot{s} = \beta u - \gamma s

    where :math:`u` and :math:`s` are the copies of unspliced and spliced RNA for a particular gene in a cell,
    :math:`\alpha`, :math:`\beta` and :math:`\gamma` are the rate constants for transcription, splicing, and degradation
    rate. Transcription rate is a function of the cell state and other variables while splicing and degradation rate are
    often regarded as constants for specific cell types in most cases. Assuming pseudo-steady state for cells with
    extreme high unspliced and spliced RNA expressions, we can perform linear regression between the spliced and
    unspliced RNA to estimate the parameters (:math:`\beta u = \gamma s`). Then the velocity can be determined by
    applying the corresponding ordinary differential equations (ODEs) in conjunction with the obtained parameter values.
    In the original paper the velocity is defined as :math:`v = u - \tilde{\gamma} s`,
    where :math:`\tilde{\gamma} = \frac{\gamma}{\beta}`.

    The parameters configuration for this method is as follows: `experiment_type` -> `conventional`, `assumption_mRNA`
    -> `ss`, `model` -> `deterministic`. Additionally, we offer multiple options for the `est_method` parameter,
    including `ols`, `rlm`, and `ransac`. For a comprehensive understanding of the est_method options and their
    respective functionalities, please refer to the detailed documentation provided within the code's docstring.

-   **stochastic splicing**: The stochastic splicing method, denoted as GMM (Generalized method of moments) in the code,
    enhances the estimation of kinetic parameters by incorporating the first, second, and mixed moments of unspliced and
    spliced RNAs for each gene across multiple cells. In contrast to the original velocity method, which only relies on
    the first moment, this approach takes advantage of the additional moments, thereby augmenting the robustness of
    parameter estimation. The parameters configuration for this method is as follows: `experiment_type` -> `conventional`,
    `model` -> `stochastic`, `assumption_mRNA` -> `ss`, `est_method` -> `gmm`.

-   **Negative binomial**: The negative binomial is an alternative procedure based on the observation that in most cases
    total RNA counts at steady state follow the NB distribution. It reformulates the moments condition with the
    NB-distribution variable and solves the formulated equation with a nonlinear least squares optimizer. The
    configuration for this method is as follows: `experiment_type` -> `conventional`, `model` -> `stochastic`,
    `assumption_mRNA` -> `ss`, `est_method` -> `negbin`.


Labeled data:
~~~~~~~~~~~~~~~~~~~~~~~~~~
The metabolic labeling dataset measures the synthesis or degradation of labeled RNA within a known period of time in an
experimentally programmable manner. This dataset offers a more direct assessment of the kinetics underlying gene
expression. Consequently, leveraging the metabolic labeling dataset presents an opportunity to address certain
challenges encountered in velocity estimation. Here we introduce multiple dynamics models tailored to different types
of experiments.

-   **One shot**: In “one-shot” experiments, there is only one labeling time point, and the splicing process is not
    explicitly considered. In this method, we substitute the spliced and unspliced RNA variables from the original
    velocity method with labeled and total RNAs, then we can incorporate the second moments using the negative binomial
    method to get a more accurate parameter. The configuration for this method is as follows: `experiment_type` ->
    `one_shot`/`one-shot`, `assumption_mRNA` -> `ss`. For one shot experiment, we can tune the model and est_method just
    like conventional methods.

-   **Kinetic**: Kinetic experiments represent a time series of 4sU or other nucleotide analog treatment to observe the
    accumulation of metabolically labeled RNA over time. Our implementation includes two distinct estimation methods:
    the two-step approach and the direct approach. The two-step approach is a generalization of one-shot method on
    multiple time points, which relies on two consecutive linear regressions to estimate the degradation rate.
    Furthermore, the two-step approach can be further generalized to the kinetic assumption. On the other hand, the
    direct approach directly employs the kinetic model to estimate the rate parameters. In terms of velocity
    calculation, the two-step approach generally outperforms the direct approach. The configuration for this method is
    as follows: `experiment_type` -> `kin`, `assumption_mRNA` -> `ss`/`kin`, `model` -> `deterministic`/`stochastic`,
    `est_method` -> `two-step`/`direct`.

-   **Degradation**: In degradation experiments, samples are chased after an extended 4sU (or other nucleotide analog)
    labeling period and the wash-out to observe the decay of the abundance of the (labeled) unspliced and spliced RNA
    decay over time. Splicing and degradation rate can be estimated from this process using the nonlinear least squares.
    The configuration for this method is as follows: `experiment_type` -> `deg`, `assumption_mRNA` -> `kinetic`, `model`
    -> `deterministic`/`stochastic`/`mixture`.

-   **Mix**: Furthermore, our framework offers support for a range of mixture experiments like mixed steady-state and
    stimulation labeling experiment (`mix_std_stm`) and mixed kinetic and degradation experiment
    (`mix_kin_deg`/`mix_pulse_chase`). More details can be found in the code.