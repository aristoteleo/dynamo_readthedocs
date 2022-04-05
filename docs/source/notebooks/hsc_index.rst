scNT-seq human hematopoiesis
===================================================
In our dynamo paper, we firstly highlight its power to overcome fundamental limitations of conventional splicing-based
RNA velocity analyses to enable accurate velocity estimations on a metabolically labeled human hematopoiesis scRNA-seq
dataset, which will be introduced first in this `labeling scRNA-seq` section. Under the `differential geometry` section,
we will show how dynamo's differential geometry analyses reveal mechanisms driving early megakaryocyte appearance, dual
origin of basophil lineage and elucidate asymmetrical regulation within the PU.1-GATA1 circuit. Lastly, under the
`vector field predictions` section, we will show we can leverage the least-action-path method to accurately predicts
drivers of numerous hematopoietic transitions. Finally, in this section, we will also show how *in silico* perturbations
predicts cell-fate diversions induced by gene perturbations. As you can, dynamo represents an important step in
advancing quantitative and predictive modeling of cell-state transitions. What is also worth to mention is that although
the demonstration here focuses on our labeling data, dynamo can be equally used to analyze conventional scRNA-seq (the
splicing data) and even used in conjunctions with other velocity toolkits. It is just that some splicing data has poor
results because splicing data has biased capture of introns and lack the real time scale because of unknown intrinsic
splicing rate. Such limitations come from the data itself cannot in general be solved with computational models.

.. toctree::
    ./tutorial_hsc_dynamo_megakaryocytes_appearance/tutorial_hsc_dynamo_megakaryocytes_appearance
    ./tutorial_hsc_dynamo_basophil_lineage/tutorial_hsc_dynamo_basophil_lineage
    ./tutorial_hsc_dynamo_cellwise_analysis/tutorial_hsc_dynamo_cellwise_analysis
