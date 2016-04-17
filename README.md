# octave-promethee
PROMETHEE II related octave/matlab functions

#### Context

PROMETHEE II is a widely used (outranking) method used in the context of multi-criteria decision aid (MCDA). 
Please refer to the [Wikipedia](https://en.wikipedia.org/wiki/Preference_ranking_organization_method_for_enrichment_evaluation) page for an introduction and some references to the method.

This repository gathers octave/matlab functions I have developed during my research work. If you're interested in this research, please have a look at my [This work has not been published as such, but is already presented in my [thesis](that.ulb.ac.be/dspace/bitstream/2013/209033/1/283816bc-c6ba-43b2-8d4c-10360bbe6909.txt).


#### Content

Currently, all functions are gathered in one folder. I intend to structure the content in a more convenient way. In the meantime, you will find below the main functions.

| File | Description |
|:---|:---|
| `phi.m` | Compute the PROMETHEE II net flow score. This function is octave/matlab "obtimized" in that is uses matrix computation as much as possible. No precise tests have been made to compare to other implementations, but performance seems reasonably good. The function's name "phi" comes from the usual greek symbol to represent the net flow score associated to an action. |
| `psi_pla.m` | An approximation (hence the symbol "psi" rather than "phi") using the Piecewise Linear Approximation (PLA) method, explained in a research paper. The corresponding technical paper can be downloaded freely [here](http://code.ulb.ac.be/dbfiles/EppDes2012atechreport.pdf). Note that this approximation method only uses the "V-type with indifference" preference function. |
| `psi_eda.m` | A second net flow score approximation using the Empirical Distribution based Approximation (EDA) method. This work has not been published as such, but is already presented in my thesis that can be downloaded freely [here](that.ulb.ac.be/dspace/bitstream/2013/209033/1/283816bc-c6ba-43b2-8d4c-10360bbe6909.txt). Note that this approximation method only uses the "V-type with indifference" preference function. |
