Fabber models for Chemical Exchange Saturation Transfer (CEST) MRI
==================================================================

These models use the Fabber_
Bayesian model fitting framework [1]_ to implement a model
for CEST-MRI.

.. note::
    The Quantiphyse_ application contains a widget for performing
    CEST analysis which uses ``FABBER_CEST`` internally. In addition
    the Baycest_ package within FSL_ can be used to run CEST analysis
    using FABBER_CEST.
    
The CEST model is a multi-pool model which is fitted to the sampled Z-spectrum. There are a large number of 
parameters within the Bloch-McConnell equations that describe each pool making model fitting prone to inaccuracy 
and increasing the risk of over fitting. A solution is to provide prior information about the parameters which 
necessitates the use of a Bayesian method. 

FABBER_CEST exploits a Bayesian non-linear fitting algorithm, which is 
essentially a probabilistic version of non-linear least squares, along with a multi-pool implementation of the 
Bloch-McConnell equations. This algorithm provides a (relatively) fast means to quantify CEST data whilst reducing 
some of the problems associated with traditional least squares fitting algorithms.

Getting FABBER_CEST
--------------------

The CEST models are included in FSL_. We
stongly recommend version 6.0.1 or later.

If you need an updated version of the model which has not yet been released to
FSL, you will either need to 
`build from source <https://fabber-core.readthedocs.io/en/latest/building.html#building-new-or-updated-model-libraries>`_ 
using an existing FSL 6.0.1 or later installation, or download 
the pre-built `Fabber bundle <https://fabber-core.readthedocs.io/en/latest/getting.html#standalone-fabber-distribution>`_ 
which contains the latest CEST release alongside other models in a standalone package.

The CEST model
--------------

The CEST model is activated by calling ``fabber_cest`` with the option ``--model=cest``.

Examples
--------

References
----------

.. [1] *Chappell, M.A., Groves, A.R., Woolrich, M.W., "Variational Bayesian
   inference for a non-linear forward model", IEEE Trans. Sig. Proc., 2009,
   57(1), 223–236.*

.. _Fabber: https://fabber-core.readthedocs.io/

.. _FSL: https://fsl.fmrib.ox.ac.uk/fsl/

.. _Quantiphyse: https://quantiphyse.readthedocs.io/

.. _Baycest https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/baycest
