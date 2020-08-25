# ot-neuro-review-supp
This repository contains supplementary code used to generate figure 2 from

> Optical tweezers exploring neuroscience,
> Isaac C. D. Lenton, Ethan K. Scott, Halina Rubinsztein-Dunlop, and Itia A. Favre-Bulle,
> Frontiers in Physics, 2020.

If you find any of this code useful, please consider citing our paper.


Requirements
------------

  * [Optical Tweezers Toolbox](https://github.com/ilent2/ott) (version 1, should work with 1.5.6)
  * [OTSLM](https://github.com/ilent2/otslm) (1.0.1)
  * Matlab (tested on R2018a)

Files
-----

  * `sim001_slm_examples.m` SLM phase patterns generated with OTSLM and simulated using OTT.
  
  * `sim002_scanned_beam.m` Dynamics simulation showing a particle being dragged by a
    focussed Gaussian beam.  Uses interpolation for force (a similar idea to our
    [artificial neural networks approach](https://doi.org/10.1088/2632-2153/abae76)).
    
  * `sim003_feedback.m` An illustrative example showing how feedback can be used to
    create a custom optical potential shape.
