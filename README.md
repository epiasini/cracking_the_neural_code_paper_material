# Overview
This MATLAB code runs the simulations used in figures 3 and S1 of the
paper *"Cracking the Neural Code for Sensory Perception by Combining
Statistics, Intervention and Behavior"*, by S. Panzeri, C. D. Harvey,
E. Piasini, P. E. Latham and T. Fellin. The main function to be used
to generate the figures is `intersection_rasters`.

To reproduce the figures in the paper, use the following settings:

- Fig 3A: `intersection_rasters('gaussian', 'difference')`
- Fig 3B: `intersection_rasters('gaussian', 'r1')`
- Fig 3C: `intersection_rasters('gaussian', 'sum')`
- Fig S1A3: `intersection_rasters('elongated_gaussian', 'difference', false, true, 0.01)`
- Fig S1A4: `intersection_rasters('elongated_gaussian', 'difference', false, true, 0.06)`
- Fig S1A5: `intersection_rasters('elongated_gaussian', 'difference', false, true, 0.2)`
- Fig S1B: `intersection_rasters('gaussian', 'r1', true)`

For more details, see the documentation of `intersection_rasters.m`.

# License
This software is licensed under version 3 of the GPL or any later version. See COPYING for details.
