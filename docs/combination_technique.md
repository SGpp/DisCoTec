# Sparse Grid Combination Technique with Time Stepping

The sparse grid combination technique (Griebel et al.
[1992](https://ins.uni-bonn.de/media/public/publication-media/griesiam.ps.gz),
Garcke [2013](https://link.springer.com/chapter/10.1007/978-3-642-31703-3_3),
Harding [2016](https://link.springer.com/chapter/10.1007/978-3-319-28262-6_4))
can be used to alleviate the curse of dimensionality encountered in
high-dimensional simulations.
Instead of using your solver on a single structured full grid (where every
dimension is finely resolved), you would use it on many different structured
full grids (each of them differently resolved).
We call these coarsely-resolved grids component grids.
Taken together, all component grids form a sparse grid approximation, which can
be explicitly obtained by a linear superposition of the individual grid
functions, with the so-called combination coefficients.

![schematic of a combination scheme in 2D](../gfx/combischeme-2d.svg)

In this two-dimensional combination scheme, all combination coefficients are 1
and -1, respectively.
Figure originally published in (Pollinger [2024](https://elib.uni-stuttgart.de/handle/11682/14229), adapted from Pollinger et al. [2023](https://dl.acm.org/doi/10.1145/3581784.3607036)).

Between time steps, the grids exchange data through an intermediate multi-scale
represenation, which is summarized as the "combination" step in DisCoTec.
Assuming a certain smoothness in the solution, this allows for a good
approximation of the finely-resolved function, while achieving drastic
reductions in compute and memory requirements.