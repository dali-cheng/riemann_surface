<h3 align="center">Experimental observation of energy-band Riemann surface</h3>

  <p align="justify">
    Scripts and data for the paper "Experimental observation of energy-band Riemann surface," available on <a href="https://arxiv.org/abs/2510.08819">arXiv</a>. If you find this repository helpful, please consider citing our work.
    <br />


## About
This repository contains scripts and experimental data for the above-mentioned manuscript. 
Download all the files into the same folder, and run the `.m` codes to generate the figures in the manuscript.
`find_closest.m`, `convert_to_0_to_a.m`, `convert_to_plus_minus_a.m`, and `lorentzian_fittype.m` are user-defined helper functions. 
Tested on MATLAB R2024b.
| File name | Comments |
|-----------|----------------|
| Fig2C.m | Generates Figure 2(C) |
| Fig3A.m | Generates Figure 3(A) |
| Fig3B4A.m | Generates Figures 3(B) and 4(A) |
| Fig3C.m | Generates Figure 3(C) |
| Fig3D4C.m | Generates Figures 3(D) and 4(C) |
| Fig3E.m | Generates Figure 3(E) |
| Fig4B.m | Generates Figure 4(B) |
| Fig4D.m | Generates Figure 4(D) |
| Fig4E.m | Generates Figure 4(E) |
| Fig4F.m | Generates Figure 4(F) |
| Fig4G.m | Generates Figure 4(G) |
| Fig4H.m | Generates Figure 4(H) |
| Fig5A.m | Generates Figure 5(A) |
| Fig5B.m | Generates Figure 5(B) |
| Fig5C.m | Generates Figure 5(C) |
| Fig5D.m | Generates Figure 5(D) |
| FigS5.m | Generates Figure S5 in the Supplementary Materials |
| spcsp.m | <a href="https://www.mathworks.com/matlabcentral/fileexchange/59463-smoothing-cubic-splines-with-periodic-conditions">Smoothing cubic splines with periodic conditions</a>, by M. Zanetti |
| sigma_X.mat | Experimental data with Im(_k_) = _σ_ = X. The variable `bandstr` contains reshaped raw data from the oscilloscope on the time-dependent transmission spectrum of the resonator, and `delta_x_mean`, `real_E`, and `imag_E` are results of Lorentzian-lineshape fitting of the spectrum. |
| gbz_theory.mat | Pre-calculated theoretical GBZ of the model. |
| OBC_GBZ.mat | Experimental results of the OBC spectrum and the GBZ of the model. Values are determined from the self-intersecting points of the complex windings as described in the manuscript. |
| OBC_GBZ_simulation.mat | Simulation results of the OBC spectrum and the GBZ of the model. Values are determined from the self-intersecting points of the complex windings as described in the manuscript. |


## License
See `LICENSE.txt` for more information.


## Contact
Questions or comments about this repository should be addressed to:
* Dali Cheng - chengdl@stanford.edu
* Shanhui Fan - shanhui@stanford.edu


## Acknowledgments
This work was supported by a MURI project from the U.S. Air Force Office of Scientific Research (Grant No. FA9550-22-1-0339). C. R.-C. is supported by a Stanford Science Fellowship.
