# Analysis of data from low-resolution NMR relaxometry

In this job, we apply advanced sparse representation tools to the problem of LR-NMR relaxometry. 

What can it do?
------
It provides a tool for processing, analysing, and inspecting Low-Resolution Nuclear Magnetic Resonance Relaxometry data, due to its capability of solving the Laplace Inverse problem using Numerical Optimization Methods in Octave.

How does it work?
------
We use a primal-dual interior method for convex objectives which can be adjusted to solve the LR-NMR relaxometry inverse problem with non-negativity constraints and an L1 regularization term that stabilizes the solution process without introducing the typical L2 peak broadening.

Since a CPMG synthetical signal with three peaks (400, 400 and 200), relaxation times of 8.73, 21.54 and 53.15 ms (respectively), and SNR of 18829. An example of the results is given in the figures below, which represent the results of simulations using L1 and L2 regularization (CO), just L2 regularization (L2), and the original signal (SG). 

![t2_pdco](https://raw.githubusercontent.com/xancandal/lr-nmr-analysis/master/t2_pdco.png)
![residuals](https://raw.githubusercontent.com/xancandal/lr-nmr-analysis/master/residuals.png)
![decay](https://raw.githubusercontent.com/xancandal/lr-nmr-analysis/master/decay.png)
![components](https://raw.githubusercontent.com/xancandal/lr-nmr-analysis/master/components.png)
