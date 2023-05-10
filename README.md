# MATLAB

MATLAB code and simulations described in the website:
https://atsushitakagi.wordpress.com/
are hosted here. 

If parts of the MATLAB code are used for publication,
please mention in the acknowledgements!

As of 15 September 2022 the following code is available:

- Kalman filter basics

    Learn how to implement a Kalman filter to obtain model-based estimates from noisy measurements, and to fuse measurements from multiple sensors.

- Finite horizion Linear Quadratic Regulator (simulation of one-dimensional reaching)

    The basics of learning how to simulate and control a dynamic system, and demonstrate how a bell-shaped velocity profile is generated.

- Human-human goal integration

    This code shows how to create a model of the partner to estimate their movement goal described in Takagi et al. (2017) (https://www.nature.com/articles/s41562-017-0054). It also contains the original data from the human-human experiment conducted on the TVINS interface at ATR, Japan.
    
- Interpersonal goal integration with compliant coupling

    The code reproduces the simulations from Takagi et al. (2018), which shows how the compliance in the physical coupling between physically interacting partners determines the ability to estimate the partner's goal or target (https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005971). In this version, the estimation of the partner's goal is simplified to improve computational speed (instead of estimating it from the force, we directly give the noisy haptic target reading. For the full implementation, see the "Human-human goal integration" MATLAB code). It contains the original data collected from the dual-robotic wrist interface at Imperial College London.


by Atsushi Takagi
