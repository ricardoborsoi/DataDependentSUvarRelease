#  A Data Dependent Multiscale Model for Hyperspectral Unmixing with Spectral Variability    #

This package contains the authors' implementation of the paper [1].

We consider a multiscale strategy using superpixels in order to address the spectral unmixing problem with endmember variability. We use the spatial regularity information about both the abundances and the endmembers (that is, that thtese variables are smooth according to the superpixel-based multiscale transformation) in order to 1) introduce a priori information in order to improve the abundance estimation quality, and 2) reformulate the optimization problem in order to significantly reduce the computational complexity of the method.

The code is implemented in MATLAB and includes:  
-  example1.m                - a demo script comparing the algorithms (DC1)  
-  example2.m                - a demo script comparing the algorithms (DC2)  
-  example3.m                - a demo script comparing the algorithms (DC3)  
-  demo_houston.m            - a demo script comparing the algorithms (Houston)  
-  demo_cuprite.m            - a demo script comparing the algorithms (Cuprite)  
-  ./MUAV/                   - contains the MATLAB files associated with the MUAV algorithm  
-  ./other_methods/          - contains the ELMM and PLMM methods  
-  ./utils/                  - useful functions  
-  ./DATA/                   - files used in the examples  
-  README                    - this file  


## IMPORTANT:
If you use this software please cite the following in any resulting
publication:

    [1] A Data Dependent Multiscale Model for Hyperspectral Unmixing with Spectral Variability
        R.A. Borsoi, T. Imbiriba, J.C.M. Bermudez.
        IEEE Transactions on Image Processing, 2020.



## INSTALLING & RUNNING:
Just start MATLAB and run one of the demo scripts (e.g. example1.m, example2.m, etc).

## NOTES:
1.  Parts of our code are based on the tooblox for the ELMM algorithm provided by Lucas Drumetz.

2.  The ELMM algorithm was provided by Lucas Drumetz.  
    Drumetz, L., Veganzones, M.-A., Henrot, S., Phlypo, R., Chanussot, J., & Jutten, C.
    Blind hyperspectral unmixing using an extended linearmixing model to address spectral variability.
    IEEE Transactions on Image Processing, 2016.

3.  The PLMM algorithm was provided by Pierre-Antoine Thouvenin.  
    Thouvenin, P.-A., Dobigeon, N., & Tourneret, J.-Y.
    Hyperspectral unmixing with spectral variability using a perturbed linear mixing model.
    IEEE Transactions on Signal Processing, 2016.

4.  The data used for example 3 was originally provided by Lucas Drumetz and his collaborators






