# TDAvec.py

`TDAvec.py` is a python interface to `TDAvec` R package, which is available on [CRAN](https://cran.r-project.org/web/packages/TDAvec/index.html)

First of all, it allows access to all implemented in the original R package vectorizations functions:

* computeAlgebraicFunctions:	Compute Algebraic Functions from a Persistence Diagram
* computeBettiCurve:	A Vector Summary of the Betti Curve
* computeComplexPolynomial:	Compute Complex Polynomial Coefficients from a Persistence Diagram
* computeEulerCharacteristic:	A Vector Summary of the Euler Characteristic Curve
* computeNormalizedLife:	A Vector Summary of the Normalized Life Curve
* computePersistenceBlock:	A Vector Summary of the Persistence Block
* computePersistenceImage:	A Vector Summary of the Persistence Surface
* computePersistenceLandscape:	Vector Summaries of the Persistence Landscape Functions
* computePersistenceSilhouette:	A Vector Summary of the Persistence Silhouette Function
* computePersistentEntropy:	A Vector Summary of the Persistent Entropy Summary Function
* computeStats:	Compute Descriptive Statistics for Births, Deaths, Midpoints, and Lifespans in a Persistence Diagram
* computeTemplateFunction:	Compute a Vectorization of a Persistence Diagram based on Tent Template Functions
* computeTropicalCoordinates:	Compute Tropical Coordinates from a Persistence Diagram

All these functions can easily be called using `tdavec.tdavec_core` package.

In addition, we provide also `sklearn`-type interface to the same functionality, which could be more familiar for python programmers.

Note that the package was tested only on python 3.12. 

# Setup

`TDAvec.py` is available on `pypi`. To install it simply type

    pip install tdavec

into your environment. 

You can also install the current verion from the GitHub with

    pip install git+https://github.com/ALuchinsky/tdavect

Alternatively, you can install it from the source. In order to do this clone mentioned above github repository and run the followin commants from the project root directory:


    pip install numpy==1.26.4 ripser==0.6.8
    python3 setup.py build_ext --inplace
    pip install .

after that you should have `tdavec` package installed in your environment. 


In order to check if the intallation process was completed, you can run python and evaluate the following lines:

    > from tdavec import test_package
    > X, D, PS = test_package()

This function will create a simple point cloud, build a persistence diagram, caclulate the Persistence Silhouette from it, and return these three objects.

