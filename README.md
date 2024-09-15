# Python Interface to TDAvec package

In this repository we describe the python interface to R package `TDAvec`, which can be installed from CRAN

## Installation instructions


First, you need to create a clean python environment. From the project root directory run:

    > python3 -m venv venv
    > source venv/bin/activate
    > pip install numpy==1.26.4 ripser==0.6.8

Now compile and install the package:

    > python3 setup.py build_ext --inplace
    > pip install .

after that you should have `tdavec` package installed in your environment. To check this run python and try the following commands:

    > from tdavec.TDAvectorizer import createEllipse, TDAvectorizer
    >  ee  = createEllipse()
    >  v = TDAvectorizer()
    >  v.fit([ee])
    >  len(v.diags) # ==> 1 since there is only one diagram
    >  len(v.diags[0]) # ==> 1 since there is two dimensions
    >  len(v.diags[0][0]) # ==> 99 since there are 99 hom0 features

or run the provided unit test:

    > python tdavec/unit_test.py

Alternatively you can install the package directly from GitHub:

    > pip install git+https://github.com/ALuchinsky/tdavect
