Files for python version of the TDAvec package

# Ways to compile and publish:

Run the following commands from project root directory:

    > cd src
    > make build_wheel
    > make upload_test
    > cd -

when password is asked, enter the password from API_token file (not included)

Be sure that you are uploading a new version of the package!

# How to install

Run

    > python3 -m pip install --index-url https://test.pypi.org/simple/ --no-deps tdavec

to install the package. You can see that it is installed by running by trying ipython:

    > ipython
    []  tdavec.TDAvectorizer import *
    []  X = createEllipse(100)
    []  vec = TDAvectorizer()
    []  vec.fit(X)
    []  D = vec.diags()

As a result, in object D we will have a persistence diagram, which is simply list of two lists: birth and death values for dimension 0 and 1.

