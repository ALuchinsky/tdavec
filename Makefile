build_wheel:
	python setup.py bdist_wheel 

clean:
	rm -rf dist/ build/ tdavec.egg-info tdavec/TDAvec.c
	python3 -m pip uninstall tdavec -y

upload_test:
	python3 -m twine upload --repository testpypi dist/*

install_package:
	python3 -m pip install --index-url https://test.pypi.org/simple/ --no-deps tdavec

