mkdir dist
rm dist/*
python setup.py sdist
python setup.py bdist_wheel
twine check dist/*
#twine upload --repository-url https://test.pypi.org/legacy/ dist/* --verbose
twine upload dist/*
