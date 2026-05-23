mkdir dist
rm dist/*
python -m build
twine check dist/*
twine upload dist/*
