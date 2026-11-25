# Check that CITATION.cff has the same version as tiskitpy/version.py
v=$(sed -n "s/^__version__ *= *['\"]\(.*\)['\"].*/\1/p" tiskitpy/version.py)
if ! grep -q "^version: $v$" CITATION.cff; then
    echo "CITATION.cff version is not $v (from tiskitpy/version.py), aborting"
    exit 1
fi

mkdir dist
rm dist/*
python -m build
twine check dist/*
twine upload dist/*
