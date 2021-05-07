import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()
    
version={}
with open("rptransient/version.py") as fp:
    exec(fp.read(),version)

setuptools.setup(
    name="rptransient",
    version=version['__version__'],
    author="Wayne Crawford (original code by E Wielandt)",
    author_email="crawford@ipgp.fr",
    description="Remove regular data transients from seismological data",
    long_description=long_description,
    long_description_content_type="text/x-rst; charset=UTF-8",
    url="https://github.com/WayneCrawford/rptransient",
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=['obspy','PyYAML'],
    entry_points={
         'console_scripts': [
         ]
    },
    python_requires='>=3.8',
    classifiers=(
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics"
    ),
    keywords='oceanography, marine, OBS'
)