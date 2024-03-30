import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()
    
version={}
with open("tiskitpy/version.py") as fp:
    exec(fp.read(),version)

setuptools.setup(
    name="tiskitpy",
    version=version['__version__'],
    author="Wayne Crawford",
    author_email="crawford@ipgp.fr",
    description="TIme Series toolKIT",
    long_description=long_description,
    long_description_content_type="text/markdown; charset=UTF-8",
    url="https://github.com/WayneCrawford/tiskitpy",
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=['obspy>=1.3.0','PyYAML', 'numpy', 'scipy', 'matplotlib',
                      'xarray'],
    entry_points={
         'console_scripts': [
            'tiskitpy_decimate_SDS=tiskitpy.scripts.decimate_SDS:main', 
            'tiskitpy_get_SDS_inventory=tiskitpy.scripts.get_SDS_inventory:main'
         ]
    },
    python_requires='>=3.8',
    classifiers=[
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
    ],
    keywords='oceanography, marine, OBS'
)
