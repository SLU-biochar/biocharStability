from setuptools import setup

import codecs
import os.path

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")

packages = [
    "numpy>=1.23.3",
    "pandas>=1.5.0",
    "geopandas>=0.10.2",
    "matplotlib>=3.6.2",
    "seaborn>=0.11.2",
    "plotly==5.8.2",
    "python-ternary==1.0.8",
    "bokeh>=2.4.3",
    "scipy>=1.9.1",
    "scikit-learn>=1.1.2",
    "stats_arrays>=0.6.5",
    "lmfit>=1.0.3",
    "adjusttext==0.7.3.1",
    "numdifftools==0.9.39",
]

setup(
    name='biocharStability',
    version=get_version("biocharStability/__init__.py"),
    packages=["biocharStability"],
    author="Elias S. Azzi",
    author_email="elias@ecoleaf.consulting",
    license="CC BY-SA 4.0", 
    install_requires=packages,
    url="https://github.com/SLU-biochar/biocharStability",
    description="Biochar incubation database and the code to analyse it!",
    long_description_content_type="text/markdown",
    long_description=open("README.md").read(),
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3 :: Only',
    ],
)