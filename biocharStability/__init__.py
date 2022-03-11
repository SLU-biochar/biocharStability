'''
module: biocharStability

The module `biocharStability` contains multiple functions used for analysing biochar stability data,  developing biochar stability models, and analysing them.

The functions are split in 4 submodules:
- analyse.py
- dashboard.py
- utils.py
- visualize.py

Explore each of them (using the documentation and its menu on the left). Also check out the demo notebooks.

'''

__version__ = "0.1.0"
__author__ = "Elias S. Azzi"

from .utils import *
from .analyse import *
from .visualize import *
from .dashboard import *



