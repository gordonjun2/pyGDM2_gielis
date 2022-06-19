# encoding: utf-8
#
#Copyright (C) 2017-2020, P. R. Wiecha
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
    pyGDM2 - a full-field electrodynamical solver python toolkit.
    Based on the Green Dyadic Method.
    
"""

from . import structures

from . import materials
from . import fields

from . import core
from . import linear
from . import tools
from . import visu
from . import visu3d

__name__ = 'pyGDM2'
__version__ = '1.0.11'
__date__ = "02/25/2020"
__author__ = 'Peter R. Wiecha'

__all__ = ["core", "propagators", "materials", "structures", "fields", 
           "linear", "linear_py", "nonlinear", 
           "tools", "visu", "EO"]