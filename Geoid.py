# authors:
# David Hernandez Lopez, david.hernandez@uclm.es

import os
import sys
import math
import json

current_path = os.path.dirname(__file__)
sys.path.append(os.path.join(current_path, '..'))

from .CRSsTools import CRSsTools
from . import CRSsDefines as cd

from pyLibGDAL.Raster import Raster

class Geoid:
    def __init__(self,
                 precision = cd.GEOID_FULL_PRECISION_CODE):
        self.precision = precision
        self.crs_tools = CRSsTools()
        self.raster = None

    def deflection(self,
                   coordinates,
                   crs_id,
                   interpolation_method = None):
        str_error = ''
        dov_n = None
        dov_e = None
        if not self.raster:
            str_error = ('Geoid is not initialized')
            return str_error, dov_n, dov_e
        if not isinstance(coordinates, list):
            str_error = ('Argument coordinates must be a list and is a: {}'.format(str(type(coordinates))))
            return str_error, dov_n, dov_e
        if len(coordinates) < 2:
            str_error = ('Argument coordinates must be a list with two values at leas')
            return str_error, dov_n, dov_e
        if not isinstance(crs_id, str):
            str_error = ('Argument crs_id must be a string and is a: {}'.format(str(type(crs_id))))
            return str_error, dov_n, dov_e
        band_position = 0
        if not interpolation_method:
            interpolation_method = cd.GEOID_DEFLECTION_INTERPOLATION_METHOD
        str_error, du_dr, du_dc = self.raster.interpolate_derivate(coordinates,
                                                                   crs_id,
                                                                   band_position,
                                                                   interpolation_method)

        return str_error, dov_n, dov_e

    def set_from_raster_file(self,
                             file_path):
        str_error = ''
        if not isinstance(file_path, str):
            str_error = ('File path must be a string and is a: {}'.format(str(type(file_path))))
            return str_error
        self.raster = None
        self.raster = Raster(self.precision)
        str_error = self.raster.set_from_file(file_path,)
        if str_error:
            str_error = ("Setting Geoid from file:\n{}\nError:\n{}".format(file_path, str_error))
            return str_error
        str_error = self.raster.load(True, None) # fully, bands
        return str_error
