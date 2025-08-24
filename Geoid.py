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
                   a, e2,
                   band_position,
                   geoid_model_accuracy,
                   interpolation_method):
        str_error = ''
        dov_n = None
        dov_e = None
        sdev_dov_n = None
        sdev_dov_e = None
        cov_dov_n_e = None
        if not self.raster:
            str_error = ('Geoid is not initialized')
            return str_error, dov_n, dov_e, sdev_dov_n, sdev_dov_e, cov_dov_n_e
        if not isinstance(coordinates, list):
            str_error = ('Argument coordinates must be a list and is a: {}'.format(str(type(coordinates))))
            return str_error, dov_n, dov_e, sdev_dov_n, sdev_dov_e, cov_dov_n_e
        if len(coordinates) < 2:
            str_error = ('Argument coordinates must be a list with two values at leas')
            return str_error, dov_n, dov_e, sdev_dov_n, sdev_dov_e, cov_dov_n_e
        if not isinstance(crs_id, str):
            str_error = ('Argument crs_id must be a string and is a: {}'.format(str(type(crs_id))))
            return str_error, dov_n, dov_e, sdev_dov_n, sdev_dov_e, cov_dov_n_e
        # band_position = 0
        if not interpolation_method:
            interpolation_method = cd.GEOID_DEFLECTION_INTERPOLATION_METHOD
        str_error, du_dr, du_dc = self.raster.interpolate_derivate(coordinates,
                                                                   crs_id,
                                                                   band_position,
                                                                   interpolation_method)
        if str_error:
            return str_error, dov_n, dov_e, sdev_dov_n, sdev_dov_e, cov_dov_n_e
        length_e = None
        length_n = None
        str_error, raster_crs_is_geographic = self.crs_tools.is_geographic(self.raster.crs_id)
        if str_error:
            return str_error, dov_n, dov_e, sdev_dov_n, sdev_dov_e, cov_dov_n_e
        if raster_crs_is_geographic:
            lat_rad = coordinates[1] * math.pi / 180
            rn = a / math.sqrt(1 - e2 * (math.sin(lat_rad) ** 2))
            grid_size_e_rad = self.raster.size_fc * math.pi / 180
            length_e = rn * math.cos(lat_rad) * grid_size_e_rad
            grid_size_n_rad = self.raster.size_sc * math.pi / 180
            rm = rn * (1 - e2) / (1 - e2 * (math.sin(lat_rad) ** 2))
            length_n = rm * grid_size_n_rad
        else:
            length_e = self.raster.size_fc
            length_n = self.raster.size_sc
        dov_n_rad = math.atan2(du_dr, length_n)
        dov_e_rad = -math.atan2(du_dc, length_e)
        dov_n = dov_n_rad * 180. / math.pi * 3600.
        dov_e = dov_e_rad * 180. / math.pi * 3600.
        if geoid_model_accuracy:
            n_geoid_points = 4
            if interpolation_method == cd.GEOID_DEFLECTION_INTERPOLATION_METHOD_BICUBIC:
                n_geoid_points = 16
            # f(x) = atan(x) -> f'(x) = 1 / (1 + x ** 2)
            # atan2(x) = tan(sin/cos)
            # std_f(x) = f'(x) * std_x
            der_dov_n_rad = 1. / (1. + (du_dr / length_n) ** 2.)
            der_dov_e_rad = 1. / (1. + (du_dc / length_e) ** 2.)
            sdev_dov_n_rad = der_dov_n_rad * du_dr * geoid_model_accuracy / (10 ** 6) * length_n / length_n
            sdev_dov_e_rad = der_dov_e_rad * du_dc * geoid_model_accuracy / (10 ** 6) * length_e / length_e
            # sdev_dov_n_rad = der_dov_n_rad * math.sqrt(n_geoid_points) * (geoid_model_accuracy / length_n)
            # sdev_dov_e_rad = der_dov_e_rad * math.sqrt(n_geoid_points) * (geoid_model_accuracy / length_e)
            sdev_dov_n = sdev_dov_n_rad * 180. / math.pi * 3600.
            sdev_dov_e = sdev_dov_e_rad * 180. / math.pi * 3600.
            correlation_dov_n_e = 1.
            # correlation_xy = covariance_xy / (sd_x * sd_y)
            cov_dov_n_e = sdev_dov_n * sdev_dov_e
        return str_error, dov_n, dov_e, sdev_dov_n, sdev_dov_e, cov_dov_n_e

    def ondulation(self,
                   coordinates,
                   crs_id,
                   a, e2,
                   band_position,
                   geoid_model_accuracy, # None
                   interpolation_method):
        str_error = ''
        ondulation = None
        sdev_ondulation = None
        if not self.raster:
            str_error = ('Geoid is not initialized')
            return str_error, ondulation, sdev_ondulation
        if not isinstance(coordinates, list):
            str_error = ('Argument coordinates must be a list and is a: {}'.format(str(type(coordinates))))
            return str_error, ondulation, sdev_ondulation
        if len(coordinates) < 2:
            str_error = ('Argument coordinates must be a list with two values at leas')
            return str_error, ondulation, sdev_ondulation
        if not isinstance(crs_id, str):
            str_error = ('Argument crs_id must be a string and is a: {}'.format(str(type(crs_id))))
            return str_error, ondulation, sdev_ondulation
        # band_position = 0
        if not interpolation_method:
            interpolation_method = cd.GEOID_DEFLECTION_INTERPOLATION_METHOD
        return self.raster.interpolate(coordinates,
                                       crs_id,
                                       band_position,
                                       interpolation_method)
        # str_error, ondulation = self.raster.interpolate_derivate(coordinates,
        #                                                          crs_id,
        #                                                          band_position,
        #                                                          interpolation_method)
        # if str_error:
        #     return str_error, ondulation, sdev_ondulation
        #
        # return str_error, ondulation, sdev_ondulation

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
