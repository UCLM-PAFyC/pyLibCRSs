# authors:
# David Hernandez Lopez, david.hernandez@uclm.es

# https://pyproj4.github.io/pyproj/stable/examples.html
# https://pyproj4.github.io/pyproj/stable/advanced_examples.html
import os
import json

from pyproj import CRS, Transformer, database, pyproj
from pyproj._crs import Datum
from pyproj.enums import TransformDirection

from . import CRSsDefines as cd
import numpy as np


class CRSsTools:
    def __init__(self):
        self.proj_version = pyproj.__version__
        self.json_data_file = None
        self.CRSs = {}
        self.CRSOperations = {}
        self.CRS_geo3d_id_by_enu_id = {}
        self.CRSOperation_by_enu_id = {}
        self.data = {}
        self.exists_crss_info = False
        # self.CRSs_compound_ids = []
        self.data["CRSs_ecef_ids"] = []
        self.data["CRSs_geo2d_ids"] = []
        self.data["CRSs_geo3d_ids"] = []
        self.data["CRSs_projected_ids"] = []
        self.data["CRSs_vertical_ids"] = []
        self.data["CRSs_compound_ids"] = []
        self.data["datum_by_crs_id"] = {}
        self.data["CRSs_geo2d_ids_by_datum"] = {}
        self.data["CRSs_geo3d_ids_by_datum"] = {}
        self.data["CRSs_ecef_ids_by_datum"] = {}
        self.data["CRSs_projected_ids_by_datum"] = {}
        self.data["base_crs_id_by_crs_id"] = {}
        self.data["CRSs_geo2d_ids_by_base_crs_id"] = {}
        self.data["CRSs_geo3d_ids_by_base_crs_id"] = {}
        self.data["CRSs_ecef_ids_by_base_crs_id"] = {}
        self.data["CRSs_projected_ids_by_base_crs_id"] = {}
        self.data["CRSs_vertical_ids_by_base_crs_id"] = {}
        self.data['CRS_info_by_id'] = {}
        self.initialize()

    def covariance_matrix_operation(self,
                                    crs_source_id,
                                    crs_target_id,
                                    source_point,
                                    target_point,
                                    source_covariance_matrix):
        str_error = ""
        target_covariance_matrix = None
        jacobian_matrix = None
        if crs_source_id == crs_target_id:
            return str_error, target_covariance_matrix, jacobian_matrix
        str_aux_error, jacobian_matrix = self.get_jacobian(crs_source_id,
                                                           crs_target_id,
                                                           source_point,
                                                           target_point)
        if str_aux_error:
            str_error = CRSsTools.__name__ + "." + self.covariance_matrix_operation.__name__
            str_error += ("\nGetting jacobian matrix, error:\n {}".format(str_aux_error))
            return str_error, target_covariance_matrix, jacobian_matrix
        jacobian_transpose = jacobian_matrix.transpose()
        target_covariance_matrix = np.matmul(source_covariance_matrix, jacobian_matrix.transpose())
        target_covariance_matrix = np.matmul(jacobian_matrix, target_covariance_matrix)
        return str_error, target_covariance_matrix, jacobian_matrix

    # def get_crs_compound_ids(self):
    #     if not self.set_crss_info:
    #         self.set_crss_info()
    #     return self.data["CRSs_ecef_ids"]
    #     # if len(self.CRSs_compound_ids) > 0:
    #     #     return self.CRSs_compound_ids
    #     # crs_info_list = database.query_crs_info(auth_name='EPSG',
    #     #                                         pj_types="COMPOUND_CRS")
    #     # crs_list = ["EPSG:" + info[1] for info in crs_info_list]
    #     # self.CRSs_compound_ids = sorted(crs_list)
    #     # return self.CRSs_compound_ids

    def get_compound_epgs_codes_from_json(self,
                                          crs_as_dict):
        str_error = ''
        epsg_code = cd.NO_EPSG_CODE
        vertical_epsg_code = cd.NO_EPSG_CODE
        if not isinstance(crs_as_dict, dict):
            str_error = CRSsTools.__name__ + "." + self.get_compound_epgs_codes_from_wkt.__name__
            str_error += ("\nArgument must be a dict and is: {}".format(str(type(crs_as_dict))))
            return str_error, epsg_code, vertical_epsg_code
        if not cd.CRS_AS_JSON_SCHEMA_TAG in crs_as_dict:
            str_error = CRSsTools.__name__ + "." + self.get_compound_epgs_codes_from_wkt.__name__
            str_error += ("\nNo: {} in json string:\n{}".format(cd.CRS_AS_JSON_SCHEMA_TAG, json_string))
            return str_error, epsg_code, vertical_epsg_code
        schema =  crs_as_dict[cd.CRS_AS_JSON_SCHEMA_TAG]
        return str_error, epsg_code, vertical_epsg_code

    def get_compound_epgs_codes_from_wkt(self,
                                         wkt):
        str_error = ''
        epsg_code = cd.NO_EPSG_CODE
        vertical_epsg_code = cd.NO_EPSG_CODE
        if not isinstance(wkt, str):
            str_error = CRSsTools.__name__ + "." + self.get_compound_epgs_codes_from_wkt.__name__
            str_error += ("\nWKT must be a string an is: {}".format(str(type(wkt))))
            return str_error, epsg_code, vertical_epsg_code
        wkt = wkt.strip()
        if not wkt.startswith(cd.CRSTOOLS_COMPOUNT_WKT_TAG):
            str_error = CRSsTools.__name__ + "." + self.get_compound_epgs_codes_from_wkt.__name__
            str_error += ("\nWKT must be start with: {}\nand is:\n{}".format(cd.CRSTOOLS_COMPOUNT_WKT_TAG), wkt)
            return str_error, epsg_code, vertical_epsg_code
        if not wkt.startswith(cd.CRSTOOLS_COMPOUNT_WKT_TAG):
            str_error = CRSsTools.__name__ + "." + self.get_compound_epgs_codes_from_wkt.__name__
            str_error += ("\nWKT must be start with: {}\nand is:\n{}".format(cd.CRSTOOLS_COMPOUNT_WKT_TAG), wkt)
            return str_error, epsg_code, vertical_epsg_code

        return str_error, epsg_code, vertical_epsg_code

    def get_crs_ecef_ids(self):
        return self.data["CRSs_ecef_ids"]

    def get_crs_geo2d_ids(self):
        return self.data["CRSs_geo2d_ids"]

    def get_crs_ecef_ids_for_crs_geo2d_id(self,
                                          crs_geo2d_id):
        datum_name = None
        if crs_geo2d_id in self.data["datum_by_crs_id"]:
            datum_name = self.data["datum_by_crs_id"][crs_geo2d_id]
            if datum_name in self.data["CRSs_ecef_ids_by_datum"]:
                return self.data["CRSs_ecef_ids_by_datum"][datum_name]
        base_crs_id = None
        if crs_geo2d_id in self.data["base_crs_id_by_crs_id"]:
            base_crs_id = self.data["base_crs_id_by_crs_id"][crs_geo2d_id]
            if base_crs_id in self.data["CRSs_ecef_ids_by_base_crs_id"]:
                return self.data["CRSs_ecef_ids_by_base_crs_id"][base_crs_id]
        # return self.data["CRSs_ecef_ids"]
        return None

    def get_crs_geo3d_ids_for_crs_geo2d_id(self,
                                           crs_geo2d_id):
        datum_name = None
        if crs_geo2d_id in self.data["datum_by_crs_id"]:
            datum_name = self.data["datum_by_crs_id"][crs_geo2d_id]
            if datum_name in self.data["CRSs_geo3d_ids_by_datum"]:
                return self.data["CRSs_geo3d_ids_by_datum"][datum_name]
        base_crs_id = None
        if crs_geo2d_id in self.data["base_crs_id_by_crs_id"]:
            base_crs_id = self.data["base_crs_id_by_crs_id"][crs_geo2d_id]
            if base_crs_id in self.data["CRSs_geo3d_ids_by_base_crs_id"]:
                return self.data["CRSs_geo3d_ids_by_base_crs_id"][base_crs_id]
        # return self.data["CRSs_geo3d_ids"]
        return None

    def get_crs_geo3d_ids(self):
        return self.data["CRSs_geo3d_ids"]

    def get_crs_projected_ids(self):
        return self.data["CRSs_projected_ids"]

    def get_crs_info_as_text(self, crs_id):
        str_error = ""
        crs_info_as_text = None
        if not crs_id in self.data['CRS_info_by_id']:
            str_error = CRSsTools.__name__ + "." + self.get_ellipsoid.__name__
            str_error += ("\nCRS: {} not found".format(crs_id))
            return str_error, crs_info_as_text
        crs_info = self.data['CRS_info_by_id'][crs_id]
        crs_info_keys = crs_info.keys()
        crs_info_as_text = ''
        for key in crs_info_keys:
            crs_info_as_text = crs_info_as_text + '- ' + key + ': '
            value = crs_info[key]
            if isinstance(value, dict):
                str_value = ''
                for key_in_value in value:
                    str_value += '\n  - '
                    str_value += key_in_value
                    str_value += ': '
                    str_value += value[key_in_value]
            else:
                str_value = value
            crs_info_as_text += str_value
            crs_info_as_text += '\n'
        return str_error, crs_info_as_text

    def get_crs_projected_ids_for_crs_geo2d_id(self,
                                               crs_geo2d_id):
        datum_name = None
        if crs_geo2d_id in self.data["datum_by_crs_id"]:
            datum_name = self.data["datum_by_crs_id"][crs_geo2d_id]
            if datum_name in self.data["CRSs_projected_ids_by_datum"]:
                return self.data["CRSs_projected_ids_by_datum"][datum_name]
        base_crs_id = None
        if crs_geo2d_id in self.data["base_crs_id_by_crs_id"]:
            base_crs_id = self.data["base_crs_id_by_crs_id"][crs_geo2d_id]
            if base_crs_id in self.data["CRSs_projected_ids_by_base_crs_id"]:
                return self.data["CRSs_projected_ids_by_base_crs_id"][base_crs_id]
        # return self.data["CRSs_projected_ids"]
        return None

    def get_crs_summary(self, crs_id):
        str_error = ""
        crs_summary = None
        if not crs_id in self.data['CRS_info_by_id']:
            str_error = CRSsTools.__name__ + "." + self.get_ellipsoid.__name__
            str_error += ("\nCRS: {} not found".format(crs_id))
            return str_error, crs_summary
        crs_summary = crs_id + ', ' + self.data['CRS_info_by_id'][crs_id]['name']
        return str_error, crs_summary

    def get_crs_vertical_ids(self):
        return self.data["CRSs_vertical_ids"]

    def get_crs_vertical_ids_for_crs_geo2d_id(self,
                                              crs_geo2d_id):
        # datum_name = None
        # if crs_geo2d_id in self.data["datum_by_crs_id"]:
        #     datum_name = self.data["datum_by_crs_id"][crs_geo2d_id]
        #     if datum_name in self.data["CRSs_projected_ids_by_datum"]:
        #         return self.data["CRSs_projected_ids_by_datum"][datum_name]
        base_crs_id = None
        if crs_geo2d_id in self.data["base_crs_id_by_crs_id"]:
            base_crs_id = self.data["base_crs_id_by_crs_id"][crs_geo2d_id]
            if base_crs_id in self.data["CRSs_vertical_ids_by_base_crs_id"]:
                return self.data["CRSs_vertical_ids_by_base_crs_id"][base_crs_id]
        # return self.data["CRSs_vertical_ids"]
        return None

    def get_ellipsoid(self, crs_id):
        str_error = ""
        ellipsoid = None
        str_aux_error, is_geographic = self.is_geographic(crs_id)
        if str_aux_error:
            str_error = CRSsTools.__name__ + "." + self.get_ellipsoid.__name__
            str_error += ("\nGetting CRS: {} is geographic, error:\n{}".format(crs_id, str_aux_error))
            return str_error, ellipsoid
        crs = self.CRSs[crs_id]
        if not is_geographic:
            if not crs.geodetic_crs:
                str_error = CRSsTools.__name__ + "." + self.get_ellipsoid.__name__
                str_error += ("\nCRS: {} must be geographic,".format(crs_id))
                return str_error, ellipsoid
        if is_geographic:
            ellipsoid = crs.get_geod()
        else:
            ellipsoid = crs.geodetic_crs.get_geod()
        return str_error, ellipsoid

    def get_jacobian(self,
                     crs_source_id,
                     crs_target_id,
                     source_point,
                     target_point):
        str_error = ""
        jacobian_matrix = None
        if crs_source_id == crs_target_id:
            return str_error, jacobian_matrix
        if not isinstance(source_point, list):
            str_error = CRSsTools.__name__ + "." + self.get_jacobian.__name__
            str_error += f"\nSource point must be a list"
            return str_error, jacobian_matrix
        source_dimension = len(source_point)
        if not isinstance(target_point, list):
            str_error = CRSsTools.__name__ + "." + self.get_jacobian.__name__
            str_error += f"\nTarget point must be a list"
            return str_error, jacobian_matrix
        target_dimension = len(target_point)
        str_aux_error, crs_source_is_geographic = self.is_geographic(crs_source_id)
        if str_aux_error:
            str_error = CRSsTools.__name__ + "." + self.get_jacobian.__name__
            str_error += ('\nGetting CRS: {} is geographic, error:\n{}'.format(crs_source_id, str_aux_error))
            return str_error, jacobian_matrix
        str_aux_error, crs_source_is_projected = self.is_projected(crs_source_id)
        if str_aux_error:
            str_error = CRSsTools.__name__ + "." + self.get_jacobian.__name__
            str_error += ('\nGetting CRS: {} is projected, error:\n{}'.format(crs_source_id, str_aux_error))
            return str_error, jacobian_matrix
        str_aux_error, crs_source_is_geocentric = self.is_geocentric(crs_source_id)
        if str_aux_error:
            str_error = CRSsTools.__name__ + "." + self.get_jacobian.__name__
            str_error += ('\nGetting CRS: {} is geocentric, error:\n{}'.format(crs_source_id, str_aux_error))
            return str_error, jacobian_matrix
        str_aux_error, crs_source_is_only_vertical = self.is_only_vertical(crs_source_id)
        if str_aux_error:
            str_error = CRSsTools.__name__ + "." + self.get_jacobian.__name__
            str_error += ('\nGetting CRS: {} is vertical, error:\n{}'.format(crs_source_id, str_aux_error))
            return str_error, jacobian_matrix
        str_aux_error, crs_target_is_geographic = self.is_geographic(crs_target_id)
        if str_aux_error:
            str_error = CRSsTools.__name__ + "." + self.get_jacobian.__name__
            str_error += ('\nGetting CRS: {} is geographic, error:\n{}'.format(crs_target_id, str_aux_error))
            return str_error, jacobian_matrix
        str_aux_error, crs_target_is_projected = self.is_projected(crs_target_id)
        if str_aux_error:
            str_error = CRSsTools.__name__ + "." + self.get_jacobian.__name__
            str_error += ('\nGetting CRS: {} is projected, error:\n{}'.format(crs_target_id, str_aux_error))
            return str_error, jacobian_matrix
        str_aux_error, crs_target_is_geocentric = self.is_geocentric(crs_target_id)
        if str_aux_error:
            str_error = CRSsTools.__name__ + "." + self.get_jacobian.__name__
            str_error += ('\nGetting CRS: {} is geocentric, error:\n{}'.format(crs_target_id, str_aux_error))
            return str_error, jacobian_matrix
        str_aux_error, crs_target_is_only_vertical = self.is_only_vertical(crs_target_id)
        if str_aux_error:
            str_error = CRSsTools.__name__ + "." + self.get_jacobian.__name__
            str_error += ('\nGetting CRS: {} is vertical, error:\n{}'.format(crs_target_id, str_aux_error))
            return str_error, jacobian_matrix
        jacobian_matrix = np.zeros((target_dimension, source_dimension))
        increment_lineal_2d_3d = cd.CRS_DIFFERENTIAL_2D_3D_LINEAL_INCREMENT
        increment_lineal_H = cd.CRS_DIFFERENTIAL_H_INCREMENT
        increment_deg = increment_lineal_2d_3d / cd.CRS_CONSTANT_JACOBI_SPHERE_RADIUS * 180. / np.pi
        # increment_rad = increment_lineal_2d_3d / cd.CRS_CONSTANT_JACOBI_SPHERE_RADIUS
        for pos_source_coordinate in range(source_dimension):
            inc_point = [source_point.copy()]
            increment_is_deg = False
            if crs_source_is_geographic or crs_source_is_only_vertical:
                if pos_source_coordinate < 2:
                    increment = increment_deg
                    increment_is_deg = True
                else:
                    increment = increment_lineal_H
            elif crs_source_is_geocentric:
                increment = increment_lineal_2d_3d
            else: # projected
                increment = increment_lineal_2d_3d
                if pos_source_coordinate == 2:
                    increment = increment_lineal_H
            inc_point[0][pos_source_coordinate] = inc_point[0][pos_source_coordinate] + increment
            str_error = self.operation(crs_source_id, crs_target_id, inc_point)
            inc_point = inc_point[0]
            if str_error:
                str_error = CRSsTools.__name__ + "." + self.get_jacobian.__name__
                str_error += ('\nIn operation from CRS: {} to CRS: {}, error:\n{}'.
                              format(crs_source_id, crs_target_id, str_error))
                return str_error, jacobian_matrix
            for pos_target_coordinate in range(target_dimension):
                coordinate_increment = inc_point[pos_target_coordinate] - target_point[pos_target_coordinate]
                partial = coordinate_increment / increment
                # if (crs_target_is_geographic or crs_target_is_only_vertical) and pos_target_coordinate < 2:
                #     if not increment_is_deg:
                #         increment_derivate = increment / cd.CRS_CONSTANT_JACOBI_SPHERE_RADIUS * 180. / np.pi
                #     else:
                #         increment_derivate = increment
                # else:
                #     if not increment_is_deg:
                #         increment_derivate = increment
                #     else:
                #         increment_derivate = increment * np.pi / 180. * cd.CRS_CONSTANT_JACOBI_SPHERE_RADIUS
                #         # if pos_target_coordinate < 2:
                #         #     increment_derivate = increment_lineal_2d_3d
                #         # else:
                #         #     increment_derivate = increment_lineal_H
                # partial = coordinate_increment / increment_derivate
                jacobian_matrix[pos_target_coordinate][pos_source_coordinate] = partial
        return str_error, jacobian_matrix

    def initialize(self):
        self.json_data_file = os.path.dirname(__file__) + "\\" + cd.JSON_DATA_FILE_BASE_NAME
        self.json_data_file += self.proj_version + ".json"
        self.json_data_file = os.path.normcase(self.json_data_file)
        if os.path.exists(self.json_data_file):
            with open(self.json_data_file, 'r') as data_file:
                self.data = json.load(data_file)
            data_is_valid = True
            if (not "CRSs_geo2d_ids_by_datum" in self.data
                    or not "CRSs_geo3d_ids_by_datum" in self.data
                    or not "CRSs_ecef_ids_by_datum" in self.data
                    or not "CRSs_projected_ids_by_datum" in self.data
                    or not "CRSs_geo2d_ids_by_base_crs_id" in self.data
                    or not "CRSs_geo3d_ids_by_base_crs_id" in self.data
                    or not "CRSs_ecef_ids_by_base_crs_id" in self.data
                    or not "CRSs_projected_ids_by_base_crs_id" in self.data
                    or not "CRSs_vertical_ids_by_base_crs_id" in self.data):
                data_is_valid = False
            if data_is_valid:
                return

        # Geo2D
        crs_info_list = database.query_crs_info(auth_name = 'EPSG',
                                                pj_types = "GEOGRAPHIC_2D_CRS")
        crs_list = []
        for info in crs_info_list:
            crs_id = "EPSG:" + info[1]
            crsInfo = {}
            crsInfo["auth_name"] = info.auth_name
            crsInfo["name"] = info.name
            crsInfo["code"] = info.code
            crsInfo["type"] = info.type.name
            crsInfo["deprecated"] = "False"
            if info.deprecated:
                crsInfo["deprecated"] = "True"
            crsInfo["area_of_use"] = {}
            crsInfo["area_of_use"]["name"] = info.area_of_use.name
            crsInfo["area_of_use"]["north"] = str(info.area_of_use.north)
            crsInfo["area_of_use"]["south"] = str(info.area_of_use.south)
            crsInfo["area_of_use"]["west"] = str(info.area_of_use.west)
            crsInfo["area_of_use"]["east"] = str(info.area_of_use.east)
            self.data['CRS_info_by_id'][crs_id] = crsInfo
            crs = CRS(crs_id)
            base_crs = crs.geodetic_crs
            if base_crs:
                base_crs_epsg = base_crs.to_authority('EPSG')
                if base_crs_epsg:
                    base_crs_id = "EPSG:" + base_crs_epsg[1]
                    self.data["base_crs_id_by_crs_id"][crs_id] = base_crs_id
                    if not base_crs_id in self.data["CRSs_geo2d_ids_by_base_crs_id"]:
                        self.data["CRSs_geo2d_ids_by_base_crs_id"][base_crs_id] = []
                    if not crs_id in self.data["CRSs_geo2d_ids_by_base_crs_id"][base_crs_id]:
                        self.data["CRSs_geo2d_ids_by_base_crs_id"][base_crs_id].append(crs_id)
            datum = crs.datum
            if datum:
                datum_name = datum.name
                if datum_name:
                    self.data["datum_by_crs_id"][crs_id] = datum_name
                    if not datum_name in self.data["CRSs_geo2d_ids_by_datum"]:
                        self.data["CRSs_geo2d_ids_by_datum"][datum_name] = []
                    if not crs_id in self.data["CRSs_geo2d_ids_by_datum"][datum_name]:
                        self.data["CRSs_geo2d_ids_by_datum"][datum_name].append(crs_id)
            crs_list.append(crs_id)
        # crs_list = ["EPSG:" + info[1] for info in crs_info_list]
        self.data["CRSs_geo2d_ids"] = sorted(crs_list)

        # # Geo3D
        crs_info_list = database.query_crs_info(auth_name = 'EPSG',
                                                pj_types = "GEOGRAPHIC_3D_CRS")
        crs_list = []
        for info in crs_info_list:
            crs_id = "EPSG:" + info[1]
            crsInfo = {}
            crsInfo["auth_name"] = info.auth_name
            crsInfo["name"] = info.name
            crsInfo["code"] = info.code
            crsInfo["type"] = info.type.name
            crsInfo["deprecated"] = "False"
            if info.deprecated:
                crsInfo["deprecated"] = "True"
            crsInfo["area_of_use"] = {}
            crsInfo["area_of_use"]["name"] = info.area_of_use.name
            crsInfo["area_of_use"]["north"] = str(info.area_of_use.north)
            crsInfo["area_of_use"]["south"] = str(info.area_of_use.south)
            crsInfo["area_of_use"]["west"] = str(info.area_of_use.west)
            crsInfo["area_of_use"]["east"] = str(info.area_of_use.east)
            self.data['CRS_info_by_id'][crs_id] = crsInfo
            crs = CRS(crs_id)
            base_crs = crs.geodetic_crs
            if base_crs:
                base_crs_epsg = base_crs.to_authority('EPSG')
                if base_crs_epsg:
                    base_crs_id = "EPSG:" + base_crs_epsg[1]
                    self.data["base_crs_id_by_crs_id"][crs_id] = base_crs_id
                    if not base_crs_id in self.data["CRSs_geo3d_ids_by_base_crs_id"]:
                        self.data["CRSs_geo3d_ids_by_base_crs_id"][base_crs_id] = []
                    if not crs_id in self.data["CRSs_geo3d_ids_by_base_crs_id"][base_crs_id]:
                        self.data["CRSs_geo3d_ids_by_base_crs_id"][base_crs_id].append(crs_id)
            datum = crs.datum
            if datum:
                datum_name = datum.name
                if datum_name:
                    self.data["datum_by_crs_id"][crs_id] = datum_name
                    if not datum_name in self.data["CRSs_geo3d_ids_by_datum"]:
                        self.data["CRSs_geo3d_ids_by_datum"][datum_name] = []
                    if not crs_id in self.data["CRSs_geo3d_ids_by_datum"][datum_name]:
                        self.data["CRSs_geo3d_ids_by_datum"][datum_name].append(crs_id)
            crs_list.append(crs_id)
        # crs_list = ["EPSG:" + info[1] for info in crs_info_list]
        self.data["CRSs_geo3d_ids"] = sorted(crs_list)

        # ECEF
        crs_info_list = database.query_crs_info(auth_name = 'EPSG',
                                                pj_types = "GEOCENTRIC_CRS")
        crs_list = []
        for info in crs_info_list:
            crs_id = "EPSG:" + info[1]
            crsInfo = {}
            crsInfo["auth_name"] = info.auth_name
            crsInfo["name"] = info.name
            crsInfo["code"] = info.code
            crsInfo["type"] = info.type.name
            crsInfo["deprecated"] = "False"
            if info.deprecated:
                crsInfo["deprecated"] = "True"
            crsInfo["area_of_use"] = {}
            crsInfo["area_of_use"]["name"] = info.area_of_use.name
            crsInfo["area_of_use"]["north"] = str(info.area_of_use.north)
            crsInfo["area_of_use"]["south"] = str(info.area_of_use.south)
            crsInfo["area_of_use"]["west"] = str(info.area_of_use.west)
            crsInfo["area_of_use"]["east"] = str(info.area_of_use.east)
            self.data['CRS_info_by_id'][crs_id] = crsInfo
            crs = CRS(crs_id)
            base_crs = crs.geodetic_crs
            if base_crs:
                base_crs_epsg = base_crs.to_authority('EPSG')
                if base_crs_epsg:
                    base_crs_id = "EPSG:" + base_crs_epsg[1]
                    self.data["base_crs_id_by_crs_id"][crs_id] = base_crs_id
                    if not base_crs_id in self.data["CRSs_ecef_ids_by_base_crs_id"]:
                        self.data["CRSs_ecef_ids_by_base_crs_id"][base_crs_id] = []
                    if not crs_id in self.data["CRSs_ecef_ids_by_base_crs_id"][base_crs_id]:
                        self.data["CRSs_ecef_ids_by_base_crs_id"][base_crs_id].append(crs_id)
            datum = crs.datum
            if datum:
                datum_name = datum.name
                if datum_name:
                    self.data["datum_by_crs_id"][crs_id] = datum_name
                    if not datum_name in self.data["CRSs_ecef_ids_by_datum"]:
                        self.data["CRSs_ecef_ids_by_datum"][datum_name] = []
                    if not crs_id in self.data["CRSs_ecef_ids_by_datum"][datum_name]:
                        self.data["CRSs_ecef_ids_by_datum"][datum_name].append(crs_id)
            crs_list.append(crs_id)
        # crs_list = ["EPSG:" + info[1] for info in crs_info_list]
        self.data["CRSs_ecef_ids"] = sorted(crs_list)

        # Projected
        crs_info_list = database.query_crs_info(auth_name = 'EPSG',
                                                pj_types = "PROJECTED_CRS")
        crs_list = []
        for info in crs_info_list:
            crs_id = "EPSG:" + info[1]
            crsInfo = {}
            crsInfo["auth_name"] = info.auth_name
            crsInfo["name"] = info.name
            crsInfo["code"] = info.code
            crsInfo["type"] = info.type.name
            crsInfo["deprecated"] = "False"
            if info.deprecated:
                crsInfo["deprecated"] = "True"
            crsInfo["area_of_use"] = {}
            crsInfo["area_of_use"]["name"] = info.area_of_use.name
            crsInfo["area_of_use"]["north"] = str(info.area_of_use.north)
            crsInfo["area_of_use"]["south"] = str(info.area_of_use.south)
            crsInfo["area_of_use"]["west"] = str(info.area_of_use.west)
            crsInfo["area_of_use"]["east"] = str(info.area_of_use.east)
            self.data['CRS_info_by_id'][crs_id] = crsInfo
            crs = CRS(crs_id)
            base_crs = crs.geodetic_crs
            if base_crs:
                base_crs_epsg = base_crs.to_authority('EPSG')
                if base_crs_epsg:
                    base_crs_id = "EPSG:" + base_crs_epsg[1]
                    self.data["base_crs_id_by_crs_id"][crs_id] = base_crs_id
                    if not base_crs_id in self.data["CRSs_projected_ids_by_base_crs_id"]:
                        self.data["CRSs_projected_ids_by_base_crs_id"][base_crs_id] = []
                    if not crs_id in self.data["CRSs_projected_ids_by_base_crs_id"][base_crs_id]:
                        self.data["CRSs_projected_ids_by_base_crs_id"][base_crs_id].append(crs_id)
            datum = crs.datum
            if datum:
                datum_name = datum.name
                if datum_name:
                    self.data["datum_by_crs_id"][crs_id] = datum_name
                    if not datum_name in self.data["CRSs_projected_ids_by_datum"]:
                        self.data["CRSs_projected_ids_by_datum"][datum_name] = []
                    if not crs_id in self.data["CRSs_projected_ids_by_datum"][datum_name]:
                        self.data["CRSs_projected_ids_by_datum"][datum_name].append(crs_id)
            crs_list.append(crs_id)
        # crs_list = ["EPSG:" + info[1] for info in crs_info_list]
        self.data["CRSs_projected_ids"] = sorted(crs_list)

        # Vertical
        crs_info_list = database.query_crs_info(auth_name = 'EPSG',
                                                pj_types = "VERTICAL_CRS")
        crs_list = []
        for info in crs_info_list:
            crs_id = "EPSG:" + info[1]
            crsInfo = {}
            crsInfo["auth_name"] = info.auth_name
            crsInfo["name"] = info.name
            crsInfo["code"] = info.code
            crsInfo["type"] = info.type.name
            crsInfo["deprecated"] = "False"
            if info.deprecated:
                crsInfo["deprecated"] = "True"
            crsInfo["area_of_use"] = {}
            crsInfo["area_of_use"]["name"] = info.area_of_use.name
            crsInfo["area_of_use"]["north"] = str(info.area_of_use.north)
            crsInfo["area_of_use"]["south"] = str(info.area_of_use.south)
            crsInfo["area_of_use"]["west"] = str(info.area_of_use.west)
            crsInfo["area_of_use"]["east"] = str(info.area_of_use.east)
            self.data['CRS_info_by_id'][crs_id] = crsInfo
            # crs = CRS(crs_id)
            crs_list.append(crs_id)
        # crs_list = ["EPSG:" + info[1] for info in crs_info_list]
        self.data["CRSs_vertical_ids"] = sorted(crs_list)

        # Vertical from Compound
        # crs_info_list = database.query_crs_info(auth_name='EPSG',
        #                                         pj_types="VERTICAL_CRS")
        crs_info_list = database.query_crs_info(auth_name = 'EPSG',
                                                pj_types = "COMPOUND_CRS")
        crs_list = []
        for info in crs_info_list:
            crs_id = "EPSG:" + info[1]
            crsInfo = {}
            crsInfo["auth_name"] = info.auth_name
            crsInfo["name"] = info.name
            crsInfo["code"] = info.code
            crsInfo["type"] = info.type.name
            crsInfo["deprecated"] = "False"
            if info.deprecated:
                crsInfo["deprecated"] = "True"
            crsInfo["area_of_use"] = {}
            crsInfo["area_of_use"]["name"] = info.area_of_use.name
            crsInfo["area_of_use"]["north"] = str(info.area_of_use.north)
            crsInfo["area_of_use"]["south"] = str(info.area_of_use.south)
            crsInfo["area_of_use"]["west"] = str(info.area_of_use.west)
            crsInfo["area_of_use"]["east"] = str(info.area_of_use.east)
            self.data['CRS_info_by_id'][crs_id] = crsInfo
            crs = CRS(crs_id)
            crs_list.append(crs_id)
            crss_list = crs.sub_crs_list
            crs_2d_id = None
            crs_v_id = None
            if crss_list:
                if len(crss_list) > 0:
                    crs_2d = crss_list[0]
                    crs_2d_epsg_code = crs_2d.to_epsg()
                    if crs_2d_epsg_code:
                        crs_2d_id = "EPSG:" + str(crs_2d_epsg_code)
                if len(crss_list) > 1:
                    crs_v = crss_list[1]
                    crs_v_epsg_code = crs_v.to_epsg()
                    if crs_v_epsg_code:
                        crs_v_id = "EPSG:" + str(crs_v_epsg_code)
            if not crs_2d_id or not crs_v_id:
                continue
            if crs_2d_id in self.data["base_crs_id_by_crs_id"]:
                base_crs_id = self.data["base_crs_id_by_crs_id"][crs_2d_id]
                self.data["base_crs_id_by_crs_id"][crs_v_id] = base_crs_id
                if not base_crs_id in self.data["CRSs_vertical_ids_by_base_crs_id"]:
                    self.data["CRSs_vertical_ids_by_base_crs_id"][base_crs_id] = []
                if not crs_v_id in self.data["CRSs_vertical_ids_by_base_crs_id"][base_crs_id]:
                    self.data["CRSs_vertical_ids_by_base_crs_id"][base_crs_id].append(crs_v_id)
        self.data["CRSs_compound_ids"] = sorted(crs_list)

        # sort list in containers
        for datum_id in self.data["CRSs_geo2d_ids_by_datum"]:
            values = self.data["CRSs_geo2d_ids_by_datum"][datum_id]
            self.data["CRSs_geo2d_ids_by_datum"][datum_id] = sorted(values)
        for datum_id in self.data["CRSs_geo3d_ids_by_datum"]:
            values = self.data["CRSs_geo3d_ids_by_datum"][datum_id]
            self.data["CRSs_geo3d_ids_by_datum"][datum_id] = sorted(values)
        for datum_id in self.data["CRSs_ecef_ids_by_datum"]:
            values = self.data["CRSs_ecef_ids_by_datum"][datum_id]
            self.data["CRSs_ecef_ids_by_datum"][datum_id] = sorted(values)
        for datum_id in self.data["CRSs_projected_ids_by_datum"]:
            values = self.data["CRSs_projected_ids_by_datum"][datum_id]
            self.data["CRSs_projected_ids_by_datum"][datum_id] = sorted(values)
        for crs_id in self.data["CRSs_geo2d_ids_by_base_crs_id"]:
            values = self.data["CRSs_geo2d_ids_by_base_crs_id"][crs_id]
            self.data["CRSs_geo2d_ids_by_base_crs_id"][crs_id] = sorted(values)
        for crs_id in self.data["CRSs_geo3d_ids_by_base_crs_id"]:
            values = self.data["CRSs_geo3d_ids_by_base_crs_id"][crs_id]
            self.data["CRSs_geo3d_ids_by_base_crs_id"][crs_id] = sorted(values)
        for crs_id in self.data["CRSs_ecef_ids_by_base_crs_id"]:
            values = self.data["CRSs_ecef_ids_by_base_crs_id"][crs_id]
            self.data["CRSs_ecef_ids_by_base_crs_id"][crs_id] = sorted(values)
        for crs_id in self.data["CRSs_projected_ids_by_base_crs_id"]:
            values = self.data["CRSs_projected_ids_by_base_crs_id"][crs_id]
            self.data["CRSs_projected_ids_by_base_crs_id"][crs_id] = sorted(values)
        for crs_id in self.data["CRSs_vertical_ids_by_base_crs_id"]:
            values = self.data["CRSs_vertical_ids_by_base_crs_id"][crs_id]
            self.data["CRSs_vertical_ids_by_base_crs_id"][crs_id] = sorted(values)

        with open(self.json_data_file, 'w') as data_file:
            data_file.write(json.dumps(self.data))

    def is_3d(self,
              crs_id):
        str_error = ""
        is_3d = False
        if crs_id in self.CRS_geo3d_id_by_enu_id:
            return str_error, is_3d
        if crs_id == cd.VERTICAL_ELLIPSOID_TAG:
            return str_error, is_3d
        if not crs_id in self.CRSs:
            try:
                crs_source = CRS(crs_id).to_3d()
            except:
                str_error = CRSsTools.__name__ + "." + self.is_3d.__name__
                str_error += f'\nFor CRS: {crs_id}'
                str_error += f"\nError making CRS"
                return str_error, is_3d
            self.CRSs[crs_id] = crs_source
        if (crs_id in self.data['CRSs_ecef_ids']
                or crs_id in self.data['CRSs_geo3d_ids']
                or crs_id in self.data['CRSs_compound_ids']):
            is_3d = True
        return str_error, is_3d

    def is_compound(self,
                         crs_id):
        str_error = ""
        is_compound = False
        if crs_id in self.CRS_geo3d_id_by_enu_id or cd.ENU_TAG in crs_id:
            return str_error, is_compound
        if crs_id == cd.VERTICAL_ELLIPSOID_TAG:
            return str_error, is_compound
        if not crs_id in self.CRSs:
            try:
                crs_source = CRS(crs_id).to_3d()
            except:
                str_error = CRSsTools.__name__ + "." + self.is_only_vertical.__name__
                str_error += f'\nFor CRS: {crs_id}'
                str_error += f"\nError making CRS"
                return str_error, is_compound
            self.CRSs[crs_id] = crs_source
        is_compound = self.CRSs[crs_id].is_compound
        return str_error, is_compound

    def is_geocentric(self,
                      crs_id):
        str_error = ""
        is_geocentric = False
        if crs_id in self.CRS_geo3d_id_by_enu_id or cd.ENU_TAG in crs_id:
            return str_error, is_geocentric
        if crs_id == cd.VERTICAL_ELLIPSOID_TAG:
            return str_error, is_geocentric
        if not crs_id in self.CRSs:
            try:
                crs_source = CRS(crs_id).to_3d()
            except:
                str_error = CRSsTools.__name__ + "." + self.is_geocentric.__name__
                str_error += f'\nFor CRS: {crs_id}'
                str_error += f"\nError making CRS"
                return str_error, is_geocentric
            self.CRSs[crs_id] = crs_source
        is_geocentric = self.CRSs[crs_id].is_geocentric
        return str_error, is_geocentric

    def is_geographic(self,
                      crs_id):
        str_error = ""
        is_geographic = False
        if crs_id in self.CRS_geo3d_id_by_enu_id or cd.ENU_TAG in crs_id:
            return str_error, is_geographic
        if crs_id == cd.VERTICAL_ELLIPSOID_TAG:
            return str_error, is_geographic
        if not crs_id in self.CRSs:
            try:
                crs_source = CRS(crs_id).to_3d()
            except:
                str_error = CRSsTools.__name__ + "." + self.is_geographic.__name__
                str_error += f'\nFor CRS: {crs_id}'
                str_error += f"\nError making CRS"
                return str_error, is_geographic
            self.CRSs[crs_id] = crs_source
        is_geographic = self.CRSs[crs_id].is_geographic
        return str_error, is_geographic

    def is_projected(self,
                     crs_id):
        str_error = ""
        is_projected = False
        if crs_id in self.CRS_geo3d_id_by_enu_id or cd.ENU_TAG in crs_id:
            return str_error, is_projected
        if crs_id == cd.VERTICAL_ELLIPSOID_TAG:
            return str_error, is_projected
        if not crs_id in self.CRSs:
            try:
                crs_source = CRS(crs_id).to_3d()
            except:
                str_error = CRSsTools.__name__ + "." + self.is_projected.__name__
                str_error += f'\nFor CRS: {crs_id}'
                str_error += f"\nError making CRS"
                return str_error, is_projected
            self.CRSs[crs_id] = crs_source
        is_projected = self.CRSs[crs_id].is_projected
        return str_error, is_projected

    def is_only_vertical(self,
                         crs_id):
        str_error = ""
        is_vertical = False
        if crs_id in self.CRS_geo3d_id_by_enu_id or cd.ENU_TAG in crs_id:
            return str_error, is_vertical
        if crs_id == cd.VERTICAL_ELLIPSOID_TAG:
            return str_error, True
        if not crs_id in self.CRSs:
            try:
                crs_source = CRS(crs_id).to_3d()
            except:
                str_error = CRSsTools.__name__ + "." + self.is_only_vertical.__name__
                str_error += f'\nFor CRS: {crs_id}'
                str_error += f"\nError making CRS"
                return str_error, is_vertical
            self.CRSs[crs_id] = crs_source
        is_vertical = self.CRSs[crs_id].is_vertical
        if self.CRSs[crs_id].is_compound:
            is_vertical = False
        return str_error, is_vertical

    def is_vertical(self,
                    crs_id):
        str_error = ""
        is_vertical = False
        if crs_id in self.CRS_geo3d_id_by_enu_id:
            return str_error, is_vertical
        if crs_id == cd.VERTICAL_ELLIPSOID_TAG:
            return str_error, True
        if not crs_id in self.CRSs:
            try:
                crs_source = CRS(crs_id).to_3d()
            except:
                str_error = CRSsTools.__name__ + "." + self.is_vertical.__name__
                str_error += f'\nFor CRS: {crs_id}'
                str_error += f"\nError making CRS"
                return str_error, is_vertical
            self.CRSs[crs_id] = crs_source
        is_vertical = self.CRSs[crs_id].is_vertical
        return str_error, is_vertical

    def operation(self,
                  crs_source_id,
                  crs_target_id,
                  points):  # {id,[fc,sc,tc]} or [[],]
        str_error = ""
        points_is_list = False
        points_is_dict = False
        points_is_only_one_point = False
        if isinstance(points, dict):
            points_is_dict = True
        elif isinstance(points, list):
            points_is_list = True
            if not isinstance(points[0], list):  # only one point as list of coordinates
                str_error = CRSsTools.__name__ + "." + self.operation.__name__
                str_error += f"\nPoints must be a dictionary or a list of points"
                return str_error
        if not points_is_dict and not points_is_list:
            str_error = CRSsTools.__name__ + "." + self.operation.__name__
            str_error += f"\nPoints must be a list or a dictionary"
            return str_error
        if crs_source_id == crs_target_id:
            return str_error
        str_aux_error, crs_source_is_only_vertical = self.is_only_vertical(crs_source_id)
        if str_aux_error:
            str_error = CRSsTools.__name__ + "." + self.operation.__name__
            str_error += ('\nGetting CRS: {} is vertical, error:\n{}'.format(crs_source_id, str_aux_error))
            return str_error
        str_aux_error, crs_target_is_only_vertical = self.is_only_vertical(crs_target_id)
        if str_aux_error:
            str_error = CRSsTools.__name__ + "." + self.operation.__name__
            str_error += ('\nGetting CRS: {} is vertical, error:\n{}'.format(crs_target_id, str_aux_error))
            return str_error
        if crs_source_is_only_vertical and cd.CRS_OPERATION_FROM_OR_TO_VERTICAL_CRS_GETS_LATITUDE_LONGITUDE_AXIS_ORDER:
            if points_is_dict:
                for point_id in points:
                    fc = points[point_id][0]
                    sc = points[point_id][1]
                    tc = points[point_id][2]
                    points[point_id] = [sc, fc, tc]
            elif points_is_list:
                for np in range(len(points)):
                    fc = points[np][0]
                    sc = points[np][1]
                    tc = points[np][2]
                    points[np] = [sc, fc, tc]
        if not cd.ENU_TAG in crs_source_id:
            if not crs_source_id in self.CRSs:
                try:
                    crs_source = CRS(crs_source_id).to_3d()
                except:
                    str_error = CRSsTools.__name__ + "." + self.operation.__name__
                    str_error += f'\nFor CRS: {crs_source_id}'
                    str_error += f"\nError making CRS"
                    return str_error
                self.CRSs[crs_source_id] = crs_source
        if not cd.ENU_TAG in crs_target_id:
            if not crs_target_id in self.CRSs:
                try:
                    crs_target = CRS(crs_target_id).to_3d()
                except:
                    str_error = CRSsTools.__name__ + "." + self.operation.__name__
                    str_error += f'\nFor CRS: {crs_target_id}'
                    str_error += f"\nError making CRS"
                    return str_error
                self.CRSs[crs_target_id] = crs_target
        crs_source_enu_id = None
        if cd.ENU_TAG in crs_source_id:
            crs_source_enu_id = crs_source_id
            if not crs_source_enu_id in self.CRS_geo3d_id_by_enu_id:
                str_values = crs_source_enu_id.split(cd.CRS_STRING_SEPARATOR)
                if len(str_values) != 4:
                    str_error = CRSsTools.__name__ + "." + self.operation.__name__
                    str_error += f'\nFor ENU CRS: {crs_source_enu_id}'
                    str_error += (f"\nThere are no two strings separated"
                                  f" by {cd.CRS_STRING_SEPARATOR} in {crs_source_enu_id}")
                    return str_error
                crs_enu_geographic_base_id = str_values[0].replace(cd.ENU_TAG, cd.EPSG_TAG)
                str_enu_longitude = str_values[1]
                str_enu_latitude = str_values[2]
                str_enu_altitude = str_values[3]
                if not crs_enu_geographic_base_id in self.CRSs:
                    try:
                        crs_enu_geographic_base = CRS(crs_enu_geographic_base_id).to_3d()
                    except:
                        str_error = CRSsTools.__name__ + "." + self.operation.__name__
                        str_error += f'\nFor ENU CRS: {crs_source_enu_id}'
                        str_error += f"\nInvalid string {crs_enu_geographic_base_id} for CRS base"
                        return str_error
                    if not crs_enu_geographic_base.is_geographic:
                        str_error = CRSsTools.__name__ + "." + self.operation.__name__
                        str_error += f'\nFor ENU CRS: {crs_source_enu_id}'
                        str_error += f"\nCRS base: {crs_enu_geographic_base_id} is not geographic"
                        return str_error
                    self.CRSs[crs_enu_geographic_base_id] = crs_enu_geographic_base
                if not crs_source_enu_id in self.CRS_geo3d_id_by_enu_id:
                    self.CRS_geo3d_id_by_enu_id[crs_source_enu_id] = crs_enu_geographic_base_id
                crs_enu_geographic_base = self.CRSs[crs_enu_geographic_base_id]
                crs_enu_ellipsoid = crs_enu_geographic_base.get_geod()
                crs_enu_ellipsoid_string = crs_enu_ellipsoid.initstring
                crs_enu_id = (f"+proj=topocentric {crs_enu_ellipsoid_string} "
                              f"+lon_0={str_enu_longitude} +lat_0={str_enu_latitude} +h_0={str_enu_altitude}")
                crs_enu_geocentric_id = f"+proj=cart {crs_enu_ellipsoid_string}"
                if not crs_source_enu_id in self.CRSOperation_by_enu_id:
                    str_pipeline = f"+proj=pipeline +step {crs_enu_geocentric_id} +step {crs_enu_id}"
                    crs_operation_enu = Transformer.from_pipeline(str_pipeline)
                    self.CRSOperation_by_enu_id[crs_source_enu_id] = crs_operation_enu
            crs_source_id = self.CRS_geo3d_id_by_enu_id[crs_source_enu_id]
            crs_operation_enu = self.CRSOperation_by_enu_id[crs_source_enu_id]
            if points_is_dict:
                for point_id in points:
                    source_fc = points[point_id][0]
                    source_sc = points[point_id][1]
                    source_tc = points[point_id][2]
                    [target_fc, target_sc, target_tc] = crs_operation_enu.transform(source_fc, source_sc, source_tc,
                                                                                    direction=TransformDirection.INVERSE)
                    points[point_id] = [target_fc, target_sc, target_tc]
            elif points_is_list:
                for np in range(len(points)):
                    source_fc = points[np][0]
                    source_sc = points[np][1]
                    source_tc = points[np][2]
                    [target_fc, target_sc, target_tc] = crs_operation_enu.transform(source_fc, source_sc, source_tc,
                                                                                    direction=TransformDirection.INVERSE)
                    points[np] = [target_fc, target_sc, target_tc]
        if crs_source_id == crs_target_id:
            return str_error
        crs_target_enu_id = None
        if cd.ENU_TAG in crs_target_id:
            crs_target_enu_id = crs_target_id
            if not crs_target_enu_id in self.CRS_geo3d_id_by_enu_id:
                str_values = crs_target_enu_id.split(cd.CRS_STRING_SEPARATOR)
                if len(str_values) != 4:
                    str_error = CRSsTools.__name__ + "." + self.operation.__name__
                    str_error += f'\nFor ENU CRS: {crs_target_enu_id}'
                    str_error += (f"\nThere are no two strings separated"
                                  f" by {cd.CRS_STRING_SEPARATOR} in {crs_source_enu_id}")
                    return str_error
                crs_enu_geographic_base_id = str_values[0].replace(cd.ENU_TAG, cd.EPSG_TAG)
                str_enu_longitude = str_values[1]
                str_enu_latitude = str_values[2]
                str_enu_altitude = str_values[3]
                if not crs_enu_geographic_base_id in self.CRSs:
                    try:
                        crs_enu_geographic_base = CRS(crs_enu_geographic_base_id).to_3d()
                    except:
                        str_error = CRSsTools.__name__ + "." + self.operation.__name__
                        str_error += f'\nFor ENU CRS: {crs_target_enu_id}'
                        str_error += f"\nInvalid string {crs_enu_geographic_base_id} for CRS base"
                        return str_error
                    if not crs_enu_geographic_base.is_geographic:
                        str_error = CRSsTools.__name__ + "." + self.operation.__name__
                        str_error += f'\nFor ENU CRS: {crs_target_enu_id}'
                        str_error += f"\nCRS base: {crs_enu_geographic_base_id} is not geographic"
                        return str_error
                    self.CRSs[crs_enu_geographic_base_id] = crs_enu_geographic_base
                if not crs_target_enu_id in self.CRS_geo3d_id_by_enu_id:
                    self.CRS_geo3d_id_by_enu_id[crs_target_enu_id] = crs_enu_geographic_base_id
                crs_enu_geographic_base = self.CRSs[crs_enu_geographic_base_id]
                crs_enu_ellipsoid = crs_enu_geographic_base.get_geod()
                crs_enu_ellipsoid_string = crs_enu_ellipsoid.initstring
                crs_enu_id = (f"+proj=topocentric {crs_enu_ellipsoid_string} "
                              f"+lon_0={str_enu_longitude} +lat_0={str_enu_latitude} +h_0={str_enu_altitude}")
                crs_enu_geocentric_id = f"+proj=cart {crs_enu_ellipsoid_string}"
                if not crs_target_enu_id in self.CRSOperation_by_enu_id:
                    str_pipeline = f"+proj=pipeline +step {crs_enu_geocentric_id} +step {crs_enu_id}"
                    crs_operation_enu = Transformer.from_pipeline(str_pipeline)
                    self.CRSOperation_by_enu_id[crs_target_enu_id] = crs_operation_enu
            crs_target_id = self.CRS_geo3d_id_by_enu_id[crs_target_enu_id]
        if crs_source_id != crs_target_id:
            crs_source = None
            if not crs_source_id in self.CRSs:
                str_error = CRSsTools.__name__ + "." + self.operation.__name__
                str_error += ('\nNot find CRS: {}'.format(crs_source_id))
                return str_error
            crs_source = self.CRSs[crs_source_id]
            crs_target = None
            if not crs_target_id in self.CRSs:
                str_error = CRSsTools.__name__ + "." + self.operation.__name__
                str_error += ('\nNot find CRS: {}'.format(crs_target_id))
                return str_error
            crs_target = self.CRSs[crs_target_id]
            if not crs_source_id in self.CRSOperations:
                # crs_operation = Transformer.from_crs(crs_source_id, crs_target_id, always_xy=True)
                crs_operation = Transformer.from_crs(crs_source, crs_target, always_xy=True)
                self.CRSOperations[crs_source_id] = {}
                self.CRSOperations[crs_source_id][crs_target_id] = crs_operation
            elif not crs_target_id in self.CRSOperations[crs_source_id]:
                # crs_operation = Transformer.from_crs(crs_source_id, crs_target_id, always_xy=True)
                crs_operation = Transformer.from_crs(crs_source, crs_target, always_xy=True)
                self.CRSOperations[crs_source_id][crs_target_id] = crs_operation
            crs_operation = self.CRSOperations[crs_source_id][crs_target_id]
            if points_is_dict:
                for point_id in points:
                    source_fc = points[point_id][0]
                    source_sc = points[point_id][1]
                    source_tc = points[point_id][2]
                    [target_fc, target_sc, target_tc] = crs_operation.transform(source_fc, source_sc, source_tc)
                    points[point_id] = [target_fc, target_sc, target_tc]
            elif points_is_list:
                for np in range(len(points)):
                    source_fc = points[np][0]
                    source_sc = points[np][1]
                    source_tc = points[np][2]
                    [target_fc, target_sc, target_tc] = crs_operation.transform(source_fc, source_sc, source_tc)
                    points[np] = [target_fc, target_sc, target_tc]
        if crs_target_enu_id:
            crs_operation_enu = self.CRSOperation_by_enu_id[crs_target_enu_id]
            if points_is_dict:
                for point_id in points:
                    source_fc = points[point_id][0]
                    source_sc = points[point_id][1]
                    source_tc = points[point_id][2]
                    [target_fc, target_sc, target_tc] = crs_operation_enu.transform(source_fc, source_sc, source_tc,
                                                                                    direction=TransformDirection.FORWARD)
                    points[point_id] = [target_fc, target_sc, target_tc]
            elif points_is_list:
                for np in range(len(points)):
                    source_fc = points[np][0]
                    source_sc = points[np][1]
                    source_tc = points[np][2]
                    [target_fc, target_sc, target_tc] = crs_operation_enu.transform(source_fc, source_sc, source_tc,
                                                                                    direction=TransformDirection.FORWARD)
                    points[np] = [target_fc, target_sc, target_tc]
        if crs_target_is_only_vertical and cd.CRS_OPERATION_FROM_OR_TO_VERTICAL_CRS_GETS_LATITUDE_LONGITUDE_AXIS_ORDER:
            if points_is_dict:
                for point_id in points:
                    fc = points[point_id][0]
                    sc = points[point_id][1]
                    tc = points[point_id][2]
                    points[point_id] = [sc, fc, tc]
            elif points_is_list:
                for np in range(len(points)):
                    fc = points[np][0]
                    sc = points[np][1]
                    tc = points[np][2]
                    points[np] = [sc, fc, tc]
        return str_error
