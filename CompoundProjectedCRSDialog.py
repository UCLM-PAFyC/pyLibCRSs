# authors:
# David Hernandez Lopez, david.hernandez@uclm.es

import os
import sys
import math
import json


current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(current_path, '..'))

from .CRSsTools import CRSsTools
from . import CRSsDefines as cd

from pyLibQtTools import Tools

from PyQt5 import QtCore, QtWidgets
from PyQt5.uic import loadUi
from PyQt5.QtWidgets import (QApplication, QMessageBox, QDialog, QTreeWidgetItem,
                             QFileDialog, QPushButton, QComboBox, QPlainTextEdit, QLineEdit,
                             QDialogButtonBox, QVBoxLayout, QTableWidget, QTableWidgetItem, QInputDialog)
from PyQt5.QtCore import QDir, QFileInfo, QFile, QSize, Qt, QDate

class CompoundProjectedCRSDialog(QDialog):
    """Employee dialog."""

    def __init__(self,
                 crs_tools,
                 crs_id,
                 parent=None):
        super().__init__(parent)
        loadUi(os.path.join(os.path.dirname(__file__), 'CompoundProjectedCRSDialog.ui'), self)
        self.crs_tools = crs_tools
        self.crs_projected_ids = []
        self.crs_vertical_ids = []
        self.active_crs_line_edit_widget = None
        self.crs_id = None
        self.is_accepted = False
        self.initialize(crs_id)

    def initialize(self, crs_id):
        crs_projected_id = ''
        crs_vertical_id = ''
        epsg_crs_prefix = cd.EPSG_TAG + ':'
        if epsg_crs_prefix in crs_id:
            crs_id = crs_id.replace(epsg_crs_prefix, '')
            if '+' in crs_id:
                epsg_codes_str_values = crs_id.split('+')
                if len(epsg_codes_str_values) == 2:
                    crs_projected_id = epsg_crs_prefix + str(int(epsg_codes_str_values[0]))
                    crs_vertical_id = epsg_crs_prefix + str(int(epsg_codes_str_values[1]))
            else:
                crs_projected_id = epsg_crs_prefix + str(int(crs_id))
                crs_vertical_id = cd.VERTICAL_ELLIPSOID_TAG
        if not crs_projected_id:
            crs_projected_id = cd.CRS_PROJECTED_DEFAULT
        if not crs_vertical_id:
            crs_vertical_id = cd.CRS_VERTICAL_DEFAULT
        self.crsProjectedLineEdit.setText(crs_projected_id)
        self.crsVerticalLineEdit.setText(crs_vertical_id)
        self.crs_projected_ids = self.crs_tools.get_crs_projected_ids()
        self.crs_vertical_ids = self.crs_tools.get_crs_vertical_ids()
        self.crsSearchTextLineEdit.cursorPositionChanged.connect(self.crsSearchTextChanged)
        self.acceptPushButton.clicked.connect(self.accept)

        self.crsTreeWidget.clear()
        self.crsTextEdit.clear()
        self.crsTreeWidget.setColumnCount(1)
        # self.crsTreeWidget.header().hide()
        item_header = self.crsTreeWidget.headerItem()
        item_header.setText(0, 'Select CRS by type')
        self.crs_projected_item = QTreeWidgetItem(self.crsTreeWidget)
        self.crs_projected_item.setText(0, cd.CRS_PROJECTED_LABEL)
        self.crs_projected_item.setFlags(self.crs_projected_item.flags() & ~QtCore.Qt.ItemIsEditable
                                         & ~QtCore.Qt.ItemIsSelectable)
        for crs_id in self.crs_projected_ids:
            str_error, crs_summary = self.crs_tools.get_crs_summary(crs_id)
            if str_error:
                str_error = ('Getting summary for CRS: {}, error:\n{}'.format(crs_id, str_error))
                Tools.error_msg(str_error)
                return
            item = QTreeWidgetItem(self.crs_projected_item)
            item.setText(0, crs_summary)
        self.crs_vertical_item = QTreeWidgetItem(self.crsTreeWidget)
        self.crs_vertical_item.setText(0, cd.CRS_VERTICAL_LABEL)
        self.crs_vertical_item.setFlags(self.crs_vertical_item.flags() & ~QtCore.Qt.ItemIsEditable
                                        & ~QtCore.Qt.ItemIsSelectable)
        item = QTreeWidgetItem(self.crs_vertical_item)
        item.setText(0, cd.VERTICAL_ELLIPSOID_TAG)
        for crs_id in self.crs_vertical_ids:
            str_error, crs_summary = self.crs_tools.get_crs_summary(crs_id)
            if str_error:
                str_error = ('Getting summary for CRS: {}, error:\n{}'.format(crs_id, str_error))
                Tools.error_msg(str_error)
                return
            item = QTreeWidgetItem(self.crs_vertical_item)
            item.setText(0, crs_summary)
        self.crsTreeWidget.itemSelectionChanged.connect(self.show_crs_details)
        self.crsTreeWidget.expandAll()

    def show_crs_details(self):
        selectedItems = self.crsTreeWidget.selectedItems()
        self.crsTextEdit.clear()
        if len(selectedItems) > 0:
            selected_item = selectedItems[0]
            selected_text = selected_item.text(0)
            if selected_text == cd.VERTICAL_ELLIPSOID_TAG:
                self.crsVerticalLineEdit.setText(selected_text)
                return
            crs_id = selected_text[0: selected_text.index(',')]
            str_error, crs_info_as_dict = self.crs_tools.get_crs_info_as_text(crs_id)
            if str_error:
                str_error = ('Getting info for CRS: {}, error:\n{}'.format(crs_id, str_error))
                Tools.error_msg(str_error)
                return
            self.crsTextEdit.setText(crs_info_as_dict)
            if crs_id in self.crs_projected_ids:
                self.crsProjectedLineEdit.setText(crs_id)
            elif crs_id in self.crs_vertical_ids:
                self.crsVerticalLineEdit.setText(crs_id)
        return

    def crsSearchTextChanged(self):
        text = self.crsSearchTextLineEdit.text()
        self.update_crs_tree(text)
        return

    def accept(self):
        crs_projected_id = self.crsProjectedLineEdit.text()
        if not crs_projected_id:
            str_error = ('Select Projected CRS')
            Tools.error_msg(str_error)
            return
        crs_vertical_id = self.crsVerticalLineEdit.text()
        if not crs_vertical_id:
            str_error = ('Select Vertical CRS')
            Tools.error_msg(str_error)
            return
        epsg_crs_prefix = cd.EPSG_TAG + ':'
        crs_2d_epsg_code = int(crs_projected_id.replace(epsg_crs_prefix, ''))
        self.crs_id = epsg_crs_prefix + str(crs_2d_epsg_code)
        if crs_vertical_id != cd.VERTICAL_ELLIPSOID_TAG:
            crs_vertical_epsg_code = int(crs_vertical_id.replace(epsg_crs_prefix, ''))
            self.crs_id += ('+' + str(crs_vertical_epsg_code))
        self.is_accepted = True
        super().accept()
