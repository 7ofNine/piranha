# -*- coding: iso-8859-1 -*-
# Copyright (C) 2007, 2008 by Francesco Biscani
# bluescarni@gmail.com
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the
# Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

import PySide2.QtCore
import PySide2.QtGui
import PySide2.QtWidgets
import pyranha
from .ui_main_window import Ui_main_window
from IPython import get_ipython

class main_window(PySide2.QtWidgets.QMainWindow, Ui_main_window):
    # The data model
    class __series_db_model(PySide2.QtCore.QAbstractItemModel):
        def __init__(self, parent):
            PySide2.QtCore.QAbstractItemModel.__init__(self, parent)
            # This is the interactive namespace of the IPython session.
            self.ip_ns = get_ipython().user_global_ns
            self.__headers = ("Id", "Name", "Type", "Length")
            self.__n_columns = len(self.__headers)
            self.__series_db = self.__build_series_db()
            
        def __build_series_db(self):
            retval = [(
                    id(self.ip_ns[x]),
                    x,
                    self.ip_ns[x].__short_type__,
                    len(self.ip_ns[x])
                ) for x in [x for x in self.ip_ns if type(self.ip_ns[x]) in pyranha.manipulators]]
            assert(not retval or self.__n_columns == len(retval[0]))
            retval.sort()
            return retval
            
        def __out_of_range(self,row,column):
            return row < 0 or column < 0 or row >= len(self.__series_db) or column >= self.__n_columns
            
        def check_update(self):
            new_db = self.__build_series_db()
            if self.__series_db != new_db:
                self.beginResetModel()
                self.__series_db = new_db
                self.endResetModel()
                
        def hasChildren(self,model_index):
            return not model_index.isValid()
            
        def rowCount(self,model_index):
            return len(self.__series_db)
            
        def columnCount(self,model_index):
            return self.__n_columns
            
        def index(self,row,column,parent=PySide2.QtCore.QModelIndex()):
            # We return a valid index only if parent is root item (i.e., it is invalid)
            # and if we are not going out of boundaries.
            if parent.isValid() or self.__out_of_range(row,column):
                return PySide2.QtCore.QModelIndex()
            else:
                return self.createIndex(row,column)
                
        def parent(self,index):
            return PySide2.QtCore.QModelIndex()
            
        def data(self,index,role):
            # Return empty data if the requested index is not valid or we are using it for
            # something else than the display role.
            if not index.isValid():
                return None
            if role == PySide2.QtCore.Qt.DisplayRole:
                return str(self.__series_db[index.row()][index.column()])
            if role == PySide2.QtCore.Qt.ToolTipRole and index.column() == 2:
                return str(self.ip_ns[self.__series_db[index.row()][1]].__doc__)
            return None
            
        def headerData(self,column,orientation,role):
            if role != PySide2.QtCore.Qt.DisplayRole:
                return None
            return str(self.__headers[column])
            
    def __init__(self):
        PySide2.QtWidgets.QMainWindow.__init__(self,None)
        self.setupUi(self)
        # Create the database.
        self.__series_db = self.__series_db_model(self)
        # Set model.
        self.__setup_model()
        # Setup various bits of the UI.
        self.__info_panel_setup(None)
        self.progress_bar = PySide2.QtWidgets.QProgressBar()
        self.progress_bar.setVisible(False)
        self.statusBar().addWidget(self.progress_bar)
        # Store pointer to QApplication's instance.
        self.__qapp = qApp       # qApp is global
        # Setup the timer.
        self.__timer = PySide2.QtCore.QTimer()
        self.__timer.start(500)
        # Connections.
        self.connect(self.__timer,PySide2.QtCore.SIGNAL("timeout()"),self.__slot_global_update)
        self.connect(self.__qapp,PySide2.QtCore.SIGNAL("focusChanged(QWidget *,QWidget *)"),self.__slot_update_on_activation)
        self.connect(self.series_tree_view,PySide2.QtCore.SIGNAL("activated(QModelIndex)"),self.__slot_series_tree_item_activated)
        self.show()
        
    def __setup_model(self):
        self.__proxy_model = PySide2.QtCore.QSortFilterProxyModel(self)
        self.__proxy_model.setSourceModel(self.__series_db)
        self.series_tree_view.setModel(self.__proxy_model)
        # We do not want to show the series' id.
        self.series_tree_view.hideColumn(0)
        
    def __info_panel_setup(self,series_name_):
        try:
            series_name = series_name_
            series = self.__series_db.ip_ns[series_name]
        except KeyError:
            series_name = None
        if series_name is None:
            self.series_info_groupbox.setEnabled(False)
            self.series_name_label.setText("")
            self.series_type_label.setText("")
            self.series_type_label.setToolTip("")
            self.series_length_label.setText("")
            self.series_atoms_label.setText("")
        else:
            self.series_info_groupbox.setEnabled(True)
            self.series_name_label.setText(series_name)
            self.series_type_label.setText(series.__short_type__)
            self.series_type_label.setToolTip(str(series.__doc__))
            self.series_length_label.setText(str(len(series)))
            self.series_atoms_label.setText(str(series.atoms))
            
    def __slot_series_tree_item_activated(self,index):
        # Fetch series' name from the index.
        name = str(index.sibling(index.row(),1).data())
        self.__info_panel_setup(name)
        
    def __slot_update_on_activation(self,old,now):
        if id(self.__qapp.activeWindow()) == id(self):
            self.__slot_global_update(force=True)
            self.__info_panel_setup(str(self.series_name_label.text()))
            
    def __slot_global_update(self,force=False):
        if force or not self.isActiveWindow():
            self.__series_db.check_update()
