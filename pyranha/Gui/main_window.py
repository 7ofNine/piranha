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

import PyQt4.QtCore, PyQt4.QtGui
from ui_main_window import Ui_main_window
import IPython.ipapi
import pyranha

class main_window(PyQt4.QtGui.QMainWindow,Ui_main_window):
	class __series_db_model(PyQt4.QtCore.QAbstractItemModel):
		def __init__(self,parent):
			PyQt4.QtCore.QAbstractItemModel.__init__(self,parent)
			# NOTE: this value must match the number of data headers and the number
			# of rows in each record of the database series.
			self.__n_columns = 4
			# This is the interactive name space of the IPython session.
			self.ip_ns = IPython.ipapi.get().user_ns
			self.__series_db = self.__build_series_db()
		def __build_series_db(self):
			retval = map(lambda x:
				(
					id(self.ip_ns[x]),
					x,
					str(type(self.ip_ns[x])).rpartition('.')[-1].strip('>\''),
					len(self.ip_ns[x])
				),filter(lambda x: type(self.ip_ns[x]) in pyranha.manipulators_type_tuple, self.ip_ns))
			assert(not retval or self.__n_columns == len(retval[0]))
			retval.sort()
			return retval
		def check_update(self):
			new_db = self.__build_series_db()
			if self.__series_db != new_db:
				self.__series_db = new_db
				self.reset()
		def hasChildren(self,model_index):
			return not model_index.isValid()
		def rowCount(self,model_index):
			return len(self.__series_db)
		def columnCount(self,model_index):
			return self.__n_columns
		def index(self,row,column,parent=PyQt4.QtCore.QModelIndex()):
			# We return a valid index only if parent is root item (i.e., it is invalid)
			# and if we are not going out of boundaries.
			if parent.isValid() or not self.hasIndex(row,column,parent):
				return PyQt4.QtCore.QModelIndex()
			else:
				return self.createIndex(row,column)
		def parent(self,index):
			return PyQt4.QtCore.QModelIndex()
		def data(self,index,role):
			# Return empty data if the requested index is not valid or we are using it for
			# something else than the display role.
			if not index.isValid() or role != PyQt4.QtCore.Qt.DisplayRole:
				return PyQt4.QtCore.QVariant()
			else:
				return PyQt4.QtCore.QVariant(self.__series_db[index.row()][index.column()])
		def headerData(self,column,orientation,role):
			if role != PyQt4.QtCore.Qt.DisplayRole:
				return PyQt4.QtCore.QVariant()
			if column == 0:
				return PyQt4.QtCore.QVariant("Id")
			if column == 1:
				return PyQt4.QtCore.QVariant("Name")
			if column == 2:
				return PyQt4.QtCore.QVariant("Type")
			if column == 3:
				return PyQt4.QtCore.QVariant("Length")
			assert(False)
	def __init__(self):
		PyQt4.QtGui.QMainWindow.__init__(self,None)
		# Setup various bits of the UI.
		self.setupUi(self)
		# Store original info panel caption.
		self.__default_info_panel_title = self.series_info_groupbox.title()
		self.series_atoms_label.setText("")
		# Store pointer to QApplication's instance.
		self.__qapp = PyQt4.QtGui.qApp
		# Setup the timer.
		self.__timer = PyQt4.QtCore.QTimer()
		self.__timer.start(500)
		# Create the database.
		self.__series_db = self.__series_db_model(self)
		# Set model.
		self.__setup_model()
		# Connections.
		self.connect(self.__timer,PyQt4.QtCore.SIGNAL("timeout()"),self.__slot_global_update)
		self.connect(self.__qapp,PyQt4.QtCore.SIGNAL("focusChanged(QWidget *,QWidget *)"),self.__slot_update_on_activation)
		self.connect(self.series_tree_view,PyQt4.QtCore.SIGNAL("doubleClicked(QModelIndex)"),self.__slot_series_tree_item_dclicked)
		self.show()
	def __setup_model(self):
		self.__proxy_model = PyQt4.QtGui.QSortFilterProxyModel(self)
		self.__proxy_model.setSourceModel(self.__series_db)
		self.series_tree_view.setModel(self.__proxy_model)
		# We do not want to show the series' id.
		#self.series_tree_view.hideColumn(0)
	def __info_panel_setup(self,series_name):
		if series_name == None:
			self.series_info_groupbox.setTitle(self.__default_info_panel_title)
			self.series_info_groupbox.setEnabled(False)
			self.series_atoms_label.setText("")
		else:
			self.series_info_groupbox.setTitle(self.__default_info_panel_title+": "+series_name)
			self.series_info_groupbox.setEnabled(True)
			self.series_atoms_label.setText(str(self.__series_db.ip_ns[series_name].atoms()))
	def __slot_series_tree_item_dclicked(self,index):
		# Fetch series' name from the index.
		name = str(index.sibling(index.row(),1).data().toString())
		self.series_graphics_view.set_series(name,self.__series_db.ip_ns[name])
		self.__info_panel_setup(name)
	def __slot_update_on_activation(self,old,now):
		if id(self.__qapp.activeWindow()) == id(self):
			self.__slot_global_update(force=True)
	def __slot_global_update(self,force=False):
		if force or not self.isActiveWindow():
			self.__series_db.check_update()
