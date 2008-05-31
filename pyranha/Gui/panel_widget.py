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
from ui_panel_widget import Ui_panel_widget
import IPython.ipapi
import pyranha

class panel_widget(PyQt4.QtGui.QWidget,Ui_panel_widget):
	class __series_list_model(PyQt4.QtCore.QAbstractListModel):
		def __init__(self,parent = None):
			PyQt4.QtCore.QAbstractListModel.__init__(self,parent)
			self.__series_list = []
			ip = IPython.ipapi.get()
			self.__series_list = sorted(filter(lambda x: type(ip.user_ns[x]) in pyranha.manipulators_type_tuple, ip.user_ns))
		def rowCount(self,model_index = PyQt4.QtCore.QModelIndex()):
			# print len(self.__series_list)
			return len(self.__series_list)
		def data(self,index,role):
			if not index.isValid() or index.row() >= len(self.__series_list):
				return PyQt4.QtCore.QVariant()
			if role == PyQt4.QtCore.Qt.DisplayRole:
				return PyQt4.QtCore.QVariant(self.__series_list[index.row()])
			else:
				return PyQt4.QtCore.QVariant()
	def __init__(self,parent = None):
		PyQt4.QtGui.QWidget.__init__(self,parent)
		self.setupUi(self)
		self.__timer = PyQt4.QtCore.QTimer()
		self.__timer.start(500)
		self.__series_list = self.__series_list_model(self)
		# Connections.
		self.connect(self.__timer,PyQt4.QtCore.SIGNAL("timeout()"),self.__global_update)
		# Set model.
		self.series_list_view.setModel(self.__series_list)
		self.show()
	def __global_update(self):
		if not self.isActiveWindow():
			del self.__series_list
			self.__series_list = self.__series_list_model(self)
			self.series_list_view.setModel(self.__series_list)
