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

class series_draw_area(PyQt4.QtGui.QGraphicsView):
	class __term_graphics_item(PyQt4.QtGui.QGraphicsRectItem):
		def __init__(self,x,y,w,h,pen,brush,parent=None):
			PyQt4.QtGui.QGraphicsRectItem.__init__(self,x,y,w,h,parent)
			self.setPen(pen)
			self.setBrush(brush)
			self.setAcceptHoverEvents(True)
		def hoverEnterEvent(self,event):
			self.__orig_pen = self.pen()
			self.setPen(PyQt4.QtGui.QPen(PyQt4.QtCore.Qt.yellow))
			self.setZValue(1)
		def hoverLeaveEvent(self,event):
			self.setPen(self.__orig_pen)
			del self.__orig_pen
			self.setZValue(0)
	def __init__(self,parent):
		PyQt4.QtGui.QGraphicsView.__init__(self,parent)
		self.__series_name = None
		self.__series = None
		self.__series_graphics_scene = PyQt4.QtGui.QGraphicsScene()
		self.setScene(self.__series_graphics_scene)
		self.setAlignment(PyQt4.QtCore.Qt.AlignLeft | PyQt4.QtCore.Qt.AlignBottom)
		# Store default pen.
		self.__pen = PyQt4.QtGui.QPen()
		# Default column width.
		self.__cwidth = 10
	def __populate_scene(self):
		if not self.__series or len(self.__series) == 0:
			return
		brush_plus = PyQt4.QtGui.QLinearGradient()
		brush_plus.setColorAt(0,PyQt4.QtCore.Qt.gray)
		brush_plus.setColorAt(1,PyQt4.QtCore.Qt.darkBlue)
		brush_minus = PyQt4.QtGui.QLinearGradient()
		brush_minus.setColorAt(0,PyQt4.QtCore.Qt.gray)
		brush_minus.setColorAt(1,PyQt4.QtCore.Qt.darkRed)
		top_value = None
		for i in self.__series.index0:
			top_value = abs(i.display_data())
			break
		brush = None
		n = 0
		index = self.__series.index0
		half_width = self.__cwidth/2.
		for i in index:
			value = i.display_data()
			abs_value = abs(value)
			if value >= 0:
				brush = brush_plus
			else:
				brush = brush_minus
			x = n*self.__cwidth
			x_brush = x + half_width
			brush.setStart(x_brush,0)
			brush.setFinalStop(x_brush,abs_value)
			self.__series_graphics_scene.addItem(self.__term_graphics_item(x,0,self.__cwidth,abs_value,self.__pen,brush))
			if abs_value > top_value:
				top_value = abs_value
			n = n + 1
		if n >= 50:
			x_scale_factor = float(self.width())/(self.__cwidth*50)
		else:
			x_scale_factor = float(self.width())/(self.__cwidth*n)
		if top_value != 0:
			y_scale_factor = self.height()/(top_value)
			y_scale_factor -= y_scale_factor/6.
			self.setMatrix(PyQt4.QtGui.QMatrix(x_scale_factor,0,0,-y_scale_factor,0,0))
	def set_series(self,name,series):
		self.__series_graphics_scene = PyQt4.QtGui.QGraphicsScene()
		self.setScene(self.__series_graphics_scene)
		self.__series = series
		self.__series_name = name
		self.__populate_scene()
	def series_name(self):
		return self.__series_name
	def series(self):
		return self.__series
