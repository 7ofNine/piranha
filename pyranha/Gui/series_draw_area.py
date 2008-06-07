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
	class __populate_scene_thread(PyQt4.QtCore.QThread):
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
		def __init__(self,scene,series,cwidth,parent):
			PyQt4.QtCore.QThread.__init__(self,parent)
			self.__cwidth = cwidth
			self.__series = series
			self.__series_graphics_scene = scene
		def run(self):
			brush_plus = PyQt4.QtGui.QLinearGradient()
			brush_plus.setColorAt(0,PyQt4.QtCore.Qt.gray)
			brush_plus.setColorAt(1,PyQt4.QtCore.Qt.darkBlue)
			brush_minus = PyQt4.QtGui.QLinearGradient()
			brush_minus.setColorAt(0,PyQt4.QtCore.Qt.gray)
			brush_minus.setColorAt(1,PyQt4.QtCore.Qt.darkRed)
			pen = PyQt4.QtGui.QPen()
			self.__n = 0
			self.__top_value = None
			brush = None
			half_width = self.__cwidth/2.
			for i in self.__series.index0:
				value = i.display_data()
				abs_value = abs(value)
				if value >= 0:
					brush = brush_plus
				else:
					brush = brush_minus
				x = self.__n * self.__cwidth
				x_brush = x + half_width
				brush.setStart(x_brush,0)
				brush.setFinalStop(x_brush,abs_value)
				self.__series_graphics_scene.addItem(self.__term_graphics_item(x,0,self.__cwidth,abs_value,pen,brush))
				if abs_value > self.__top_value:
					self.__top_value = abs_value
				# Report progress every 50 items added.
				if not self.__n % 50:
					PyQt4.QtCore.QObject.emit(self,PyQt4.QtCore.SIGNAL("step(int)"),self.__n)
				self.__n += 1
			self.__series_graphics_scene.moveToThread(PyQt4.QtCore.QCoreApplication.instance().thread())
		def n(self):
			return self.__n
		def top_value(self):
			return self.__top_value
	def __init__(self,parent):
		PyQt4.QtGui.QGraphicsView.__init__(self,parent)
		self.__series_name = None
		self.__series = None
		self.__series_graphics_scene = None
		self.__progress_bar = None
		self.__mutex = PyQt4.QtCore.QMutex()
		self.setAlignment(PyQt4.QtCore.Qt.AlignLeft | PyQt4.QtCore.Qt.AlignBottom)
		# Default column width.
		self.__cwidth = 10
		self.__populator_thread = None
		self.__set_busy(False)
	def __slot_population_finished(self):
		n = self.__populator_thread.n()
		top_value = self.__populator_thread.top_value()
		if n >= 50:
			x_scale_factor = float(self.width())/(self.__cwidth*50)
		else:
			x_scale_factor = float(self.width())/(self.__cwidth*n)
		if top_value != 0:
			y_scale_factor = self.height()/(top_value)
			y_scale_factor -= y_scale_factor/6.
			self.setMatrix(PyQt4.QtGui.QMatrix(x_scale_factor,0,0,-y_scale_factor,0,0))
			self.setScene(self.__series_graphics_scene)
		self.disconnect(self.__populator_thread,PyQt4.QtCore.SIGNAL("finished()"),self.__slot_population_finished)
		self.disconnect(self.__populator_thread,PyQt4.QtCore.SIGNAL("step(int)"),self.__progress_bar.setValue)
		self.__progress_bar.setVisible(False)
		self.__set_busy(False)
	def __populate_scene(self):
		if not self.__series or len(self.__series) == 0:
			return
		self.setScene(None)
		self.__series_graphics_scene = PyQt4.QtGui.QGraphicsScene()
		self.__populator_thread = self.__populate_scene_thread(self.__series_graphics_scene,self.__series,
			self.__cwidth,self)
		self.__series_graphics_scene.moveToThread(self.__populator_thread)
		main_window = self.__get_main_window()
		self.__progress_bar = main_window.progress_bar
		self.__progress_bar.setMinimum(0)
		self.__progress_bar.setMaximum(len(self.__series))
		self.__progress_bar.setVisible(True)
		self.connect(self.__populator_thread,PyQt4.QtCore.SIGNAL("finished()"),self.__slot_population_finished)
		self.connect(self.__populator_thread,PyQt4.QtCore.SIGNAL("step(int)"),self.__progress_bar.setValue)
		self.__set_busy(True)
		self.__populator_thread.start()
	def set_series(self,name,series):
		if self.__is_busy():
			print "thread is running, won't do"
			return
		self.__series = series
		self.__series_name = name
		self.__populate_scene()
	def series_name(self):
		return self.__series_name
	def series(self):
		return self.__series
	def __is_busy(self):
		PyQt4.QtCore.QMutexLocker(self.__mutex)
		return self.__busy
	def __set_busy(self,value):
		PyQt4.QtCore.QMutexLocker(self.__mutex)
		self.__busy = value
	def __get_main_window(self):
		main_window = self.parent()
		while main_window.parent() : main_window = main_window.parent()
		return main_window
