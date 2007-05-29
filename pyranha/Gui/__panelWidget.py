# Copyright (C) 2007 by Francesco Biscani
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

import sys
from PyQt4 import QtGui,QtCore
from ui___panelWidget import Ui_panelWidget
from pyranha.Core import *

class panelWidget(QtGui.QWidget):
  def __init__(self):
    QtGui.QWidget.__init__(self)
    self.timer = QtCore.QTimer()
    self.timer.start(500)
    self.ui = Ui_panelWidget()
    self.ui.setupUi(self)
    self.ui.digitsSlider.setMinimum(stream_manager.min_digits())
    self.ui.digitsSlider.setMaximum(stream_manager.max_digits())
    self.ui.digitsSlider.setValue(stream_manager.digits())
    self.ui.digitsLabel.setText(QtCore.QString.number(stream_manager.digits()))
    self.setWindowTitle("Pyranha control panel")
    self.ui.theoriesPathLineEdit.setText(settings_manager.default_theories_path())
    self.ui.fpComboBox.addItem("Scientific")
    self.ui.fpComboBox.addItem("Decimal")
    # Connections
    self.connect(self.ui.digitsSlider,QtCore.SIGNAL("valueChanged(int)"),self.__setDigits)
    self.connect(self.ui.theoriesPathButton,QtCore.SIGNAL("clicked()"),self.__setTheoriesPathDialog)
    self.connect(self.ui.resetPathButton,QtCore.SIGNAL("clicked()"),self.__resetTheoriesPath)
    self.connect(self.timer,QtCore.SIGNAL("timeout()"),self.__updateSettings)
    self.connect(self.ui.theoriesPathLineEdit,QtCore.SIGNAL("textChanged(const QString &)"),
      self.__setTheoriesPath)
    self.connect(self.ui.fpComboBox,QtCore.SIGNAL("currentIndexChanged(int)"),self.__setFpRep)
  def __setTheoriesPathDialog(self):
    dialog=QtGui.QFileDialog()
    dialog.setFileMode(QtGui.QFileDialog.DirectoryOnly)
    dialog.setViewMode(QtGui.QFileDialog.Detail)
    dialog.setWindowTitle("Select path for theories of motion")
    dialog.setDirectory(settings_manager.theories_path())
    retval=dialog.exec_()
    if retval==0:
      return
    dir=QtCore.QString(dialog.directory().absolutePath())
    self.ui.theoriesPathLineEdit.setText(dir)
  def __setDigits(self,n):
    stream_manager.set_digits(n)
    self.ui.digitsLabel.setText(QtCore.QString.number(stream_manager.digits()))
  def __resetTheoriesPath(self):
    self.ui.theoriesPathLineEdit.setText(settings_manager.default_theories_path())
  def __updateSettings(self):
    self.ui.digitsSlider.setValue(stream_manager.digits())
    if not self.ui.theoriesPathLineEdit.hasFocus():
      self.ui.theoriesPathLineEdit.setText(settings_manager.theories_path())
    if not self.ui.fpComboBox.hasFocus():
      self.ui.fpComboBox.setCurrentIndex(int(stream_manager.fp_rep()))
  def __setTheoriesPath(self,path):
    settings_manager.set_theories_path(path.__str__().__str__())
  def __setFpRep(self,n):
    stream_manager.set_fp_rep(fp_representation(n))


global panel
panel = panelWidget()
