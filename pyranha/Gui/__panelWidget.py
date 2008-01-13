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

import sys,time
from PyQt4 import QtGui,QtCore
from ui___panelWidget import Ui_panelWidget
from pyranha.Core import *

class latexRenderThread(QtCore.QThread):
  def __init__(self,str_):
    QtCore.QThread.__init__(self)
    self.str=str_
  def run(self):
    self.returnValue=latexRender(self.str)

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
    self.ui.latexRenderingProgressBar.hide()
    self.symCache=dict()
    self.renderFlag=False
    # Connections
    self.connect(self.ui.digitsSlider,QtCore.SIGNAL("valueChanged(int)"),self.__setDigits)
    self.connect(self.ui.theoriesPathButton,QtCore.SIGNAL("clicked()"),self.__setTheoriesPathDialog)
    self.connect(self.ui.resetPathButton,QtCore.SIGNAL("clicked()"),self.__resetTheoriesPath)
    self.connect(self.timer,QtCore.SIGNAL("timeout()"),self.__updateSettings)
    self.connect(self.ui.theoriesPathLineEdit,QtCore.SIGNAL("textChanged(const QString &)"),
      self.__setTheoriesPath)
    self.connect(self.ui.fpComboBox,QtCore.SIGNAL("currentIndexChanged(int)"),self.__setFpRep)
    self.show()
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
    self.__updatePsymbolList()
    self.__updateSymbolsIcons()
  def __setTheoriesPath(self,path):
    settings_manager.set_theories_path(path.__str__().__str__())
  def __setFpRep(self,n):
    stream_manager.set_fp_rep(fp_representation(n))
  def __updatePsymbolList(self):
    n_psym=psymbol_manager.__len__() 
    n_twi=self.ui.treeWidget.topLevelItemCount()
    if n_psym != n_twi:
      for i in psymbol_manager():
        if self.ui.treeWidget.findItems(i.name(),QtCore.Qt.MatchExactly).__len__()==0:
          newSym1=QtGui.QTreeWidgetItem(self.ui.treeWidget)
          newSym1.setText(0,i.name())
          newSym2=QtGui.QTreeWidgetItem(newSym1)
          newSym2.setText(0,i.powers_string().replace(stream_manager.data_separator(),'\n'))
          self.ui.treeWidget.addTopLevelItem(newSym1)
          self.ui.treeWidget.addTopLevelItem(newSym2)
  def __updateSymbolsIcons(self):
    if self.renderFlag:
        print "I want to render, but slot is busy. Waiting for next opportunity."
        return
    self.renderFlag=True
    if self.ui.latexRenderCheckBox.checkState():
      # Build a list of tree widget items that have no associated icon.
      twiList=[i for i in self.ui.treeWidget.findItems("*",QtCore.Qt.MatchWildcard) if i.icon(0).isNull()]
      if twiList:
        self.ui.latexRenderingProgressBar.setMaximum(twiList.__len__())
        self.ui.latexRenderingProgressBar.setValue(0)
        self.ui.latexRenderingProgressBar.show()
      for i in twiList:
        i.setIcon(0,QtGui.QIcon(self.__latexRender(i.text(0))))
        self.ui.latexRenderingProgressBar.setValue(self.ui.latexRenderingProgressBar.value()+1)
      self.ui.latexRenderingProgressBar.hide()
    else:
      for i in [j for j in self.ui.treeWidget.findItems("*",QtCore.Qt.MatchWildcard) if not j.icon(0).isNull()]:
        i.setIcon(0,QtGui.QIcon())
    self.renderFlag=False
  def __latexRender(self,str):
    if self.symCache.has_key(str):
      return self.symCache[str]
    print "no cache for " + str
    retval=latexRender(str)
    self.symCache[str]=retval
    return retval


def __latexCleanup(baseName):
  assert QtCore.QFile(baseName+".log").remove()
  assert QtCore.QFile(baseName+".aux").remove()
  assert QtCore.QFile(baseName+".dvi").remove()


def latexRender(str_):
  str=QtCore.QString(r"\documentclass{minimal}\begin{document}$")+str_+r"$\end{document}"
  tmpFileTex=QtCore.QTemporaryFile(QtCore.QDir.tempPath()+str_.__str__().replace("\\","backslash_")+"_pyranha_tmp_XXXXXX.tex")
  if tmpFileTex.open():
    print "Opened file " + tmpFileTex.fileName()
    out=QtCore.QTextStream(tmpFileTex)
    out << str << "\n"
    out.flush()
    print "writing " + str + " to file"
  else:
    print "Error opening file " + tmpFileTex.fileName()
    print "Aborting"
    return QtGui.QPixmap(":/images/symbol_broken.png")
  # Basename of latex files, including directory.
  baseName=QtCore.QDir.tempPath()+QtCore.QFileInfo(tmpFileTex).baseName()
  latexProcess=QtCore.QProcess()
  latexProcess.setWorkingDirectory(QtCore.QDir.tempPath())
  arguments=QtCore.QStringList()
  arguments << "-interaction=nonstopmode" << QtCore.QFileInfo(tmpFileTex).fileName()
  latexProcess.start("latex",arguments)
  while latexProcess.state()!=QtCore.QProcess.NotRunning:
    QtCore.QCoreApplication.processEvents()
    time.sleep(.05)
  latexErr=QtCore.QString(latexProcess.readAllStandardOutput())
  if latexProcess.exitCode():
    print "Latex exited with an error."
    print latexErr
    __latexCleanup(baseName)
    return QtGui.QPixmap(":/images/symbol_broken.png")
  assert tmpFileTex.remove()
  tmpFilePng=QtCore.QTemporaryFile(QtCore.QDir.tempPath()+str_.__str__().replace("\\","backslash_")+"_pyranha_tmp_XXXXXX.png")
  if tmpFilePng.open():
    print "Opened file " + tmpFilePng.fileName()
    pass
  else:
    print "Error opening file " + tmpFilePng.fileName()
    print "Aborting."
    __latexCleanup(baseName)
    return QtGui.QPixmap(":/images/symbol_broken.png")
  dvipngProcess=QtCore.QProcess()
  dvipngProcess.setWorkingDirectory(QtCore.QDir.tempPath())
  arguments.clear()
  arguments << "--dvinum*" << "-T" << "tight" << "-D" << "120" << "-bg" << "Transparent" << (baseName+".dvi") << "-o" <<  QtCore.QFileInfo(tmpFilePng).fileName()
  dvipngProcess.start("dvipng",arguments)
  while dvipngProcess.state()!=QtCore.QProcess.NotRunning:
    QtCore.QCoreApplication.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)
    time.sleep(.05)
  dvipngErr=QtCore.QString(dvipngProcess.readAllStandardOutput())
  if dvipngProcess.exitCode():
    print "dvipng exited with an error."
    print dvipngErr
    __latexCleanup(baseName)
    return QtGui.QPixmap()
  print "Png file complete path is: " + tmpFilePng.fileName()
  retval=QtGui.QPixmap(tmpFilePng.fileName())
  assert tmpFilePng.remove()
  # Clean up temp files.
  __latexCleanup(baseName)
  return retval

global panel
panel = panelWidget()
