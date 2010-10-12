#!/usr/bin/env python

"""PySide port of the opengl/hellogl example from Qt v4.x"""

import sys
import math
import numpy as np
from PyQt4 import QtCore, QtGui, QtOpenGL

try:
    from OpenGL import GL
    from OpenGL import GLU
except ImportError:
    app = QtGui.QApplication(sys.argv)
    QtGui.QMessageBox.critical(None, "OpenGL hellogl",
                            "PyOpenGL must be installed to run this example.",
                            QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default,
                            QtGui.QMessageBox.NoButton)
    sys.exit(1)


class Window(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)

        self.glWidget = GLWidget()

        self.xSlider = self.createSlider(QtCore.SIGNAL("xRotationChanged(int)"),
                                         self.glWidget.setXRotation)
        self.ySlider = self.createSlider(QtCore.SIGNAL("yRotationChanged(int)"),
                                         self.glWidget.setYRotation)
        self.zSlider = self.createSlider(QtCore.SIGNAL("zRotationChanged(int)"),
                                         self.glWidget.setZRotation)
        self.zoomSlider = self.createSlider(QtCore.SIGNAL("zoomChanged(int)"),
                                         self.glWidget.setZoom)
        self.xposSlider = self.createSlider(QtCore.SIGNAL("xposChanged(int)"),
                                       self.glWidget.setXPos)
        self.yposSlider = self.createSlider(QtCore.SIGNAL("yposChanged(int)"),
                                         self.glWidget.setYPos)
        self.zposSlider = self.createSlider(QtCore.SIGNAL("zosChanged(int)"),
                                         self.glWidget.setZPos)

        
        mainLayout = QtGui.QHBoxLayout()
        mainLayout.addWidget(self.glWidget)
        mainLayout.addWidget(self.xSlider)
        mainLayout.addWidget(self.ySlider)
        mainLayout.addWidget(self.zSlider)
        mainLayout.addWidget(self.zoomSlider)
        
        mainLayout.addWidget(self.xposSlider)
        mainLayout.addWidget(self.yposSlider)
        mainLayout.addWidget(self.zposSlider)
        self.setLayout(mainLayout)

        self.xSlider.setValue(15 * 16)
        self.ySlider.setValue(345 * 16)
        self.zSlider.setValue(0 * 16)
        self.xposSlider.setValue(180*16)
        self.yposSlider.setValue(180*16)
        self.zposSlider.setValue(180*16)
        self.setWindowTitle(self.tr("Hello GL"))

    def createSlider(self, changedSignal, setterSlot):
        slider = QtGui.QSlider(QtCore.Qt.Vertical)

        slider.setRange(0, 360 * 16)
        slider.setSingleStep(16)
        slider.setPageStep(15 * 16)
        slider.setTickInterval(15 * 16)
        slider.setTickPosition(QtGui.QSlider.TicksRight)

        self.glWidget.connect(slider, QtCore.SIGNAL("valueChanged(int)"), setterSlot)
        self.connect(self.glWidget, changedSignal, slider, QtCore.SLOT("setValue(int)"))

        return slider


class GLWidget(QtOpenGL.QGLWidget):
    def __init__(self, parent=None):
        QtOpenGL.QGLWidget.__init__(self, parent)

        self.object = 0
        self.xRot = 0
        self.yRot = 0
        self.zRot = 0
        self.xpos = 0
        self.ypos = 0
        self.zpos = 10
        self.v_range = 2
        self.lastPos = QtCore.QPoint()

        self.trolltechGreen = QtGui.QColor.fromCmykF(0.40, 0.0, 1.0, 0.0)
        self.trolltechPurple = QtGui.QColor.fromCmykF(0.39, 0.39, 0.0, 0.0)

    def xRotation(self):
        return self.xRot

    def yRotation(self):
        return self.yRot

    def zRotation(self):
        return self.zRot

    def minimumSizeHint(self):
        return QtCore.QSize(50, 50)

    def sizeHint(self):
        return QtCore.QSize(400, 400)

    def setXRotation(self, angle):
        angle = self.normalizeAngle(angle)
        if angle != self.xRot:
            self.xRot = angle
            self.emit(QtCore.SIGNAL("xRotationChanged(int)"), angle)
            self.updateGL()

    def setYRotation(self, angle):
        angle = self.normalizeAngle(angle)
        if angle != self.yRot:
            self.yRot = angle
            self.emit(QtCore.SIGNAL("yRotationChanged(int)"), angle)
            self.updateGL()

    def setZRotation(self, angle):
        angle = self.normalizeAngle(angle)
        if angle != self.zRot:
            self.zRot = angle
            self.emit(QtCore.SIGNAL("zRotationChanged(int)"), angle)
            self.updateGL()

    
    def setZPos(self, factor):
        factor = self.normalizeAngle(factor)
        if factor != self.zpos-10:
            self.zpos = factor/360.0-8
            self.emit(QtCore.SIGNAL("zposChanged(int)"), factor)
            self.updateGL()
            
            
    def setXPos(self, factor):
        factor = self.normalizeAngle(factor)
        if factor != self.xpos:
            self.xpos = factor/360.0-8
            self.emit(QtCore.SIGNAL("xposChanged(int)"), factor)
            self.updateGL()

            
    def setYPos(self, factor):
        factor = self.normalizeAngle(factor)
        if factor != self.ypos:
            self.ypos = factor/360.0 -8
            self.emit(QtCore.SIGNAL("yposChanged(int)"), factor)
            self.updateGL()

            
    def setZoom(self, factor):
        factor = self.normalizeAngle(factor)
        if factor != self.v_range:
            self.v_range = factor/360.0 +2
            self.emit(QtCore.SIGNAL("zoomChanged(int)"), factor)
            self.updateGL()

    
            
    def initializeGL(self):
        self.qglClearColor(self.trolltechPurple.darker())
        self.object = self.makeObject()
        GL.glShadeModel(GL.GL_FLAT)
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_CULL_FACE)

    def paintGL(self):
        v_range = self.v_range
        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glLoadIdentity()
        GL.glOrtho(-v_range,v_range,-v_range,v_range, 0, 150.0)
        GL.glMatrixMode(GL.GL_MODELVIEW)

        
        
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        GL.glLoadIdentity()
             
        GL.glTranslated(0,0,-1)
        GL.glRotated(self.xRot / 16.0, 1.0, 0.0, 0.0)
        GL.glRotated(self.yRot / 16.0, 0.0, 1.0, 0.0)
        GL.glRotated(self.zRot / 16.0, 0.0, 0.0, 1.0)
        GL.glTranslated(-self.xpos, -self.ypos, -self.zpos)
        GL.glCallList(self.object)
        
        

    def resizeGL(self, width, height):
        side = min(width, height)
        GL.glViewport((width - side) / 2, (height - side) / 2, side, side)
        v_range = self.v_range
        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glLoadIdentity()
        GL.glOrtho(-v_range,v_range,-v_range,v_range, 0, 150.0)
        GL.glMatrixMode(GL.GL_MODELVIEW)
        
    def mousePressEvent(self, event):
        self.lastPos = QtCore.QPoint(event.pos())

    def mouseMoveEvent(self, event):
        dx = event.x() - self.lastPos.x()
        dy = event.y() - self.lastPos.y()

        if event.buttons() & QtCore.Qt.LeftButton:
            self.setXRotation(self.xRot + 8 * dy)
            self.setYRotation(self.yRot + 8 * dx)
        elif event.buttons() & QtCore.Qt.RightButton:
            self.setXRotation(self.xRot + 8 * dy)
            self.setZRotation(self.zRot + 8 * dx)

        self.lastPos = QtCore.QPoint(event.pos())

    def makeObject(self):
        genList = GL.glGenLists(1)
        GL.glNewList(genList, GL.GL_COMPILE)

        s3 = np.sqrt(3)
        self.sphere(0.5,0,0,0,QtGui.QColor.fromCmykF(1,1,0,.5))
        self.sphere(0.5,1,0,0,QtGui.QColor.fromCmykF(1,1,0,.5))
        self.sphere(0.5,.5,np.sqrt(3)/2,0,QtGui.QColor.fromCmykF(1,1,0,.5))
        #self.sphere(0.5,.5,-np.sqrt(3)/2,0)
        self.sphere(0.5,.5,1/(2*s3),np.sqrt(6)/3,QtGui.QColor.fromCmykF(.25,0,.75,.44))
        self.sphere(0.5,.5,1/(2*s3),-np.sqrt(6)/3,QtGui.QColor.fromCmykF(.25,0,.75,.44))
        #self.sphere(0.5,.5,.5,0)
        GL.glEndList()

        return genList
    def sphere(self,R,x1,y1,z1,scolor = None):
        if scolor is None:
            self.qglColor(QtGui.QColor.fromCmykF(0, 1.0, 1.0, 0.6))
        else:
            self.qglColor(scolor)
        panel_count = 50
        plane_count = 50
        def gen_x(lat,long):
            return R * np.cos(2*np.pi/panel_count * long) * np.sin(np.pi/plane_count * lat) +x1
        def gen_y(lat,long):
            return R * np.sin(2*np.pi/panel_count * long) * np.sin(np.pi/plane_count * lat)+y1
        def gen_z(lat):
            return R * np.cos(np.pi/plane_count * lat)+z1
        for j in range(0,plane_count+1):
            GL.glBegin(GL.GL_QUAD_STRIP)
            for k in range(0,panel_count+1):
                               
                GL.glVertex3d(gen_x(j+1,k), gen_y(j+1,k),gen_z(j+1) )
                GL.glVertex3d(gen_x(j,k), gen_y(j,k),gen_z(j) )
            GL.glEnd()    
                
                
                
                
    def quad(self, x1, y1, x2, y2, x3, y3, x4, y4):
        self.qglColor(self.trolltechGreen)

        GL.glVertex3d(x1, y1, -0.05)
        GL.glVertex3d(x2, y2, -0.05)
        GL.glVertex3d(x3, y3, -0.05)
        GL.glVertex3d(x4, y4, -0.05)

        GL.glVertex3d(x4*0.8, y4*0.8, +0.05)
        GL.glVertex3d(x3*0.8, y3*0.8, +0.05)
        GL.glVertex3d(x2*0.8, y2*0.8, +0.05)
        GL.glVertex3d(x1*0.8, y1*0.8, +0.05)

    def extrude(self, x1, y1, x2, y2):
        self.qglColor(self.trolltechGreen.darker(250 + int(100 * x1)))

        GL.glVertex3d(x1*0.8, y1*0.8, +0.05)
        GL.glVertex3d(x2*0.8, y2*0.8, +0.05)
        GL.glVertex3d(x2, y2, -0.05)
        GL.glVertex3d(x1, y1, -0.05)

    def normalizeAngle(self, angle):
        while angle < 0:
            angle += 360 * 16
        while angle > 360 * 16:
            angle -= 360 * 16
        return angle


if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())
