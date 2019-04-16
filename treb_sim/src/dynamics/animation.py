'''animation:  animate a simulation in a PyQt5 widget'''

from PyQt5 import QtGui, QtCore, QtWidgets
from PyQt5.QtCore import Qt
from dynamics.constants import meter2foot
import numpy as np

class EnergyBar(QtWidgets.QWidget):
    '''widget to display potential, kinetic, and dissipated energies'''
    def __init__(self, sim, frame):
        super().__init__()
        self.sim, self.frame = sim, frame
    def sizeHint(self): #IGNORE:invalid-name
        '''return minimum size for widget'''
        return QtCore.QSize(400,20)
    def paintEvent(self, _): #IGNORE:invalid-name
        '''paint widget'''
        painter = QtGui.QPainter()
        painter.begin(self)
        # clear background
        width,height = self.size().width(), self.size().height()
        painter.setBrush(QtGui.QBrush(QtGui.QColor(190,190,190)))
        painter.drawRect(0, 0, width, height)

        # draw potential energy bar
        potential_fraction = float((self.frame.PE() - self.frame.PEmin) /
                                     self.sim.total_energy)
        painter.setBrush(QtGui.QBrush(QtGui.QColor(25,255,25)))
        painter.drawRect(0, 0, int(width*potential_fraction), height)

        # draw kinetic energy bar
        kinetic_fraction = float(self.frame.KE() / self.sim.total_energy)
        painter.setBrush(QtGui.QBrush(QtGui.QColor(255,25,25)))
        painter.drawRect(width*potential_fraction, 0,
                         int(width*kinetic_fraction), height)
        painter.end()

class EnergiesBox(QtWidgets.QGroupBox):
    '''box of energy bars, one per frame and spring'''
    def __init__(self, sim):
        super().__init__('Energies')
        grid = QtWidgets.QGridLayout()
        row = 0
        for frame in sim.frames:
            label = QtWidgets.QLabel(frame.name)
            grid.addWidget(label, row, 0)
            label.show()
            frame.energy_bar = EnergyBar(sim, frame)
            grid.addWidget(frame.energy_bar, row, 1)
            row=row+1
        self.setLayout(grid)

class RangeBox(QtWidgets.QGroupBox):
    '''display scalar value'''
    def __init__(self, max_range):
        super().__init__('Range')
        self.range_bar = QtWidgets.QProgressBar()
        self.range_bar.setRange(0, int(max_range))
        hbox = QtWidgets.QHBoxLayout()
        self.text_box = QtWidgets.QLabel()
        hbox.addWidget(self.text_box)
        hbox.addWidget(self.range_bar)
        self.setLayout(hbox)
    def set_range(self, range_):
        '''set value'''
        self.range_bar.setValue(int(range_))
        self.text_box.setText("%g feet" % (meter2foot(range_)))

class TimeSlider(QtWidgets.QHBoxLayout):
    timeChanged = QtCore.pyqtSignal()
    def __init__(self, sim):
        self.sim = sim
        super().__init__()
        self._text = QtWidgets.QLabel()
        self.addWidget(self._text)
        self._slider = QtWidgets.QSlider(Qt.Horizontal)
        self.addWidget(self._slider)
        self._slider.setRange(0, sim.Y.shape[0]-1)
        self._slider.setValue(0)
        self._slider.setPageStep(sim.Y.shape[0]/100)
        self._slider.setSingleStep(1)
        self._slider.valueChanged.connect(self._value_changed)
    def _value_changed(self, _):
        self._text.setText("{:.3f} sec".format(self.time))
        self.timeChanged.emit()
    @property
    def time(self):
        return self.sim.time_step * self._slider.value()

class Animation(QtWidgets.QWidget):
    def __init__(self, sim, dist):
        super().__init__()
        self.sim, self.dist = sim, dist
        self.set_sim_state(time_idx=0)
        self.sim.total_energy = sum([frame.KE() + frame.PE() - frame.PEmin
                                      for frame in sim.frames])

        vbox = QtWidgets.QVBoxLayout()
        self.time_slider = TimeSlider(sim)
        draw_forces_button = QtWidgets.QCheckBox("Draw Force Vectors")
        self.drawing=Drawing(self.time_slider, draw_forces_button, sim)
        energies_box = EnergiesBox(sim)
        self.range_box = RangeBox(np.max(self.dist(self.sim, self.sim.Y)))

        vbox.addWidget(self.drawing)
        vbox.addLayout(self.time_slider)
        vbox.addWidget(energies_box)
        vbox.addWidget(self.range_box)
        vbox.addWidget(draw_forces_button)
        self.setLayout(vbox)

        self.time_slider.timeChanged.connect(self.update)
        draw_forces_button.stateChanged.connect(self.update)

        self.update()
        self.show()

    def set_sim_state(self, time_idx):
        c_idx=0
        for constraint in self.sim.constraints:
            constraint.enabled = self.sim.constraints_enabled[time_idx,c_idx]
            c_idx=c_idx+1
        self.sim.time_idx = time_idx
        self.sim.deriv(time_idx * self.sim.time_step, self.sim.Y[time_idx])

    def update(self):
        '''set the simulation state and update the widget'''
        time = self.time_slider.time
        time_idx = round(time/self.sim.time_step)
        time_idx = min(self.sim.num_steps-1, time_idx)
        self.set_sim_state(time_idx)
        self.range_box.set_range(self.dist(self.sim, self.sim.Y[time_idx]))
        self.drawing.update()
        super().update()

class Drawing(QtWidgets.QGraphicsView):
    def __init__(self, time_slider, draw_forces_button, sim):
        self.time_slicer = time_slider
        self.draw_forces_button, self.sim = draw_forces_button, sim
        self.scene = QtWidgets.QGraphicsScene()

        super().__init__(self.scene)
        self.setDragMode(self.ScrollHandDrag)
        self.setRenderHint(QtGui.QPainter.Antialiasing)
        self.setRenderHint(QtGui.QPainter.HighQualityAntialiasing)
        self.setTransformationAnchor(self.AnchorUnderMouse)
        self.setTransform(QtGui.QTransform.fromScale(1.0, -1.0))
        self.sim.time_idx = 0
        self.scene.setSceneRect(-10.0, -5.0,
                                20.0, 20.0) # in drawing units
        self.scale(32.0,32.0)
    def sizeHint(self):
        return QtCore.QSize(400,400)
    def wheelEvent(self, event):
        scale_factor = 1.15
        if event.delta() > 0:
            self.scale(scale_factor, scale_factor)
        else:
            self.scale(1/scale_factor, 1/scale_factor)
    def update(self):
        self.scene.clear()
        # draw the objects and constraints
        for frame in self.sim.frames:
            frame.draw(self.scene)
        for constraint in self.sim.constraints:
            constraint.draw(self.scene,
                            self.draw_forces_button.isChecked())
        for spring in self.sim.springs:
            spring.draw(self.scene,
                        self.draw_forces_button.isChecked())
        super().update()
