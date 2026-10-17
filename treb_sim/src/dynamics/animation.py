'''animation:  animate a simulation in a PyQt5 widget'''

from PyQt5 import QtGui, QtCore, QtWidgets
from PyQt5.QtCore import Qt
from dynamics.constants import meter2foot
from dynamics.sketchup_frame import SketchupFrame
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
        painter.drawRect(int(width*potential_fraction), 0,
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
        range_ = np.asarray(range_).item()
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
        self._slider.setPageStep(max(1, sim.Y.shape[0]//100))
        self._slider.setSingleStep(1)
        self._slider.valueChanged.connect(self._value_changed)
    def _value_changed(self, _):
        self._text.setText("{:.3f} sec".format(self.time))
        self.timeChanged.emit()
    @property
    def time(self):
        return self.sim.t[self._slider.value()]

class Animation(QtWidgets.QWidget):
    def __init__(self, sim, dist, *, view_mode='3d'):
        super().__init__()
        self.sim, self.dist = sim, dist
        self.set_sim_state(time_idx=0)
        self.sim.total_energy = sum([frame.KE() + frame.PE() - frame.PEmin
                                      for frame in sim.frames])

        vbox = QtWidgets.QVBoxLayout()
        self.time_slider = TimeSlider(sim)
        draw_forces_button = QtWidgets.QCheckBox("Draw Force Vectors")
        self.frame_button = QtWidgets.QPushButton("Show Frame")
        self.frame_button.setObjectName('frame_visibility')
        self.frame_button.setCheckable(True)
        self.frame_button.setChecked(True)
        self.frame_button.setEnabled(hasattr(sim, 'machineFrame'))
        self.frame_button.setMaximumHeight(30)
        if view_mode == '3d':
            from dynamics.animation_3d import Drawing3D
            drawing_class = Drawing3D
        elif view_mode == '2d':
            drawing_class = Drawing
        else:
            raise ValueError('view_mode must be 2d or 3d')
        self.drawing=drawing_class(self.time_slider, draw_forces_button, sim, self.frame_button)
        energies_box = EnergiesBox(sim)
        self.range_box = RangeBox(np.max(self.dist(self.sim, self.sim.Y)))

        vbox.addWidget(self.drawing)
        if view_mode == '3d':
            controls = QtWidgets.QHBoxLayout()
            controls.addWidget(QtWidgets.QLabel('Drag to rotate; Ctrl+drag to pan; wheel to zoom'))
            controls.addStretch()
            side_button = QtWidgets.QPushButton('Side View')
            reset_button = QtWidgets.QPushButton('Reset View')
            side_button.setMaximumHeight(28)
            reset_button.setMaximumHeight(28)
            side_button.clicked.connect(self.drawing.side_camera)
            reset_button.clicked.connect(self.drawing.reset_camera)
            controls.addWidget(side_button)
            controls.addWidget(reset_button)
            vbox.addLayout(controls)
        vbox.addLayout(self.time_slider)
        vbox.addWidget(energies_box)
        vbox.addWidget(self.range_box)
        vbox.addWidget(draw_forces_button)
        vbox.addWidget(self.frame_button)
        self.setLayout(vbox)

        self.time_slider.timeChanged.connect(self.update)
        draw_forces_button.stateChanged.connect(self.update)
        self.frame_button.toggled.connect(self.update)

        self.update()
        self.show()

    def set_sim_state(self, time_idx):
        c_idx=0
        for constraint in self.sim.constraints:
            constraint.enabled = self.sim.constraints_enabled[time_idx,c_idx]
            c_idx=c_idx+1
        self.sim.time_idx = time_idx
        self.sim.deriv(self.sim.t[time_idx], self.sim.Y[time_idx])

    def update(self):
        '''set the simulation state and update the widget'''
        time_idx = self.time_slider._slider.value()
        time_idx = min(self.sim.num_steps-1, time_idx)
        self.set_sim_state(time_idx)
        self.range_box.set_range(self.dist(self.sim, self.sim.Y[time_idx]))
        if hasattr(self.drawing, 'update_scene'):
            self.drawing.update_scene()
        self.drawing.update()
        super().update()

class Drawing(QtWidgets.QGraphicsView):
    def __init__(self, time_slider, draw_forces_button, sim, frame_button):
        self.time_slicer = time_slider
        self.draw_forces_button, self.sim = draw_forces_button, sim
        self.frame_button = frame_button
        self.frame_art = SketchupFrame() if hasattr(sim, 'machineFrame') else None
        self.initial_bounds = None
        if self.frame_art:
            # Include the throwing arm's full sweep in the initial camera view.
            arm_axle = sim.axle.xframe0.A1
            radius = 0.0
            for obj in sim.armFrame.objects:
                points = obj.obj2frame(getattr(obj, 'corners', obj.cg)).A.T
                extent = np.max(np.linalg.norm(points - arm_axle, axis=1))
                radius = max(radius, extent + getattr(obj, 'radius', 0.0))
            self.initial_bounds = self.frame_art.bounds.united(
                QtCore.QRectF(-radius, -radius, 2 * radius, 2 * radius))
        self._user_zoomed = False
        self.scene = QtWidgets.QGraphicsScene()

        super().__init__(self.scene)
        self.setDragMode(self.ScrollHandDrag)
        self.setRenderHint(QtGui.QPainter.Antialiasing)
        self.setRenderHint(QtGui.QPainter.HighQualityAntialiasing)
        self.setTransformationAnchor(self.AnchorUnderMouse)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setTransform(QtGui.QTransform.fromScale(1.0, -1.0))
        self.sim.time_idx = 0
        self.scene.setSceneRect(-10.0, -5.0,
                                20.0, 20.0) # in drawing units
        self.scale(32.0,32.0)
        if self.frame_art:
            self.scene.setSceneRect(self.initial_bounds.adjusted(-1, -1, 1, 1))
    def sizeHint(self):
        return QtCore.QSize(600,450)
    def wheelEvent(self, event):
        self._user_zoomed = True
        scale_factor = 1.15
        if event.angleDelta().y() > 0:
            self.scale(scale_factor, scale_factor)
        else:
            self.scale(1/scale_factor, 1/scale_factor)
    def resizeEvent(self, event):
        super().resizeEvent(event)
        if self.frame_art and not self._user_zoomed:
            self.fitInView(self.initial_bounds.adjusted(-0.5, -0.5, 0.5, 0.5),
                           Qt.KeepAspectRatio)
    def update(self):
        self.scene.clear()
        if self.frame_art and self.frame_button.isChecked():
            self.frame_art.draw(self.scene, self.sim)
        # draw the objects and constraints
        for frame in self.sim.frames:
            if frame is getattr(self.sim, 'machineFrame', None):
                continue
            frame.draw(self.scene)
        for constraint in self.sim.constraints:
            constraint.draw(self.scene,
                            self.draw_forces_button.isChecked())
        for spring in self.sim.springs:
            spring.draw(self.scene,
                        self.draw_forces_button.isChecked())
        super().update()
