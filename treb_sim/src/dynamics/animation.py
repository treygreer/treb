from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import Qt
from dynamics.constants import meter2foot
import numpy as np
#import cv  # for movie creation

class EnergyBar(QtGui.QWidget):
    def __init__(self, sim, frame):
        super().__init__()
        self.sim, self.frame = sim, frame
        #self.set_size_request(400,20)

    def draw(self):
        return
        # Create the cairo context
        self.cr = self.window.cairo_create()
        w,h = self.window.get_size()
        #print("w=%d h=%d" % (w, h))

        # clear background
        self.cr.set_source_rgb(0.75,0.75,0.75)
        self.cr.rectangle(0, 0, w, h)
        self.cr.fill()

        # draw potential energy bar
        potential_fraction = (self.frame.PE() - 
                              self.frame.PEmin) / self.sim.total_energy
        self.cr.set_source_rgb(0.1,1.0,0.1)
        self.cr.rectangle(0, 0, w*potential_fraction, h)
        self.cr.fill()

        # draw kinetic energy bar
        kinetic_fraction = self.frame.KE() / self.sim.total_energy
        self.cr.set_source_rgb(1.0,0.1,0.1)
        self.cr.rectangle(w*potential_fraction, 0, w*kinetic_fraction, h)
        self.cr.fill()

class EnergiesBox(QtGui.QGroupBox):
    def __init__(self, sim):
        super().__init__('Energies')
        grid = QtGui.QGridLayout()
        row = 0
        for frame in sim.frames:
            label = QtGui.QLabel(frame.name)
            grid.addWidget(label, 0, row)
            label.show()
            frame.energy_bar = EnergyBar(sim, frame)
            grid.addWidget(frame.energy_bar, 1, row)
            row=row+1
        self.setLayout(grid)

class RangeBox(QtGui.QGroupBox):
    def __init__(self, max_range):
        super().__init__('Range')
        self.range_bar = QtGui.QProgressBar()
        self.range_bar.setRange(0, int(max_range))
        hbox = QtGui.QHBoxLayout()
        self.text_box = QtGui.QLabel()
        hbox.addWidget(self.text_box)
        hbox.addWidget(self.range_bar)
        self.setLayout(hbox)
    def set_range(self, range_):
        self.range_bar.setValue(int(range_))
        self.text_box.setText("%g feet" % (meter2foot(range_)))

class TimeSlider(QtGui.QSlider):
    def __init__(self, sim):
        super().__init__(Qt.Horizontal)
        self.setRange(0, sim.Y.shape[0])
        self.setValue(0)
        self.setPageStep(sim.Y.shape[0]/100)
        self.setSingleStep(1)

class Animation(QtGui.QWidget):
    def __init__(self, sim, dist, movie_flag=False):
        super().__init__()
        self.sim, self.dist, self.movie_flag = sim, dist, movie_flag
        self.sim.total_energy = sum([frame.KE() + frame.PE() - frame.PEmin
                                      for frame in sim.frames])

        vbox = QtGui.QVBoxLayout()
        time_slider = TimeSlider(sim)
        draw_forces_button = QtGui.QCheckBox("Draw Force Vectors")
        self.drawing=Drawing(time_slider, draw_forces_button, sim)
        energies_box = EnergiesBox(sim)
        self.range_box = RangeBox(np.max(self.dist(self.sim.Y)))

        vbox.addWidget(self.drawing)
        vbox.addWidget(time_slider)
        vbox.addWidget(energies_box)
        vbox.addWidget(self.range_box)
        vbox.addWidget(draw_forces_button)
        self.setLayout(vbox)

        time_slider.valueChanged.connect(self.draw)
        draw_forces_button.stateChanged.connect(self.draw)

        self.update_range()

        self.show()

    def draw(self):
        self.drawing.draw()
        for frame in self.sim.frames:
            frame.energy_bar.draw()
        self.update_range()

    def update_range(self):
        self.range_box.set_range(self.dist(self.sim.Y[self.sim.time_idx]))
    '''
    def movie(self):
        #create a video writer
        width, height = self.vbox.window.get_size()
        writer = cv.CreateVideoWriter('movie.avi', -1, 60, (width, height), is_color=1)
        gtk.gdk.window_process_all_updates()
        for time in np.arange(0.0, 1.0, 0.001):
            self.time_adj.set_value(time)
            gtk.gdk.window_process_all_updates()
            pb = gtk.gdk.Pixbuf(gtk.gdk.COLORSPACE_RGB, False, 8, width, height)
            pb = pb.get_from_drawable(self.vbox.window, self.vbox.window.get_colormap(), 0, 0, 0, 0, width, height)
            pixel_str = pb.get_pixels()
            pixel_array = np.fromstring(pixel_str, dtype=np.uint8)
            pixel_array = pixel_array.reshape((height, width, 3), order='C')
            cv_mat = cv.fromarray(pixel_array)
            #cv_frame = cv.DecodeImage(cv_mat, iscolor=1)
            cv.WriteFrame(writer, cv_mat)
    '''


class Drawing(QtGui.QWidget):
    def __init__(self, time_adj, draw_forces, sim):
        super().__init__()
        self.time_adj, self.draw_forces, self.sim = time_adj,draw_forces,sim
        #self.set_size_request(400,400)
        self.sim.time_idx = 0

    def draw(self):
        return
        # Create the cairo context
        self.cr = self.window.cairo_create()

        # make window 15x15 meters in drawing units, with origin at center
        self.draw_w, self.draw_h = 15,15  # in drawing units
        w,h = self.window.get_size()
        #self.cr.translate(w*0.5, h*0.65)
        self.cr.translate(w*0.5, h*0.5)
        scale = max(w, h) / max(self.draw_w,self.draw_h)
        self.draw_w = w/scale; self.draw_h = h / scale
        self.cr.scale(scale,-scale)

        # clear screen
        self.cr.set_source_rgb(1.0,1.0,1.0)
        self.cr.rectangle(-self.draw_w,-self.draw_h,
                          2.0*self.draw_w,  2.0*self.draw_h)
        self.cr.fill()

        # set the simulation state
        time = self.time_adj.value
        drawForces = self.draw_forces.get_active()
        self.sim.time_idx = round(time/self.sim.time_step)
        self.sim.time_idx = min(self.sim.num_steps-1, self.sim.time_idx)
        #print "time=", time, "time_idx=", self.sim.time_idx
        c=0
        for constraint in self.sim.constraints:
            constraint.enabled = self.sim.constraints_enabled \
                [self.sim.time_idx,c]
            c=c+1
        self.sim._deriv(time, self.sim.Y[self.sim.time_idx])

        # draw the objects and constraints
        for f in self.sim.frames:
            f.draw(self.cr)
        for c in self.sim.constraints:
            c.draw(self.cr, drawForces)

