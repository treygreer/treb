from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import Qt
from dynamics.constants import meter2foot
import numpy as np
#import cv  # for movie creation

class EnergiesBox(QtGui.QGroupBox):
    def __init(self):
        super().__init__('Energies')

        energies_table = gtk.Table(rows=len(sim.frames), columns=2,
                                   homogeneous=False)
        energies_frame.add(energies_table)

        row = 0
        for frame in sim.frames:
            label = gtk.Label(frame.name)
            energies_table.attach(label, 0, 1, row, row+1, xoptions=gtk.SHRINK,
                                  xpadding=2)
            label.show()
            frame.energy_bar = EnergyBar(sim, frame)
            energies_table.attach(frame.energy_bar, 1, 2, row, row+1,
                                  xoptions=gtk.EXPAND,
                                  xpadding=2, ypadding=2)
            frame.energy_bar.show()
            row=row+1
        energies_table.show()
        self.sim.total_energy = sum([frame.KE() + frame.PE() - frame.PEmin
                                      for frame in sim.frames])

class RangeBox(QtGui.QGroupBox):
    def __init(self):
        super().__init__('Range')
        self.range_bar = gtk.ProgressBar()
        self.range_bar.show()
        range_frame.add(self.range_bar)

class Animation(QtGui.QVBoxLayout):
    def __init__(self, sim, dist, movie_flag=False):
        super().__init__()
        self.sim, self.dist, self.movie_flag = sim, dist, movie_flag

        time_slider = QtGui.QSlider(Qt.Horizontal)
        time_slider.setRange(0, sim.Y.shape[0])
        time_slider.setValue(0)
        time_slider.setPageStep(sim.Y.shape[0]/100)
        time_slider.setSingleStep(1)

        draw_forces_button = QtGui.QCheckBox("Draw Force Vectors")

        energies_frame = EnergiesBox()
        range_frame = RangeBox()

        self.drawing=Drawing(self.time_adj, draw_forces, sim)
        self.vbox.pack_start(self.drawing, expand=True)
        self.drawing.show()
        self.vbox.pack_end(time_slider, expand=False)
        time_slider.show()
        self.vbox.pack_end(energies_frame, expand=False)
        energies_frame.show()
        self.vbox.pack_end(range_frame, expand=False)
        range_frame.show()
        self.vbox.pack_end(draw_forces, expand=False)
        draw_forces.show()
        self.time_adj.connect("value_changed", lambda w: self.draw())
        draw_forces.connect("toggled", lambda w: self.draw())
        self.update_range()

        self.window.add(self.vbox)
        self.vbox.show()
        self.window.present()
        gtk.main()

    def draw(self):
        self.drawing.draw()
        for frame in self.sim.frames:
            frame.energy_bar.draw()
        self.update_range()


    def update_range(self):
        range_ = self.dist(self.sim.Y[self.sim.time_idx])
        maxrange = np.max(self.dist(self.sim.Y))
        fraction = range_/maxrange
        if fraction < 0:
            fraction = 0
        self.range_bar.set_fraction(fraction)
        self.range_bar.set_text("%g feet" % (meter2foot(range_)))
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


class EnergyBar(QtGui.QWidget):
    def __init__(self, sim, frame):
        super().__init__()
        self.sim, self.frame = sim, frame
        self.set_size_request(400,20)

    # Handle the expose-event by drawing
    def do_expose_event(self, event):
        self.draw()

    def draw(self):
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

class Drawing(QtGui.QWidget):
    def __init__(self, time_adj, draw_forces, sim):
        super().__init__()
        gtk.DrawingArea.__init__(self)
        self.time_adj, self.draw_forces, self.sim = time_adj,draw_forces,sim
        self.set_size_request(400,400)
        self.sim.time_idx = 0

    # Handle the expose-event by drawing
    def do_expose_event(self, event):
        self.draw()

    def draw(self):
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

