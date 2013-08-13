#import pdb
import pygtk
pygtk.require('2.0')
import gtk
#import cairo
from dynamics.constants import meter2foot
#from numpy import arange
import numpy
#import cv

# Create a GTK+ widget on which we will draw using Cairo
class Animation:
    def __init__(self, sim, dist, movie_flag=False):
        self.sim, self.dist, self.movie_flag = sim, dist, movie_flag

        self.window = gtk.Window()
        self.vbox = gtk.VBox(homogeneous=False, spacing=1)

        self.time_adj = gtk.Adjustment(value=0, lower=0, upper=(sim.Y.shape[0]-1)*sim.time_step,
                                       step_incr=sim.time_step, page_incr=10*sim.time_step)
        time_scale = gtk.HScale(self.time_adj)
        time_scale.set_digits(3)
        time_scale.set_value_pos(gtk.POS_LEFT)

        draw_forces = gtk.CheckButton(label="Draw Force Vectors")

        energies_frame = gtk.Frame(label="Energies")
        energies_table = gtk.Table(rows=len(sim.frames), columns=2, homogeneous=False)
        energies_frame.add(energies_table)

        row = 0
        for frame in sim.frames:
            label = gtk.Label(frame.name)
            energies_table.attach(label, 0, 1, row, row+1, xoptions=gtk.SHRINK, xpadding=2)
            label.show()
            frame.energy_bar = EnergyBar(sim, frame)
            energies_table.attach(frame.energy_bar, 1, 2, row, row+1, xoptions=gtk.EXPAND,
                                  xpadding=2, ypadding=2)
            frame.energy_bar.show()
            row=row+1
        energies_table.show()
        self.sim.total_energy = sum([frame.KE() + frame.PE() - frame.PEmin for frame in sim.frames])

        range_frame = gtk.Frame(label="Range")
        self.range_bar = gtk.ProgressBar()
        self.range_bar.show()
        range_frame.add(self.range_bar)

        self.drawing=Drawing(self.time_adj, draw_forces, sim)
        self.vbox.pack_start(self.drawing, expand=True)
        self.drawing.show()
        self.vbox.pack_end(time_scale, expand=False)
        time_scale.show()
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
        range = self.dist(self.sim.Y[self.sim.time_idx])
        maxrange = numpy.max(self.dist(self.sim.Y))
        fraction = range/maxrange
        if fraction < 0:
            fraction = 0
        self.range_bar.set_fraction(fraction)
        self.range_bar.set_text("%g feet" % (meter2foot(range)))

#     def movie(self):
#         #create a video writer
#         width, height = self.vbox.window.get_size()
#         writer = cv.CreateVideoWriter('movie.avi', -1, 60, (width, height), is_color=1)
#         gtk.gdk.window_process_all_updates()
#         for time in arange(0.0, 1.0, 0.001):
#             self.time_adj.set_value(time)
#             gtk.gdk.window_process_all_updates()
#             pb = gtk.gdk.Pixbuf(gtk.gdk.COLORSPACE_RGB, False, 8, width, height)
#             pb = pb.get_from_drawable(self.vbox.window, self.vbox.window.get_colormap(), 0, 0, 0, 0, width, height)
#             pixel_str = pb.get_pixels()
#             pixel_array = numpy.fromstring(pixel_str, dtype=numpy.uint8)
#             pixel_array = pixel_array.reshape((height, width, 3), order='C')
#             cv_mat = cv.fromarray(pixel_array)
#             #cv_frame = cv.DecodeImage(cv_mat, iscolor=1)
#             cv.WriteFrame(writer, cv_mat)


class EnergyBar(gtk.DrawingArea):
    # Draw in response to an expose-event
    __gsignals__ = { "expose-event": "override" }

    def __init__(self, sim, frame):
        gtk.DrawingArea.__init__(self)
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
        potential_fraction = (self.frame.PE() - self.frame.PEmin) / self.sim.total_energy
        self.cr.set_source_rgb(0.1,1.0,0.1)
        self.cr.rectangle(0, 0, w*potential_fraction, h)
        self.cr.fill()

        # draw kinetic energy bar
        kinetic_fraction = self.frame.KE() / self.sim.total_energy
        self.cr.set_source_rgb(1.0,0.1,0.1)
        self.cr.rectangle(w*potential_fraction, 0, w*kinetic_fraction, h)
        self.cr.fill()

class Drawing(gtk.DrawingArea):
    # Draw in response to an expose-event
    __gsignals__ = { "expose-event": "override" }

    def __init__(self, time_adj, draw_forces, sim):
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
            constraint.enabled = self.sim.constraints_enabled[self.sim.time_idx,c]
            c=c+1
        self.sim._deriv(time, self.sim.Y[self.sim.time_idx])

        # draw the objects and constraints
        for f in self.sim.frames:
            f.draw(self.cr)
        for c in self.sim.constraints:
            c.draw(self.cr, drawForces)

