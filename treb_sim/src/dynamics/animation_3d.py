"""3D rigid-assembly rendering of the existing planar simulation."""

import json
from pathlib import Path

import numpy as np
from PyQt5 import QtCore, QtGui, QtWidgets
import pyqtgraph.opengl as gl
from pyqtgraph.opengl import shaders


ASSETS = Path(__file__).resolve().parents[2] / 'assets'

# SketchUp contains open and mirrored surfaces. Illuminate both orientations
# with sufficient ambient light to keep their material colors recognizable.
shaders.ShaderProgram('assembly', [
    shaders.VertexShader('''
        uniform mat4 u_mvp;
        uniform mat3 u_normal;
        attribute vec4 a_position;
        attribute vec3 a_normal;
        attribute vec4 a_color;
        varying vec4 v_color;
        varying vec3 v_normal;
        void main() {
            v_normal = normalize(u_normal * a_normal);
            v_color = a_color;
            gl_Position = u_mvp * a_position;
        }
    '''),
    shaders.FragmentShader('''
        #ifdef GL_ES
        precision mediump float;
        #endif
        varying vec4 v_color;
        varying vec3 v_normal;
        void main() {
            float light = abs(dot(normalize(v_normal), normalize(vec3(0.5, -0.7, 1.0))));
            gl_FragColor = vec4(v_color.rgb * (0.58 + 0.42 * light), v_color.a);
        }
    ''')
])


def plane_matrix(angle, origin):
    """Embed a simulator X/Y rotation in the viewer's X/Z plane."""
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, 0, -s, origin[0]], [0, 1, 0, 0],
                     [s, 0, c, origin[1]], [0, 0, 0, 1]], dtype=float)


def body_matrix(sim, binding):
    constraint = getattr(sim, binding['constraint'])
    frame = getattr(constraint, f"frame{binding['side']}")
    anchor = getattr(constraint, f"xframe{binding['side']}")
    return plane_matrix(frame.theta - binding['reference_angle'],
                        frame.frame2world(anchor).A1)


def object_matrix(obj):
    origin = obj.frame.frame2world(obj.origin).A1
    return plane_matrix(obj.frame.theta + obj.theta, origin)


def square_tube(length, width, wall=0.0047625):
    """Hollow square display geometry; wall thickness does not affect physics."""
    outer, inner = width / 2, width / 2 - min(wall, width / 4)
    points = []
    for x in (-length / 2, length / 2):
        for radius in (outer, inner):
            points.extend([(x, -radius, -radius), (x, radius, -radius),
                           (x, radius, radius), (x, -radius, radius)])
    faces = []
    for i in range(4):
        j = (i + 1) % 4
        for quad in ((i, j, j + 8, i + 8), (i + 4, i + 12, j + 12, j + 4),
                     (i, i + 4, j + 4, j), (i + 8, j + 8, j + 12, i + 12)):
            a, b, c, d = quad
            faces.extend(((a, b, c), (a, c, d)))
    return np.array(points, dtype=np.float32), np.array(faces, dtype=np.uint32)


def cylinder(length, radius, sides=24):
    angle = np.arange(sides) * 2 * np.pi / sides
    ring = np.column_stack((np.zeros(sides), radius * np.cos(angle), radius * np.sin(angle)))
    points = np.vstack((ring + [-length / 2, 0, 0], ring + [length / 2, 0, 0],
                        [-length / 2, 0, 0], [length / 2, 0, 0]))
    faces = []
    for i in range(sides):
        j = (i + 1) % sides
        faces.extend(((i, j, j + sides), (i, j + sides, i + sides),
                      (2 * sides, j, i), (2 * sides + 1, i + sides, j + sides)))
    return points.astype(np.float32), np.array(faces, dtype=np.uint32)


class Drawing3D(gl.GLViewWidget):
    def __init__(self, time_slider, draw_forces_button, sim, frame_button):
        super().__init__()
        self.sim = sim
        self.draw_forces_button = draw_forces_button
        self.frame_button = frame_button
        self.setBackgroundColor((244, 246, 249))
        self.setMinimumHeight(320)
        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.metadata = json.loads((ASSETS / 'sketchup_3d.json').read_text(encoding='utf-8'))
        self.body_items = {}
        with np.load(ASSETS / 'sketchup_3d.npz', allow_pickle=False) as meshes:
            for body in self.metadata['bindings']:
                if not hasattr(sim, 'machineFrame'):
                    break
                data = gl.MeshData(vertexes=meshes[f'{body}_vertices'],
                                   faces=meshes[f'{body}_faces'],
                                   faceColors=meshes[f'{body}_colors'])
                item = gl.GLMeshItem(meshdata=data, smooth=False, shader='assembly', glOptions='opaque')
                self.body_items[body] = item
                self.addItem(item)
            ground_z = float(meshes['machine_vertices'][:, 2].min())
            low = meshes['machine_vertices'].min(axis=0).astype(float)
            high = meshes['machine_vertices'].max(axis=0).astype(float)
            # Include the long arm's complete sweep, not just the cocked pose.
            low = np.minimum(low, [-3.2, -0.2, -3.2])
            high = np.maximum(high, [3.2, 0.2, 3.2])
            self.camera_center = (low + high) / 2
            self.camera_corners = np.array([[x, y, z] for x in (low[0], high[0])
                                           for y in (low[1], high[1])
                                           for z in (low[2], high[2])])
            self._manual_camera = False

        self.object_items = []
        if hasattr(sim, 'connectorFrame'):
            # The SketchUp export lacks the long rod's swept surface. Preserve
            # the simulator rod and the later square reinforcement tube.
            rod = sim.connectorFrame.rod
            tube = sim.connectorFrame.stiffener
            for obj, geometry, color in ((rod, cylinder(rod.l, rod.w / 2), (0.3, 0.34, 0.4, 1)),
                                         (tube, square_tube(tube.l, tube.w), (0.35, 0.4, 0.46, 1))):
                vertexes, faces = geometry
                item = gl.GLMeshItem(vertexes=vertexes, faces=faces, color=color,
                                     smooth=False, shader='assembly', glOptions='opaque')
                self.object_items.append((obj, item))
                self.addItem(item)

        self.pumpkin = None
        if hasattr(sim, 'pumpkinFrame'):
            sphere = gl.MeshData.sphere(rows=12, cols=24, radius=sim.pumpkinFrame.pumpkin.radius)
            self.pumpkin = gl.GLMeshItem(meshdata=sphere, color=(1, 0.42, 0.06, 1),
                                         smooth=True, shader='assembly', glOptions='opaque')
            self.addItem(self.pumpkin)
        self.sling = gl.GLLinePlotItem(pos=np.zeros((2, 3)), color=(0.25, 0.28, 0.33, 1), width=2)
        self.forces = gl.GLLinePlotItem(pos=np.empty((0, 3)), color=(0.85, 0.1, 0.12, 1),
                                        width=2, mode='lines')
        self.trajectory = gl.GLLinePlotItem(pos=np.empty((0, 3)), color=(0.9, 0.4, 0.05, 0.65), width=1)
        for item in (self.sling, self.forces, self.trajectory):
            self.addItem(item)
        grid = gl.GLGridItem(color=(100, 110, 125, 65))
        grid.setSize(18, 12)
        grid.setSpacing(1, 1)
        grid.translate(2, 0, ground_z)
        self.addItem(grid)
        self.reset_camera()
        self.update_scene()

    def sizeHint(self):
        return QtCore.QSize(850, 520)

    def reset_camera(self):
        self._manual_camera = False
        self.opts['center'] = QtGui.QVector3D(*self.camera_center.tolist())
        self.setCameraPosition(elevation=22, azimuth=-65)
        self.setCameraPosition(distance=self.fit_distance())

    def side_camera(self):
        self._manual_camera = False
        self.opts['center'] = QtGui.QVector3D(*self.camera_center.tolist())
        self.setCameraPosition(elevation=0, azimuth=-90)
        self.setCameraPosition(distance=self.fit_distance())

    def fit_distance(self):
        aspect = max(0.1, self.width() / max(1, self.height()))
        horizontal = np.tan(np.deg2rad(self.opts['fov'] / 2))
        vertical = horizontal / aspect
        transform = self.viewMatrix()
        required = []
        for corner in self.camera_corners:
            point = transform.map(QtGui.QVector3D(*corner.tolist()))
            depth = point.z() + self.opts['distance']
            required.append(depth + max(abs(point.x()) / horizontal, abs(point.y()) / vertical))
        return max(required) * 1.08

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if hasattr(self, 'camera_corners') and not self._manual_camera:
            self.setCameraPosition(distance=self.fit_distance())

    def mouseMoveEvent(self, event):
        if event.buttons():
            self._manual_camera = True
        super().mouseMoveEvent(event)

    def wheelEvent(self, event):
        self._manual_camera = True
        super().wheelEvent(event)

    def update_scene(self):
        for body, item in self.body_items.items():
            matrix = body_matrix(self.sim, self.metadata['bindings'][body])
            item.setTransform(QtGui.QMatrix4x4(*matrix.ravel().tolist()))
        if 'machine' in self.body_items:
            self.body_items['machine'].setVisible(self.frame_button.isChecked())
        for obj, item in self.object_items:
            item.setTransform(QtGui.QMatrix4x4(*object_matrix(obj).ravel().tolist()))
        if self.pumpkin is not None:
            pos = self.sim.pumpkinFrame.pumpkin.obj2world(self.sim.pumpkinFrame.pumpkin.cg).A1
            self.pumpkin.resetTransform()
            self.pumpkin.translate(pos[0], 0, pos[1])
            index = 6 * self.sim.pumpkinFrame.idx
            start = getattr(self.sim, 'release_index', self.sim.time_idx)
            states = self.sim.Y[start:self.sim.time_idx + 1]
            self.trajectory.setData(pos=np.column_stack((states[:, index], np.zeros(len(states)),
                                                         states[:, index + 1])))
        if hasattr(self.sim, 'sling'):
            sling = self.sim.sling
            a, b = (sling.frame0.frame2world(sling.xframe0).A1,
                    sling.frame1.frame2world(sling.xframe1).A1)
            self.sling.setData(pos=np.array([[a[0], 0, a[1]], [b[0], 0, b[1]]]))
            self.sling.setVisible(bool(sling.enabled))
        else:
            self.sling.setVisible(False)
        self.forces.setVisible(self.draw_forces_button.isChecked())
        segments = []
        if self.draw_forces_button.isChecked():
            def add_force(location, force):
                start = np.array([location[0], 0, location[1]])
                vector = np.asarray(force).reshape(-1)[:2] / 400e3
                end = start + [vector[0], 0, vector[1]]
                segments.extend((start, end))
                norm = np.linalg.norm(vector)
                if norm > 1e-8:
                    direction = (end - start) / norm
                    cross = np.array([-direction[2], 0, direction[0]])
                    size = min(0.08, norm / 3)
                    for sign in (-1, 1):
                        segments.extend((end, end - size * direction + sign * size * cross / 2))

            for constraint in self.sim.constraints:
                if not constraint.enabled or constraint.forces is None:
                    continue
                if hasattr(constraint, 'frame0'):
                    locations = [constraint.frame0.frame2world(constraint.xframe0).A1,
                                 constraint.frame1.frame2world(constraint.xframe1).A1]
                elif hasattr(constraint, 'xframe'):
                    locations = [constraint.frame.frame2world(constraint.xframe).A1]
                else:
                    continue
                for location, force in zip(locations, constraint.forces):
                    add_force(location, force.A1)
            for spring in self.sim.springs:
                if hasattr(spring, 'x_frame'):
                    add_force(spring.frame.frame2world(spring.x_frame).A1, spring.force().A1)
            self.forces.setData(pos=np.array(segments, dtype=float).reshape(-1, 3))
