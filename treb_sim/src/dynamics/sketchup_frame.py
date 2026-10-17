"""Drawing-only SketchUp side elevation; it does not change physical objects."""

import json
from pathlib import Path

from PyQt5 import QtCore, QtGui, QtWidgets


ASSET = Path(__file__).resolve().parents[2] / 'assets' / 'sketchup_frame.json'


class SketchupFrame:
    def __init__(self):
        data = json.loads(ASSET.read_text(encoding='utf-8'))
        left, bottom, right, top = data['bounds']
        self.bounds = QtCore.QRectF(left, bottom, right - left, top - bottom)
        self.paths = []
        for shape in data['shapes']:
            path = QtGui.QPainterPath()
            path.setFillRule(QtCore.Qt.OddEvenFill)
            for ring in (shape['exterior'], *shape['holes']):
                path.moveTo(*ring[0])
                for point in ring[1:]:
                    path.lineTo(*point)
                path.closeSubpath()
            self.paths.append((path, QtGui.QColor.fromRgbF(*shape['color'])))

    def draw(self, scene, sim):
        frame = sim.machineFrame
        axle = frame.frame2world(sim.axle.xframe1).A1
        rotation = frame.R.A
        transform = QtGui.QTransform(rotation[0, 0], rotation[1, 0],
                                     rotation[0, 1], rotation[1, 1], axle[0], axle[1])
        group = QtWidgets.QGraphicsItemGroup()
        group.setData(0, 'sketchup_frame')
        scene.addItem(group)
        pen = QtGui.QPen(QtGui.QColor(65, 57, 49))
        pen.setCosmetic(True)
        pen.setWidthF(0.7)
        for path, color in self.paths:
            item = QtWidgets.QGraphicsPathItem(path)
            item.setPen(pen)
            item.setBrush(QtGui.QBrush(color))
            group.addToGroup(item)
        group.setTransform(transform)
        group.setZValue(-10)
        return group
