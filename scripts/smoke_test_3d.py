"""Validate 3D assembly motion and native OpenGL rendering (requires a desktop)."""

import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'treb_sim' / 'src'))

import numpy as np
from PyQt5.QtCore import QTimer
from PyQt5.QtWidgets import QApplication, QCheckBox

import first_in_fright_2026 as model
from dynamics.animation_3d import ASSETS, body_matrix, square_tube


def check_anchors(sim, metadata):
    anchors = metadata['anchors']
    # Each source segment must remain attached at both physical joints as it
    # moves, rather than merely looking aligned in one camera orientation.
    segments = [('machine', 'axle', 'hinge', 'hinge', 1),
                ('weight', 'hinge', 'counterweight_joint', 'rampPin', 0),
                ('upper_link', 'hanger', 'link_joint', 'linkPin', 0),
                ('lower_link', 'counterweight_joint', 'link_joint', 'linkPin', 1),
                ('arm', 'axle', 'short_arm_joint', 'connectorPin1', 0),
                ('connector', 'short_arm_joint', 'link_joint', 'connectorPin2', 1)]
    for body, pivot, end, joint_name, side in segments:
        local = np.array(anchors[end]) - np.array(anchors[pivot])
        mapped = body_matrix(sim, metadata['bindings'][body]) @ np.r_[local, 1]
        constraint = getattr(sim, joint_name)
        frame = getattr(constraint, f'frame{side}')
        actual = frame.frame2world(getattr(constraint, f'xframe{side}')).A1
        np.testing.assert_allclose(mapped[[0, 2]], actual, atol=0.0005)
        np.testing.assert_allclose(mapped[1], 0, atol=1e-8)


def main():
    app = QApplication([])
    sim = model.treb(time_step=0.001)
    metadata = json.loads((ASSETS / 'sketchup_3d.json').read_text(encoding='utf-8'))
    assert not any('wheel' in p['path'].lower() or 'tire' in p['path'].lower()
                   for p in metadata['included'])
    points, faces = square_tube(sim.connectorFrame.stiffener.l, sim.connectorFrame.stiffener.w)
    np.testing.assert_allclose(np.ptp(points, axis=0)[1], np.ptp(points, axis=0)[2])
    assert len(np.unique(np.abs(points[:, 1]))) == 2, 'Tube must have a hollow square section'
    view = model.Animation(sim, model.display_range)
    view.resize(1100, 950)
    state_copy = sim.Y.copy()
    output = Path(__file__).resolve().parents[1] / 'models' / 'extracted'
    output.mkdir(parents=True, exist_ok=True)

    def verify():
        try:
            assert view.drawing.context().isValid(), 'OpenGL context did not initialize'
            for index in (0, sim.release_index - 1, sim.release_index, len(sim.t) - 1):
                view.time_slider._slider.setValue(index)
                app.processEvents()
                check_anchors(sim, metadata)
                assert view.drawing.sling.visible() == bool(sim.sling.enabled)
                assert not view.drawing.grabFramebuffer().isNull()
            view.time_slider._slider.setValue(0)
            app.processEvents()
            view.grab().save(str(output / 'viewer_3d_cocked.png'))
            view.time_slider._slider.setValue(sim.release_index)
            view.findChild(QCheckBox).setChecked(True)
            app.processEvents()
            view.grab().save(str(output / 'viewer_3d_release.png'))
            view.frame_button.click()
            assert not view.drawing.body_items['machine'].visible()
            assert view.drawing.object_items[1][1].visible(), 'Square tube must remain visible'
            view.drawing.orbit(90, 15)
            app.processEvents()
            view.grab().save(str(output / 'viewer_3d_frame_hidden.png'))
            np.testing.assert_array_equal(sim.Y, state_copy)
            view.frame_button.click()
            view.drawing.side_camera()
            app.processEvents()
            assert not view.drawing.grabFramebuffer().isNull()
            print('3D: joint alignment, native GL rendering, orbit, frame toggle, tube and release OK')
            view.close()
            app.exit(0)
        except Exception:
            import traceback
            traceback.print_exc()
            view.close()
            app.exit(1)

    QTimer.singleShot(1500, verify)
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
