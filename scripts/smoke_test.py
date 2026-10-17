"""Exercise both historical models and their Qt views without a display."""

import os
import json
from pathlib import Path
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "treb_sim" / "src"))

import numpy as np
from PyQt5.QtCore import QTimer
from PyQt5.QtCore import QPointF
from PyQt5.QtWidgets import QApplication, QCheckBox

import first_in_fright_2012
import first_in_fright_2013
import first_in_fright_2026
from flight import Flight
from dynamics.sketchup_frame import ASSET


def check_frame_view(view, sim):
    asset = json.loads(ASSET.read_text(encoding='utf-8'))
    assert not any('wheel' in part['path'].lower() or 'tire' in part['path'].lower()
                   for part in asset['included'])
    assert {part['layer'] for part in asset['included']} <= {'frame', 'trailer', 'Layer0'}
    assert 'cocked' not in {part['layer'] for part in asset['included']}
    group = next(item for item in view.drawing.scene.items() if item.data(0) == 'sketchup_frame')
    for name in ('hinge', 'hanger'):
        mapped = group.sceneTransform().map(QPointF(*asset['anchors'][name]))
        constraint = getattr(sim, name)
        physical = sim.machineFrame.frame2world(constraint.xframe1).A1
        np.testing.assert_allclose([mapped.x(), mapped.y()], physical, atol=1e-4)
    states = sim.Y.copy()
    # The added square tube remains a physical and visible connector object.
    tube = sim.connectorFrame.stiffener
    assert tube in sim.connectorFrame.objects and tube.mass > 0
    view.frame_button.click()
    assert not any(item.data(0) == 'sketchup_frame' for item in view.drawing.scene.items())
    np.testing.assert_array_equal(sim.Y, states)
    view.frame_button.click()
    assert any(item.data(0) == 'sketchup_frame' for item in view.drawing.scene.items())


def main():
    app = QApplication([])
    for model in (first_in_fright_2012, first_in_fright_2013):
        sim = model.treb(debug=False)
        # Historical 3.516e6 lb*in^2 must produce about 1028.92 kg*m^2,
        # not the 31.98 kg*m^2 obtained by an extra division by 32.174.
        np.testing.assert_allclose(sim.weightFrame.cw.moment, 1028.921021, rtol=1e-5)
        assert np.isfinite(sim.Y).all(), "Non-finite simulation state"
        assert sim.maxrange > 0, "No positive launch range"
        assert 0 < sim.tmax < sim.t[-1], "Launch optimum outside simulation"
        view = model.Animation(sim, model.Y2range, view_mode='2d')
        view.findChild(QCheckBox).setChecked(True)
        for index in (0, len(sim.t) // 2, len(sim.t) - 1):
            view.time_slider._slider.setValue(index)
            app.processEvents()
            assert not view.grab().isNull(), "Qt rendering failed"
        view.close()
        print(f"{model.__name__}: simulation, range, slider and force rendering OK")

    sim = first_in_fright_2026.treb(time_step=0.001)
    np.testing.assert_allclose(sim.pumpkinFrame.pumpkin.mass, 8 * 0.45359)
    np.testing.assert_allclose(sim.pumpkinFrame.sling.mass, 0.45359 * 11 / 16)
    np.testing.assert_allclose(sim.armFrame.tip_hardware.mass, 0.45359 * 9 / 16)
    np.testing.assert_allclose(sim.armFrame.tip_hardware.origin, sim.armFrame.release_pin.origin)
    assert sim.armFrame.tip_hardware.moment == 0.0
    baseline = first_in_fright_2026.treb(time_step=0.001, release_at_max_range=False,
                                        arm_tip_hardware_weight=0.0)
    added_mass = sim.armFrame.tip_hardware.mass
    np.testing.assert_allclose(sim.armFrame.mass - baseline.armFrame.mass, added_mass)
    # Verify the physical inertia change about the axle, independent of the
    # shifted CG used internally by the rigid-body solver.
    def axle_inertia(run):
        return run.armFrame.moment + run.armFrame.mass * np.dot(run.axle.xframe0.A1, run.axle.xframe0.A1)
    distance = np.linalg.norm(sim.armFrame.tip_hardware.origin.A1 - sim.axle.xframe0.A1)
    np.testing.assert_allclose(axle_inertia(sim) - axle_inertia(baseline),
                               added_mass * distance**2, rtol=1e-10)
    np.testing.assert_allclose(sim.weightFrame.cw.mass, 4581 * 0.45359)
    np.testing.assert_allclose(sim.weightFrame.cw.moment, 1028.921021, rtol=1e-5)
    release = sim.release_index
    np.testing.assert_allclose(sim.launchDegrees, 30.0, atol=1e-7)
    assert np.all(np.diff(sim.t) > 0)
    assert sim.constraints_enabled[:release, sim.sling.idx].all()
    assert not sim.constraints_enabled[release:, sim.sling.idx].any()
    assert np.isfinite(sim.Y).all()
    # Removing an ideal sling must not impart an impulse. The subsequent free
    # trajectory must agree with ballistic motion from the saved release state.
    pumpkin = sim.Y[:, 6 * sim.pumpkinFrame.idx:6 * sim.pumpkinFrame.idx + 6]
    elapsed = sim.t[release:] - sim.release_time
    launch = pumpkin[release]
    np.testing.assert_allclose(pumpkin[release:, 3], launch[3], atol=1e-6)
    np.testing.assert_allclose(pumpkin[release:, 4], launch[4] - 9.80665 * elapsed,
                               atol=1e-6)
    np.testing.assert_allclose(pumpkin[release:, 1],
                               launch[1] + launch[4] * elapsed - 0.5 * 9.80665 * elapsed**2,
                               atol=1e-5)
    np.testing.assert_allclose(sim.range[release:], sim.maxrange)
    view = first_in_fright_2026.Animation(sim, first_in_fright_2026.display_range, view_mode='2d')
    view.findChild(QCheckBox).setChecked(True)
    for index in (release - 1, release, release + 1, len(sim.t) - 1):
        view.time_slider._slider.setValue(index)
        app.processEvents()
        np.testing.assert_allclose(view.time_slider.time, sim.t[index], atol=1e-12)
        assert sim.time_idx == index
        assert not view.grab().isNull()
    check_frame_view(view, sim)
    view.close()
    print("2026: requested masses, impulse-free release, fixed range and viewer OK")

    # With drag disabled, flight should match the analytic ballistic trajectory.
    flight = Flight(Cd=0, max_time=10)
    flight.run([0, 0], [10, 10])
    np.testing.assert_allclose(flight.Y[:, 0], 10 * flight.t, atol=1e-5)
    np.testing.assert_allclose(
        flight.Y[:, 1], 10 * flight.t - 0.5 * 9.80665 * flight.t**2, atol=1e-5
    )
    QTimer.singleShot(0, app.quit)
    assert app.exec_() == 0, "Qt event loop failed"
    print("Ballistic flight and Qt event loop OK")


if __name__ == "__main__":
    main()
