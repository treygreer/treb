"""2026: 2013 machine, 2012 counterweight, 8 lb pumpkin, 11 oz sling, 9 oz tip hardware."""

import sys

import numpy as np
from scipy.integrate import ode
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from PyQt5.QtWidgets import QApplication

import first_in_fright_2013 as machine
from dynamics.animation import Animation
from dynamics.constants import meter2foot


# Maximum range at 30 degrees; see scripts/optimize_2026.py and the saved search.
OPTIMUM_SLING_LENGTH = 9.54  # feet; near the broad optimum, rounded for use
LAUNCH_ELEVATION = 30.0  # degrees above horizontal
PUMPKIN_WEIGHT_LB = 8.0
SLING_WEIGHT_LB = 11.0 / 16.0
ARM_TIP_HARDWARE_WEIGHT_LB = 9.0 / 16.0


def treb(sling_length=OPTIMUM_SLING_LENGTH, *, release_at_max_range=True,
         launch_elevation=LAUNCH_ELEVATION, time_step=0.0001,
         sim_duration=2.0, debug=False, **kwargs):
    parameters = dict(cw_weight=4581.0, cw_moment_arm=10.41,
                      cw_moment=3.516e6, pumpkin_weight=PUMPKIN_WEIGHT_LB, sling_weight=SLING_WEIGHT_LB,
                      arm_tip_hardware_weight=ARM_TIP_HARDWARE_WEIGHT_LB)
    parameters.update(kwargs)
    # First run with the sling attached to find the best release state.
    sim = machine.treb(sling_length=sling_length, release_time=0.0,
                       time_step=time_step, sim_duration=sim_duration,
                       debug=debug, **parameters)
    sim.sling_length_feet = sling_length
    sim.pumpkin_weight_lb = parameters['pumpkin_weight']
    sim.sling_weight_lb = parameters['sling_weight']
    sim.arm_tip_hardware_weight_lb = parameters['arm_tip_hardware_weight']
    if release_at_max_range:
        _release_at_best_range(sim, launch_elevation)
    return sim


def release_metrics(sim, time, state):
    """Evaluate a candidate using the pumpkin velocity and rigid arm's full KE."""
    pumpkin = state[6 * sim.pumpkinFrame.idx:6 * sim.pumpkinFrame.idx + 6]
    arm = state[6 * sim.armFrame.idx:6 * sim.armFrame.idx + 6]
    rotational = 0.5 * sim.armFrame.moment * arm[5]**2
    translational = 0.5 * sim.armFrame.mass * np.dot(arm[3:5], arm[3:5])
    return dict(time_s=float(time), state=state.copy(),
                angle_deg=float(np.rad2deg(np.arctan2(pumpkin[4], pumpkin[3]))),
                range_m=float(machine.Y2range(sim, state).item()),
                arm_kinetic_energy_j=float(rotational + translational),
                arm_rotational_energy_j=float(rotational),
                pumpkin_speed_m_s=float(np.linalg.norm(pumpkin[3:5])))


def find_release(sim, launch_elevation=None):
    """Find maximum range at the requested elevation, or unrestricted if None.

    Search every descending crossing of the specified firing angle. Interpolate
    the saved attached trajectory so the release is at the actual angle, rather
    than at the first output sample after it.
    """
    if launch_elevation is None:
        index = int(np.argmax(sim.range))
        return release_metrics(sim, sim.t[index], sim.Y[index])
    if not 0 < launch_elevation < 90:
        raise ValueError('Launch elevation must be between 0 and 90 degrees')
    index = 6 * sim.pumpkinFrame.idx
    vx, vy = sim.Y[:, index + 3], sim.Y[:, index + 4]
    angles = np.rad2deg(np.arctan2(vy, vx))
    forward = (vx > 0) & (vy > 0)
    crossings = np.flatnonzero(forward[:-1] & forward[1:] &
                              (angles[:-1] >= launch_elevation) &
                              (angles[1:] <= launch_elevation))
    if not len(crossings):
        raise ValueError('No descending forward launch at the requested elevation')
    trajectory = CubicSpline(sim.t, sim.Y, axis=0)

    def angle_error(time):
        state = trajectory(time)
        return np.rad2deg(np.arctan2(state[index + 4], state[index + 3])) - launch_elevation

    candidates = []
    for crossing in crossings:
        time = brentq(angle_error, sim.t[crossing], sim.t[crossing + 1], xtol=1e-13)
        candidates.append(release_metrics(sim, time, trajectory(time)))
    return max(candidates, key=lambda candidate: (candidate['range_m'],
                                                 -candidate['arm_kinetic_energy_j']))


def _release_at_best_range(sim, launch_elevation):
    """Keep the attached trajectory through the selected launch, then release.

    A new integrator starts at the saved launch state. This avoids continuing
    VODE with derivative history from the attached sling after release.
    The historical range calculation (including its landing plane) is retained.
    """
    launch = find_release(sim, launch_elevation)
    release_index = int(np.searchsorted(sim.t, launch['time_s']))
    if release_index == 0 or release_index == len(sim.t) - 1:
        raise ValueError("Best release is at the trajectory boundary; extend the run")
    sim.release_index = release_index
    pre_release_range = sim.range[:release_index].copy()
    sim.release_time = launch['time_s']
    sim.tmax = sim.release_time
    sim.maxrange = launch['range_m']
    sim.launchDegrees = launch['angle_deg']
    sim.arm_energy_at_release = launch['arm_kinetic_energy_j']

    states = np.empty((sim.num_steps, sim.Y.shape[1]))
    states[:release_index + 1] = sim.Y[:release_index + 1]
    states[release_index] = launch['state']
    times = np.arange(sim.num_steps) * sim.time_step
    times[release_index] = sim.release_time
    enabled = np.zeros((sim.num_steps, len(sim.constraints)))
    enabled[:release_index + 1] = sim.constraints_enabled[:release_index + 1]
    for constraint, active in zip(sim.constraints, enabled[release_index]):
        constraint.enabled = bool(active)
    sim.sling.enabled = False
    enabled[release_index, sim.sling.idx] = False

    integrator = ode(sim.deriv).set_integrator('vode', rtol=1e-9, atol=1e-6)
    integrator.set_initial_value(states[release_index], sim.release_time)
    for index in range(release_index + 1, sim.num_steps):
        integrator.integrate(times[index])
        if not integrator.successful():
            raise RuntimeError(f"Integration failed at {times[index]:g} seconds")
        states[index] = integrator.y
        # Evaluate contact decisions at the accepted state, rather than at the
        # last internal derivative evaluation left behind by the integrator.
        sim.deriv(times[index], states[index])
        previous = [constraint.enabled for constraint in sim.constraints]
        machine.continue_sim(sim, times[index], states[index])
        current = [constraint.enabled for constraint in sim.constraints]
        enabled[index] = current
        if current != previous:
            integrator.set_initial_value(states[index], times[index])

    sim.Y, sim.t, sim.constraints_enabled = states, times, enabled
    # Recompute diagnostic forces at each stored state for this new model.
    sim.F = np.zeros((sim.num_steps, len(sim.constraints), 2, 2))
    for index, time in enumerate(times):
        for constraint, active in zip(sim.constraints, enabled[index]):
            constraint.enabled = bool(active)
        sim.deriv(time, states[index])
        for constraint in sim.constraints:
            if constraint.enabled:
                for frame_index, force in enumerate(constraint.forces):
                    sim.F[index, constraint.idx, frame_index] = force.A1
    sim.Fmax = float(np.max(sim.hanger.Fvec()))
    # After release the viewer reports the prediction at release, rather than
    # predicting a new launch from a later, gravity-only flight state.
    sim.range = np.empty(sim.num_steps)
    sim.range[:release_index] = pre_release_range
    sim.range[release_index:] = sim.maxrange
    for frame in sim.frames:
        frame.PEmin = float(np.min(frame.PEvec()))
    print(f"2026: sling={sim.sling_length_feet:.4f} ft, "
          f"release={sim.release_time:.6f} s, "
          f"range={meter2foot(sim.maxrange):.2f} ft, "
          f"launch angle={sim.launchDegrees:.2f} degrees, "
          f"arm KE={sim.arm_energy_at_release:.1f} J")


def display_range(sim, states):
    """Range adapter for the viewer's full trajectory and current-state calls."""
    if states.ndim == 2:
        return sim.range
    return sim.range[sim.time_idx]


if __name__ == '__main__':
    sim = treb(debug=True)
    app = QApplication(sys.argv)
    app.setApplicationName("First in Fright 2026")
    anim = Animation(sim, display_range)
    anim.setWindowTitle(f"First in Fright 2026 - {sim.pumpkin_weight_lb:g} lb pumpkin - {sim.launchDegrees:.0f} degree launch")
    sys.exit(app.exec_())
