"""Search sling length for the 2026 configuration, retaining historical physics."""

import argparse
import contextlib
import io
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "treb_sim" / "src"))

import numpy as np
from scipy.optimize import minimize_scalar

from first_in_fright_2026 import treb, find_release, PUMPKIN_WEIGHT_LB, SLING_WEIGHT_LB, ARM_TIP_HARDWARE_WEIGHT_LB
from dynamics.constants import meter2foot


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--time-step', type=float, default=0.0001)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--elevation', type=float, default=30.0)
    parser.add_argument('--unrestricted', action='store_true',
                        help='Optimize range without a firing-angle constraint')
    parser.add_argument('--refine-bounds', nargs=2, type=float, metavar=('LOW', 'HIGH'),
                        help='Refine a previously screened interval without repeating the broad search')
    args = parser.parse_args()
    trials = []
    elevation = None if args.unrestricted else args.elevation

    def objective(length):
        with contextlib.redirect_stdout(io.StringIO()):
            sim = treb(sling_length=float(length), release_at_max_range=False,
                       time_step=args.time_step)
        launch = find_release(sim, elevation)
        trial = dict(sling_length_ft=float(length),
                     range_ft=float(meter2foot(launch['range_m'])),
                     release_time_s=launch['time_s'],
                     launch_angle_deg=launch['angle_deg'],
                     arm_kinetic_energy_j=launch['arm_kinetic_energy_j'],
                     arm_rotational_energy_j=launch['arm_rotational_energy_j'],
                     pumpkin_speed_m_s=launch['pumpkin_speed_m_s'])
        trials.append(trial)
        if args.output:
            args.output.with_suffix('.progress.json').write_text(json.dumps(trials, indent=2) + '\n', encoding='utf-8')
        print(f"{length:.6f} ft: {trial['range_ft']:.3f} ft "
              f"at {launch['time_s']:.6f} s, "
              f"{launch['angle_deg']:.3f} degrees, "
              f"arm KE {launch['arm_kinetic_energy_j']:.1f} J", flush=True)
        return -trial['range_ft']

    # Broad search first, then refine around the best interior coarse sample.
    if args.refine_bounds:
        bounds = tuple(args.refine_bounds)
        if not 0 < bounds[0] < bounds[1]:
            parser.error('Refinement bounds must be positive and increasing')
        search_bounds = list(bounds)
    else:
        lengths = np.arange(4.0, 16.01, 2.0)
        values = [objective(length) for length in lengths]
        index = int(np.argmin(values))
        if index in (0, len(lengths) - 1):
            raise RuntimeError("Optimum is at search boundary; widen the search")
        bounds = (lengths[index - 1], lengths[index + 1])
        search_bounds = [4.0, 16.0]
    result = minimize_scalar(objective,
                             bounds=bounds,
                             method='bounded', options={'xatol': 0.002})
    if not result.success:
        raise RuntimeError(result.message)
    best = max(trials, key=lambda trial: trial['range_ft'])
    report = dict(configuration="2013 machine with 2012 counterweight",
                  pumpkin_weight_lb=PUMPKIN_WEIGHT_LB, sling_weight_lb=SLING_WEIGHT_LB,
                  arm_tip_hardware_weight_lb=ARM_TIP_HARDWARE_WEIGHT_LB,
                  counterweight_weight_lb=4581.0, time_step_s=args.time_step,
                  launch_elevation_deg=elevation, objective='Maximum predicted range at the specified elevation',
                  spin_lift_modeled=False,
                  search_bounds_ft=search_bounds, length_tolerance_ft=0.002,
                  range_definition="Historical drag prediction to axle-height landing plane",
                  best=best, trials=trials)
    print(json.dumps(report, indent=2))
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
