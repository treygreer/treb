"""Export the native 2026 run and 3D assets for the offline browser viewer."""
import base64
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'treb_sim' / 'src'))
import first_in_fright_2026 as model
from dynamics.animation_3d import ASSETS, body_matrix, object_matrix, cylinder, square_tube
from dynamics.constants import meter2foot


def packed(array, dtype='<f4'):
    array = np.asarray(array, dtype=dtype)
    assert np.isfinite(array).all()
    return dict(shape=list(array.shape), dtype=dtype,
                base64=base64.b64encode(array.tobytes()).decode('ascii'))


def main():
    sim = model.treb()
    metadata = json.loads((ASSETS / 'sketchup_3d.json').read_text())
    meshes = {}
    with np.load(ASSETS / 'sketchup_3d.npz', allow_pickle=False) as source:
        ground = float(source['machine_vertices'][:, 2].min())
        for name in metadata['bindings']:
            meshes[name] = {key: packed(source[f'{name}_{key}'], '<u4' if key == 'faces' else '<f4')
                            for key in ('vertices', 'faces', 'colors')}
    objects = [('rod', sim.connectorFrame.rod), ('tube', sim.connectorFrame.stiffener)]
    for name, obj in objects:
        vertices, faces = cylinder(obj.l, obj.w / 2) if name == 'rod' else square_tube(obj.l, obj.w)
        meshes[name] = dict(vertices=packed(vertices), faces=packed(faces, '<u4'),
                            color=[.35, .4, .46, 1])
    # Preserve the exact release and the full 0.1 ms resolution near it.
    indices = np.unique(np.r_[np.arange(0, sim.num_steps, 20), sim.num_steps - 1,
                              np.arange(max(0, sim.release_index - 100), sim.release_index + 101)])
    poses, pumpkins, slings, energies, forces = [], [], [], [], []
    for index in indices:
        for constraint, active in zip(sim.constraints, sim.constraints_enabled[index]):
            constraint.enabled = bool(active)
        sim.deriv(sim.t[index], sim.Y[index])
        matrices = [body_matrix(sim, binding) for binding in metadata['bindings'].values()]
        matrices += [object_matrix(obj) for _, obj in objects]
        poses.append(matrices)
        pos = sim.pumpkinFrame.pumpkin.obj2world(sim.pumpkinFrame.pumpkin.cg).A1
        pumpkins.append([pos[0], 0, pos[1]])
        ends = [frame.frame2world(anchor).A1 for frame, anchor in
                [(sim.sling.frame0, sim.sling.xframe0), (sim.sling.frame1, sim.sling.xframe1)]]
        slings.append([[p[0], 0, p[1]] for p in ends])
        energies.append([[float(frame.PE() - frame.PEmin), float(frame.KE())] for frame in sim.frames])
        segments = []

        def add_force(location, force):
            start = np.array([location[0], 0, location[1]])
            vector = np.asarray(force).reshape(-1)[:2] / 400e3
            end = start + [vector[0], 0, vector[1]]
            segments.extend((start, end))
            norm = np.linalg.norm(vector)
            if norm > 1e-8:
                direction = (end - start) / norm
                cross = np.array([-direction[2], 0, direction[0]])
                size = min(.08, norm / 3)
                for sign in (-1, 1):
                    segments.extend((end, end - size * direction + sign * size * cross / 2))

        for c in sim.constraints:
            if not c.enabled or c.forces is None:
                continue
            if hasattr(c, 'frame0'):
                locations = [c.frame0.frame2world(c.xframe0).A1, c.frame1.frame2world(c.xframe1).A1]
            elif hasattr(c, 'xframe'):
                locations = [c.frame.frame2world(c.xframe).A1]
            else:
                continue
            for location, force in zip(locations, c.forces):
                add_force(location, force.A1)
        for spring in sim.springs:
            if hasattr(spring, 'x_frame'):
                add_force(spring.frame.frame2world(spring.x_frame).A1, spring.force().A1)
        forces.append(packed(np.array(segments).reshape(-1, 3)))
    release = int(np.flatnonzero(indices == sim.release_index)[0])
    result = dict(title='First in Fright 2026', meshes=meshes, names=list(meshes),
                  frameNames=[f.name for f in sim.frames], ground=ground,
                  time=sim.t[indices].tolist(), poses=packed(poses), pumpkin=packed(pumpkins),
                  sling=packed(slings), energy=packed(energies), forces=forces,
                  rangeFeet=(meter2foot(sim.range[indices])).tolist(), releaseIndex=release,
                  releaseTime=sim.release_time, launchAngle=sim.launchDegrees,
                  maxRangeFeet=float(meter2foot(sim.maxrange)), slingLengthFeet=sim.sling_length_feet,
                  pumpkinWeightLb=sim.pumpkin_weight_lb,
                  slingWeightLb=sim.sling_weight_lb, armTipHardwareWeightLb=sim.arm_tip_hardware_weight_lb,
                  pumpkinRadius=sim.pumpkinFrame.pumpkin.radius,
                  totalEnergy=float(np.sum(energies[0])), timeStep=sim.time_step,
                  armEnergyAtRelease=sim.arm_energy_at_release, sourceSha256=metadata['source_sha256'])
    output = ROOT / 'web' / 'data' / 'viewer.json'
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, separators=(',', ':'), allow_nan=False), encoding='utf-8')
    print(f'Exported {len(indices)} samples, {output.stat().st_size:,} bytes to {output}')


if __name__ == '__main__':
    main()
