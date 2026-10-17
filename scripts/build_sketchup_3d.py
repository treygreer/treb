"""Bake rigid SketchUp assemblies: uv run --with openskp python scripts/build_sketchup_3d.py."""

import hashlib
import json
from pathlib import Path

import numpy as np
from openskp import SkpFile


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'models' / 'treb.skp'
OUT = ROOT / 'treb_sim' / 'assets'


def main():
    source = SkpFile.open(str(SOURCE))
    source.parse()
    scene = source.build_scene()
    primitives = {scene.mesh_index[p.geom_name].path: p for p in scene.glb_primitives}

    def native(primitive):
        v = np.asarray(primitive.positions, dtype=float).reshape(-1, 3)
        return np.column_stack((v[:, 0], -v[:, 2], v[:, 1]))

    def center(path):
        v = native(primitives[path])
        return (v.min(0) + v.max(0)) / 2

    axle = center('ROOT / Frame#1 / Component_4446')
    hinge = center('ROOT / Frame#1 / Component_4594')
    hanger = center('ROOT / Frame#1 / Component_4742')
    joint = center('ROOT / Component_8984 / Upper Link#1 / Component_9754')
    ramp = center('ROOT / Component_8984 / Lower Link#1 / Component_26948')
    short_arm = center('ROOT / Component_13205 / arm#1 / Component_13482')
    # All rigid-body pivots lie on the machine's center plane. The selected
    # source meshes are one side of each axle/bearing; preserve lateral widths.
    for point in (axle, hinge, hanger, joint, ramp, short_arm):
        point[0] = 0.0

    def angle(a, b):
        v = b - a
        return float(np.arctan2(v[2], v[1]))

    bindings = {
        'machine': dict(pivot=axle, reference_angle=0.0, constraint='axle', side=1),
        'arm': dict(pivot=axle, reference_angle=np.pi / 2, constraint='axle', side=0),
        'weight': dict(pivot=hinge, reference_angle=angle(hinge, ramp), constraint='hinge', side=0),
        'upper_link': dict(pivot=hanger, reference_angle=angle(hanger, joint), constraint='hanger', side=0),
        'lower_link': dict(pivot=ramp, reference_angle=angle(ramp, joint), constraint='rampPin', side=1),
        'connector': dict(pivot=short_arm, reference_angle=angle(short_arm, joint), constraint='connectorPin1', side=1),
    }
    parts = {name: [] for name in bindings}
    included, excluded = [], []
    for primitive in scene.glb_primitives:
        meta = scene.mesh_index[primitive.geom_name]
        path = meta.path.lower()
        body = None
        if any(word in path for word in ('wheel', 'tire')):
            excluded.append(meta.path)
            continue
        if meta.layer in ('frame', 'trailer', 'Layer0'):
            body = 'machine'
        elif meta.layer == 'cocked':
            for label, candidate in (('upper link', 'upper_link'), ('lower link', 'lower_link'),
                                     ('counterweight', 'weight'), ('connector half', 'connector'),
                                     (' / arm', 'arm')):
                if label in path:
                    body = candidate
                    break
        if body is None:
            excluded.append(meta.path)
            continue
        included.append(dict(path=meta.path, layer=meta.layer, body=body))
        pivot = bindings[body]['pivot']
        points = native(primitive) - pivot
        # Render coordinates: simulator horizontal X, lateral Y, simulator height Z.
        points = np.column_stack((points[:, 1], points[:, 0], points[:, 2]))
        faces = np.asarray(primitive.indices).reshape(-1, 3)
        triangles = points[faces]
        area = np.linalg.norm(np.cross(triangles[:, 1] - triangles[:, 0],
                                       triangles[:, 2] - triangles[:, 0]), axis=1)
        faces = faces[area > 1e-12]
        color = scene.gltf_materials[primitive.material_index].get(
            'pbrMetallicRoughness', {}).get('baseColorFactor', [0.65, 0.65, 0.65, 1.0])
        parts[body].append((points, faces, np.tile(color, (len(faces), 1))))

    arrays = {}
    for body, meshes in parts.items():
        assert meshes, f'Missing SketchUp assembly: {body}'
        vertexes, faces, colors, offset = [], [], [], 0
        for v, f, c in meshes:
            vertexes.append(v)
            faces.append(f + offset)
            colors.append(c)
            offset += len(v)
        arrays[f'{body}_vertices'] = np.concatenate(vertexes).astype(np.float32)
        arrays[f'{body}_faces'] = np.concatenate(faces).astype(np.uint32)
        arrays[f'{body}_colors'] = np.concatenate(colors).astype(np.float32)
    anchors = {name: ((point - axle)[[1, 0, 2]]).tolist() for name, point in
               dict(axle=axle, hinge=hinge, hanger=hanger, link_joint=joint,
                    counterweight_joint=ramp, short_arm_joint=short_arm).items()}
    metadata = dict(units='metres', source='models/treb.skp',
                    source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                    axes='X horizontal, Y across machine, Z vertical', anchors=anchors,
                    bindings={name: {k: v for k, v in b.items() if k != 'pivot'}
                              for name, b in bindings.items()}, included=included, excluded=excluded)
    OUT.mkdir(exist_ok=True)
    np.savez_compressed(OUT / 'sketchup_3d.npz', **arrays)
    (OUT / 'sketchup_3d.json').write_text(json.dumps(metadata, indent=2) + '\n', encoding='utf-8')
    print('Bound mesh parts:', {body: len(meshes) for body, meshes in parts.items()})
    print('Wheels, tires, fired duplicates and helper geometry omitted')


if __name__ == '__main__':
    main()
