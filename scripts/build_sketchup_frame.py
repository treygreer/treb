"""Bake stationary SketchUp geometry for the Qt viewer.

Run: uv run --with openskp python scripts/build_sketchup_frame.py
The generated asset uses metres, relative to the throwing-arm axle.
"""

import hashlib
import json
from pathlib import Path

import numpy as np
from openskp import SkpFile
from shapely import Polygon, unary_union


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'models' / 'treb.skp'
DESTINATION = ROOT / 'treb_sim' / 'assets' / 'sketchup_frame.json'


def main():
    source = SkpFile.open(str(SOURCE))
    source.parse()
    scene = source.build_scene()
    primitives = {scene.mesh_index[p.geom_name].path: p for p in scene.glb_primitives}

    def native_points(primitive):
        points = np.asarray(primitive.positions, dtype=float).reshape(-1, 3)
        # OpenSKP's render geometry is metres with Y up and Z = -SketchUp Y.
        return np.column_stack((points[:, 0], -points[:, 2], points[:, 1]))

    def center(path):
        points = native_points(primitives[path])
        return (points.min(axis=0) + points.max(axis=0)) / 2

    axle = center('ROOT / Frame#1 / Component_4446')
    hinge = center('ROOT / Frame#1 / Component_4594')
    hanger = center('ROOT / Frame#1 / Component_4742')
    anchors = dict(axle=[0.0, 0.0], hinge=(hinge[1:] - axle[1:]).tolist(),
                   hanger=(hanger[1:] - axle[1:]).tolist())
    # Three independent pivot centers agree with the simulator without scaling.
    np.testing.assert_allclose(anchors['hinge'], [8 * 0.0254, -48 * 0.0254], atol=1e-4)
    np.testing.assert_allclose(anchors['hanger'], [11.38508 * 0.3048, -24 * 0.0254], atol=1e-4)

    polygons_by_color = {}
    included, excluded = [], []
    for primitive in scene.glb_primitives:
        meta = scene.mesh_index[primitive.geom_name]
        if meta.layer not in ('frame', 'trailer', 'Layer0'):
            excluded.append(dict(path=meta.path, layer=meta.layer, reason='moving assembly or helper'))
            continue
        if any(word in meta.path.lower() for word in ('wheel', 'tire')):
            excluded.append(dict(path=meta.path, layer=meta.layer, reason='wheel or tire assembly'))
            continue
        included.append(dict(path=meta.path, layer=meta.layer))
        points = native_points(primitive)[:, 1:] - axle[1:]
        faces = np.asarray(primitive.indices).reshape(-1, 3)
        rgba = scene.gltf_materials[primitive.material_index].get(
            'pbrMetallicRoughness', {}).get('baseColorFactor', [0.6, 0.6, 0.6, 1])
        color = tuple(round(float(channel), 5) for channel in rgba[:3])
        for triangle in points[faces]:
            polygon = Polygon(triangle)
            if polygon.is_valid and polygon.area > 1e-10:
                polygons_by_color.setdefault(color, []).append(polygon)

    shapes = []
    for color, polygons in polygons_by_color.items():
        merged = unary_union(polygons).simplify(0.00005, preserve_topology=True)
        pieces = list(merged.geoms) if merged.geom_type == 'MultiPolygon' else [merged]
        for polygon in pieces:
            shapes.append(dict(color=list(color), exterior=np.round(polygon.exterior.coords, 6).tolist(),
                               holes=[np.round(ring.coords, 6).tolist() for ring in polygon.interiors]))
    bounds = unary_union([p for parts in polygons_by_color.values() for p in parts]).bounds
    asset = dict(source='models/treb.skp', source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                 units='metres', projection='SketchUp Y/Z side elevation, relative to main arm axle',
                 source_axle_inches=(axle / 0.0254).tolist(), anchors=anchors,
                 bounds=list(bounds), included=included, excluded=excluded, shapes=shapes)
    DESTINATION.parent.mkdir(exist_ok=True)
    DESTINATION.write_text(json.dumps(asset, indent=2) + '\n', encoding='utf-8')
    print(f"Saved {len(shapes)} shapes from {len(included)} stationary mesh parts to {DESTINATION}")
    print('Anchors (m):', anchors)


if __name__ == '__main__':
    main()
