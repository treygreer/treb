"""Extract the SketchUp model: uv run --with openskp python scripts/extract_sketchup.py."""

from dataclasses import asdict
import importlib.metadata
import json
from pathlib import Path
import zipfile

import numpy as np
import openskp
import trimesh
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'models' / 'treb.skp'
OUT = ROOT / 'models' / 'extracted'


def definition_summary(definition):
    points = np.array([(v.x, v.y, v.z) for v in definition.vertices.values()])
    return dict(id=definition.id, name=definition.name,
                vertex_count=len(definition.vertices), face_count=len(definition.faces),
                bounds_inches=dict(min=points.min(axis=0).tolist(),
                                   max=points.max(axis=0).tolist()) if len(points) else None,
                instances=[asdict(i) for i in definition.instances],
                texts=[asdict(t) for t in definition.texts],
                dimensions=[asdict(d) for d in definition.dimensions])


def main():
    OUT.mkdir(exist_ok=True)
    source = openskp.SkpFile.open(str(SOURCE))
    model = source.parse()
    scene = source.build_scene()
    report = dict(source='models/treb.skp', parser='openskp',
                  parser_version=importlib.metadata.version('openskp'),
                  sketchup_version=model.version, source_units=model.units,
                  source_geometry_coordinates='SketchUp inches; Z up',
                  export_geometry_coordinates='metres; Y up',
                  definition_count=len(model.definitions),
                  prototype_vertex_count=sum(len(d.vertices) for d in model.definitions.values()),
                  prototype_face_count=sum(len(d.faces) for d in model.definitions.values()),
                  layers=[asdict(l) for l in model.layers],
                  pages=[asdict(p) for p in model.pages],
                  materials=[dict(id=k, name=v.name) for k, v in model.materials_by_id.items()],
                  dimensions=[dict(**asdict(d), endpoint_distance_inches=float(
                      np.linalg.norm(np.asarray(d.a) - np.asarray(d.b)))) for d in model.dimensions],
                  root=definition_summary(model.root),
                  definitions=[definition_summary(d) for d in model.definitions.values()])

    with zipfile.ZipFile(SOURCE) as archive:
        for name in archive.namelist():
            if name == 'meta/model_thumbnail.png':
                (OUT / 'thumbnail.png').write_bytes(archive.read(name))
            elif name.startswith(('matched_photos/', 'scene_thumbnails/')) and name.endswith(('.jpg', '.png')):
                target = OUT / name
                target.parent.mkdir(exist_ok=True)
                target.write_bytes(archive.read(name))

    visible = trimesh.Scene()
    counterweight = trimesh.Scene()
    concrete_meshes = []
    cw_triangles = []
    for primitive in scene.glb_primitives:
        meta = scene.mesh_index[primitive.geom_name]
        if scene.layer_hidden.get(meta.layer, False):
            continue
        points = np.asarray(primitive.positions).reshape(-1, 3)
        faces = np.asarray(primitive.indices).reshape(-1, 3)
        mesh = trimesh.Trimesh(vertices=points, faces=faces, process=True)
        material = scene.gltf_materials[primitive.material_index]
        rgba = material.get('pbrMetallicRoughness', {}).get('baseColorFactor', [0.6, 0.6, 0.6, 1])
        mesh.visual.face_colors = np.array(rgba) * 255
        visible.add_geometry(mesh, geom_name=primitive.geom_name)
        if 'counterweight' in meta.path and meta.layer == 'cocked':
            counterweight.add_geometry(mesh.copy(), geom_name=primitive.geom_name)
            # These two instance paths identify the concrete halves in this model.
            is_concrete = ('Component_11421' in meta.path or 'Component_17264' in meta.path)
            if is_concrete:
                concrete_meshes.append(mesh)
            # Return from GLB coordinates to SketchUp Z-up, now in feet.
            feet = np.column_stack((points[:, 0], -points[:, 2], points[:, 1])) / 0.3048
            cw_triangles.append((feet[faces], '#a6a6a6' if is_concrete else '#356c94'))

    visible.export(OUT / 'treb_visible.glb')
    counterweight.export(OUT / 'counterweight.glb')
    concrete = trimesh.util.concatenate(concrete_meshes)
    extents = sorted((concrete.bounding_box_oriented.primitive.extents / 0.0254).tolist())
    report['counterweight_concrete'] = dict(
        oriented_envelope_inches=extents, watertight=bool(concrete.is_watertight),
        interpretation='Two mirrored concrete half-components in the cocked configuration',
        limitation='Envelope measurements only; no reliable solid volume or stored physical density')
    (OUT / 'model_inventory.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(projection='3d')
    for triangles, color in cw_triangles:
        ax.add_collection3d(Poly3DCollection(triangles, facecolor=color, edgecolor='#333333',
                                            linewidth=0.15, alpha=0.9))
    points = np.concatenate([t.reshape(-1, 3) for t, _ in cw_triangles])
    low, high = points.min(axis=0), points.max(axis=0)
    ax.set_xlim(low[0] - 0.3, high[0] + 0.3)
    ax.set_ylim(low[1] - 0.3, high[1] + 0.3)
    ax.set_zlim(low[2] - 0.3, high[2] + 0.3)
    ax.set_box_aspect(high - low)
    ax.view_init(elev=22, azim=-45)
    ax.set(xlabel='X (ft)', ylabel='Y (ft)', zlabel='Z (ft)',
           title='Counterweight assembly extracted from the 2010 SketchUp model')
    ax.legend(handles=[Patch(color='#a6a6a6', label='Concrete halves'),
                       Patch(color='#356c94', label='Supporting beam and hardware')])
    fig.tight_layout()
    fig.savefig(OUT / 'counterweight.png', dpi=160)
    plt.close(fig)
    print(f"Extracted {len(model.definitions)} definitions and {len(model.dimensions)} dimension records")
    print(f"Concrete oriented envelope: {extents} inches; watertight={concrete.is_watertight}")
    print(f"Saved inventory, photos, GLB geometry and counterweight figure to {OUT}")


if __name__ == '__main__':
    main()
