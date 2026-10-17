# treb

Trebuchet dynamics and flight simulator, with the historical 2012 and 2013
First in Fright models and a PyQt5 viewer.

## SketchUp model

The [First in Fright 2010 Trebuchet model on SketchUp 3D Warehouse](https://3dwarehouse.sketchup.com/model/2a951363504f4de367b420905e4d0003/First-in-Fright-2010-Trebuchet)
is available locally as [models/treb.skp](models/treb.skp).

The Qt viewer now renders the SketchUp geometry in **3D**, aligned at the
throwing-arm axle, counterweight hinge and linkage pins. The frame, arm,
counterweight, upper and lower links, and connecting rod assemblies follow
their simulated rigid-body motion. Wheels and tires are omitted. The
**Show Frame** toggle hides or shows the stationary structure. The simulator's
added square reinforcement tube remains on the connecting rod, rendered with
a hollow square section; its displayed wall thickness does not affect physics.
The model's title says 2010, but its aluminum throwing-arm tip
indicates that it includes later updates.

Drag to orbit, Ctrl+drag or middle-drag to pan, and use the wheel to zoom.
**Side View** returns to a side elevation; **Reset View** restores the 3D
camera. The time slider, forces and energy bars retain their previous behavior.
Only rendering is three-dimensional; the underlying dynamics remain planar.

Geometry is prebuilt in `treb_sim/assets/sketchup_3d.npz` with body bindings in
`sketchup_3d.json`, so SketchUp and its parser are not required to run the viewer.
The original SketchUp file, runtime geometry and saved web simulation data are
versioned. Extracted previews, screenshots, generated `web/site/` files and ZIP
downloads are local build output; the published viewer files live on `gh-pages`.
The older 2D rendering remains available programmatically with
`Animation(sim, range_function, view_mode='2d')`.
To regenerate the 3D geometry or the older side elevation:

```powershell
uv run --with openskp python scripts/build_sketchup_3d.py
uv run --with openskp python scripts/build_sketchup_frame.py
```

## Run

Install [uv](https://docs.astral.sh/uv/), then run from the repository root:

```powershell
uv sync
uv run python treb_sim/src/first_in_fright_2013.py
```

uv manages Python 3.13 and the dependencies in `.venv`. The 2013 example simulates
two seconds with release at 0.661 seconds, then opens the viewer. Use the slider
to inspect the motion, the checkbox to display forces, and the mouse wheel to
zoom.

To run the 2012 model:

```powershell
uv run python treb_sim/src/first_in_fright_2012.py
```

To run the 2026 model:

```powershell
uv run python treb_sim/src/first_in_fright_2026.py
```

The 2026 configuration uses the 2013 machine with the 2012 counterweight
assembly (4,581 lb, historical CG and corrected inertia), an 8 lb pumpkin,
an **11 oz sling**, and **9 oz of additional hardware at the throwing-arm tip**.
The added hardware is represented as a point mass at the sling attachment;
the rigid-arm calculation includes its mass, CG shift and parallel-axis inertia.
The historical 9 lb release-pin assembly remains, with the new 9 oz added to it.
Other geometry, including the 8-inch pumpkin diameter and drag
coefficient, follows the historical model. The default sling is 9 ft 6.5 in
(9.54 ft in the model), near
the maximum found by a sling-only search across 4–16 ft at a **30-degree firing
elevation**. Predicted range is approximately 2,177 ft, with release at
0.68313 seconds and approximately 1.63 kJ remaining in the arm. The viewer computes the maximum-range release state at that
elevation for the chosen configuration and then detaches the sling, restarting the
integrator at release to avoid an artificial impulse. Its range display stays
at the launch prediction after release; the displayed free trajectory retains
the historical gravity-only dynamics.

The 8 lb search is saved in `treb_sim/2026_optimization_8lb_11oz_9oz_tip.json`.
The completed broad screening and partial refinement are recorded in
`treb_sim/2026_8lb_coarse_search.txt`; the narrow refinement found approximately
9.538 ft. The broad command exceeded its execution time limit after locating
the peak, so the final refinement used `--refine-bounds 9.50 9.58`.
The earlier 6 lb results remain in `treb_sim/2026_optimization_11oz_9oz_tip.json`.
Reproduce it with:

```powershell
uv run python scripts/optimize_2026.py --elevation 30 --output treb_sim/2026_optimization_8lb_11oz_9oz_tip.json
```

This maximizes predicted range over sling length and the descending release
crossings of the specified elevation, with other machine parameters fixed.
It retains the historical drag calculation and axle-height landing plane.
The release state is interpolated between simulation samples to hit the
requested elevation, then the free trajectory starts at that precise time.
The time slider uses the saved timestamps, including the exact release event.
The flight sampling introduces small numerical variations near the maximum,
so the optimum should be read as about 9.5 ft.
This command saves search results; the default length is set explicitly in the
2026 source file. Historical launchers keep their existing behavior.

The earlier 1 lb sling configuration, without the added 9 oz, remains recorded
in `treb_sim/2026_optimization_30deg.json`. For that earlier configuration,
compared with the previous approximately 35-degree setting at the same rounded
sling length, the 30-degree launch leaves approximately **3.04 kJ** of kinetic
energy in the arm instead of **3.61 kJ**, a reduction of about **15.6%**, at a
predicted range cost of about **1.2%**. This includes the arm's translational and
rotational kinetic energy. The comparison and timestep checks are saved in
`treb_sim/2026_launch_comparison.json`.

The pumpkin's observed backspin lift is not modeled; the lower elevation is an
explicit launch constraint, not a calibrated prediction of lift. The original
unrestricted-angle search remains in `treb_sim/2026_optimization.json`; use
`--unrestricted` to run that objective again. For programmatic runs,
`treb(launch_elevation=None)` releases at unrestricted maximum range.

## Offline web viewer

The generated **2026 web viewer** is in `web/site/`. It uses the same 3D
assemblies, including the square connecting-rod reinforcement and omitted
wheels, with touch orbit/pan/zoom, frame and force toggles, time slider, playback,
energy bars, range display and a jump-to-release button.
The **2D / 3D** button switches between a flat side elevation and perspective
orbit view. Both show the same selected state and retain force/frame settings.
Drag to pan in 2D, and pinch or scroll to zoom; returning to 3D restores its
camera position.
If WebGL cannot start, the viewer uses software SVG rendering of the same 3D
assemblies. This retains camera controls but may be slower. Startup errors are
shown on the page rather than leaving the loading message indefinitely.

Open `web/site/treb-2026-offline.html` directly in a browser for the standalone
download. All JavaScript, geometry and simulation data are embedded in this
approximately 4 MB file; it makes no requests for companion files or libraries.
Android's handling of local HTML files varies by file manager, so the installable
version below is the preferred phone workflow.

To install in Chrome on Android, host the contents of `web/site/` together at
an **HTTPS** URL, open it once online, wait for **Offline ready**, then use
**Install for Offline Use** or Chrome's **Add to Home screen / Install** menu.
The service worker caches the viewer and download for later launches without
network access. The installed app needs neither the computer nor Python running.
Clearing Chrome's site storage also removes its offline cache; revisit online
to restore it.

### Website deployment handoff

The website and viewer are published separately:

- **Team homepage:** https://www.firstinfright.org/ from
  https://github.com/treygreer/firstinfright, root of `main`.
- **Simulation viewer:** https://treygreer.github.io/treb/ from this simulator
  repository's generated `gh-pages` branch, published from `/` with no custom domain.

The homepage links to the viewer, and the viewer links back to the homepage.
The simulator source checkout generates the viewer in `web/site/`; it is
published separately from the team website. The website checkout is in
`.site-deploy/firstinfright`; the viewer publishing checkout is in
`.site-deploy/treb-pages`. These local deployment checkouts are ignored by the
simulator repository. The viewer was initially published from commit `1df0ada`;
the team homepage was initially published from commit `0c2ffb5`. Git uses the existing GitHub
Credential Manager login for publishing; credentials are not stored in this
project.

Domainmonger manages the domain account and currently hosts its authoritative
DNS through `ns5.domainmonger.com` through `ns8.domainmonger.com`.
GitHub Pages hosts the website. The website repository's `CNAME` contains
`www.firstinfright.org`. The DNS configuration is:

| Host | Type | Value |
| --- | --- | --- |
| `@` | A | `185.199.108.153` |
| `@` | A | `185.199.109.153` |
| `@` | A | `185.199.110.153` |
| `@` | A | `185.199.111.153` |
| `www` | CNAME | `treygreer.github.io` |
| `_github-pages-challenge-treygreer` | TXT | Preserve the existing GitHub verification token |

Checked October 14, 2026: those DNS records resolve and the domain redirects
to `https://www.firstinfright.org/` with valid HTTPS. Before the publishing split,
this address served the older viewer with cache `treb-viewer-5cc2b3ecd637`.
The team homepage retires that root-scoped worker/cache so returning
visitors see the homepage. Install the viewer again from its new `/treb/` URL
to use its separate offline cache.
The updated viewer is verified over HTTPS at `https://treygreer.github.io/treb/`,
including the 8 lb pumpkin, 11 oz sling, 9 oz tip hardware, favicon, and offline reopening.
After both deployments completed, the live homepage's simulation link and
offline viewer reopening passed the end-to-end check on October 14, 2026.
The October 17 update to an 8 lb pumpkin and the feet-and-inches sling label
also passed the live end-to-end check, including the constant-gain time slider.

For a future viewer update, rebuild here, then copy the **contents** of
`web/site/` into the simulator's `gh-pages` publishing checkout and push that
branch. Team homepage changes belong in `firstinfright` on `main`; preserve its
`CNAME` and repository documentation.
Do not upload `web/src/template.html`; it is an unbuilt source template.
GitHub's repository **Settings → Pages** controls publishing, the custom domain
and HTTPS. Domainmonger's DNS editor controls the records above; DNS does not
need to change when updating the viewer.

The team website favicon has been renamed to `favicon.ico` and returns HTTP 200.
Its source is now `web/src/assets/favicon.ico`. The build copies it to
`web/site/favicon.ico`, adds an explicit favicon link with a content-based
version in the URL, and precaches it for offline use. The standalone download
embeds the icon to remain self-contained. Publish those viewer files together
on the simulator's `gh-pages` branch; the team homepage has its own explicit
favicon link.

Recheck the live viewer and offline reopening from `web/` with:

```powershell
node check-live.mjs
node check-split-live.mjs
```

This is a saved simulation viewer, not a browser physics solver: changing
machine parameters or optimizing a new throw still requires Python and a new
export. The solver uses its normal 0.1 ms step. Display samples are 2 ms apart,
with full 0.1 ms sampling around release and the exact interpolated release
included. Energy and force diagnostics are evaluated at those saved states.
The web time slider maps its position to elapsed seconds, with a constant time
gain across the entire run. It selects the nearest saved state; the denser
release samples no longer occupy a disproportionate portion of the slider.
Playback advances by wall-clock time at the selected speed.

Regenerate from the repository root (requires Node.js and npm as well as uv):

```powershell
uv run python scripts/export_web.py
uv run python scripts/web_icons.py
Set-Location web
npm ci --cache .npm-cache
npm run build
npm test
```

For a desktop preview, run `uv run python -m http.server 8000 --directory web/site`
from the repository root and open `http://localhost:8000`. Localhost supports
service workers for development. A phone accessing the computer's plain HTTP
LAN address can view it online but needs HTTPS for offline installation.
The browser check uses installed Windows Chrome and verifies a phone-sized
layout and touch rotation, control behavior, cached offline reopening and the
self-contained file with networking disabled. A physical Pixel has not been
tested.

## Verify

```powershell
uv run python treb_sim/src/test.py
uv run python scripts/smoke_test.py
uv run python scripts/smoke_test_3d.py
```

The smoke check runs both models, checks finite states and launch ranges, renders
their Qt viewers offscreen with forces enabled at several slider positions, and
checks drag-free flight against the analytic ballistic trajectory.
The separate 3D check requires a Windows desktop and OpenGL context; it verifies
the moving assemblies' joint alignment, native rendering, camera orbit, frame
visibility, the square tube and sling release. The offscreen smoke check uses
the retained 2D rendering.

## Counterweight inertia

`cw_moment` is the moving counterweight assembly's inertia about its own center
of mass, supplied in lb·in² and converted directly to kg·m². These input units
are inferred from the historical conversion and plausible assembly dimensions;
the original measurement record is not in the repository. The default
`3.516e6` gives approximately 1,028.92 kg·m². The hinge offset is represented
separately by the constraints and must not be added to this input.

Before the October 20, 2013 commit `f7136e2` ("updated masses and moments to
conform with 2011 measurements"), the model had a 4,000 lb counterweight and a
separate 400 lb supporting beam. That commit replaced them with a 4,581 lb
assembly, a massless beam, and the explicit inertia above. This suggests the
assembly properties may include the beam, but the original record is absent.
The 2013 model adds 2,000 lb of tractor weights while retaining the historical
center of mass and inertia; their detailed mass distribution is not documented.

The corrected conversion slightly changes the best attainable range and shifts
the best release time. A fixed release time can change the predicted range much
more substantially; release timing has not been retuned.
