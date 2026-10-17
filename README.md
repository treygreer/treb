# First in Fright simulation viewer

Live viewer: https://treygreer.github.io/treb/

Team website: https://www.firstinfright.org/

This `gh-pages` branch contains the generated offline web viewer. The simulator
source is in the repository's main development branch. GitHub Pages should
publish `gh-pages` from `/`, with no custom domain on this repository.

Current configuration: 2013 machine with the 2012 counterweight, 8 lb pumpkin,
11 oz sling and 9 oz added hardware at the throwing-arm tip. The optimized
rounded sling length is 9.54 ft, at a 30-degree launch elevation, predicting
approximately 2,177 ft.

The viewer includes 2D/3D switching, touch camera controls, force and energy
displays and playback. Install from Chrome once it reports **Offline ready**.
The slider position is linear in elapsed seconds, even around the densely
sampled release event. Playback uses elapsed wall-clock time at the chosen speed.
`treb-2026-offline.html` is a self-contained download.

To regenerate, run `scripts/export_web.py` in the simulator checkout, then
`npm run build` in `web/`. Publish the contents of `web/site/` to this branch.
All paths are relative so the viewer's cache and manifest work under `/treb/`.
The team homepage is maintained independently in `treygreer/firstinfright`.
