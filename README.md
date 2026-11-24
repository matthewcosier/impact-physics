# Impact Physics

Impact Physics is a fixed-step, deformable-body solver written in JavaScript. Model a panel as a lattice of masses and yielding beams, attach it to neighbouring parts with breakable joints, then use the nodal displacements to deform your own render mesh. The solver runs in Node or a browser without a physics or rendering dependency.

## Try the synthetic impact

Use Node 20.19 or newer. The example drives a finite-mass sphere into a moulded panel held by four breakable mounts. It prints peak contact force, displacement, residual beam extension and the energy balance.

```sh
git clone https://github.com/matthewcosier/impact-physics.git
cd impact-physics
npm test
npm run example
```

To use the public entry point from another project, install from GitHub. Pin a commit in reproducible builds.

```sh
npm install github:matthewcosier/impact-physics
```

## Fit a model

1. Choose a reference frame in metres. Use kilograms, seconds, newtons and pascals throughout. Each node has a reference position, a mass and an optional initial velocity.
2. Create each part's mechanical domain. Supply reference nodes and beam properties directly, or call `buildLattice(model, {id, cage, mass, volumeOfMaterial, material, sectionFactor})`. A cage supplies `origin`, three cell dimensions, three cell counts and occupied integer cell coordinates. Include the cells traversed by your surface. Mass includes fittings. Material volume is the amount of solid material in the part.
3. Register contacts and mounts. A point returned by `bind(part, [x,y,z])` interpolates motion from its cell and distributes a force back to the same nodes. Contacts describe exposed surfaces. A joint supplies coincident sites `a` and `b`, stiffness `k` in N/m, damage onset and full separation in metres, and a damping ratio. Omitting `b` attaches the site to the fixture. For an assembly initially supported in static equilibrium, `model.createReleases = () => new ReleasedParts(model)` enables release by connectivity. Supply `model.parts` with owner IDs and node IDs. Intact joints to the support mark a component as retained. Detached islands acquire gravity and a ground response. Alternative support, gravity and ground conditions can be supplied by the application.
4. Advance the mechanical world at a fixed timestep and read `world.q`. Its three values per node are displacements from the reference pose. `world.v` is velocity. `world.telemetry()` reports force, impulse, damage, plastic work and numerical balance. `peakForce` is an instantaneous numerical peak. `peakMeanForce` is the largest absolute force averaged over a moving 1 ms window; it allows comparison across timesteps when a hard contact delivers a short impulse. An optional array passed to `world.step([{point, force: [fx, fy, fz]}])` applies external loads for one step.

The example contains a complete model definition and run. [Synthetic panel example](examples/panel.js)

## Draw a mesh

Bind each reference vertex once. Keep the original positions when updating your geometry. This example assumes the mechanical and rendering coordinates use the same basis and units.

```js
import { bind, skinPoint } from 'impact-physics';

const attribute = geometry.getAttribute('position'); // Three.js BufferGeometry
const rest = Array.from({ length: attribute.count }, (_, i) => [
  attribute.getX(i), attribute.getY(i), attribute.getZ(i)
]);
const sites = rest.map(p => bind(panel, p));

function display(world) {
  for (let i = 0; i < sites.length; i++) {
    const p = sites[i];
    attribute.setXYZ(i,
      rest[i][0] + skinPoint(p, world.q, 0),
      rest[i][1] + skinPoint(p, world.q, 1),
      rest[i][2] + skinPoint(p, world.q, 2));
  }
  attribute.needsUpdate = true;
  geometry.computeVertexNormals();
}
```

Use a simulation clock or worker to advance the solver independently of display refresh. Increasing the timestep changes numerical stability. Refine it and compare impact impulse, maximum force, deformation and energy balance for each new model.

## Contact and deformation

- **Guided impact.** A sphere, vertical cylinder or rectangular block travels along a guide in the XY plane. A block turns with its guide. Set `angle` in degrees, `center` in metres, `speed` in m/s, and the dimensions, stiffness, damping and mass in `createImpactor`. Register representative exposed material points. The registered points carry compliant contact. Supplying `model.surfaceCells` also engages normal impulse constraints on the fitted vertices and backing points. The guided constraint checks the current pose and a velocity prediction. Set `surfaceMargin` (metres, default 0.002) and `approachTime` (seconds, default 0.003) to suit your mesh spacing and impact speed. A zero value disables the corresponding margin or prediction. Independent triangle interiors, very fast crossings between steps and arbitrary mesh self-collision require additional coverage and tests.
- **Fixed wall.** The plane is at `center[0]` with the vehicle side at greater X. The model and its backing sled approach along −X. The wall applies an impulse to prevent crossing at registered points. Supply `model.boundaryPoints` for custom witnesses. For a rendered lattice mesh, give each part `objects: [{positions: byteOffset, vertices: count}]` into an ArrayBuffer of XYZ float32 positions and set `model.surfaceCells = surfaceCells(parts, buffer)`. The wall checks the deformed render vertices in candidate cells. For straight-edged triangles, keeping their vertices on the permitted side also keeps their interiors there. With no buffer, `surfaceCells(parts)` encloses complete occupied cells, which may reach the plane earlier than a fine skin.
- **Constitutive models.** The catalogue contains 30 starting material cards. Axial links, hyperelastic or elastoplastic tetrahedra, crushable foam and direction-dependent fabric coupons support different loading histories. `MaterialPoint` advances a uniaxial specimen independently. Parameter values require calibration to the specimen and loading rate.
- **Structural damage.** Beam rest lengths change after yield. Cohesive joints soften irreversibly, transfer less load on unloading and stop transmitting force after rupture. Optional deformation damping removes velocity relative to the group's instantaneous rigid motion. Ordinary lattice beams retain connectivity. Crack propagation through a newly torn render mesh, general triangle self-contact and validated crashworthiness prediction are outside this release.

## Optional WebAssembly acceleration

The JavaScript beam loop is the reference implementation. A Rust toolchain with the WebAssembly target can compile the same axial-beam equations locally.

```sh
rustup target add wasm32-unknown-unknown
rustc --target wasm32-unknown-unknown --crate-type cdylib -C opt-level=3 -C lto=yes -C panic=abort -C strip=symbols native/beams.rs -o beams.wasm
```

Pass the compiled module to `useBeamKernel(world, module)` before stepping that world. Material, joint and contact calculations continue in JavaScript.

## Licence

MIT. [Licence text](LICENSE)
