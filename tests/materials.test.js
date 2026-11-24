import { MetalPoint } from '../src/metal.js';
import test from 'node:test';
import assert from 'node:assert/strict';
import { materials } from '../src/materials.js';
import { MaterialPoint } from '../src/material-point.js';
import { ImpactWorld } from '../src/physics.js';
import { rubberResponse, SolidTetrahedron } from '../src/solid.js';

const close = (a, b, tol = 1e-9) =>
  assert.ok(Math.abs(a - b) <= tol, `${a} differs from ${b} by more than ${tol}`);
const pos = (position, mass = 0) => ({ position, mass });
const load = (n, force) => ({ point: { nodes: [n], weights: [1] }, force });
test('every sheet can run an isolated material point from rest', () => {
  for (const m of Object.values(materials)) {
    const p = new MaterialPoint(m);
    close(p.load(0), 0);
    assert.ok(m.density > 0 && m.typicalThickness > 0);
    assert.ok(Math.abs(p.load(m.law === 'foam' ? -0.01 : 0.0001)) > 0);
  }
});
test('foams follow the pressure data and dissipate work over a full quasistatic cycle', () => {
  for (const id of ['seatFoam', 'bolsterFoam', 'dashFoam', 'absorber']) {
    const m = materials[id],
      p = new MaterialPoint(m);
    for (const [e, stress] of m.compressionCurve) close(-p.load(-e), stress, 1e-8);
    const peak = m.compressionCurve.at(-1)[0];
    const permanent = m.crushFraction * (peak - m.crushOnset);
    close(p.release(m.recoverySeconds), -permanent + (-peak + permanent) / Math.E);
    close(p.release(100), -permanent);
    const prior = p.strain;
    assert.throws(() => p.load(-0.99), RangeError);
    assert.throws(() => p.load(0.01), RangeError);
    assert.equal(p.strain, prior);
    const cycle = new MaterialPoint(m);
    let work = 0;
    for (let i = 1; i <= 200; i++) {
      const e = -peak * (1 - Math.abs(i - 100) / 100);
      const old = cycle.strain,
        stress = cycle.stress;
      work += 0.5 * (stress + cycle.load(e)) * (e - old);
    }
    assert.ok(work > 0);
  }
});
test('seat, bolster and dash differ in pressure and residual thickness', () => {
  const pads = ['seatFoam', 'bolsterFoam', 'dashFoam'].map((k) => new MaterialPoint(materials[k]));
  for (const [i, p] of pads.entries()) close(-p.load(-0.5), [8500, 17000, 53000][i]);
  const strains = pads.map((p) => p.release(1));
  close(strains[0], -0.5 * Math.exp(-1 / 0.65));
  close(strains[2], -0.0348 - 0.4652 * Math.exp(-1 / 0.26));
  assert.ok(strains[0] < strains[1]);
});
test('leather has a compliant toe, directional tension, slack compression and irreversible tear', () => {
  const along = new MaterialPoint(materials.leather),
    across = new MaterialPoint(materials.leather, true);
  assert.ok(along.load(0.01) / 0.01 < along.load(0.12) / 0.12);
  close(across.load(0.12) / along.stress, materials.leather.transverseRatio);
  assert.ok(-along.load(-0.1) < along.load(0.1) / 50);
  assert.equal(along.load(0.5), 0);
  assert.ok(along.failed);
  assert.equal(along.load(0.03), 0);
});
test('metal unloads at Young’s modulus and retains set, glass ruptures elastically, rubber follows uniaxial nominal stress', () => {
  const m = materials.mildSteel,
    p = new MaterialPoint(m);
  close(p.load(0.0001), m.E * 0.0001, 1e-6);
  const sigma = p.load(0.04);
  const set = 0.04 - sigma / m.E;
  close(p.load(0.04 - 0.0001), sigma - m.E * 0.0001, 1e-4);
  close(p.release(1), set);
  const brittle = new MaterialPoint(materials.glass);
  close(brittle.load(0.0001), materials.glass.E * 0.0001);
  close(brittle.load(0.002), 0);
  assert.ok(brittle.failed);
  const rubber = new MaterialPoint(materials.seal);
  close(rubber.load(0.2), materials.seal.G * (1.2 - 1 / 1.2 ** 2));
});
test('generic material links carry equal and opposite axial traction, with reference mass independent of pose', () => {
  for (const id of ['seatFoam', 'dashFoam', 'leather']) {
    const m = materials[id],
      A = 0.01,
      L = 0.1;
    const world = new ImpactWorld({
      nodes: [pos([0, 0, 0]), pos([0, L, 0])],
      materialLinks: [{ a: 0, b: 1, area: A, material: m }],
    });
    close(world.mass[0] + world.mass[1], m.density * A * L);
    const e = id === 'leather' ? 0.1 : -0.5;
    world.q.set([0.2, 0.1, 0.1, 0.2 + L * (1 + e), 0.1 - L, 0.1]);
    const link = world.materialLinks[0];
    link.sample(world.q);
    link.scatter(world.force, link.axialForce);
    const stress = new MaterialPoint(m).load(e);
    close(world.force[0], stress * A, 1e-7);
    close(world.force[0] + world.force[3], 0);
    close(world.force[1], 0);
  }
});
test('free material under opposed loads has finite compression, balanced momentum and step-converged recovery', () => {
  function run(dt, m) {
    const world = new ImpactWorld(
      {
        nodes: [pos([0, 0, 0], 0.1), pos([0.1, 0, 0], 0.1)],
        materialLinks: [{ a: 0, b: 1, area: 0.01, material: m }],
      },
      null,
      dt,
    );
    const press = [load(0, [10, 0, 0]), load(1, [-10, 0, 0])];
    for (let i = 0; i < Math.round(0.2 / dt); i++) world.step(press);
    const compressed = world.materialLinks[0].strain;
    for (let i = 0; i < Math.round(0.1 / dt); i++) world.step();
    const t = world.telemetry();
    close(t.momentumError, 0);
    assert.ok(t.damping > 0);
    assert.ok(Math.abs(t.energyError) < 0.01 * Math.abs(t.externalWork));
    return [compressed, world.materialLinks[0].strain];
  }
  const a = run(1 / 8000, materials.seatFoam),
    b = run(1 / 16000, materials.seatFoam),
    c = run(1 / 8000, materials.dashFoam);
  assert.ok(a[0] < a[1] && a[1] < 0 && a[0] < c[0]);
  close(a[0], b[0], 0.0001);
  close(a[1], b[1], 0.0001);
});
test('a pair of point forces delivers a pure couple with the correct angular impulse and work', () => {
  const w = new ImpactWorld({ nodes: [pos([0, -0.5, 0], 2), pos([0, 0.5, 0], 2)] }, null, 0.001);
  const loads = [load(0, [4, 0, 0]), load(1, [-4, 0, 0])];
  for (let i = 0; i < 100; i++) w.step(loads);
  close(w.telemetry().momentumError, 0);
  close(w.v[0] - w.v[3], 0.4);
  close(w.externalWork, 0.08);
  close(w.telemetry().energyError, 0);
});
test('3D rubber has objective strain energy and the analytical shear stress, pressure and stress-work gradient', () => {
  const m = materials.seal;
  const rotated = [0, -1, 0, 1, 0, 0, 0, 0, 1];
  const rest = rubberResponse(rotated, m);
  close(rest.energyDensity, 0);
  for (const x of rest.P) close(x, 0);
  const gamma = 0.3,
    F = [1, gamma, 0, 0, 1, 0, 0, 0, 1],
    r = rubberResponse(F, m);
  close(r.energyDensity, 0.5 * m.G * gamma ** 2, 1e-7);
  close(r.P[1], m.G * gamma);
  const volumetric = rubberResponse([1.1, 0, 0, 0, 1.1, 0, 0, 0, 1.1], m);
  close(volumetric.P[0], (m.bulkModulus * Math.log(1.1 ** 3)) / 1.1, 1e-6);
  const combined = [0.97, 0.15, 0.03, -0.09, 1.08, 0.12, 0.1, -0.03, 1.04];
  const general = rubberResponse(combined, m);
  for (let i = 0; i < 9; i++) {
    const lo = combined.slice(),
      hi = combined.slice(),
      h = 1e-6;
    hi[i] += h;
    lo[i] -= h;
    close(
      general.P[i],
      (rubberResponse(hi, m).energyDensity - rubberResponse(lo, m).energyDensity) / (2 * h),
      0.02,
    );
  }
  assert.throws(() => rubberResponse([-1, 0, 0, 0, 1, 0, 0, 0, 1], m), RangeError);
});
test('a rubber tetrahedron transmits shear, conserves force and moment, and its forces are the negative energy gradient', () => {
  const nodes = [
    [0, 0, 0],
    [0.1, 0.01, 0],
    [0, 0.1, 0.02],
    [0.01, 0, 0.1],
  ].map((p) => pos(p));
  const spec = { nodes: [0, 1, 2, 3], material: materials.seal };
  const solid = new SolidTetrahedron(spec, nodes),
    w = new ImpactWorld({ nodes, solids: [spec] }, null, 1e-6);
  close(
    w.mass.reduce((a, b) => a + b),
    solid.volume * materials.seal.density,
  );
  for (let i = 0; i < 4; i++) w.q[3 * i] = 0.2 * nodes[i].position[1];
  const force = new Float64Array(12),
    q = w.q;
  solid.evaluate(q, force);
  for (let i = 0; i < q.length; i++) {
    const h = 1e-7;
    q[i] += h;
    const hi = solid.evaluate(q);
    q[i] -= 2 * h;
    const lo = solid.evaluate(q);
    q[i] += h;
    close(force[i], -(hi - lo) / (2 * h), 1e-4);
  }
  const torque = [0, 0, 0];
  for (let d = 0; d < 3; d++) {
    close(force[d] + force[d + 3] + force[d + 6] + force[d + 9], 0, 1e-9);
    for (let n = 0; n < 4; n++) {
      const a = (d + 1) % 3,
        b = (d + 2) % 3;
      torque[d] +=
        (nodes[n].position[a] + q[3 * n + a]) * force[3 * n + b] -
        (nodes[n].position[b] + q[3 * n + b]) * force[3 * n + a];
    }
    close(torque[d], 0, 1e-8);
  }
  const energy = solid.evaluate(q);
  w.step();
  assert.ok(w.telemetry().kinetic > 0);
  assert.ok(Math.abs(w.telemetry().elastic + w.telemetry().kinetic - energy) / energy < 0.0001);
});
test('pads in series settle to the same traction with different compressions regardless of element order', () => {
  function run(reverse) {
    const pads = [
      { a: 0, b: 1, material: materials.seatFoam, area: 0.01 },
      { a: 1, b: 2, material: materials.dashFoam, area: 0.01 },
    ];
    const w = new ImpactWorld(
      {
        nodes: [pos([0, 0, 0], 0.4), pos([0.1, 0, 0], 0.4), pos([0.2, 0, 0], 0.4)],
        materialLinks: reverse ? pads.reverse() : pads,
      },
      null,
      1 / 2000,
    );
    const loads = [load(0, [10, 0, 0]), load(2, [-10, 0, 0])];
    for (let i = 0; i < 12000; i++) w.step(loads);
    const [seat, dash] = w.materialLinks.toSorted((a, b) => a.a - b.a);
    close(seat.axialForce, -10, 0.01);
    close(dash.axialForce, -10, 0.01);
    assert.ok(seat.strain < dash.strain);
    return seat.strain;
  }
  close(run(false), run(true), 1e-8);
});

test('metal under hydrostatic pressure stays elastic and shear returns to the expanding von Mises surface', () => {
  const m = materials.mildSteel,
    p = new MetalPoint(m);
  const hydro = p.response([0.001, 0, 0, 0, 0.001, 0, 0, 0, 0.001], true);
  close(hydro.accumulated, 0);
  close(hydro.stress[4], 3 * p.K * 0.001, 1e-6);
  const strain = [0, 0.004, 0, 0.004, 0, 0, 0, 0, 0];
  const s = p.response(strain, true);
  assert.ok(s.accumulated > 0);
  close(Math.sqrt(3) * s.stress[1], m.yieldStress + p.H * s.accumulated, 1e-6);
  close(s.plasticWorkDensity, m.yieldStress * s.accumulated);
  close(p.plastic[0] + p.plastic[4] + p.plastic[8], 0);
  strain[1] -= 0.0001;
  strain[3] = strain[1];
  const unload = p.response(strain, true);
  close(unload.stress[1], s.stress[1] - 2 * m.G * 0.0001, 1e-6);
  close(unload.accumulated, s.accumulated);
});
test('a metal solid under simultaneous stretch and shear retains the same energy and rotates its nodal forces', () => {
  const nodes = [
    [0, 0, 0],
    [0.1, 0, 0],
    [0, 0.1, 0],
    [0, 0, 0.1],
  ].map((p) => pos(p));
  const spec = { nodes: [0, 1, 2, 3], material: materials.mildSteel };
  const w = new ImpactWorld({ nodes, solids: [spec] }, null, 1e-7);
  const solid = w.solids[0],
    q = w.q;
  q[3] = 0.0001;
  q[6] = 0.0004;
  const f = new Float64Array(12),
    u = solid.evaluate(q, f, true);
  assert.ok(solid.plasticWork > 0);
  const qr = new Float64Array(12);
  for (let n = 0; n < 4; n++) {
    const [x, y] = nodes[n].position;
    qr[3 * n] = -(y + q[3 * n + 1]) - x + 0.3;
    qr[3 * n + 1] = x + q[3 * n] - y - 0.2;
    qr[3 * n + 2] = q[3 * n + 2];
  }
  const g = new Float64Array(12);
  close(solid.evaluate(qr, g), u, 1e-8);
  for (let i = 0; i < 12; i += 3) {
    close(g[i], -f[i + 1], 1e-6);
    close(g[i + 1], f[i], 1e-6);
    close(g[i + 2], f[i + 2], 1e-6);
  }
});

test('material mass can supply the entire inertia and rejected forces leave the run untouched', () => {
  const nodes = [{ position: [0, 0, 0] }, { position: [0.1, 0, 0] }];
  const w = new ImpactWorld({
    nodes,
    materialLinks: [{ a: 0, b: 1, area: 0.01, material: materials.seatFoam }],
  });
  close(w.mass[0] + w.mass[1], 0.048);
  assert.throws(
    () => w.step([load(0, [4, 0, 0]), { point: { nodes: [5], weights: [1] }, force: [3, 0, 0] }]),
    RangeError,
  );
  close(w.externalImpulseX, 0);
  close(w.time, 0);
  close(w.q[0], 0);
});
