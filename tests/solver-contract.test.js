import assert from 'node:assert/strict';
import test from 'node:test';
import { ImpactWorld } from '../src/physics.js';
import { createImpactor, contactDepth } from '../src/impact-shapes.js';
import { FlexDamping } from '../src/flex-damping.js';

const near = (a, b, e = 1e-10) => assert.ok(Math.abs(a - b) < e, `${a} vs ${b}`);
test('internal velocity damping preserves total linear and angular momentum and accounts for its energy loss', () => {
  const positions = [
    [1, 2, 3],
    [2, 2, 3],
    [1, 2.6, 3],
    [1, 2, 4.5],
  ];
  const masses = [1, 4, 2, 3];
  const q = Float64Array.from(positions.flatMap((p, i) => [0.1 * i, 0, 0.2 * i]));
  const moved = positions.map((p, i) => p.map((x, d) => x + q[3 * i + d]));
  const rigid = moved.flatMap(([x, y, z]) => [
    3 + 0.2 * z + 0.5 * y,
    2 - 0.5 * x - 0.7 * z,
    -1 + 0.7 * y - 0.2 * x,
  ]);
  const damper = new FlexDamping(
    { nodes: [0, 1, 2, 3], rate: 80, maxAcceleration: 18 },
    positions.map((position) => ({ position })),
    masses,
    0.004,
  );
  const v = Float64Array.from(rigid);
  near(damper.apply(q, v, masses), 0);
  rigid.forEach((value, i) => near(v[i], value));
  const invariants = () => {
    const p = [0, 0, 0],
      l = [0, 0, 0];
    let energy = 0;
    moved.forEach((r, i) => {
      for (let d = 0; d < 3; d++) {
        const a = (d + 1) % 3,
          b = (d + 2) % 3;
        p[d] += masses[i] * v[3 * i + d];
        l[d] += masses[i] * (r[a] * v[3 * i + b] - r[b] * v[3 * i + a]);
        energy += 0.5 * masses[i] * v[3 * i + d] ** 2;
      }
    });
    return { momentum: [...p, ...l], energy };
  };
  v[4] += 6;
  const before = invariants(),
    loss = damper.apply(q, v, masses),
    after = invariants();
  assert.ok(loss > 0);
  near(before.energy - after.energy, loss);
  before.momentum.forEach((x, i) => near(after.momentum[i], x));
});

test('a hard wall catches a point which would cross the plane within a single step', () => {
  const position = [0.01, 0, 0],
    mass = 2;
  const model = { nodes: [{ mass, position }], wallPoints: [{ nodes: [0], weights: [1], position }] };
  const wall = createImpactor({
    shape: 'wall',
    center: [0, 0, 0],
    speed: 8,
    mass: 3,
    stiffness: 100,
    damping: 0,
  });
  const w = new ImpactWorld(model, wall, 0.004);
  for (let i = 0; i < 3; i++) {
    w.step();
    assert.ok(position[0] + w.q[0] >= -1e-10);
  }
  near(w.v[0], 0);
  near(w.contactImpulse, 16);
  near(w.telemetry().momentumError, 0);
  near(w.telemetry().energyError, 0);
});
test('an oblique flat face reacts along its rotated normal and the guide carries transverse impulse', () => {
  const angle = 25,
    cs = Math.cos((angle * Math.PI) / 180),
    sn = Math.sin((angle * Math.PI) / 180);
  const actor = createImpactor({
    shape: 'block',
    center: [0, 0, 0],
    angle,
    speed: 4,
    mass: 1,
    depth: 0.2,
    width: 0.8,
    height: 0.6,
    stiffness: 3000,
    damping: 0.2,
  });
  const n = [];
  near(contactDepth(actor, 0.099 * cs, 0.099 * sn, 0, n), 0.001);
  near(n[0], cs);
  near(n[1], sn);
  const position = [0.13 * cs, 0.13 * sn, 0];
  const w = new ImpactWorld(
    { nodes: [{ mass: 1, position }], contacts: [{ nodes: [0], weights: [1], position }] },
    actor,
  );
  for (let i = 0; i < 900; i++) w.step();
  assert.ok(w.contactImpulse > 1 && w.v[1] > 0);
  near(w.impactor.y / w.impactor.x, sn / cs);
  near(w.telemetry().momentumError, 0);
  assert.ok(Math.abs(w.telemetry().energyError) < 0.01 * w.initialEnergy);
});
test('a hidden render vertex loads each finite impactor even with no outer contact witnesses', async () => {
  const { buildLattice, bind } = await import('../src/lattice.js');
  const { surfaceCells } = await import('../src/wall-contact.js');
  const { materials } = await import('../src/materials.js');
  const { skinPoint } = await import('../src/physics.js');
  for (const shape of ['sphere', 'block', 'cylinder']) {
    const model = { nodes: [], beams: [] };
    const xyz = new Float32Array([0.002, 0.03, 0.02]);
    const cage = { origin: [0, 0, 0], cell: [0.05, 0.06, 0.04], counts: [1, 1, 1], occupied: [[0, 0, 0]] };
    const part = buildLattice(model, {
      id: 'backing',
      material: materials.panel,
      mass: 0.5,
      volumeOfMaterial: 1e-6,
      cage,
    });
    part.objects = [{ positions: 0, vertices: 1 }];
    model.surfaceCells = surfaceCells([part], xyz.buffer);
    const a = bind(part, Array.from(xyz));
    const actor = createImpactor({
      shape,
      center: [-0.26, 0.03, 0.02],
      mass: 1,
      speed: 5,
      radius: 0.05,
      height: 0.1,
      depth: 0.1,
      width: 0.1,
      stiffness: 5000,
      damping: 0.2,
    });
    const w = new ImpactWorld(model, actor);
    let peak = 0;
    for (let i = 0; i < 1250; i++) {
      w.step();
      const v = a.position.map((x, d) => x + skinPoint(a, w.q, d));
      peak = Math.max(
        peak,
        contactDepth(w.impactor, v[0] - w.impactor.x, v[1] - w.impactor.y, v[2] - w.impactor.z, []),
      );
    }
    assert.ok(w.contactImpulse > 0 && w.impactor.speed < 5);
    assert.ok(peak < 1e-6);
    near(w.telemetry().momentumError, 0, 1e-8);
  }
});

test('the one millisecond load reading integrates force over the same physical time at both timesteps', () => {
  for (const dt of [1 / 16000, 1 / 32000, 0.0003]) {
    const actor = createImpactor({
      shape: 'sphere',
      center: [0, 0, 0],
      radius: 0.05,
      mass: 0.9,
      speed: 4,
      stiffness: 2000,
      damping: 0.2,
    });
    const position = [0.08, 0, 0],
      point = { position, nodes: [0], weights: [1] };
    const w = new ImpactWorld({ nodes: [{ mass: 2, position }], contacts: [point] }, actor, dt);
    const history = [];
    let peak = 0;
    for (let i = 0; i < Math.ceil(0.02 / dt); i++) {
      w.step();
      history.push(w.contactForce);
      let timeLeft = 0.001,
        impulse = 0;
      for (let j = history.length - 1; j >= 0 && timeLeft > 1e-15; j--) {
        const span = Math.min(timeLeft, dt);
        impulse += history[j] * span;
        timeLeft -= span;
      }
      peak = Math.max(peak, Math.abs(impulse / 0.001));
      near(w.meanForce, impulse / 0.001, 1e-8);
    }
    assert.ok(peak > 0);
    near(w.peakMeanForce, peak, 1e-8);
  }
});

test('small-scale models can select a physical contact margin and use a transverse guide', async () => {
  const { buildLattice } = await import('../src/lattice.js');
  const { surfaceCells } = await import('../src/wall-contact.js');
  const { materials } = await import('../src/materials.js');
  const model = { nodes: [], beams: [] };
  const cage = { origin: [0, 0, 0], cell: [0.002, 0.002, 0.002], counts: [1, 1, 1], occupied: [[0, 0, 0]] };
  const part = buildLattice(model, {
    id: 'coupon',
    material: materials.panel,
    mass: 0.1,
    volumeOfMaterial: 1e-10,
    cage,
  });
  part.objects = [{ positions: 0, vertices: 1 }];
  model.surfaceCells = surfaceCells([part], new Float32Array([0.001, 0, 0.001]).buffer);
  const spec = {
    shape: 'sphere',
    center: [0.001, -0.0013, 0.001],
    radius: 0.001,
    angle: 90,
    mass: 1,
    speed: 0.01,
    stiffness: 100,
    damping: 0.2,
    surfaceMargin: 0,
    approachTime: 0,
  };
  const w = new ImpactWorld(model, createImpactor(spec));
  for (let i = 0; i < 100; i++) w.step();
  assert.equal(w.contactImpulse, 0);
  for (let i = 0; i < 500; i++) w.step();
  assert.ok(w.contactImpulse > 0 && w.impactor.y > spec.center[1]);
  near(w.impactor.x, spec.center[0]);
  assert.throws(() => createImpactor({ ...spec, surfaceMargin: -0.001 }), RangeError);
});
