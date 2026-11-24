import test from 'node:test';
import assert from 'node:assert/strict';
import { plateImpact } from '../examples/panel.js';
import { bind, skinPoint, MaterialPoint, materials } from '../index.js';

test('a fitted synthetic plate keeps its zero-load pose and absorbs energy during impact', () => {
  for (const speed of [0, 5]) {
    const { world, centre } = plateImpact(speed);
    let travel = 0;
    for (let i = 0; i < Math.round(0.15 / world.dt); i++) {
      world.step();
      travel = Math.max(travel, skinPoint(centre, world.q, 0));
    }
    const t = world.telemetry();
    if (!speed) {
      assert.equal(travel, 0);
      assert.equal(t.damaged, 0);
      continue;
    }
    assert.ok(t.peakForce > 0 && travel > 0.001 && t.plastic > 0);
    assert.ok(Math.abs(t.momentumError) < 1e-8 && Math.abs(t.energyError) < world.initialEnergy * 0.025);
  }
});
test('a display vertex follows a rigidly rotated and translated cage', () => {
  const { world, panel } = plateImpact();
  const x = [0.014, -0.032, 0.05],
    p = bind(panel, x);
  const q = Float64Array.from(
    world.model.nodes.flatMap(({ position: [a, b] }) => [-b - a + 0.2, a - b - 0.1, 0.05]),
  );
  assert.ok(Math.abs(x[0] + skinPoint(p, q, 0) - (0.2 - x[1])) < 1e-12);
  assert.ok(Math.abs(x[1] + skinPoint(p, q, 1) - (x[0] - 0.1)) < 1e-12);
});
test('every supplied material card starts in equilibrium', () => {
  for (const card of Object.values(materials)) assert.ok(Math.abs(new MaterialPoint(card).load(0)) < 1e-12);
});
