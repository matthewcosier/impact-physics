import { pathToFileURL } from 'node:url';
import { bind, buildLattice, createImpactor, ImpactWorld, materials, skinPoint } from '../index.js';

// A 300 × 240 mm moulding on four breakable mounting lands. +X is inward.
export function plateImpact(speed = 5) {
  const model = { nodes: [], beams: [], joints: [], contacts: [] };
  const cage = { origin: [0, -0.15, -0.12], cell: [0.02, 0.075, 0.06], counts: [1, 4, 4], occupied: [] };
  for (let y = 0; y < 4; y++) for (let z = 0; z < 4; z++) cage.occupied.push([0, y, z]);
  const material = materials.panel,
    volume = 0.3 * 0.24 * 0.002;
  const panel = buildLattice(model, {
    id: 'plate',
    cage,
    material,
    volumeOfMaterial: volume,
    mass: volume * material.density,
    sectionFactor: 0.12,
  });
  const centre = bind(panel, [0, 0, 0]);
  for (const y of [-0.15, 0.15])
    for (const z of [-0.12, 0.12]) {
      model.joints.push({
        a: bind(panel, [0.02, y, z]),
        k: 18000,
        onset: 0.004,
        failure: 0.022,
        damping: 0.25,
      });
    }
  for (let y = 0; y <= 8; y++)
    for (let z = 0; z <= 8; z++)
      model.contacts.push(bind(panel, [0, -0.15 + (0.3 * y) / 8, -0.12 + (0.24 * z) / 8]));
  const striker = createImpactor({
    shape: 'sphere',
    center: [-0.065, 0, 0],
    radius: 0.045,
    speed,
    mass: 0.6,
    stiffness: 6000,
    damping: 0.18,
    pad: 0.001,
  });
  const world = new ImpactWorld(model, striker, 1 / 128000);
  return { world, centre, panel };
}

if (process.argv[1] && pathToFileURL(process.argv[1]).href === import.meta.url) {
  const { world, centre } = plateImpact();
  let peakTravel = 0;
  for (let step = 0; step < Math.round(0.15 / world.dt); step++) {
    world.step();
    peakTravel = Math.max(peakTravel, skinPoint(centre, world.q, 0));
  }
  const t = world.telemetry();
  console.log(
    JSON.stringify(
      {
        seconds: t.time,
        peakNewtons: t.peakForce,
        peakTravelMm: 1000 * peakTravel,
        permanentBeamSetMm: 1000 * t.plasticSet,
        damagedMounts: t.damaged,
        failedMounts: t.failed,
        momentumErrorNs: t.momentumError,
        energyErrorJ: t.energyError,
      },
      null,
      2,
    ),
  );
}
