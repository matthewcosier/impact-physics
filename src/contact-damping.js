// Relative-velocity contact forces are solved against a common predicted
// velocity field. Several surface samples can share the same light node.
class ContactStencil {
  constructor(nodes, weights, masses) {
    this.nodes = nodes;
    this.weights = weights;
    this.inverse = weights.map((w, j) => w / masses[nodes[j]]);
    this.mobility = weights.reduce((sum, w, j) => sum + w * this.inverse[j], 0);
  }
  read(field) {
    let x = 0,
      y = 0,
      z = 0;
    for (let j = 0; j < this.nodes.length; j++) {
      const i = 3 * this.nodes[j],
        w = this.weights[j];
      x += w * field[i];
      y += w * field[i + 1];
      z += w * field[i + 2];
    }
    return [x, y, z];
  }
  add(field, weights, x, y, z) {
    for (let j = 0; j < this.nodes.length; j++) {
      const i = 3 * this.nodes[j],
        w = weights[j];
      field[i] += w * x;
      field[i + 1] += w * y;
      field[i + 2] += w * z;
    }
  }
}

export function contactStencil(p, body, carrierWeights, mass) {
  return new ContactStencil(
    [...p.nodes, ...body.nodes],
    [...p.weights, ...carrierWeights.map((w) => -w)],
    mass,
  );
}

export function dampContacts(world) {
  const active = world.contactBodies.flatMap((b) => b.contacts);
  if (!active.length) return;
  const { dt, v, force, mass } = world;
  const predicted = (world.contactVelocity ??= new Float64Array(v.length));
  for (let n = 0; n < mass.length; n++)
    for (let d = 0; d < 3; d++) {
      const i = 3 * n + d;
      predicted[i] = v[i] + (dt * force[i]) / mass[n];
    }
  for (let iteration = 0; iteration < 160; iteration++) {
    let maxChange = 0;
    for (const c of active) {
      const [nx, ny, nz] = c.normal,
        N = c.normalStencil;
      const [vx, vy, vz] = N.read(predicted),
        vn = nx * vx + ny * vy + nz * vz;
      const next = Math.max(
        0,
        c.normalForce - (c.normalForce + c.damping * vn) / (1 + c.damping * dt * N.mobility),
      );
      const dn = next - c.normalForce;
      N.add(predicted, N.inverse, dt * dn * nx, dt * dn * ny, dt * dn * nz);
      c.normalForce = next;
      maxChange = Math.max(maxChange, Math.abs(dn) / (1 + next));
      if (!c.frictionLimit) continue;
      const T = c.tangentStencil,
        coef = c.frictionLimit / 0.2;
      const [ux, uy, uz] = T.read(predicted),
        un = ux * nx + uy * ny + uz * nz;
      const scale = 1 / (1 + coef * dt * T.mobility);
      let x = c.fx - (c.fx + coef * (ux - un * nx)) * scale;
      let y = c.fy - (c.fy + coef * (uy - un * ny)) * scale;
      let z = c.fz - (c.fz + coef * (uz - un * nz)) * scale;
      const length = Math.hypot(x, y, z);
      if (length > c.frictionLimit) {
        const cap = c.frictionLimit / length;
        x *= cap;
        y *= cap;
        z *= cap;
      }
      const dx = x - c.fx,
        dy = y - c.fy,
        dz = z - c.fz;
      T.add(predicted, T.inverse, dt * dx, dt * dy, dt * dz);
      maxChange = Math.max(maxChange, Math.hypot(dx, dy, dz) / (1 + Math.min(length, c.frictionLimit)));
      c.fx = x;
      c.fy = y;
      c.fz = z;
    }
    if (maxChange < 0.0005) break;
    if (iteration === 159) throw new Error('Contact damping failed to converge; reduce the timestep');
  }
  for (const c of active) {
    const [nx, ny, nz] = c.normal,
      N = c.normalStencil;
    const f = c.normalForce,
      [vx, vy, vz] = N.read(v);
    N.add(force, N.weights, f * nx, f * ny, f * nz);
    c.startNormalSpeed = vx * nx + vy * ny + vz * nz;
    c.body.impulses[c.owner] += dt * f;
    if (c.frictionLimit) {
      const T = c.tangentStencil;
      c.startTangential = T.read(v);
      T.add(force, T.weights, c.fx, c.fy, c.fz);
    }
  }
}

// Record work against the two actual endpoint velocities. The provisional
// unconstrained velocity used to solve a dashpot can be changed by a wall,
// guide or a rigid fitting in the very same step.
export function contactWork(world) {
  for (const body of world.contactBodies)
    for (const c of body.contacts) {
      const n = c.normal,
        v = c.normalStencil.read(world.v);
      world.loss.damping -=
        0.5 * world.dt * c.normalForce * (c.startNormalSpeed + n[0] * v[0] + n[1] * v[1] + n[2] * v[2]);
      if (c.frictionLimit) {
        const t = c.tangentStencil.read(world.v),
          t0 = c.startTangential;
        world.loss.friction -=
          0.5 * world.dt * (c.fx * (t[0] + t0[0]) + c.fy * (t[1] + t0[1]) + c.fz * (t[2] + t0[2]));
      }
    }
}
