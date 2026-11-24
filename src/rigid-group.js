// The nodes supply the mass distribution and all contact/attachment forces.
// An intact rigid fitting can translate and rotate, but has no plastic
// extension. Equal and opposite loads still reach its supporting joint.
export function product(m, x, y, z) {
  return [m[0] * x + m[1] * y + m[2] * z, m[3] * x + m[4] * y + m[5] * z, m[6] * x + m[7] * y + m[8] * z];
}
export function inverseInertia(xx, yy, zz, xy, xz, yz) {
  const a = yy * zz - yz * yz,
    b = xz * yz - xy * zz,
    c = xy * yz - xz * yy;
  const d = xx * zz - xz * xz,
    e = xy * xz - xx * yz,
    f = xx * yy - xy * xy;
  const det = xx * a + xy * b + xz * c;
  if (!(det > 0)) throw new RangeError('Singular inertia tensor');
  return [a, b, c, b, d, e, c, e, f].map((cofactor) => cofactor / det);
}
function rotation(omega, dt) {
  const [x, y, z] = omega,
    speed = Math.hypot(x, y, z),
    angle = speed * dt;
  // sin(a)/w and (1-cos(a))/w², including their limits at rest.
  const h = speed ? Math.sin(angle) / speed : dt;
  const h2 = speed ? 2 * (Math.sin(angle / 2) / speed) ** 2 : 0.5 * dt * dt;
  return [
    1 - h2 * (y * y + z * z),
    h2 * x * y - h * z,
    h2 * x * z + h * y,
    h2 * x * y + h * z,
    1 - h2 * (x * x + z * z),
    h2 * y * z - h * x,
    h2 * x * z - h * y,
    h2 * y * z + h * x,
    1 - h2 * (x * x + y * y),
  ];
}
export class RigidGroup {
  constructor(indices, model, masses) {
    if (indices.length < 3) throw new RangeError('A rigid group needs three non-collinear nodes');
    this.indices = indices;
    this.rest = indices.map((n) => model.nodes[n].position);
    this.mass = indices.reduce((sum, n) => sum + masses[n], 0);
    this.arms = new Float64Array(3 * indices.length);
  }
  advance({ q, v, force, mass, dt }) {
    const ns = this.indices,
      arms = this.arms;
    const centre = [0, 0, 0],
      linear = [0, 0, 0];
    for (const [j, n] of ns.entries())
      for (let d = 0; d < 3; d++) {
        const i = 3 * n + d,
          position = this.rest[j][d] + q[i];
        arms[3 * j + d] = position;
        centre[d] += mass[n] * position;
        linear[d] += mass[n] * v[i] + dt * force[i];
      }
    for (let d = 0; d < 3; d++) {
      centre[d] /= this.mass;
      linear[d] /= this.mass;
    }
    let xx = 0,
      yy = 0,
      zz = 0,
      xy = 0,
      xz = 0,
      yz = 0;
    const angularMomentum = [0, 0, 0];
    for (const [j, n] of ns.entries()) {
      const i = 3 * n,
        k = 3 * j,
        m = mass[n];
      const x = (arms[k] -= centre[0]),
        y = (arms[k + 1] -= centre[1]),
        z = (arms[k + 2] -= centre[2]);
      xx += m * (y * y + z * z);
      yy += m * (x * x + z * z);
      zz += m * (x * x + y * y);
      xy -= m * x * y;
      xz -= m * x * z;
      yz -= m * y * z;
      const px = m * v[i] + dt * force[i],
        py = m * v[i + 1] + dt * force[i + 1],
        pz = m * v[i + 2] + dt * force[i + 2];
      angularMomentum[0] += y * pz - z * py;
      angularMomentum[1] += z * px - x * pz;
      angularMomentum[2] += x * py - y * px;
    }
    const inverse = inverseInertia(xx, yy, zz, xy, xz, yz);
    const R = rotation(product(inverse, ...angularMomentum), dt);
    // I(end)^-1 L = R I(start)^-1 R^T L. The endpoint velocities carry
    // exactly the integrated angular impulse even for a non-cubic group.
    const [lx, ly, lz] = angularMomentum;
    const localOmega = product(
      inverse,
      R[0] * lx + R[3] * ly + R[6] * lz,
      R[1] * lx + R[4] * ly + R[7] * lz,
      R[2] * lx + R[5] * ly + R[8] * lz,
    );
    const [wx, wy, wz] = product(R, ...localOmega);
    for (const [j, n] of ns.entries()) {
      const k = 3 * j,
        i = 3 * n;
      const [x, y, z] = product(R, arms[k], arms[k + 1], arms[k + 2]);
      q[i] = centre[0] + dt * linear[0] + x - this.rest[j][0];
      q[i + 1] = centre[1] + dt * linear[1] + y - this.rest[j][1];
      q[i + 2] = centre[2] + dt * linear[2] + z - this.rest[j][2];
      v[i] = linear[0] + wy * z - wz * y;
      v[i + 1] = linear[1] + wz * x - wx * z;
      v[i + 2] = linear[2] + wx * y - wy * x;
    }
  }
}
