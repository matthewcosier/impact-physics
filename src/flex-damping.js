import { inverseInertia, product } from './rigid-group.js';

// Project velocity onto instantaneous rigid translation and rotation in the
// mass metric. Relax only the orthogonal part, at an effective rate in 1/s.
export class FlexDamping {
  constructor({ nodes, rate, maxAcceleration = Infinity }, reference, masses, dt) {
    if (
      !(maxAcceleration > 0) ||
      !Number.isFinite(rate) ||
      rate < 0 ||
      new Set(nodes).size !== nodes.length ||
      nodes.length < 3 ||
      nodes.some((n) => !Number.isInteger(n) || !reference[n])
    )
      throw new RangeError('Invalid flexural damping group');
    this.nodes = nodes;
    this.rest = Float64Array.from(nodes.flatMap((n) => reference[n].position));
    this.arms = new Float64Array(this.rest.length);
    this.mass = nodes.reduce((sum, n) => sum + masses[n], 0);
    this.maxStep = maxAcceleration * dt;
    this.fraction = -Math.expm1(-rate * dt);
  }

  apply(q, v, mass) {
    const { nodes, arms, rest } = this;
    let cx = 0,
      cy = 0,
      cz = 0,
      vx = 0,
      vy = 0,
      vz = 0;
    for (let j = 0; j < nodes.length; j++) {
      const n = nodes[j],
        m = mass[n],
        i = 3 * n,
        k = 3 * j;
      cx += m * (arms[k] = rest[k] + q[i]);
      cy += m * (arms[k + 1] = rest[k + 1] + q[i + 1]);
      cz += m * (arms[k + 2] = rest[k + 2] + q[i + 2]);
      vx += m * v[i];
      vy += m * v[i + 1];
      vz += m * v[i + 2];
    }
    cx /= this.mass;
    cy /= this.mass;
    cz /= this.mass;
    vx /= this.mass;
    vy /= this.mass;
    vz /= this.mass;
    let xx = 0,
      yy = 0,
      zz = 0,
      xy = 0,
      xz = 0,
      yz = 0,
      lx = 0,
      ly = 0,
      lz = 0;
    for (let j = 0; j < nodes.length; j++) {
      const n = nodes[j],
        m = mass[n],
        i = 3 * n,
        k = 3 * j;
      const x = (arms[k] -= cx),
        y = (arms[k + 1] -= cy),
        z = (arms[k + 2] -= cz);
      xx += m * (y * y + z * z);
      yy += m * (x * x + z * z);
      zz += m * (x * x + y * y);
      xy -= m * x * y;
      xz -= m * x * z;
      yz -= m * y * z;
      lx += m * (y * (v[i + 2] - vz) - z * (v[i + 1] - vy));
      ly += m * (z * (v[i] - vx) - x * (v[i + 2] - vz));
      lz += m * (x * (v[i + 1] - vy) - y * (v[i] - vx));
    }
    // A collinear mass distribution has two equal transverse moments and
    // cannot carry spin about the line. Use its transverse pseudoinverse.
    const trace = xx + yy + zz;
    const det = xx * (yy * zz - yz * yz) + xy * (xz * yz - xy * zz) + xz * (xy * yz - xz * yy);
    const [wx, wy, wz] =
      det > 1e-15 * trace ** 3
        ? product(inverseInertia(xx, yy, zz, xy, xz, yz), lx, ly, lz)
        : trace
          ? [(2 * lx) / trace, (2 * ly) / trace, (2 * lz) / trace]
          : [0, 0, 0];
    let twiceModalEnergy = 0;
    for (let j = 0; j < nodes.length; j++) {
      const n = nodes[j],
        i = 3 * n,
        k = 3 * j,
        x = arms[k],
        y = arms[k + 1],
        z = arms[k + 2];
      const dx = v[i] - vx - (wy * z - wz * y),
        dy = v[i + 1] - vy - (wz * x - wx * z),
        dz = v[i + 2] - vz - (wx * y - wy * x);
      arms[k] = dx;
      arms[k + 1] = dy;
      arms[k + 2] = dz;
      twiceModalEnergy += mass[n] * (dx * dx + dy * dy + dz * dz);
    }
    if (!twiceModalEnergy) return 0;
    const fraction = Math.min(this.fraction, this.maxStep / Math.sqrt(twiceModalEnergy / this.mass));
    for (let j = 0; j < nodes.length; j++) {
      const i = 3 * nodes[j],
        k = 3 * j;
      v[i] -= fraction * arms[k];
      v[i + 1] -= fraction * arms[k + 1];
      v[i + 2] -= fraction * arms[k + 2];
    }
    return 0.5 * fraction * (2 - fraction) * twiceModalEnergy;
  }
}
