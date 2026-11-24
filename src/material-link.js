import { MaterialPoint } from './material-point.js';

// Two nodes carry the ends of a straight reference specimen. Area and length
// are undeformed dimensions. Nominal stress × initial area is axial force.
export class MaterialLink {
  constructor(spec, nodes) {
    this.dampingRatio = spec.dampingRatio ?? 0;
    if (!Number.isFinite(this.dampingRatio) || this.dampingRatio < 0)
      throw new RangeError('Invalid link damping ratio');
    this.a = spec.a;
    this.b = spec.b;
    this.area = spec.area;
    if (!Number.isFinite(this.area) || this.area <= 0 || !nodes[this.a] || !nodes[this.b])
      throw new RangeError('Invalid material link endpoints or cross-section');
    this.delta = nodes[this.b].position.map((x, d) => x - nodes[this.a].position[d]);
    this.length = Math.hypot(...this.delta);
    if (!(this.length > 0) || !Number.isFinite(this.length)) throw new RangeError('Invalid reference length');
    this.volume = this.area * this.length;
    this.mass = this.volume * spec.material.density;
    this.point = new MaterialPoint(spec.material, spec.across);
    this.direction = this.delta.map((x) => x / this.length);
    this.axialForce = this.materialWork = this.strain = 0;
  }
  sample(q) {
    const i = 3 * this.a,
      j = 3 * this.b;
    for (let d = 0; d < 3; d++) this.direction[d] = this.delta[d] + q[j + d] - q[i + d];
    const length = Math.hypot(...this.direction);
    if (length < 1e-12) throw new RangeError('Material link collapsed to a point');
    for (let d = 0; d < 3; d++) this.direction[d] /= length;
    const strain = (length - this.length) / this.length;
    // Unbonded foam pads lose contact on extension.
    const stress = this.point.load(this.point.material.law === 'foam' ? Math.min(0, strain) : strain);
    const force = stress * this.area;
    this.materialWork += 0.5 * (force + this.axialForce) * (strain - this.strain) * this.length;
    this.strain = strain;
    this.axialForce = force;
  }
  scatter(force, magnitude) {
    for (let d = 0; d < 3; d++) {
      const f = this.direction[d] * magnitude;
      force[3 * this.a + d] += f;
      force[3 * this.b + d] -= f;
    }
  }
  prepareViscosity(world) {
    this.viscousForce = 0;
    this.viscosity = (this.point.recoveryViscosity() * this.area) / this.length;
    if (this.dampingRatio && !this.point.failed) {
      const mobility = 1 / world.mass[this.a] + 1 / world.mass[this.b];
      this.viscosity +=
        2 * this.dampingRatio * Math.sqrt((this.point.material.E * this.area) / (this.length * mobility));
    }
    return this.viscosity > 0;
  }
  predictedRate({ dt, v, force, mass }) {
    let rate = 0;
    for (let d = 0; d < 3; d++) {
      const i = 3 * this.a + d,
        j = 3 * this.b + d;
      rate += this.direction[d] * (v[j] - v[i] + dt * (force[j] / mass[this.b] - force[i] / mass[this.a]));
    }
    return rate;
  }
  relaxViscosity(world) {
    const c = this.viscosity;
    const mobility = 1 / world.mass[this.a] + 1 / world.mass[this.b];
    let next =
      this.viscousForce + (c * this.predictedRate(world) - this.viscousForce) / (1 + world.dt * c * mobility);
    if (this.point.material.law === 'foam') next = Math.min(next, -this.axialForce);
    const change = next - this.viscousForce;
    this.scatter(world.force, change);
    this.viscousForce = next;
    return Math.abs(change) / (1 + Math.abs(next));
  }
}

// Backward Euler solves all participating dashpots against the same nodal
// velocities. Projected iterations also enforce zero tensile traction on pads.
export function dampMaterialLinks(world) {
  const active = world.materialLinks.filter((link) => link.prepareViscosity(world));
  if (!active.length) return 0;
  for (let iteration = 0; iteration < 256; iteration++) {
    let change = 0;
    for (const link of active) change = Math.max(change, link.relaxViscosity(world));
    if (change < 1e-9) {
      return active.reduce((sum, link) => sum + world.dt * link.viscousForce * link.predictedRate(world), 0);
    }
  }
  throw new Error('Material damping did not converge. Reduce the time step.');
}
