// The supported clip starts in static equilibrium. An island that has lost
// every load path to the frame loses its support reaction and falls under gravity.
export class ReleasedParts {
  constructor(model) {
    const packages = (model.contactBodies ?? []).map((body) => ({ id: body.name, nodeIDs: body.nodes }));
    this.parts = [...model.parts, ...packages];
    this.released = new Set();
    this.nodes = [];
    this.referenceHeight = [];
    this.groundLoss = 0;
    this.groundImpulseX = 0;
    this.model = model;
    this.pickups = (model.gravityPickups ?? []).map((p) => ({ ...p, active: false, startHeight: 0 }));
    this.transferredPotential = 0;
  }

  pickupHeight(p, q) {
    return p.nodes.reduce((sum, n, i) => sum + p.weights[i] * q[3 * n + 2], 0);
  }

  reconnect(joints, q) {
    const neighbours = new Map(this.parts.map((p) => [p.id, []]));
    const supported = new Set();
    for (const joint of joints) {
      if (joint.damage >= 1) continue;
      const a = joint.a.owner,
        b = joint.b?.owner;
      if (neighbours.has(a) && neighbours.has(b)) {
        neighbours.get(a).push(b);
        neighbours.get(b).push(a);
      } else if (neighbours.has(a)) supported.add(a);
      else if (neighbours.has(b)) supported.add(b);
    }
    for (const id of supported) {
      for (const other of neighbours.get(id)) supported.add(other);
    }
    for (const part of this.parts) {
      if (supported.has(part.id) || this.released.has(part.id)) continue;
      this.released.add(part.id);
      for (const n of part.nodeIDs) {
        this.nodes.push(n);
        this.referenceHeight.push(q[n * 3 + 2]);
      }
    }
    const alive = new Set(joints.filter((j) => j.damage < 1).map((j) => j.source));
    for (const p of this.pickups) {
      if (this.released.has(p.owner)) {
        if (p.active)
          this.transferredPotential += 9.81 * p.mass * (this.pickupHeight(p.point, q) - p.startHeight);
        p.active = false;
      } else if (!p.active && p.stages.some((stage) => !stage.some((j) => alive.has(j)))) {
        p.active = true;
        p.startHeight = this.pickupHeight(p.point, q);
      }
    }
  }

  apply(world, measureOnly = false) {
    let potential = this.transferredPotential;
    const { q, v, force, mass, dt } = world;
    for (const p of this.pickups)
      if (p.active) {
        const weight = 9.81 * p.mass;
        potential += weight * (this.pickupHeight(p.point, q) - p.startHeight);
        if (!measureOnly)
          for (let j = 0; j < p.point.nodes.length; j++)
            force[p.point.nodes[j] * 3 + 2] -= weight * p.point.weights[j];
      }
    for (let k = 0; k < this.nodes.length; k++) {
      const n = this.nodes[k],
        i = n * 3;
      potential += mass[n] * 9.81 * (q[i + 2] - this.referenceHeight[k]);
      if (!measureOnly) force[i + 2] -= mass[n] * 9.81;
      const depth = -(this.model.nodes[n].position[2] + q[i + 2]);
      if (depth <= 0) continue;
      const stiffness = mass[n] * 18000;
      potential += 0.5 * stiffness * depth * depth;
      if (measureOnly) continue;
      const closing = Math.min(0, v[i + 2]),
        damping = mass[n] * 175;
      const normal = stiffness * depth - damping * closing;
      force[i + 2] += normal;
      const speed = Math.hypot(v[i], v[i + 1]);
      const drag = Math.min(mass[n] / dt, (0.46 * normal) / Math.max(0.001, speed));
      force[i] -= drag * v[i];
      force[i + 1] -= drag * v[i + 1];
      this.groundImpulseX += dt * drag * v[i];
      this.groundLoss += dt * (drag * speed * speed + damping * closing * closing);
    }
    return potential;
  }
}
