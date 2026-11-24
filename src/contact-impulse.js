// End-of-step nonpenetration impulses, distributed by the same shape functions
// that move the render skin. The finite impactor retains one guide coordinate.
export class ContactImpulses {
  constructor(world) {
    this.startV = new Float64Array(world.v.length);
    this.impulses = new Float64Array(world.v.length);
    this.groups = new Map();
    for (const group of world.rigidGroups) {
      const state = {
        group,
        q: new Float64Array(3 * group.indices.length),
        v: new Float64Array(3 * group.indices.length),
      };
      for (const n of group.indices) this.groups.set(n, state);
    }
    this.rigid = [...new Set(this.groups.values())];
  }
  begin(w) {
    this.startV.set(w.v);
    this.impulses.fill(0);
    this.axial = 0;
    this.total = 0;
    this.startSpeed = w.impactor.speed;
    for (const { group, q, v } of this.rigid)
      for (const [j, n] of group.indices.entries())
        for (let d = 0; d < 3; d++) {
          q[3 * j + d] = w.q[3 * n + d];
          v[3 * j + d] = w.v[3 * n + d];
        }
  }
  push(w, nodes, weights, normal, depth) {
    const c = w.impactor,
      guided = c.shape !== 'wall',
      axis = normal[0] * c.axisX + normal[1] * c.axisY;
    const mobilityGuide = guided ? (axis * axis) / c.mass : 0;
    let mobility = mobilityGuide;
    for (let j = 0; j < nodes.length; j++) mobility += weights[j] ** 2 / w.mass[nodes[j]];
    let J = depth / (w.dt * mobility);
    const state = this.groups.get(nodes[0]);
    if (state) {
      const { group, q, v } = state;
      // Reintegrate with the contact force from the same pre-integration pose.
      // A trial impulse measures the compliance of the finite rigid rotation.
      const replay = (impulse) => {
        for (let j = 0; j < nodes.length; j++)
          for (let d = 0; d < 3; d++) w.force[3 * nodes[j] + d] += (impulse * weights[j] * normal[d]) / w.dt;
        for (let j = 0; j < group.indices.length; j++)
          for (let d = 0; d < 3; d++) {
            const i = 3 * group.indices[j] + d;
            w.q[i] = q[3 * j + d];
            w.v[i] = v[3 * j + d];
          }
        group.advance(w);
      };
      const position = () => {
        let s = 0;
        for (let j = 0; j < nodes.length; j++)
          for (let d = 0; d < 3; d++) s += weights[j] * normal[d] * w.q[3 * nodes[j] + d];
        return s;
      };
      const start = position();
      replay(J);
      const moved = position() - start + w.dt * J * mobilityGuide;
      if (moved > 1e-12) {
        const refined = (J * depth) / moved;
        replay(refined - J);
        J = refined;
      }
    } else {
      for (let j = 0; j < nodes.length; j++)
        for (let d = 0; d < 3; d++) {
          const n = nodes[j],
            i = 3 * n + d,
            dv = (J * weights[j] * normal[d]) / w.mass[n];
          w.v[i] += dv;
          w.q[i] += w.dt * dv;
        }
    }
    for (let j = 0; j < nodes.length; j++)
      for (let d = 0; d < 3; d++) this.impulses[3 * nodes[j] + d] += J * weights[j] * normal[d];
    if (guided) {
      const dv = (-J * axis) / c.mass;
      c.speed += dv;
      c.x += w.dt * dv * c.axisX;
      c.y += w.dt * dv * c.axisY;
      w.externalImpulseX += J * (normal[0] - axis * c.axisX);
    }
    this.axial += J * axis;
    this.total += J;
  }
  finish(w) {
    if (!this.total) return;
    // Midpoint velocities cover the physical step, including the simultaneous
    // elastic and support forces. A static normal reaction does zero work.
    let work = 0;
    for (let i = 0; i < w.v.length; i++) work += this.impulses[i] * (this.startV[i] + w.v[i]);
    if (w.impactor.shape !== 'wall') work -= this.axial * (this.startSpeed + w.impactor.speed);
    w.loss.damping -= 0.5 * work;
    w.contactImpulse += this.axial;
    w.contactForce += this.axial / w.dt;
    w.peakForce = Math.max(w.peakForce, Math.abs(w.contactForce));
  }
}
