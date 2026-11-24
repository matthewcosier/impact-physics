import { contactStencil } from './contact-damping.js';
// A braced eight-node cage transports a compound of box volumes. Its affine
// motion is the least-squares fit to the eight corners. The same interpolation
// and virtual-work transpose drive the collision surface, reaction and display.
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const sign = (i, d) => (i & (1 << d) ? 1 : -1);
export function cageWeights(body, local, out = new Float64Array(8)) {
  const x = local[0] / body.half[0],
    y = local[1] / body.half[1],
    z = local[2] / body.half[2];
  for (let i = 0; i < 8; i++) out[i] = (1 + sign(i, 0) * x + sign(i, 1) * y + sign(i, 2) * z) / 8;
  return out;
}
export function cageFrame(body, q) {
  const c = body.center.slice(),
    axis = body.referenceAxes?.map((a) => [...a]) ?? [
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1],
    ];
  for (let i = 0; i < 8; i++) {
    const n = 3 * body.nodes[i];
    for (let d = 0; d < 3; d++) {
      const v = q[n + d] / 8;
      c[d] += v;
      for (let j = 0; j < 3; j++) axis[j][d] += (sign(i, j) * v) / body.half[j];
    }
  }
  const inverse = [cross(axis[1], axis[2]), cross(axis[2], axis[0]), cross(axis[0], axis[1])];
  const det = dot(axis[0], inverse[0]);
  if (!(det > 1e-6)) throw new RangeError('Collapsed contact cage: ' + body.name);
  for (const row of inverse) for (let d = 0; d < 3; d++) row[d] /= det;
  return { c, axis, inverse };
}
export class BoxContact {
  constructor(body, mass) {
    this.body = body;
    this.witnesses = body.witnesses.map((p) => ({
      ...p,
      entryBox: null,
      entryAxis: 0,
      entrySign: 0,
      mobility: p.nodes.reduce((a, n, j) => a + p.weights[j] ** 2 / mass[n], 0),
    }));
    this.weightBounds = this.witnesses.map((p) => p.weights.reduce((a, w) => a + Math.abs(w), 0));
    this.volumeWeightBound = 0;
    for (const box of body.volumes)
      for (let i = 0; i < 8; i++) {
        const local = box.center.map((x, d) => x + sign(i, d) * box.half[d]);
        const bound = cageWeights(body, local).reduce((a, w) => a + Math.abs(w), 0);
        this.volumeWeightBound = Math.max(this.volumeWeightBound, bound);
      }
    this.contacts = [];
    this.impulses = {};
    this.maxDepth = 0;
    this.depthByOwner = {};
    this.peakTravel = 0;
  }
  candidates(frame, q, margin) {
    const low = [Infinity, Infinity, Infinity],
      high = [-Infinity, -Infinity, -Infinity];
    for (const box of this.body.volumes)
      for (let d = 0; d < 3; d++) {
        let middle = frame.c[d],
          width = 0;
        for (let j = 0; j < 3; j++) {
          middle += frame.axis[j][d] * box.center[j];
          width += Math.abs(frame.axis[j][d]) * box.half[j];
        }
        low[d] = Math.min(low[d], middle - width);
        high[d] = Math.max(high[d], middle + width);
      }
    // Both the interpolated witness and the transported box can move
    // before the next refresh. Extrapolating weights need their L1 norm.
    return this.witnesses.filter((p, index) => {
      if (p.entryBox) return true;
      const pad = margin * (this.weightBounds[index] + this.volumeWeightBound);
      for (let d = 0; d < 3; d++) {
        let x = p.position[d];
        for (let j = 0; j < p.nodes.length; j++) x += p.weights[j] * q[3 * p.nodes[j] + d];
        if (x < low[d] - pad || x > high[d] + pad) return false;
      }
      return true;
    });
  }
  apply(world, measureOnly = false, rebuild = true, margin = 0) {
    const body = this.body,
      { q, v, force, mass, dt } = world,
      frame = cageFrame(body, q),
      { c, inverse } = frame;
    if (rebuild || !this.active) this.active = this.candidates(frame, q, margin);
    const lens = inverse.map((r) => Math.hypot(...r));
    const local = [0, 0, 0],
      projected = [0, 0, 0],
      weights = new Float64Array(8);
    let energy = 0;
    if (!measureOnly) this.contacts.length = 0;
    if (!measureOnly)
      this.peakTravel = Math.max(
        this.peakTravel,
        Math.hypot(
          c[0] - body.center[0] - (world.sled?.x ?? 0),
          c[1] - body.center[1],
          c[2] - body.center[2],
        ),
      );
    for (const p of this.active) {
      let x = p.position[0] - c[0],
        y = p.position[1] - c[1],
        z = p.position[2] - c[2];
      for (let j = 0; j < p.nodes.length; j++) {
        const i = p.nodes[j] * 3,
          w = p.weights[j];
        x += w * q[i];
        y += w * q[i + 1];
        z += w * q[i + 2];
      }
      for (let d = 0; d < 3; d++) {
        const r = inverse[d];
        local[d] = r[0] * x + r[1] * y + r[2] * z;
      }
      let depth = 0,
        face = 0,
        side = 1,
        selected = null;
      for (const box of body.volumes) {
        let nearest = Infinity,
          axis = 0,
          normal = 1;
        for (let d = 0; d < 3; d++) {
          const diff = local[d] - box.center[d],
            normalSide = p.entryBox === box && d === p.entryAxis ? p.entrySign : diff < 0 ? -1 : 1,
            gap = (box.half[d] - normalSide * diff) / lens[d];
          if (gap <= depth) {
            nearest = gap;
            break;
          }
          if (gap < nearest) {
            nearest = gap;
            axis = d;
            normal = normalSide;
          }
        }
        if (nearest > depth) {
          selected = box;
          depth = nearest;
          face = axis;
          side = normal;
        }
      }
      if (!depth) {
        if (!measureOnly) p.entryBox = null;
        continue;
      }
      // Keep the incident face on this axis until separation. A soft
      // penalty can compress past the mid-plane of a thin plate.
      if (!measureOnly && p.entryBox !== selected) {
        p.entryBox = selected;
        p.entryAxis = face;
        p.entrySign = side;
      }
      const k = p.contactStiffness ?? body.stiffness;
      const stop = p.backstop,
        travel = stop ? Math.max(0, depth - stop.travel) : 0;
      const extra = travel ? stop.stiffness : 0;
      energy += 0.5 * (k * depth ** 2 + extra * travel ** 2);
      if (measureOnly) continue;
      const row = inverse[face],
        scale = side / lens[face],
        nx = row[0] * scale,
        ny = row[1] * scale,
        nz = row[2] * scale;
      for (let d = 0; d < 3; d++) {
        const r = inverse[d];
        projected[d] = local[d] + depth * (r[0] * nx + r[1] * ny + r[2] * nz);
      }
      cageWeights(body, projected, weights);
      const normalStencil = contactStencil(p, body, weights, mass);
      const magnitude = k * depth + extra * travel;
      normalStencil.add(force, normalStencil.weights, magnitude * nx, magnitude * ny, magnitude * nz);
      const frictionLimit = (body.friction ?? 0) * magnitude;
      // For tangential forces apply the reaction at the witness itself.
      // The elastic normal force acts through the projected face point.
      const tangentStencil = frictionLimit
        ? contactStencil(p, body, cageWeights(body, local, weights), mass)
        : null;
      this.contacts.push({
        body: this,
        owner: p.owner,
        normal: [nx, ny, nz],
        normalStencil,
        tangentStencil,
        damping: 2 * body.damping * Math.sqrt((k + extra) / normalStencil.mobility),
        frictionLimit,
        normalForce: 0,
        fx: 0,
        fy: 0,
        fz: 0,
      });
      this.impulses[p.owner] = (this.impulses[p.owner] ?? 0) + dt * magnitude;
      this.maxDepth = Math.max(this.maxDepth, depth);
      this.depthByOwner[p.owner] = Math.max(this.depthByOwner[p.owner] ?? 0, depth);
    }
    return energy;
  }
}
