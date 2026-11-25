import { bind as bindPoint } from './lattice.js';
// A static plane supporting the complete skin. Inside a trilinear cell every
// vertex is a convex combination of its eight corners. A corner bound rejects
// whole cells; the remaining cells use the actual rendered vertices.
export function surfaceCells(parts, binary, bind = bindPoint) {
  const cells = new Map();
  for (const part of parts) {
    for (const nodes of part.cells.values()) cells.set(nodes, []);
    if (!binary) continue;
    for (const entry of part.objects) {
      const verts = new Float32Array(binary, entry.positions, 3 * entry.vertices);
      for (let i = 0; i < verts.length; i += 3) {
        const p = bind(part, [verts[i], verts[i + 1], verts[i + 2]]);
        if (p.uvw.some((v) => v < -1e-6 || v > 1 + 1e-6))
          throw new RangeError('Surface vertex outside its occupied lattice');
        cells.get(p.nodes).push(...p.uvw);
      }
    }
  }
  const corners = Array.from({ length: 8 }, (_, i) => [i & 1, (i >> 1) & 1, (i >> 2) & 1]).flat();
  return Array.from(cells, ([nodes, uvw]) => ({ nodes, uvw: Float32Array.from(binary ? uvw : corners) }));
}

export const weights = (u, v, w) => [
  (1 - u) * (1 - v) * (1 - w),
  u * (1 - v) * (1 - w),
  (1 - u) * v * (1 - w),
  u * v * (1 - w),
  (1 - u) * (1 - v) * w,
  u * (1 - v) * w,
  (1 - u) * v * w,
  u * v * w,
];

export class WallContact {
  constructor(world, cells, points) {
    this.cells = cells;
    this.points = points;
    this.weightBounds = points.map((p) => p.weights.reduce((sum, v) => sum + Math.abs(v), 0));
    this.generation = -1;
    this.restX = Float64Array.from(world.model.nodes, (n) => n.position[0]);
    this.x = new Float64Array(this.restX.length);
    this.corner = new Float64Array(8);
  }
  solve(w, impulses) {
    const plane = w.impactor.x,
      x = this.x,
      h = this.corner;
    for (let i = 0; i < x.length; i++) x[i] = this.restX[i] + w.q[3 * i];
    const normal = [1, 0, 0];
    for (let pass = 0; pass < 96; pass++) {
      // A chassis skin can contain thousands of affine material points. An
      // L1 weight bound also covers points extrapolated beyond their nodes.
      const generation = w.refreshContactBounds();
      if (generation !== this.generation) {
        this.activePoints = this.points.filter((p, i) => {
          let x = p.position[0];
          for (let j = 0; j < p.nodes.length; j++) x += p.weights[j] * w.q[3 * p.nodes[j]];
          return x <= plane + 2.1 * w.contactMargin * this.weightBounds[i];
        });
        this.generation = generation;
      }
      let applied = false;
      for (const { nodes, uvw } of this.cells) {
        let min = Infinity;
        for (let k = 0; k < 8; k++) {
          h[k] = x[nodes[k]];
          min = Math.min(min, h[k]);
        }
        if (min >= plane) continue;
        let worst = plane - 1e-8,
          idx = -1;
        for (let j = 0; j < uvw.length; j += 3) {
          const u = uvw[j],
            v = uvw[j + 1],
            z = uvw[j + 2];
          const a = h[0] + u * (h[1] - h[0]),
            b = h[2] + u * (h[3] - h[2]);
          const c = h[4] + u * (h[5] - h[4]),
            d = h[6] + u * (h[7] - h[6]);
          const low = a + v * (b - a),
            high = c + v * (d - c),
            s = low + z * (high - low);
          if (s < worst) {
            worst = s;
            idx = j;
          }
        }
        if (idx >= 0) {
          impulses.push(w, nodes, weights(uvw[idx], uvw[idx + 1], uvw[idx + 2]), normal, plane - worst);
          for (const n of impulses.groups.get(nodes[0])?.group.indices ?? nodes)
            x[n] = this.restX[n] + w.q[3 * n];
          applied = true;
        }
      }
      for (const p of this.activePoints) {
        let s = p.position[0];
        for (let j = 0; j < p.nodes.length; j++) s += p.weights[j] * w.q[3 * p.nodes[j]];
        if (s < plane - 1e-8) {
          impulses.push(w, p.nodes, p.weights, normal, plane - s);
          for (const n of impulses.groups.get(p.nodes[0])?.group.indices ?? p.nodes)
            x[n] = this.restX[n] + w.q[3 * n];
          applied = true;
        }
      }
      if (!applied) return;
    }
    throw new Error('Wall contact did not converge at ' + w.time + ' s');
  }
}
