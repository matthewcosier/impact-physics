import { contactDepth, horizontalExtent } from './impact-shapes.js';
import { weights } from './wall-contact.js';

// The compliant skin usually takes the initial load. This constraint also
// catches exposed backing, detached pieces and vertices between its witnesses.
export class GuidedContact {
  constructor(w, cells, points) {
    this.cells = cells;
    this.points = points;
    this.h = new Float64Array(24);
    this.rest = Float64Array.from(w.model.nodes.flatMap((n) => n.position));
    // A 2 mm contact skin covers the chord sag between nearby curved contacts.
    this.shape = { ...w.impactor, pad: 0.002 };
    this.extents = [
      horizontalExtent(this.shape, 0),
      horizontalExtent(this.shape, 1),
      this.shape.shape === 'sphere' ? this.shape.radius : this.shape.height / 2,
    ];
    // Resolve normal approach over a short physical rise time as skin and
    // backing engage. This avoids assigning a finite momentum jump to one
    // arbitrarily short integration step. The present pose is checked too.
    this.horizon = 0.003;
    this.extents = this.extents.map((v) => v + this.shape.pad);
    this.cache = -1;
    this.normal = [0, 0, 0];
  }
  possible(w, nodes, margin, time = 0) {
    const c = w.impactor,
      origin = [c.x + time * c.speed * c.axisX, c.y + time * c.speed * c.axisY, c.z];
    for (let d = 0; d < 3; d++) {
      let lo = Infinity,
        hi = -Infinity;
      for (let k = 0; k < 8; k++) {
        const i = 3 * nodes[k] + d,
          v = this.rest[i] + w.q[i] + time * w.v[i];
        lo = Math.min(lo, v);
        hi = Math.max(hi, v);
      }
      if (lo > origin[d] + this.extents[d] + margin || hi < origin[d] - this.extents[d] - margin)
        return false;
    }
    return true;
  }
  worst(w, cell, time) {
    const { nodes, uvw } = cell,
      c = w.impactor,
      h = this.h,
      n = this.normal;
    // Coefficients of the trilinear field, relative to the moving guide.
    for (let d = 0; d < 3; d++) {
      const k = 8 * d;
      for (let j = 0; j < 8; j++) {
        const i = 3 * nodes[j] + d;
        h[k + j] = this.rest[i] + w.q[i] + time * w.v[i];
      }
      h[k + 7] -= h[k + 6] + h[k + 5] + h[k + 3] - h[k + 4] - h[k + 2] - h[k + 1] + h[k];
      h[k + 6] -= h[k + 4] + h[k + 2] - h[k];
      h[k + 5] -= h[k + 4] + h[k + 1] - h[k];
      h[k + 3] -= h[k + 1] + h[k + 2] - h[k];
      h[k + 4] -= h[k];
      h[k + 2] -= h[k];
      h[k + 1] -= h[k];
      h[k] -= d === 0 ? c.x + time * c.speed * c.axisX : d === 1 ? c.y + time * c.speed * c.axisY : c.z;
    }
    let best = 1e-7,
      idx = -1,
      bx = 0,
      by = 0,
      bz = 0;
    for (let i = 0; i < uvw.length; i += 3) {
      const u = uvw[i],
        v = uvw[i + 1],
        z = uvw[i + 2];
      const x = h[0] + u * h[1] + v * (h[2] + u * h[3]) + z * (h[4] + u * h[5] + v * (h[6] + u * h[7]));
      if (Math.abs(x) > this.extents[0]) continue;
      const y = h[8] + u * h[9] + v * (h[10] + u * h[11]) + z * (h[12] + u * h[13] + v * (h[14] + u * h[15]));
      if (Math.abs(y) > this.extents[1]) continue;
      const zz =
        h[16] + u * h[17] + v * (h[18] + u * h[19]) + z * (h[20] + u * h[21] + v * (h[22] + u * h[23]));
      if (Math.abs(zz) > this.extents[2]) continue;
      const depth = contactDepth(this.shape, x, y, zz, n);
      if (depth > best) {
        best = depth;
        idx = i;
        bx = n[0];
        by = n[1];
        bz = n[2];
      }
    }
    if (idx >= 0)
      return { weights: weights(uvw[idx], uvw[idx + 1], uvw[idx + 2]), depth: best, normal: [bx, by, bz] };
  }
  solve(w, impulses) {
    const gen = w.refreshContactBounds();
    if (this.cache !== gen) {
      this.active = this.cells.filter((cell) => this.possible(w, cell.nodes, 2.1 * w.contactMargin + 0.04));
      this.cache = gen;
      this.activePoints = this.points.filter((p) => {
        const margin = 0.04 + w.contactMargin * (1 + p.weights.reduce((a, b) => a + Math.abs(b), 0));
        return p.position.every((x, d) => {
          for (let j = 0; j < p.nodes.length; j++) x += p.weights[j] * w.q[3 * p.nodes[j] + d];
          return (
            Math.abs(x - (d === 0 ? w.impactor.x : d === 1 ? w.impactor.y : w.impactor.z)) <=
            this.extents[d] + margin
          );
        });
      });
    }
    for (let pass = 0; pass < 64; pass++) {
      let changed = false;
      for (const horizon of [0, this.horizon]) {
        const fraction = w.dt / (w.dt + horizon);
        for (const cell of this.active) {
          if (!this.possible(w, cell.nodes, 0, horizon)) continue;
          const hit = this.worst(w, cell, horizon);
          if (hit) {
            impulses.push(w, cell.nodes, hit.weights, hit.normal, (hit.depth + 1e-8) * fraction);
            changed = true;
          }
        }
        for (const p of this.activePoints) {
          const c = w.impactor,
            position = [
              p.position[0] - c.x - horizon * c.speed * c.axisX,
              p.position[1] - c.y - horizon * c.speed * c.axisY,
              p.position[2] - c.z,
            ];
          for (let j = 0; j < p.nodes.length; j++)
            for (let d = 0; d < 3; d++) {
              const i = 3 * p.nodes[j] + d;
              position[d] += p.weights[j] * (w.q[i] + horizon * w.v[i]);
            }
          const depth = contactDepth(this.shape, ...position, this.normal);
          if (depth > 1e-7) {
            impulses.push(w, p.nodes, p.weights, this.normal, (depth + 1e-8) * fraction);
            changed = true;
          }
        }
      }
      if (!changed) return;
    }
    throw new Error('Guided contact did not converge at ' + w.time + ' s');
  }
}
