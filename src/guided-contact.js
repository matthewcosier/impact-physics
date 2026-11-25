import { contactDepth, horizontalExtent } from './impact-shapes.js';
import { weights } from './wall-contact.js';

// Preserve contact order while caching only witnesses close enough to reach
// the solid. For a trilinear cell, each point's change is bounded by the
// largest change at a corner. Signed affine weights use their L1 norm.
const cacheTravel = 0.004;
function nearShape(c, x, y, z, margin) {
  if (c.shape === 'sphere') return x * x + y * y + z * z <= (c.radius + c.pad + margin) ** 2;
  const band = c.pad + margin;
  if (Math.abs(z) > c.height / 2 + band) return false;
  if (c.shape === 'cylinder') return x * x + y * y <= (c.radius + band) ** 2;
  return (
    Math.abs(x * c.axisX + y * c.axisY) <= c.depth / 2 + band &&
    Math.abs(y * c.axisX - x * c.axisY) <= c.width / 2 + band
  );
}
function trilinear(h, k, u, v, z) {
  return (
    h[k] +
    u * h[k + 1] +
    v * (h[k + 2] + u * h[k + 3]) +
    z * (h[k + 4] + u * h[k + 5] + v * (h[k + 6] + u * h[k + 7]))
  );
}
function pointOffset(w, p, time, c, out) {
  let x = p.position[0] - c.x - time * c.speed * c.axisX,
    y = p.position[1] - c.y - time * c.speed * c.axisY,
    z = p.position[2] - c.z;
  for (let j = 0; j < p.nodes.length; j++) {
    const i = 3 * p.nodes[j],
      a = p.weights[j];
    x += a * (w.q[i] + time * w.v[i]);
    y += a * (w.q[i + 1] + time * w.v[i + 1]);
    z += a * (w.q[i + 2] + time * w.v[i + 2]);
  }
  out[0] = x;
  out[1] = y;
  out[2] = z;
}
function pointGroups(points) {
  const groups = [];
  for (let i = 0; i < points.length; i += 64) {
    const group = points.slice(i, i + 64);
    groups.push({
      points: group,
      nodes: [...new Set(group.flatMap((p) => Array.from(p.nodes)))],
      norm: Math.max(...group.map((p) => p.weights.reduce((a, w) => a + Math.abs(w), 0))),
      caches: [],
    });
  }
  return groups;
}

// The compliant skin usually takes the initial load. This constraint also
// catches exposed backing, detached pieces and vertices between its witnesses.
export class GuidedContact {
  constructor(w, cells, points) {
    this.cells = cells.map((c) => ({ ...c, caches: [] }));
    this.offset = [0, 0, 0];
    this.points = points;
    this.h = new Float64Array(24);
    this.rest = Float64Array.from(w.model.nodes.flatMap((n) => n.position));
    // A 2 mm contact skin covers the chord sag between nearby curved contacts.
    this.shape = { ...w.impactor, pad: w.impactor.surfaceMargin ?? 0.002 };
    this.extents = [
      horizontalExtent(this.shape, 0),
      horizontalExtent(this.shape, 1),
      this.shape.shape === 'sphere' ? this.shape.radius : this.shape.height / 2,
    ];
    // Resolve normal approach over a short physical rise time as skin and
    // backing engage. This avoids assigning a finite momentum jump to one
    // arbitrarily short integration step. The present pose is checked too.
    this.horizon = w.impactor.approachTime ?? 0.003;
    if (![this.shape.pad, this.horizon].every((x) => Number.isFinite(x) && x >= 0))
      throw new RangeError('Expected nonnegative contact skin thickness and prediction time');
    this.extents = this.extents.map((v) => v + this.shape.pad);
    this.horizons = this.horizon ? [0, this.horizon] : [0];
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
    const phase = time === 0 ? 0 : 1;
    const cached = (cell.caches[phase] ??= {
      last: new Float64Array(24),
      indices: new Int32Array(uvw.length / 3),
      count: -1,
    });
    let rebuild = cached.count < 0;
    // Coefficients of the trilinear field, relative to the moving guide.
    for (let d = 0; d < 3; d++) {
      const k = 8 * d;
      for (let j = 0; j < 8; j++) {
        const i = 3 * nodes[j] + d;
        const origin =
          d === 0 ? c.x + time * c.speed * c.axisX : d === 1 ? c.y + time * c.speed * c.axisY : c.z;
        h[k + j] = this.rest[i] + w.q[i] + time * w.v[i] - origin;
        if (Math.abs(h[k + j] - cached.last[k + j]) >= cacheTravel) rebuild = true;
      }
      h[k + 7] -= h[k + 6] + h[k + 5] + h[k + 3] - h[k + 4] - h[k + 2] - h[k + 1] + h[k];
      h[k + 6] -= h[k + 4] + h[k + 2] - h[k];
      h[k + 5] -= h[k + 4] + h[k + 1] - h[k];
      h[k + 3] -= h[k + 1] + h[k + 2] - h[k];
      h[k + 4] -= h[k];
      h[k + 2] -= h[k];
      h[k + 1] -= h[k];
    }
    let best = 1e-7,
      idx = -1,
      bx = 0,
      by = 0,
      bz = 0;
    if (rebuild) {
      for (let k = 0; k < 8; k++)
        for (let d = 0; d < 3; d++) {
          const j = 3 * nodes[k] + d,
            origin =
              d === 0 ? c.x + time * c.speed * c.axisX : d === 1 ? c.y + time * c.speed * c.axisY : c.z;
          cached.last[8 * d + k] = this.rest[j] + w.q[j] + time * w.v[j] - origin;
        }
      cached.count = 0;
      for (let i = 0; i < uvw.length; i += 3) {
        const u = uvw[i],
          v = uvw[i + 1],
          z = uvw[i + 2];
        if (
          nearShape(
            this.shape,
            trilinear(h, 0, u, v, z),
            trilinear(h, 8, u, v, z),
            trilinear(h, 16, u, v, z),
            Math.sqrt(3) * cacheTravel,
          )
        )
          cached.indices[cached.count++] = i;
      }
    }
    for (let j = 0; j < cached.count; j++) {
      const i = cached.indices[j],
        u = uvw[i],
        v = uvw[i + 1],
        z = uvw[i + 2],
        x = trilinear(h, 0, u, v, z);
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
  pointCandidates(w, group, time) {
    const phase = time === 0 ? 0 : 1,
      c = w.impactor;
    const cached = (group.caches[phase] ??= { last: new Float64Array(3 * group.nodes.length), near: null });
    let dirty = !cached.near;
    const origin = [c.x + time * c.speed * c.axisX, c.y + time * c.speed * c.axisY, c.z];
    for (let j = 0; !dirty && j < group.nodes.length; j++)
      for (let d = 0; d < 3; d++) {
        const i = 3 * group.nodes[j] + d;
        if (Math.abs(w.q[i] + time * w.v[i] - origin[d] - cached.last[3 * j + d]) >= cacheTravel) {
          dirty = true;
          break;
        }
      }
    if (dirty) {
      for (let j = 0; j < group.nodes.length; j++)
        for (let d = 0; d < 3; d++) {
          const i = 3 * group.nodes[j] + d;
          cached.last[3 * j + d] = w.q[i] + time * w.v[i] - origin[d];
        }
      cached.near = group.points.filter((p) => {
        pointOffset(w, p, time, c, this.offset);
        return nearShape(this.shape, ...this.offset, Math.sqrt(3) * group.norm * cacheTravel);
      });
    }
    return cached.near;
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
      this.groups = pointGroups(this.activePoints);
    }
    for (let pass = 0; pass < 64; pass++) {
      let changed = false;
      for (const horizon of this.horizons) {
        const fraction = w.dt / (w.dt + horizon);
        for (const cell of this.active) {
          if (!this.possible(w, cell.nodes, 0, horizon)) continue;
          const hit = this.worst(w, cell, horizon);
          if (hit) {
            impulses.push(w, cell.nodes, hit.weights, hit.normal, (hit.depth + 1e-8) * fraction);
            changed = true;
          }
        }
        for (const group of this.groups)
          for (const p of this.pointCandidates(w, group, horizon)) {
            pointOffset(w, p, horizon, w.impactor, this.offset);
            const depth = contactDepth(this.shape, ...this.offset, this.normal);
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
