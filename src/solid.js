import { MetalPoint } from './metal.js';
// Finite-strain compressible neo-Hookean solid. F and P are row-major.
// P is nominal (first Piola) stress, work-conjugate to F per reference volume.
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const dot = (a, b) => a.reduce((sum, x, i) => sum + x * b[i], 0);

export function rubberResponse(F, material) {
  if (F.length !== 9 || !F.every(Number.isFinite)) throw new RangeError('Invalid deformation gradient');
  const { G, bulkModulus: K } = material;
  if (!(G > 0 && K > 0 && Number.isFinite(G) && Number.isFinite(K)))
    throw new RangeError('Invalid elastic moduli');
  const rows = [F.slice(0, 3), F.slice(3, 6), F.slice(6, 9)];
  const cof = [...cross(rows[1], rows[2]), ...cross(rows[2], rows[0]), ...cross(rows[0], rows[1])];
  const J = cof[0] * F[0] + cof[1] * F[1] + cof[2] * F[2];
  if (!(J > 0)) throw new RangeError('Inverted or collapsed solid element');
  const i1 = dot(F, F),
    isochoric = J ** (-2 / 3),
    distortion = isochoric * i1 - 3,
    logJ = Math.log(J);
  const gent = material.solidLaw === 'Gent';
  const limit = material.chainLimit;
  if (gent && !(Number.isFinite(limit) && limit > 0 && distortion < limit))
    throw new RangeError('Invalid or exceeded Gent chain extensibility');
  const scale = (G * isochoric) / (gent ? 1 - distortion / limit : 1);
  const shearEnergy = gent ? (-G * limit * Math.log1p(-distortion / limit)) / 2 : (G * distortion) / 2;
  const pressure = (K * logJ - (scale * i1) / 3) / J;
  const P = F.map((f, i) => scale * f + pressure * cof[i]);
  return { P, J, energyDensity: shearEnergy + 0.5 * K * logJ ** 2 };
}

export class SolidTetrahedron {
  constructor(spec, nodes) {
    this.indices = spec.nodes;
    const law = spec.material?.solidLaw;
    if (law !== 'J2' && law !== 'neo-Hookean' && law !== 'Gent')
      throw new RangeError('No 3D law specified for this material');
    this.metal = law === 'J2' ? new MetalPoint(spec.material) : null;
    this.plasticWork = 0;
    if (spec.nodes.length !== 4 || new Set(spec.nodes).size !== 4 || !spec.nodes.every((i) => nodes[i]))
      throw new RangeError('Expected four distinct solid nodes');
    this.material = spec.material;
    const [base, ...others] = spec.nodes.map((i) => nodes[i].position);
    this.edges = others.map((p) => p.map((x, d) => x - base[d]));
    const [a, b, c] = this.edges;
    const det = dot(a, cross(b, c));
    if (!Number.isFinite(det) || det === 0) throw new RangeError('Degenerate reference tetrahedron');
    this.volume = Math.abs(det) / 6;
    this.mass = this.volume * this.material.density;
    this.gradients = [cross(b, c), cross(c, a), cross(a, b)].map((row) => row.map((x) => x / det));
  }
  evaluate(q, force = null, commit = false) {
    const base = 3 * this.indices[0],
      F = Array(9).fill(0);
    for (let k = 0; k < 3; k++) {
      const i = 3 * this.indices[k + 1],
        edge = this.edges[k],
        grad = this.gradients[k];
      for (let row = 0; row < 3; row++) {
        const dx = edge[row] + q[i + row] - q[base + row];
        for (let col = 0; col < 3; col++) F[3 * row + col] += dx * grad[col];
      }
    }
    const { energyDensity, P } = this.metal
      ? this.metalResponse(F, commit)
      : rubberResponse(F, this.material);
    if (force)
      for (let k = 0; k < 3; k++) {
        const grad = this.gradients[k],
          i = 3 * this.indices[k + 1];
        for (let row = 0; row < 3; row++) {
          const f =
            -this.volume * (P[3 * row] * grad[0] + P[3 * row + 1] * grad[1] + P[3 * row + 2] * grad[2]);
          force[i + row] += f;
          force[base + row] -= f;
        }
      }
    return this.volume * energyDensity;
  }
  // Objective total-Lagrangian extension of the small-strain metal law.
  // Use for small material strains with finite rigid rotations.
  metalResponse(F, commit) {
    if (!(dot(F.slice(0, 3), cross(F.slice(3, 6), F.slice(6))) > 0))
      throw new RangeError('Inverted or collapsed solid element');
    const strain = Array(9),
      P = Array(9);
    for (let i = 0; i < 3; i++)
      for (let j = i; j < 3; j++) {
        const c = F[i] * F[j] + F[3 + i] * F[3 + j] + F[6 + i] * F[6 + j];
        strain[3 * i + j] = strain[3 * j + i] = 0.5 * (c - Number(i === j));
      }
    const { stress: S, energyDensity, plasticWorkDensity } = this.metal.response(strain, commit);
    if (commit) this.plasticWork = plasticWorkDensity * this.volume;
    for (let i = 0; i < 3; i++)
      for (let j = 0; j < 3; j++)
        P[3 * i + j] = F[3 * i] * S[j] + F[3 * i + 1] * S[3 + j] + F[3 * i + 2] * S[6 + j];
    return { P, energyDensity };
  }
}
