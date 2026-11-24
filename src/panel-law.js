import { RestrainedTorsion } from './panel-torsion.js';
// Rate-independent bilinear section response under opposed, prescribed grip couples.
// EI and GJ act on the actual section dimensions. Opening the unloaded grips after
// fracture is a separate bench action; it does not inject energy into the test.
export class PanelBench {
  constructor(specimen, material, mode = 'bend') {
    this.specimen = specimen;
    this.material = material;
    this.mode = mode;
    this.stationLength =
      mode === 'twist' && specimen.gripAxis ? specimen.ds / specimen.gripAxis[0] : specimen.ds;
    const twist = mode === 'twist',
      root3 = Math.sqrt(3);
    this.laws = specimen.sections.map((s) => {
      const modulus = twist ? material.G : material.E;
      const geometric = twist ? s.J : s.I;
      const lever = twist ? s.torsionLever : s.bendLever;
      const yieldStress = twist
        ? material.yieldStress / root3
        : Math.min(material.yieldStress, material.compressionYield);
      const brittle = material.law === 'brittle' || material.law === 'laminate';
      const strength = brittle
        ? Math.min(material.ultimate, material.E * material.fractureStrain, material.compressionYield)
        : material.ultimate;
      const ultimate = twist ? strength / root3 : strength;
      const fractureStrain = twist ? material.fractureStrain * root3 : material.fractureStrain;
      const hardening = (ultimate - yieldStress) / Math.max(fractureStrain - yieldStress / modulus, 1e-6);
      const law = {
        brittle,
        stiffness: modulus * geometric,
        hardening: Math.min(0.95 * modulus, hardening) * geometric,
        yield: ((brittle ? ultimate : yieldStress) * geometric) / lever,
        limit: (ultimate * geometric) / lever,
        plastic: 0,
        k: 0,
      };
      if (twist && s.fibres) law.torsion = new RestrainedTorsion(s, material, law);
      return law;
    });
    this.firstYield = Math.min(...this.laws.map((s) => s.torsion?.firstYield ?? s.yield));
    this.capacity = Math.min(...this.laws.map((s) => s.torsion?.capacity ?? s.limit));
    this.breakIndex = -1;
    this.gap = 0;
    this.setMoment(0);
  }
  setMoment(request) {
    if (!Number.isFinite(request)) throw new RangeError('Expected finite applied moment');
    const intact = this.breakIndex < 0;
    const sign = Math.sign(request) || 1;
    const limits = this.laws.map((s) =>
      s.torsion ? sign * s.torsion.torque(sign * s.torsion.limitRate) : s.limit,
    );
    const capacity = Math.min(...limits);
    const applied = intact ? sign * Math.min(Math.abs(request), capacity) : 0;
    for (const s of this.laws) {
      if (s.torsion) {
        s.k = s.torsion.rate(applied, true);
        s.freeRate = s.torsion.rate(0);
        continue;
      }
      // Linear kinematic hardening; unloading preserves set until reverse yield.
      const H = (s.stiffness * s.hardening) / (s.stiffness - s.hardening);
      const trial = applied - H * s.plastic;
      if (!s.brittle && Math.abs(trial) > s.yield && H > 0)
        s.plastic += (trial - Math.sign(trial) * s.yield) / H;
      s.k = s.plastic + applied / s.stiffness;
    }
    if (intact && Math.abs(request) >= capacity * (1 - 1e-12)) {
      this.breakIndex = limits.indexOf(capacity);
      this.failureMoment = applied;
      return this.setMoment(0);
    }
    this.moment = applied;
    this.angle = this.laws.reduce((a, s) => a + s.k * this.stationLength, 0);
    this.setAngle = this.laws.reduce((a, s) => a + (s.freeRate ?? s.plastic) * this.stationLength, 0);
    this.buildFrames();
  }
  openGrips() {
    if (this.breakIndex >= 0) this.gap = this.specimen.length * 0.11;
  }
  buildFrames() {
    const { sections, ds, start } = this.specimen;
    const mid = (i, key) =>
      (sections[Math.max(0, i - 1)][key] + sections[Math.min(sections.length - 1, i)][key]) / 2;
    this.frames = [];
    let angle = 0,
      Y = mid(0, 'cy'),
      Z = mid(0, 'cz'),
      X = start;
    for (let i = 0; i <= sections.length; i++) {
      const cy = mid(i, 'cy'),
        cz = mid(i, 'cz');
      this.frames.push({ X, Y, Z, cy, cz, angle });
      if (i === sections.length) break;
      const da = this.laws[i].k * this.stationLength,
        phi = angle + da / 2;
      const dy = mid(i + 1, 'cy') - cy,
        dz = mid(i + 1, 'cz') - cz,
        c = Math.cos(phi),
        s = Math.sin(phi);
      if (this.mode === 'bend') {
        X += c * ds - s * dz;
        Y += dy;
        Z += s * ds + c * dz;
      } else {
        X += ds;
        Y += c * dy - s * dz;
        Z += s * dy + c * dz;
      }
      angle += da;
    }
  }
  map(x, y, z, distal = false) {
    const { ds, start, sections } = this.specimen;
    const r = (x - start) / ds,
      t = Math.max(0, Math.min(sections.length, r));
    const i = Math.min(sections.length - 1, Math.floor(t)),
      a = t - i;
    const f = this.frames[i],
      g = this.frames[i + 1];
    const mix = (key) => f[key] + a * (g[key] - f[key]);
    const angle = mix('angle'),
      c = Math.cos(angle),
      s = Math.sin(angle),
      dx = (r - t) * ds,
      dy = y - mix('cy'),
      dz = z - mix('cz');
    const gap = distal ? this.gap : 0;
    if (this.mode === 'twist' && this.specimen.gripAxis) {
      const [ax, ay, az] = this.specimen.gripAxis,
        [ox, oy, oz] = this.specimen.gripOrigin;
      const X = x - ox,
        Y = y - oy,
        Z = z - oz,
        along = (ax * X + ay * Y + az * Z) * (1 - c);
      return [
        ox + c * X + s * (ay * Z - az * Y) + ax * along + gap,
        oy + c * Y + s * (az * X - ax * Z) + ay * along,
        oz + c * Z + s * (ax * Y - ay * X) + az * along,
      ];
    }
    if (this.mode === 'bend')
      return [mix('X') + c * dx - s * dz + gap, mix('Y') + dy, mix('Z') + s * dx + c * dz];
    return [mix('X') + dx + gap, mix('Y') + c * dy - s * dz, mix('Z') + s * dy + c * dz];
  }
}
