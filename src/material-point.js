// Engineering (nominal) stress and strain, with tension positive. Each
// specimen holds its own history. All tabulated data and outputs use SI.
function interpolate(curve, x) {
  let i = 1;
  while (i < curve.length - 1 && curve[i][0] < x) i++;
  const [x0, y0] = curve[i - 1],
    [x1, y1] = curve[i];
  return y0 + ((x - x0) / (x1 - x0)) * (y1 - y0);
}
export class MaterialPoint {
  constructor(material, across = false) {
    if (!['foam', 'membrane', 'ductile', 'elastomer', 'brittle', 'laminate'].includes(material?.law))
      throw new RangeError('Unsupported material law');
    const positive = (key) => Number.isFinite(material[key]) && material[key] > 0;
    if (!positive('density') || !positive('E')) throw new RangeError('Invalid density or modulus');
    if (material.law === 'foam') {
      const curve = material.compressionCurve;
      if (
        !Array.isArray(curve) ||
        curve.length < 2 ||
        curve[0][0] !== 0 ||
        curve[0][1] !== 0 ||
        !curve.every(
          ([x, y], i) =>
            Number.isFinite(x) &&
            Number.isFinite(y) &&
            x < 1 &&
            (i === 0 || (x > curve[i - 1][0] && y >= curve[i - 1][1])),
        ) ||
        !(
          material.crushFraction >= 0 &&
          material.crushFraction <= 1 &&
          material.crushOnset >= 0 &&
          material.crushOnset < 1 &&
          positive('recoverySeconds')
        )
      )
        throw new RangeError('Invalid compression or recovery data');
    } else {
      if (!positive('fractureStrain') || !positive('ultimate'))
        throw new RangeError('Invalid failure envelope');
      if (
        material.law === 'ductile' &&
        (!positive('yieldStress') ||
          !positive('compressionYield') ||
          material.ultimate < material.yieldStress ||
          material.ultimate >= material.E * material.fractureStrain)
      )
        throw new RangeError('Invalid ductile hardening data');
      if (['laminate', 'brittle'].includes(material.law) && !positive('compressionYield'))
        throw new RangeError('Invalid compressive strength');
      if (
        material.law === 'membrane' &&
        (!positive('toeStrain') || !(material.compressionRatio >= 0 && material.compressionRatio <= 1))
      )
        throw new RangeError('Invalid membrane take-up or slack data');
      if (
        material.law === 'elastomer' &&
        (!positive('G') || !(material.compressionLimit > 0 && material.compressionLimit < 1))
      )
        throw new RangeError('Invalid rubber compression data');
      if (material.transverseRatio !== undefined && !positive('transverseRatio'))
        throw new RangeError('Invalid transverse ratio');
    }
    this.material = material;
    this.across = across;
    this.maxCompression = this.set = this.strain = this.stress = this.plasticTravel = 0;

    this.failed = false;
  }
  load(e) {
    const m = this.material;
    if (!Number.isFinite(e) || e <= -1) throw new RangeError('Invalid engineering strain');
    if (m.law === 'foam' && (e > 0 || -e > m.compressionCurve.at(-1)[0]))
      throw new RangeError('Outside tabulated foam compression');
    if (m.law === 'elastomer' && e < -m.compressionLimit)
      throw new RangeError('Outside rubber compression range');
    const previousStrain = this.strain;
    this.strain = e;
    if (this.failed) return (this.stress = 0);
    const ratio = this.across ? (m.transverseRatio ?? 1) : 1,
      E = m.E * ratio;
    if (m.law === 'foam') {
      const c = -e,
        curve = m.compressionCurve;
      const advancing = c > this.maxCompression;
      this.maxCompression = Math.max(this.maxCompression, c);
      const peak = this.maxCompression;
      const crush = m.crushFraction * Math.max(0, peak - m.crushOnset);
      this.set = -crush;
      if (advancing || (c === peak && previousStrain === e && this.stress < 0))
        return (this.stress = -interpolate(curve, c));
      const rest = crush;
      if (c <= rest) return (this.stress = 0);
      const unload = (c - rest) / (peak - rest);
      return (this.stress = -Math.min(interpolate(curve, c), interpolate(curve, peak) * unload ** 1.65));
    }
    if (e >= m.fractureStrain) return this.rupture();
    if (m.law === 'elastomer') {
      const stretch = 1 + e,
        stress = m.G * (stretch - 1 / stretch ** 2);
      return stress >= m.ultimate ? this.rupture() : (this.stress = stress);
    }
    if (m.law === 'membrane') {
      const stress =
        e <= 0 ? E * m.compressionRatio * e : E * (e + m.toeStrain * Math.expm1(-e / m.toeStrain));
      return stress >= m.ultimate * ratio ? this.rupture() : (this.stress = stress);
    }
    const trial = (e - this.set) * E;
    if (m.law === 'brittle' || m.law === 'laminate') {
      if (trial >= m.ultimate * ratio || trial <= -m.compressionYield * ratio) return this.rupture();
      return (this.stress = trial);
    }
    const limit = (trial >= 0 ? m.yieldStress : m.compressionYield) * ratio;
    const tangent = Math.max(0, (m.ultimate * ratio - limit) / Math.max(1e-8, m.fractureStrain - limit / E));
    const H = (E * tangent) / (E - tangent);
    const over = Math.abs(trial) - limit - H * this.plasticTravel;
    if (over > 0) {
      const increment = over / (E + H);
      this.plasticTravel += increment;
      this.set += Math.sign(trial) * increment;
    }
    const stress = (e - this.set) * E;
    return m.ultimate > m.yieldStress && stress >= m.ultimate * ratio
      ? this.rupture()
      : (this.stress = stress);
  }
  rupture() {
    this.failed = true;
    return (this.stress = 0);
  }
  // Parallel viscosity (Pa·s) for a free-recovery segment. The secant
  // viscosity gives de/dt = -(e-set)/tau at zero total traction. The envelope
  // supplies the quasistatic stress, so rate and equilibrium can be tested apart.
  recoveryViscosity() {
    if (!this.material.recoverySeconds || this.failed) return 0;
    const elasticStrain = this.strain - this.set;
    if (Math.abs(elasticStrain) < 1e-12) return 0;
    return Math.max(0, (this.material.recoverySeconds * this.stress) / elasticStrain);
  }
  release(seconds) {
    if (!(Number.isFinite(seconds) && seconds >= 0)) throw new RangeError('Invalid free recovery time');
    if (this.failed) throw new RangeError('Reset ruptured specimen before testing recovery');
    const tau = this.material.recoverySeconds;
    this.strain = tau ? this.set + (this.strain - this.set) * Math.exp(-seconds / tau) : this.set;

    this.stress = 0;
    return this.strain;
  }
}
