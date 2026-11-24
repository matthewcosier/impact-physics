// Isotropic J2 radial return with linear isotropic hardening. Strain and
// stress use symmetric row-major tensors (off-diagonals are tensor shear).
export class MetalPoint {
  constructor(sheet) {
    if (sheet.solidLaw !== 'J2' || sheet.yieldStress !== sheet.compressionYield)
      throw new RangeError('Expected a pressure-independent metal sheet');
    const { E, G, ultimate, yieldStress, fractureStrain } = sheet;
    this.G = G;
    this.K = (E * G) / (3 * (3 * G - E));
    this.yieldStress = yieldStress;
    this.maxPlastic = fractureStrain - ultimate / E;
    this.H = (ultimate - yieldStress) / this.maxPlastic;
    if (
      ![G, this.K, this.maxPlastic, yieldStress].every((x) => Number.isFinite(x) && x > 0) ||
      !Number.isFinite(this.H) ||
      this.H < 0
    )
      throw new RangeError('Invalid metal moduli or hardening data');
    this.plastic = Array(9).fill(0);
    this.accumulated = 0;
  }
  response(strain, commit = false) {
    if (
      strain.length !== 9 ||
      !strain.every(Number.isFinite) ||
      strain[1] !== strain[3] ||
      strain[2] !== strain[6] ||
      strain[5] !== strain[7]
    )
      throw new RangeError('Expected a finite symmetric strain tensor');
    const trial = strain.map((v, i) => v - this.plastic[i]);
    const trace = trial[0] + trial[4] + trial[8];
    const s = trial.map((e, i) => 2 * this.G * (e - (i % 4 === 0 ? trace / 3 : 0)));
    const equivalent = Math.sqrt(1.5 * s.reduce((sum, v) => sum + v * v, 0));
    const over = equivalent - this.yieldStress - this.H * this.accumulated;
    const dp = Math.max(0, over) / (3 * this.G + this.H);
    const accumulated = this.accumulated + dp;
    if (accumulated > this.maxPlastic) throw new RangeError('Outside metal hardening range');
    const ratio = dp ? (3 * this.G * dp) / equivalent : 0;
    const stress = s.map((v, i) => (1 - ratio) * v + (i % 4 === 0 ? this.K * trace : 0));
    let energyDensity = 0.5 * this.H * accumulated ** 2;
    for (let i = 0; i < 9; i++) {
      const increment = dp ? (1.5 * dp * s[i]) / equivalent : 0;
      energyDensity += 0.5 * stress[i] * (trial[i] - increment);
      if (commit) this.plastic[i] += increment;
    }
    if (commit) this.accumulated = accumulated;
    return { stress, energyDensity, accumulated, plasticWorkDensity: this.yieldStress * accumulated };
  }
}
