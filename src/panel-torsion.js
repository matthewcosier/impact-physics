// Fixed-span, unwrinkled-section approximation. A longitudinal fibre at
// radius r follows a helix. Nominal stress times d(stretch)/dk supplies its
// torque by virtual work. The open-section shear response acts in parallel.
export class RestrainedTorsion {
  constructor(section, material, shear) {
    this.shear = shear;
    this.GJ = shear.stiffness;
    this.backModulus = shear.brittle ? 0 : (this.GJ * shear.hardening) / (this.GJ - shear.hardening);
    this.E = material.E;
    this.tensileYield = shear.brittle ? material.ultimate : material.yieldStress;
    this.compressionYield = material.compressionYield;
    const ultimateStrain = Math.min(
      material.fractureStrain,
      shear.brittle ? material.ultimate / material.E : Infinity,
    );
    this.H = shear.brittle
      ? 0
      : (material.ultimate - this.tensileYield) /
        Math.max(1e-10, material.fractureStrain - material.ultimate / material.E);
    this.fibres = section.fibres.map((f) => ({ ...f, set: 0, accumulated: 0 }));
    this.lastRate = 0;
    this.plasticWork = 0; // N·m of work per metre of grip span
    const shearYieldRate = shear.yield / this.GJ;
    // Use the material strain limit directly, including a perfectly plastic plateau.
    const shearLimitRate = shear.brittle
      ? shear.limit / this.GJ
      : (Math.sqrt(3) * material.fractureStrain) / section.torsionLever;
    this.limitRate = Math.min(
      shearLimitRate,
      Math.sqrt(ultimateStrain * (2 + ultimateStrain)) / section.outerRadius,
    );
    const yieldStrain = this.tensileYield / this.E;
    this.firstYield = this.torque(
      Math.min(
        this.limitRate,
        shearYieldRate,
        Math.sqrt(yieldStrain * (2 + yieldStrain)) / section.outerRadius,
      ),
    );
    this.capacity = this.torque(this.limitRate);
  }
  torque(k, commit = false) {
    const shear = this.shear;
    let plastic = shear.plastic;
    const trial = this.GJ * (k - plastic) - this.backModulus * plastic;
    if (!shear.brittle && Math.abs(trial) > shear.yield)
      plastic += (trial - Math.sign(trial) * shear.yield) / (this.GJ + this.backModulus);
    let moment = this.GJ * (k - plastic);
    if (commit) {
      this.plasticWork += shear.yield * Math.abs(plastic - shear.plastic);
      shear.plastic = plastic;
    }
    for (const f of this.fibres) {
      const rr = f.r2 * k,
        stretch = Math.sqrt(1 + rr * k),
        strain = (rr * k) / (1 + stretch);
      const trialStress = this.E * (strain - f.set),
        direction = Math.sign(trialStress);
      const yieldStress = direction < 0 ? this.compressionYield : this.tensileYield;
      const dp = shear.brittle
        ? 0
        : Math.max(0, Math.abs(trialStress) - yieldStress - this.H * f.accumulated) / (this.E + this.H);
      const set = f.set + direction * dp;
      moment += (f.area * this.E * (strain - set) * rr) / stretch;
      if (commit) {
        f.set = set;
        f.accumulated += dp;
        this.plasticWork += f.area * yieldStress * dp;
      }
    }
    return moment;
  }
  // After yield, compressed fibres can make the straight equilibrium
  // unstable. Follow the stable branch connected to the previous grip pose.
  rate(moment, commit = false) {
    const side = Math.sign(this.lastRate) || Math.sign(moment) || 1;
    const target = side * moment,
      T = (k) => side * this.torque(side * k);
    let lo = 0,
      hi = this.limitRate;
    if (target === T(0) && (!this.lastRate || T(hi * 1e-8) >= target)) hi = 0;
    else if (target < T(0)) {
      let a = 0,
        b = Math.abs(this.lastRate);
      for (let i = 0; i < 40; i++) {
        const x = (2 * a + b) / 3,
          y = (a + 2 * b) / 3;
        if (T(x) < T(y)) b = y;
        else a = x;
      }
      lo = (a + b) / 2;
      if (T(lo) > target) {
        // Reverse twist passes through the straight state. Reverse plastic
        // flow there must enter the history before solving the other branch.
        if (commit && this.lastRate) this.torque(0, true);
        lo = -this.limitRate;
        hi = 0;
      }
    }
    for (let i = 0; i < 45 && lo !== hi; i++) {
      const mid = (lo + hi) / 2;
      if (T(mid) < target) lo = mid;
      else hi = mid;
    }
    const k = (side * (lo + hi)) / 2;
    if (commit) {
      this.torque(k, true);
      this.lastRate = k;
    }
    return k;
  }
}
