import { ContactImpulses } from './contact-impulse.js';
import { GuidedContact } from './guided-contact.js';
import { contactDepth, horizontalExtent } from './impact-shapes.js';
import { dampContacts, contactWork } from './contact-damping.js';
import { WallContact } from './wall-contact.js';
import { FlexDamping } from './flex-damping.js';
import { RigidGroup } from './rigid-group.js';
import { BoxContact } from './box-contact.js';
import { SolidTetrahedron } from './solid.js';
import { MaterialLink, dampMaterialLinks } from './material-link.js';
// SI units. Each node carries a displacement from its reference position and a velocity.
// A skin point or attachment distributes its force with the same weights used to read its motion.
export const FIXED_DT = 1 / 16000;

export function skinPoint(point, field, axis) {
  let result = 0;
  for (let j = 0; j < point.nodes.length; j++) result += point.weights[j] * field[3 * point.nodes[j] + axis];
  return result;
}

function scatter(point, field, x, y, z) {
  for (let j = 0; j < point.nodes.length; j++) {
    const i = 3 * point.nodes[j],
      w = point.weights[j];
    field[i] += x * w;
    field[i + 1] += y * w;
    field[i + 2] += z * w;
  }
}

function inverseMass(point, mass) {
  return point.nodes.reduce((sum, n, j) => sum + point.weights[j] ** 2 / mass[n], 0);
}

export class ImpactWorld {
  constructor(model, striker = null, dt = (striker && model.timestepFor?.(striker)) || FIXED_DT) {
    if (!(dt > 0 && Number.isFinite(dt))) throw new RangeError('Expected a positive finite time step');
    if (striker) {
      if (
        ![striker.x, striker.y, striker.z, striker.speed, striker.damping, striker.pad].every(
          Number.isFinite,
        ) ||
        striker.damping < 0 ||
        striker.pad < 0
      )
        throw new RangeError('Invalid impact object pose or contact parameters');
      const axisX = striker.axisX ?? 1,
        axisY = striker.axisY ?? 0;
      if (!Number.isFinite(axisX + axisY) || Math.abs(axisX * axisX + axisY * axisY - 1) > 1e-10)
        throw new RangeError('The impact guide must have a unit direction');
      const shape = striker.shape ?? 'cylinder';
      if (!['wall', 'block', 'sphere', 'cylinder'].includes(shape))
        throw new RangeError('Unknown impact object');
      const dimensions =
        shape === 'block'
          ? [striker.depth, striker.width, striker.height]
          : shape === 'cylinder'
            ? [striker.radius, striker.height]
            : shape === 'sphere'
              ? [striker.radius]
              : [];
      if (![striker.mass, striker.k, ...dimensions].every((n) => Number.isFinite(n) && n > 0))
        throw new RangeError('Expected positive mass, stiffness and impact object dimensions');
    }
    this.model = model;
    this.dt = dt;
    this.steps = 0;
    for (const n of model.nodes) {
      if (
        n.position?.length !== 3 ||
        !n.position.every(Number.isFinite) ||
        (n.mass !== undefined && (!Number.isFinite(n.mass) || n.mass < 0))
      )
        throw new RangeError('Invalid node reference position or additional mass');
    }
    this.mass = Float64Array.from(model.nodes, (n) => n.mass ?? 0);
    this.solids = (model.solids ?? []).map((spec) => new SolidTetrahedron(spec, model.nodes));
    for (const solid of this.solids) for (const n of solid.indices) this.mass[n] += solid.mass / 4;
    this.materialLinks = (model.materialLinks ?? []).map((spec) => new MaterialLink(spec, model.nodes));
    for (const link of this.materialLinks) {
      this.mass[link.a] += link.mass / 2;
      this.mass[link.b] += link.mass / 2;
    }
    if (!this.mass.every((m) => Number.isFinite(m) && m > 0)) throw new RangeError('Invalid nodal mass');
    this.q = new Float64Array(3 * this.mass.length);
    this.v = Float64Array.from(model.nodes.flatMap((n) => n.velocity ?? [0, 0, 0]));
    if (this.v.length !== this.q.length || !this.v.every(Number.isFinite))
      throw new RangeError('Expected three finite velocities per node');
    this.force = new Float64Array(this.q.length);
    const grouped = new Set();
    this.rigidGroups = (model.rigidGroups ?? []).map((ns) => {
      for (const n of ns) {
        if (!Number.isInteger(n) || !model.nodes[n] || grouped.has(n))
          throw new RangeError('Invalid rigid group membership');
        grouped.add(n);
      }
      return new RigidGroup(ns, model, this.mass);
    });
    this.flexDamping = (model.flexDamping ?? []).map(
      (group) => new FlexDamping(group, model.nodes, this.mass, dt),
    );
    this.freeNodes = Int32Array.from(model.nodes.map((_, n) => n).filter((n) => !grouped.has(n)));
    this.impactor = striker ? { axisX: 1, axisY: 0, ...striker } : null;
    const cells = model.surfaceCells ?? model.wallCells,
      points = model.boundaryPoints ?? model.wallPoints;
    this.boundary =
      striker?.shape === 'wall'
        ? new WallContact(this, cells ?? [], points ?? model.contacts ?? [])
        : striker && cells
          ? new GuidedContact(this, cells, points ?? [])
          : null;
    this.constraint = this.boundary ? new ContactImpulses(this) : null;
    // Wall runs use the lab frame. The cut-plane pickups ride a free
    // backing carriage, and the complete clip approaches the stationary wall.
    this.sled = striker?.shape === 'wall' ? { x: 0, speed: -striker.speed, mass: striker.mass } : null;
    if (this.sled) {
      this.impactor.speed = 0;
      for (let i = 0; i < this.v.length; i += 3) this.v[i] += this.sled.speed;
    }
    const mobile = this.sled ?? this.impactor;
    this.initialEnergy = mobile ? 0.5 * mobile.mass * mobile.speed ** 2 : 0;
    this.initialMomentum = mobile ? mobile.mass * mobile.speed * (mobile.axisX ?? 1) : 0;
    for (let n = 0; n < this.mass.length; n++) {
      const i = 3 * n;
      this.initialEnergy += 0.5 * this.mass[n] * (this.v[i] ** 2 + this.v[i + 1] ** 2 + this.v[i + 2] ** 2);
      this.initialMomentum += this.mass[n] * this.v[i];
    }
    this.beams = (model.beams ?? []).map((b) => {
      const delta = model.nodes[b.b].position.map((p, d) => p - model.nodes[b.a].position[d]);
      const length = Math.sqrt(delta[0] ** 2 + delta[1] ** 2 + delta[2] ** 2);
      if (!(length > 0 && b.k > 0 && b.yieldStrain > 0 && b.damping >= 0))
        throw new RangeError('Invalid beam geometry or material');
      if (
        b.dampingCoefficient !== undefined &&
        !(Number.isFinite(b.dampingCoefficient) && b.dampingCoefficient >= 0)
      )
        throw new RangeError('Invalid beam damping coefficient');
      const reducedMass = 1 / (1 / this.mass[b.a] + 1 / this.mass[b.b]);
      return {
        owner: b.owner,
        k: b.k,
        i: b.a * 3,
        j: b.b * 3,
        dx: delta[0],
        dy: delta[1],
        dz: delta[2],
        length,
        rest: length,
        yieldLength: length * b.yieldStrain,
        compressiveYieldLength: length * (b.compressiveYield ?? b.yieldStrain),
        c: b.dampingCoefficient ?? 2 * b.damping * Math.sqrt(b.k * reducedMass),
      };
    });
    this.joints = (model.joints ?? []).map((j) => {
      if (!(j.k > 0 && j.onset > 0 && j.failure > j.onset && j.damping >= 0))
        throw new RangeError('Invalid joint law');
      if (j.b && j.a.position.some((x, i) => Math.abs(x - j.b.position[i]) > 1e-7))
        throw new RangeError('Joint sites must coincide in the reference pose');
      if (
        j.dampingCoefficient !== undefined &&
        !(Number.isFinite(j.dampingCoefficient) && j.dampingCoefficient >= 0)
      )
        throw new RangeError('Invalid joint damping coefficient');
      const mobility = inverseMass(j.a, this.mass) + (j.b ? inverseMass(j.b, this.mass) : 0);
      return {
        a: j.a,
        b: j.b,
        type: j.type,
        source: j,
        k: j.k,
        onset: j.onset,
        failure: j.failure,
        opening: 0,
        damage: 0,
        impulseX: 0,
        peak: 0,
        c: j.dampingCoefficient ?? 2 * j.damping * Math.sqrt(j.k / mobility),
      };
    });
    this.contacts = (striker ? (model.contacts ?? []) : []).map((p) => ({
      nodes: p.nodes,
      weights: p.weights,
      position: p.position,
      owner: p.owner,
      c:
        2 *
        striker.damping *
        Math.sqrt(striker.k / (inverseMass(p, this.mass) + (this.sled ? 0 : 1 / striker.mass))),
    }));
    this.interfaces = (model.interfaces ?? []).map((s) => {
      const mobility = inverseMass(s.a, this.mass) + inverseMass(s.b, this.mass);
      return {
        a: s.a,
        b: s.b,
        delta: s.a.position.map((v, d) => v - s.b.position[d]),
        radius: s.radius,
        minRadius: s.radius - (s.crushTravel ?? 0),
        yieldTravel: s.yieldTravel ?? 0,
        k: s.k,
        c: 2 * s.damping * Math.sqrt(s.k / mobility),
        impulse: 0,
      };
    });
    this.contactBodies = [...(model.contactBodies ?? []), ...(model.internalContactBodies ?? [])].map(
      (b) => new BoxContact(b, this.mass),
    );
    this.antiRollBar = model.antiRollBar;
    this.roadSupports = (model.roadSupports ?? []).map((s) => ({ ...s, slip: 0 }));
    for (const support of this.roadSupports) {
      if (
        ![support.lateralStiffness ?? 0, support.lateralDamping ?? 0, support.friction ?? 0].every(
          (n) => Number.isFinite(n) && n >= 0,
        ) ||
        !Number.isInteger(support.node) ||
        !model.nodes[support.node] ||
        ![support.stiffness, support.preload].every((n) => Number.isFinite(n) && n > 0) ||
        !Number.isFinite(support.damping) ||
        support.damping < 0
      )
        throw new RangeError('Invalid preloaded road support');
    }
    this.releases = model.createReleases?.() ?? null;
    this.connectionsChanged = false;
    this.contactReference = new Float64Array(this.q.length);
    this.activeContacts = null;
    this.contactGeneration = 0;
    this.contactMargin = 0.012;
    for (const p of this.contacts) p.weightNorm = p.weights.reduce((sum, w) => sum + Math.abs(w), 0);
    this.loss = { plastic: 0, fracture: 0, damping: 0, friction: 0 };
    this.externalWork = 0;
    this.externalImpulseX = 0;
    this.sledReaction = 0;
    this.fixtureImpulse = 0;
    this.contactImpulse = 0;
    this.peakForce = 0;
    this.peakMeanForce = 0;
    this.meanForce = 0;
    this.forceWindowTime = 0.001;
    this.forceHistory = new Float64Array(Math.ceil(this.forceWindowTime / dt));
    this.windowImpulse = 0;
    this.contactForce = 0;
  }

  // Existing cylinder fixtures can still access the guided body here.
  get cylinder() {
    return this.impactor;
  }

  get time() {
    return this.steps * this.dt;
  }

  solveBeams() {
    const { q, v, force, dt } = this;
    for (const b of this.beams) {
      const { i, j } = b;
      const x = b.dx + (q[j] - q[i]),
        y = b.dy + (q[j + 1] - q[i + 1]),
        z = b.dz + (q[j + 2] - q[i + 2]);
      const length = Math.sqrt(x * x + y * y + z * z);
      if (length < 1e-10) continue;
      const strainLength = length - b.rest;
      const yieldLength = strainLength < 0 ? b.compressiveYieldLength : b.yieldLength;
      if (Math.abs(strainLength) > yieldLength) {
        const old = b.rest;
        // Densification stops additional compressive flow once 55% of the reference length is crushed.
        b.rest = Math.max(0.45 * b.length, length - Math.sign(strainLength) * yieldLength);
        this.loss.plastic += Math.abs(b.rest - old) * b.k * yieldLength;
      }
      const rate = ((v[j] - v[i]) * x + (v[j + 1] - v[i + 1]) * y + (v[j + 2] - v[i + 2]) * z) / length;
      const magnitude = (b.k * (length - b.rest) + b.c * rate) / length;
      const fx = magnitude * x,
        fy = magnitude * y,
        fz = magnitude * z;
      force[i] += fx;
      force[i + 1] += fy;
      force[i + 2] += fz;
      force[j] -= fx;
      force[j + 1] -= fy;
      force[j + 2] -= fz;
      this.loss.damping += b.c * rate * rate * dt;
    }
  }

  step(loads = null) {
    const { q, v, force, dt } = this;
    if (loads)
      for (const { point, force: f } of loads) {
        if (f?.length !== 3 || !f.every(Number.isFinite))
          throw new RangeError('Expected a finite world force');
        if (
          !point?.nodes?.length ||
          point.nodes.length !== point.weights?.length ||
          !point.nodes.every((n) => Number.isInteger(n) && n >= 0 && n < this.mass.length) ||
          !point.weights.every(Number.isFinite) ||
          Math.abs(point.weights.reduce((a, b) => a + b, 0) - 1) > 1e-9
        )
          throw new RangeError('Invalid material-point force weights or nodes');
      }
    for (const group of this.flexDamping) this.loss.damping += group.apply(q, v, this.mass);
    force.fill(0);
    this.solveBeams();
    for (const solid of this.solids) solid.evaluate(q, force);
    for (const link of this.materialLinks) {
      link.sample(q);
      link.scatter(force, link.axialForce);
    }
    this.sledReaction = 0;
    for (const j of this.joints) this.applyJoint(j);
    if (this.connectionsChanged) {
      this.releases?.reconnect(this.joints, q);
      this.connectionsChanged = false;
    }
    this.releases?.apply(this);
    this.applyRoadSupports();
    this.applyBearings();
    const generation = this.refreshContactBounds();
    this.applyBoxContacts(false, generation);
    this.applyContacts(false, generation);
    // A force is expressed in world newtons at any weighted material point.
    // Equal/opposite forces at separated points apply a pure couple.
    if (loads)
      for (const { point, force: f } of loads) {
        scatter(point, force, ...f);
        this.externalImpulseX += f[0] * dt;
        for (let d = 0; d < 3; d++) this.externalWork += 0.5 * dt * f[d] * skinPoint(point, v, d);
      }
    if (this.materialLinks.length) this.loss.damping += dampMaterialLinks(this);
    if (this.contactBodies.length) dampContacts(this);
    this.constraint?.begin(this);
    if (this.sled) {
      this.sled.speed += (dt * this.sledReaction) / this.sled.mass;
      this.sled.x += dt * this.sled.speed;
    } else if (this.impactor) {
      this.impactor.speed -= (dt * this.contactForce) / this.impactor.mass;
      this.impactor.x += dt * this.impactor.speed * this.impactor.axisX;
      this.impactor.y += dt * this.impactor.speed * this.impactor.axisY;
    }
    for (const group of this.rigidGroups) group.advance(this);
    for (let j = 0; j < this.freeNodes.length; j++) {
      const n = this.freeNodes[j],
        factor = dt / this.mass[n];
      for (let i = 3 * n; i < 3 * n + 3; i++) {
        v[i] += force[i] * factor;
        q[i] += dt * v[i];
      }
    }
    if (this.boundary) {
      this.boundary.solve(this, this.constraint);
      this.constraint.finish(this);
    }
    if (this.contactBodies.length) contactWork(this);
    if (loads)
      for (const { point, force: f } of loads)
        for (let axis = 0; axis < 3; axis++)
          this.externalWork += 0.5 * dt * f[axis] * skinPoint(point, v, axis);
    for (const link of this.materialLinks) link.sample(q);
    for (const solid of this.solids)
      if (solid.metal) {
        const previousWork = solid.plasticWork;
        solid.evaluate(q, null, true);
        this.loss.plastic += solid.plasticWork - previousWork;
      }
    if (this.impactor) {
      const count = this.forceHistory.length,
        j = this.steps % count;
      this.windowImpulse += dt * (this.contactForce - this.forceHistory[j]);
      this.forceHistory[j] = this.contactForce;
      // Include the correct fraction of the oldest sample even for a step
      // which does not divide one millisecond.
      const excess = count * dt - this.forceWindowTime;
      this.meanForce =
        (this.windowImpulse - excess * this.forceHistory[(j + 1) % count]) / this.forceWindowTime;
      this.peakMeanForce = Math.max(this.peakMeanForce, Math.abs(this.meanForce));
    }
    this.steps++;
  }

  // Incremental vertical force about the laden equilibrium. On lift-off
  // the tyre normal force is zero and the balanced static load remains.
  applyRoadSupports(measureOnly = false) {
    let potential = 0;
    if (this.antiRollBar) {
      const { towers, hubs, rate } = this.antiRollBar;
      const a = 3 * towers[0] + 2,
        b = 3 * hubs[0] + 2,
        c = 3 * towers[1] + 2,
        d = 3 * hubs[1] + 2;
      const diff = this.q[a] - this.q[b] - (this.q[c] - this.q[d]),
        load = rate * diff;
      potential += 0.5 * load * diff;
      if (!measureOnly) {
        this.force[a] -= load;
        this.force[b] += load;
        this.force[c] += load;
        this.force[d] -= load;
      }
    }
    for (const support of this.roadSupports) {
      const {
        node,
        stiffness: k,
        damping,
        preload,
        lateralStiffness = 0,
        lateralDamping = 0,
        friction = 0,
      } = support;
      const i = node * 3 + 2,
        z = this.q[i],
        lift = preload / k,
        vz = this.v[i];
      potential += z < lift ? 0.5 * k * z * z : preload * (z - 0.5 * lift);
      const normal = Math.max(0, preload - k * z);
      const reaction = normal > 0 ? Math.max(0, normal - damping * vz) : 0;
      let extension = this.q[i - 1] - support.slip;
      if (!measureOnly) {
        this.force[i] += reaction - preload;
        this.loss.damping += this.dt * (normal - reaction) * vz;
        if (lateralStiffness) {
          // The lateral contact patch sticks up to the available tyre
          // friction. It slides beyond that force and releases on lift-off.
          const cap = friction * reaction,
            limit = cap / lateralStiffness;
          const retained = Math.max(-limit, Math.min(limit, extension));
          this.loss.damping += 0.5 * lateralStiffness * (extension * extension - retained * retained);
          support.slip += extension - retained;
          extension = retained;
          const elasticForce = lateralStiffness * retained,
            vy = this.v[i - 1];
          const viscous = Math.max(-cap - elasticForce, Math.min(cap - elasticForce, lateralDamping * vy));
          this.force[i - 1] -= elasticForce + viscous;
          this.loss.damping += this.dt * viscous * vy;
        }
      }
      potential += 0.5 * lateralStiffness * extension * extension;
    }
    return potential;
  }

  applyJoint(j) {
    if (j.damage === 1) return;
    const { q, v, dt } = this;
    // Both weighted points start at the same world-space site. A missing B is a rear fixture pickup.
    const x = skinPoint(j.a, q, 0) - (j.b ? skinPoint(j.b, q, 0) : (this.sled?.x ?? 0));
    const y = skinPoint(j.a, q, 1) - (j.b ? skinPoint(j.b, q, 1) : 0);
    const z = skinPoint(j.a, q, 2) - (j.b ? skinPoint(j.b, q, 2) : 0);
    const r = Math.sqrt(x * x + y * y + z * z);
    if (r > j.opening) {
      const previous = j.opening;
      j.opening = r;
      if (r > j.onset) {
        j.damage = r >= j.failure ? 1 : ((j.failure / r) * (r - j.onset)) / (j.failure - j.onset);
        if (j.damage === 1) this.connectionsChanged = true;
        const growth = Math.min(r, j.failure) - Math.max(previous, j.onset);
        this.loss.fracture += (0.5 * j.k * j.onset * j.failure * growth) / (j.failure - j.onset);
      }
    }
    const vx = skinPoint(j.a, v, 0) - (j.b ? skinPoint(j.b, v, 0) : (this.sled?.speed ?? 0));
    const vy = skinPoint(j.a, v, 1) - (j.b ? skinPoint(j.b, v, 1) : 0);
    const vz = skinPoint(j.a, v, 2) - (j.b ? skinPoint(j.b, v, 2) : 0);
    const intact = 1 - j.damage;
    const fx = intact * (j.k * x + j.c * vx),
      fy = intact * (j.k * y + j.c * vy),
      fz = intact * (j.k * z + j.c * vz);
    scatter(j.a, this.force, -fx, -fy, -fz);
    if (j.b) scatter(j.b, this.force, fx, fy, fz);
    else {
      this.fixtureImpulse += fx * dt;
      this.sledReaction += fx;
    }
    j.impulseX += fx * dt;
    j.peak = Math.max(j.peak, Math.sqrt(fx * fx + fy * fy + fz * fz));
    this.loss.damping += intact * j.c * (vx * vx + vy * vy + vz * vz) * dt;
  }

  applyBearings(measureOnly = false) {
    let elastic = 0;
    for (const pair of this.interfaces) {
      const { a, b, delta, radius } = pair;
      const x = delta[0] + skinPoint(a, this.q, 0) - skinPoint(b, this.q, 0);
      if (Math.abs(x) >= radius) continue;
      const y = delta[1] + skinPoint(a, this.q, 1) - skinPoint(b, this.q, 1);
      const z = delta[2] + skinPoint(a, this.q, 2) - skinPoint(b, this.q, 2);
      const len = Math.sqrt(x * x + y * y + z * z);
      if (len >= radius || len < 1e-10) continue;
      if (!measureOnly && pair.yieldTravel && radius - len > pair.yieldTravel) {
        pair.radius = Math.max(pair.minRadius, len + pair.yieldTravel);
        this.loss.plastic += pair.k * pair.yieldTravel * (radius - pair.radius);
      }
      const overlap = pair.radius - len;
      elastic += 0.5 * pair.k * overlap * overlap;
      if (measureOnly) continue;
      const rate = Math.min(
        0,
        ((skinPoint(a, this.v, 0) - skinPoint(b, this.v, 0)) * x +
          (skinPoint(a, this.v, 1) - skinPoint(b, this.v, 1)) * y +
          (skinPoint(a, this.v, 2) - skinPoint(b, this.v, 2)) * z) /
          len,
      );
      const load = pair.k * overlap - pair.c * rate,
        scale = load / len;
      scatter(a, this.force, x * scale, y * scale, z * scale);
      scatter(b, this.force, -x * scale, -y * scale, -z * scale);
      pair.impulse += load * this.dt;
      this.loss.damping += pair.c * rate * rate * this.dt;
    }
    return elastic;
  }

  refreshContactBounds() {
    if (!this.contactBodies.length && !this.contacts.length && !this.boundary) return 0;
    const q = this.q,
      reference = this.contactReference,
      limit = this.contactMargin,
      objectX = this.impactor?.x ?? 0;
    let rebuild =
      !this.contactGeneration ||
      Math.abs(objectX - this.contactX) >= limit ||
      Math.abs((this.impactor?.y ?? 0) - this.contactY) >= limit;
    if (!rebuild)
      for (let i = 0; i < reference.length; i++) {
        if (Math.abs(q[i] - reference[i]) >= limit) {
          rebuild = true;
          break;
        }
      }
    if (rebuild) {
      reference.set(q);
      this.contactX = objectX;
      this.contactY = this.impactor?.y ?? 0;
      this.contactGeneration++;
    }
    return this.contactGeneration;
  }
  applyBoxContacts(measureOnly = false, generation = this.refreshContactBounds()) {
    if (!this.contactBodies.length) return 0;
    const rebuild = generation !== this.boxGeneration;
    this.boxGeneration = generation;
    let energy = 0;
    for (const body of this.contactBodies)
      energy += body.apply(this, measureOnly, rebuild, this.contactMargin);
    return energy;
  }
  contactCandidates(generation = this.refreshContactBounds()) {
    if (generation !== this.impactorGeneration) {
      this.impactorGeneration = generation;
      const c = this.impactor,
        limit = this.contactMargin,
        q = this.q;
      this.activeContacts = this.contacts.filter((p) => {
        const margin = c.pad + limit * (1 + p.weightNorm);
        const x = p.position[0] + skinPoint(p, q, 0) - c.x;
        if (c.shape === 'wall') return x <= margin;
        return (
          Math.abs(x) <= horizontalExtent(c, 0) + margin &&
          Math.abs(p.position[1] + skinPoint(p, q, 1) - c.y) <= horizontalExtent(c, 1) + margin &&
          Math.abs(p.position[2] + skinPoint(p, q, 2) - c.z) <=
            (c.shape === 'sphere' ? c.radius : c.height * 0.5) + margin
        );
      });
    }
    return this.activeContacts;
  }

  applyContacts(measureOnly = false, generation) {
    if (!this.impactor) return 0;
    const { impactor: c, q, v, dt } = this,
      normal = [0, 0, 0];
    let thrust = 0,
      elastic = 0;
    for (const p of measureOnly ? this.contacts : this.contactCandidates(generation)) {
      const x = p.position[0] + skinPoint(p, q, 0) - c.x;
      if (c.shape === 'wall' ? x >= c.pad : Math.abs(x) > horizontalExtent(c, 0) + c.pad) continue;
      const y = p.position[1] + skinPoint(p, q, 1) - c.y,
        z = p.position[2] + skinPoint(p, q, 2) - c.z;
      const depth = contactDepth(c, x, y, z, normal);
      if (depth <= 0) continue;
      const [nx, ny, nz] = normal;
      elastic += 0.5 * c.k * depth * depth;
      if (measureOnly) continue;
      const rate = Math.min(
        0,
        (skinPoint(p, v, 0) - c.speed * c.axisX) * nx +
          (skinPoint(p, v, 1) - c.speed * c.axisY) * ny +
          skinPoint(p, v, 2) * nz,
      );
      const load = c.k * depth - p.c * rate;
      scatter(p, this.force, nx * load, ny * load, nz * load);
      const axial = (nx * c.axisX + ny * c.axisY) * load;
      thrust += axial;
      this.externalImpulseX += (nx * load - axial * c.axisX) * dt;
      this.loss.damping += p.c * rate * rate * dt;
    }
    if (measureOnly) return elastic;
    this.contactForce = thrust;
    this.contactImpulse += thrust * dt;
    this.peakForce = Math.max(this.peakForce, Math.abs(thrust));
  }

  telemetry() {
    this.syncBeamState?.();
    const mobile = this.sled ?? this.impactor;
    let kinetic = mobile ? 0.5 * mobile.mass * mobile.speed ** 2 : 0;
    let px = mobile ? mobile.mass * mobile.speed * (mobile.axisX ?? 1) : 0,
      elastic = this.applyContacts(true) + this.applyBearings(true),
      plasticSet = 0;
    for (let n = 0; n < this.mass.length; n++) {
      const i = n * 3;
      kinetic += 0.5 * this.mass[n] * (this.v[i] ** 2 + this.v[i + 1] ** 2 + this.v[i + 2] ** 2);
      px += this.mass[n] * this.v[i];
    }
    for (const b of this.beams) {
      const { i, j } = b,
        q = this.q;
      const length = Math.hypot(b.dx + q[j] - q[i], b.dy + q[j + 1] - q[i + 1], b.dz + q[j + 2] - q[i + 2]);
      elastic += 0.5 * b.k * (length - b.rest) ** 2;
      plasticSet = Math.max(plasticSet, Math.abs(b.rest - b.length));
    }
    elastic += this.applyBoxContacts(true) + this.applyRoadSupports(true);
    for (const solid of this.solids) elastic += solid.evaluate(this.q);
    for (const joint of this.joints) {
      if (joint.damage === 1) continue;
      let separation2 = 0;
      for (let d = 0; d < 3; d++) {
        const delta =
          skinPoint(joint.a, this.q, d) -
          (joint.b ? skinPoint(joint.b, this.q, d) : d === 0 ? (this.sled?.x ?? 0) : 0);
        separation2 += delta * delta;
      }
      elastic += 0.5 * joint.k * (1 - joint.damage) * separation2;
    }
    const materialWork = this.materialLinks.reduce((sum, link) => sum + link.materialWork, 0);
    const groundLoss = this.releases?.groundLoss ?? 0;
    const gravityAndGround = this.releases?.apply(this, true) ?? 0;
    const dissipated =
      this.loss.plastic + this.loss.fracture + this.loss.damping + this.loss.friction + groundLoss;
    return {
      time: this.time,
      speed: mobile?.speed ?? 0,
      force: this.contactForce,
      peakForce: this.peakForce,
      peakMeanForce: this.peakMeanForce,
      impulse: this.contactImpulse,
      fixtureImpulse: this.fixtureImpulse,
      plasticSet,
      kinetic,
      elastic,
      dissipated,
      ...this.loss,
      releasedParts: this.releases?.released.size ?? 0,
      groundLoss,
      gravityAndGround,
      materialWork,
      externalWork: this.externalWork,
      energyError:
        kinetic +
        elastic +
        materialWork +
        gravityAndGround +
        dissipated -
        this.initialEnergy -
        this.externalWork,
      momentumError:
        px +
        (this.sled ? -this.contactImpulse : this.fixtureImpulse) +
        (this.releases?.groundImpulseX ?? 0) -
        this.initialMomentum -
        this.externalImpulseX,
      damaged: this.joints.filter((j) => j.damage > 0).length,
      failed: this.joints.filter((j) => j.damage === 1).length,
    };
  }
}
