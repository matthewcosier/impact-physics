// SI units. A guide constrains a finite mass to a straight line in XY.
export function createImpactor({
  shape,
  center,
  speed,
  mass,
  angle = 0,
  stiffness,
  damping,
  pad = 0,
  radius,
  height,
  width,
  depth,
}) {
  if (!Array.isArray(center) || center.length !== 3 || !center.every(Number.isFinite))
    throw new RangeError('Expected a finite three-dimensional impactor centre');
  if (!Number.isFinite(angle) || Math.abs(angle) > 45 || (shape === 'wall' && angle !== 0))
    throw new RangeError('Expected a guide angle within 45 degrees or a head-on wall');
  const yaw = (angle * Math.PI) / 180;
  return {
    shape,
    x: center[0],
    y: center[1],
    z: center[2],
    speed,
    mass,
    axisX: Math.cos(yaw),
    axisY: Math.sin(yaw),
    yaw,
    k: stiffness,
    damping,
    pad,
    radius,
    height,
    width,
    depth,
  };
}

// Return engagement depth and write the outward world normal. For a box,
// the Euclidean outside distance handles rounded contact pads at edges and
// corners. Inside, the nearest face supplies the separating direction.
export function contactDepth(c, x, y, z, n) {
  if (c.shape === 'wall') {
    n[0] = 1;
    n[1] = 0;
    n[2] = 0;
    return c.pad - x;
  }
  if (c.shape === 'sphere') {
    const distance = Math.hypot(x, y, z);
    n[0] = distance ? x / distance : 1;
    n[1] = distance ? y / distance : 0;
    n[2] = distance ? z / distance : 0;
    return c.pad + c.radius - distance;
  }
  if (c.shape === 'block') {
    const cs = c.axisX ?? 1,
      sn = c.axisY ?? 0,
      localX = cs * x + sn * y;
    y = cs * y - sn * x;
    x = localX;
    const dx = Math.abs(x) - c.depth / 2,
      dy = Math.abs(y) - c.width / 2,
      dz = Math.abs(z) - c.height / 2;
    if (Math.max(dx, dy, dz) >= c.pad) return 0;
    const ox = Math.max(dx, 0),
      oy = Math.max(dy, 0),
      oz = Math.max(dz, 0);
    const distance = Math.sqrt(ox * ox + oy * oy + oz * oz);
    if (distance > 0) {
      n[0] = (Math.sign(x) * ox) / distance;
      n[1] = (Math.sign(y) * oy) / distance;
      n[2] = (Math.sign(z) * oz) / distance;
      rotateNormal(n, cs, sn);
      return c.pad - distance;
    }
    const nearest = Math.max(dx, dy, dz);
    const axis = dx === nearest ? 0 : dy === nearest ? 1 : 2;
    n[0] = axis === 0 ? Math.sign(x || 1) : 0;
    n[1] = axis === 1 ? Math.sign(y || 1) : 0;
    n[2] = axis === 2 ? Math.sign(z || 1) : 0;
    rotateNormal(n, cs, sn);
    return c.pad - nearest;
  }
  const radius = Math.sqrt(x * x + y * y),
    side = radius - c.radius,
    cap = Math.abs(z) - c.height / 2;
  if (side >= c.pad || cap >= c.pad) return 0;
  const os = Math.max(0, side),
    oc = Math.max(0, cap),
    corner = Math.sqrt(os * os + oc * oc);
  const radial = corner > 0 ? os / corner : Number(side >= cap);
  n[0] = radial * (radius > 1e-12 ? x / radius : 1);
  n[1] = radial * (radius > 1e-12 ? y / radius : 0);
  n[2] = Math.sign(z || 1) * (corner > 0 ? oc / corner : 1 - radial);
  return c.pad - (corner + Math.min(0, Math.max(side, cap)));
}

function rotateNormal(n, cs, sn) {
  const x = n[0];
  n[0] = cs * x - sn * n[1];
  n[1] = sn * x + cs * n[1];
}
export function horizontalExtent(c, axis) {
  if (c.shape !== 'block') return c.radius;
  const cs = Math.abs(c.axisX ?? 1),
    sn = Math.abs(c.axisY ?? 0);
  return 0.5 * (axis === 0 ? cs * c.depth + sn * c.width : sn * c.depth + cs * c.width);
}
