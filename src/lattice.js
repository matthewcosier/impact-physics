const offsets = Array.from({ length: 8 }, (_, n) => [n & 1, (n >> 1) & 1, (n >> 2) & 1]);
const dist2 = (a, b) => (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2;

// Append an occupied hexahedral lattice. Volume is the amount of load-bearing
// material (m³); mass is the full moving mass (kg), including supported fittings.
export function buildLattice(
  model,
  { id, cage: grid, mass, volumeOfMaterial, material, sectionFactor = 1, rigid = false },
) {
  if (
    ![
      mass,
      volumeOfMaterial,
      sectionFactor,
      material.E,
      material.G,
      material.yieldStress,
      material.compressionYield,
    ].every((n) => Number.isFinite(n) && n > 0) ||
    !Number.isFinite(material.damping) ||
    material.damping < 0 ||
    !grid.occupied?.length ||
    grid.cell?.length !== 3 ||
    !grid.cell.every((n) => Number.isFinite(n) && n > 0)
  )
    throw new RangeError('Invalid lattice mass, material or grid');
  if (
    grid.origin?.length !== 3 ||
    !grid.origin.every(Number.isFinite) ||
    grid.counts?.length !== 3 ||
    !grid.counts.every((n) => Number.isInteger(n) && n > 0) ||
    grid.occupied.some(
      (c) => c.length !== 3 || c.some((v, d) => !Number.isInteger(v) || v < 0 || v >= grid.counts[d]),
    ) ||
    new Set(grid.occupied.map((c) => c.join())).size !== grid.occupied.length
  )
    throw new RangeError('Invalid occupied lattice cells');
  const nodes = new Map(),
    cells = new Map(),
    beams = new Map();
  const volume = volumeOfMaterial / grid.occupied.length;
  for (const coord of grid.occupied) {
    const corners = offsets.map((offset) => {
      const ijk = offset.map((v, d) => v + coord[d]),
        key = ijk.join();
      if (!nodes.has(key)) {
        nodes.set(key, model.nodes.length);
        model.nodes.push({
          position: ijk.map((v, d) => grid.origin[d] + v * grid.cell[d]),
          mass: 0,
          owner: id,
        });
      }
      const index = nodes.get(key);
      model.nodes[index].mass += mass / (8 * grid.occupied.length);
      return index;
    });
    cells.set(coord.join(), corners);
    // The 28 links in each occupied hexahedron include edge, face and body diagonals.
    for (let a = 0; a < 8; a++)
      for (let b = a + 1; b < 8; b++) {
        const i = corners[a],
          j = corners[b],
          key = [Math.min(i, j), Math.max(i, j)].join();
        if (!beams.has(key))
          beams.set(key, {
            a: i,
            b: j,
            k: 0,
            yieldStrain: material.yieldStress / material.E,
            compressiveYield: material.compressionYield / material.E,
            damping: material.damping,
            owner: id,
          });
        const l2 = dist2(model.nodes[i].position, model.nodes[j].position);
        const diagonal = offsets[a].filter((v, d) => v !== offsets[b][d]).length > 1;
        const modulus = diagonal ? 2.8 * material.G : material.E;
        beams.get(key).k += (modulus * sectionFactor * volume) / (7 * l2);
      }
  }
  const nodeIDs = [...nodes.values()];
  if (rigid) (model.rigidGroups ??= []).push(nodeIDs);
  else model.beams.push(...beams.values());
  return { id, cage: grid, cells, nodeIDs, mass };
}

// Use the surface witness to select a cell. A bonded site may extrapolate a short
// distance from that cell to put the two physical attachment points at a common location.
export function bind(part, position, witness = position) {
  const { origin, cell, counts } = part.cage;
  const coord = witness.map((v, d) =>
    Math.max(0, Math.min(counts[d] - 1, Math.floor((v - origin[d]) / cell[d]))),
  );
  const nodes = part.cells.get(coord.join());
  if (!nodes) throw new Error('Surface witness outside occupied cells for ' + part.id);
  const uvw = position.map((v, d) => (v - origin[d]) / cell[d] - coord[d]);
  const weights = offsets.map((offset) => offset.reduce((w, v, d) => w * (v ? uvw[d] : 1 - uvw[d]), 1));
  return { nodes, weights, position, uvw, owner: part.id };
}
