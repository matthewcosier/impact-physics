// The module has no imports. This optional hot loop shares state with the JS
// solver; JS remains the reference path and runs if WebAssembly is unavailable.
export function useBeamKernel(world, module) {
  if (!module) return;
  const { exports: wasm } = new WebAssembly.Instance(module);
  const n = world.q.length,
    count = world.beams.length,
    stride = 11;
  const ptr = wasm.arena(3 * n + stride * count + 2);
  const buf = wasm.memory.buffer;
  const data = (offset, length) => new Float64Array(buf, ptr + offset * 8, length);
  for (const [i, key] of ['q', 'v', 'force'].entries()) {
    const field = data(i * n, n);
    field.set(world[key]);
    world[key] = field;
  }
  const records = data(3 * n, stride * count),
    losses = data(3 * n + stride * count, 2);
  world.beams.forEach((b, j) => {
    const at = j * stride;
    records.set(
      [b.i, b.j, b.dx, b.dy, b.dz, b.length, b.rest, b.yieldLength, b.compressiveYieldLength, b.k, b.c],
      at,
    );
  });
  world.solveBeams = () => {
    wasm.beams(ptr, n, count, world.dt);
    world.loss.plastic += losses[0];
    world.loss.damping += losses[1];
  };
  world.syncBeamState = () => {
    for (let j = 0; j < count; j++) world.beams[j].rest = records[j * stride + 6];
  };
}
