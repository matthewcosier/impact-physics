// Compiled arithmetic for the same axial return mapping as physics.js.
// The caller owns one contiguous arena per simulation; no per-step allocation.
#[no_mangle]
pub extern "C" fn arena(words: usize) -> *mut f64 {
    let mut data = vec![0.0f64; words].into_boxed_slice();
    let ptr = data.as_mut_ptr();
    std::mem::forget(data);
    ptr
}

#[no_mangle]
pub unsafe extern "C" fn beams(base: *mut f64, dofs: usize, count: usize, dt: f64) {
    let q = std::slice::from_raw_parts(base, dofs);
    let v = std::slice::from_raw_parts(base.add(dofs), dofs);
    let f = std::slice::from_raw_parts_mut(base.add(2*dofs), dofs);
    let records = std::slice::from_raw_parts_mut(base.add(3*dofs), count*11);
    let losses = base.add(3*dofs+count*11);
    let mut plastic=0.0;
    let mut damping=0.0;
    for b in records.chunks_exact_mut(11) {
        let i=b[0] as usize;
        let j=b[1] as usize;
        let x=b[2]+(q[j]-q[i]); let y=b[3]+(q[j+1]-q[i+1]); let z=b[4]+(q[j+2]-q[i+2]);
        let len=(x*x+y*y+z*z).sqrt();
        if len<1.0e-10 {continue}
        let old=b[6];
        let elong=len-old;
        let limit=if elong<0.0 {b[8]} else {b[7]};
        if elong.abs()>limit {
            b[6]=(0.45*b[5]).max(len-limit.copysign(elong));
            plastic+=(b[6]-old).abs()*b[9]*limit;
        }
        let rate=((v[j]-v[i])*x+(v[j+1]-v[i+1])*y+(v[j+2]-v[i+2])*z)/len;
        let scale=(b[9]*(len-b[6])+b[10]*rate)/len;
        for (axis,h) in [x,y,z].iter().enumerate() {
            let force=scale*h;
            f[i+axis]+=force;
            f[j+axis]-=force;
        }
        damping+=b[10]*rate*rate*dt;
    }
    *losses=plastic;
    *losses.add(1)=damping;
}
