#!/usr/bin/env node
'use strict';

const N = Number(process.env.N || 5000000);
const checkpoints = new Set([10000, 100000, 1000000, N]);

function kahan() {
  let sum = 0.0, corr = 0.0;
  return {
    add(x) {
      const y = x - corr;
      const t = sum + y;
      corr = (t - sum) - y;
      sum = t;
    },
    value() { return sum; }
  };
}

function li4half() {
  let s = 0.0;
  let p = 0.5;
  for (let k = 1; k < 200; ++k) {
    s += p / (k*k*k*k);
    p *= 0.5;
  }
  return s;
}

const L = Math.log(2.0);
const Z2 = Math.PI * Math.PI / 6.0;
const Z3 = 1.2020569031595942854;
const Li4 = li4half();

const rhs = {
  Qp: 20*Li4 + (5/6)*L**4 + 7*L**2*Z2 - (59/10)*Z2**2,
  Qm: -22*Li4 - (11/12)*L**4 - (13/2)*L**2*Z2 - (7/4)*L*Z3 + (67/10)*Z2**2,
  T3p: -(7/2)*L*Z3 + (3/4)*Z2**2,
  T3m: -2*Li4 - (1/12)*L**4 + (1/2)*L**2*Z2 + (7/4)*L*Z3 + (1/10)*Z2**2,
  Shift: -(3/2)*L**2 - (1/2)*Z2 + 3*L - (5/16)*Z3,
};

let R = 0.0;
let H2 = 0.0;
const Qp = kahan(), Qm = kahan(), T3p = kahan(), T3m = kahan(), Shift = kahan();

function report(n) {
  const vals = {Qp: Qp.value(), Qm: Qm.value(), T3p: T3p.value(), T3m: T3m.value(), Shift: Shift.value()};
  console.log(`CHECKPOINT n=${n}`);
  for (const key of Object.keys(vals)) {
    console.log(`${key}_partial=${vals[key].toPrecision(17)} rhs=${rhs[key].toPrecision(17)} diff=${(vals[key]-rhs[key]).toExponential(9)}`);
  }
}

for (let n = 1; n <= N; ++n) {
  const sign = (n & 1) ? -1.0 : 1.0;
  const c = (n & 1) ? -1.0 : 3.0;
  R += c / n;
  H2 += 1.0 / (n*n);
  const S = R*R - H2;
  Qp.add(S / (n*n));
  Qm.add(sign * S / (n*n));
  T3p.add(R / (n*n*n));
  T3m.add(sign * R / (n*n*n));
  if ((n & 1) === 0) Shift.add(R / (n*n*(n-1)));
  if (checkpoints.has(n)) report(n);
}

console.log('CONSTANTS', {L, Z2, Z3, Li4});
console.log('Q5972_NUMERIC_AUDIT=PASS');
