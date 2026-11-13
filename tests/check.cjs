const fs = require('node:fs');
const assert = require('node:assert/strict');
const path = require('node:path');
const {execFileSync} = require('node:child_process');
const root = path.resolve(__dirname, '..');
const one = {clearance: 1.2, volume: 8, absorption: 0.9, bioavailability: 0.8};
const two = {clearance: 1.2, central_volume: 8, peripheral_volume: 15, intercompartmental_clearance: 2.5, absorption: 0.9, bioavailability: 0.8};
const dose = 100, rate = 12, duration = 3;
let checks = 0, worst = 0;
function evaluate(name, t, p) {
  const record = `record(${Object.entries(p).map(([k,v]) => `${k}=${v}`).join(',')})`;
  const args = name.endsWith('infusion') ? `${t},${rate},${duration},${record}` : `${t},${dose},${record}`;
  return JSON.parse(execFileSync(process.execPath, [path.join(__dirname, 'evaluate.cjs')], {
    input: JSON.stringify([`a.${name}(${args})`]), encoding: 'utf8', env: process.env,
  }))[0];
}
function near(actual, expected, tolerance, label) {
  const error = Math.abs(actual - expected) / Math.max(1, Math.abs(expected));
  worst = Math.max(worst, error);
  assert.ok(error <= tolerance, `${label}: actual=${actual}, expected=${expected}, error=${error}`);
  checks++;
}
// Independent oracle integrates depot/central/peripheral ODEs and absorption
// sensitivities by RK4. Infusion segments have constant input, so no RK stage
// straddles the discontinuous forcing. No eigenvalue formula enters this route.
function oracle(kind, t, p) {
  if (t < 0) return {value: 0, absorptionDerivative: 0};
  const v = p.volume ?? p.central_volume;
  const k10 = p.clearance / v;
  const k12 = p.volume ? 0 : p.intercompartmental_clearance / v;
  const k21 = p.volume ? 0 : p.intercompartmental_clearance / p.peripheral_volume;
  let y = [kind === 'oral' ? dose * p.bioavailability : 0, kind === 'bolus' ? dose : 0, 0, 0, 0, 0];
  const rhs = (z, input) => [
    -p.absorption*z[0], p.absorption*z[0]-(k10+k12)*z[1]+k21*z[2]+input, k12*z[1]-k21*z[2],
    -z[0]-p.absorption*z[3], z[0]+p.absorption*z[3]-(k10+k12)*z[4]+k21*z[5], k12*z[4]-k21*z[5],
  ];
  const advance = (span, input) => {
    const n = Math.max(1, Math.ceil(span/0.001)), h = span/n;
    const shifted = (a,b,s) => a.map((v,i) => v+s*b[i]);
    for (let i=0; i<n; i++) {
      const a=rhs(y,input), b=rhs(shifted(y,a,h/2),input), c=rhs(shifted(y,b,h/2),input), d=rhs(shifted(y,c,h),input);
      y=y.map((v,j)=>v+h*(a[j]+2*b[j]+2*c[j]+d[j])/6);
    }
  };
  if (kind === 'infusion') {advance(Math.min(t,duration),rate); advance(Math.max(0,t-duration),0);}
  else advance(t,0);
  return {value:y[1]/v, absorptionDerivative:y[4]/v};
}
for (const [prefix,p] of [['one',one],['two',two]]) {
  for (const kind of ['bolus','oral','infusion']) {
    for (const t of [-2,0,0.2,3,3.0001,8,24]) {
      near(evaluate(`${prefix}_compartment_${kind}`,t,p),oracle(kind,t,p).value,3e-11,`${prefix}/${kind}/${t}`);
    }
  }
}
const k10=two.clearance/two.central_volume, k12=two.intercompartmental_clearance/two.central_volume, k21=two.intercompartmental_clearance/two.peripheral_volume;
const gap=Math.sqrt((k10+k12-k21)**2+4*k12*k21), alpha=(k10+k12+k21+gap)/2, beta=k10*k21/alpha;
for (const [prefix,base,k] of [['one',one,one.clearance/one.volume],['two',two,alpha],['two',two,beta]]) {
  for (const offset of [-1e-9,0,1e-9]) {
    const p={...base,absorption:k+offset}, t=2.3, expected=oracle('oral',t,p);
    near(evaluate(`${prefix}_compartment_oral`,t,p),expected.value,3e-11,`${prefix}/equal-rate/value/${offset}`);
    const h=1e-5;
    const derivative=(evaluate(`${prefix}_compartment_oral`,t,{...p,absorption:p.absorption+h})-evaluate(`${prefix}_compartment_oral`,t,{...p,absorption:p.absorption-h}))/(2*h);
    near(derivative,expected.absorptionDerivative,1e-9,`${prefix}/equal-rate/derivative/${offset}`);
  }
}
// Exercise both sides of the series transitions and weak exchange at similar
// micro-rates, while keeping all parameters strictly inside their domains.
for (const x of [0.019999,0.02,0.020001]) {
  const p={...one,absorption:one.clearance/one.volume+x/2.3};
  near(evaluate('one_compartment_oral',2.3,p),oracle('oral',2.3,p).value,3e-11,`series-transition/${x}`);
}
for (const q of [1e-6,1e-10]) {
  const p={...two,clearance:1,central_volume:10,peripheral_volume:q/0.1,intercompartmental_clearance:q};
  for (const kind of ['bolus','oral','infusion']) near(evaluate(`two_compartment_${kind}`,2.3,p),oracle(kind,2.3,p).value,3e-11,`close-eigenvalues/${q}/${kind}`);
}
// These checks use the exported module, including dose and observation helpers.
function probes(expressions) {
  return JSON.parse(execFileSync(process.execPath, [path.join(__dirname, 'evaluate.cjs')], {
    input: JSON.stringify(expressions), encoding: 'utf8', env: process.env,
  }));
}
const parameters = 'record(clearance=1.2,volume=8.0,absorption=0.9,bioavailability=0.8)';
const vector = probes([`a.one_compartment_oral([0.0,1.0,2.0],100.0,${parameters})`])[0];
Object.values(vector.data).forEach((value, i) => near(value, oracle('oral', i, one).value, 3e-11, 'vector concentration'));
const normal = (y, mean, sd) => -0.5*((y-mean)/sd)**2-Math.log(sd)-0.5*Math.log(2*Math.PI);
const values = probes([
  `a.dose_superposition(a.one_compartment_bolus,2.0,[0.0,1.0],[10.0,20.0],${parameters})`,
  `a.infusion_superposition(a.one_compartment_infusion,4.0,[0.0,2.0],[12.0,6.0],[3.0,1.0],${parameters})`,
  'a.one_compartment_rhs(0.0,[10.0],record(input_rate=2.0,clearance=1.2,volume=8.0))',
  'sum(a.two_compartment_rhs(0.0,[10.0,5.0],record(input_rate=2.0,clearance=1.2,central_volume=8.0,peripheral_volume=15.0,intercompartmental_clearance=2.5)))',
  'a.michaelis_menten_rhs(0.0,[20.0],record(input_rate=1.0,vmax=4.0,km=1.0,volume=10.0))',
  'logdensityof(a.additive_error([1.0,2.0],0.1),[1.1,1.8])',
  'logdensityof(a.combined_error([1.0,2.0],0.1,0.2),[1.1,1.8])',
  'logdensityof(a.lognormal_error([1.0,2.0],0.2),[1.1,1.8])',
]);
near(values[0], 10/8*Math.exp(-0.15*2)+20/8*Math.exp(-0.15), 1e-12, 'repeated bolus');
near(values[1], oracle('infusion',4,one).value+0.5*(12/1.2)*(1-Math.exp(-0.15))*Math.exp(-0.15), 3e-11, 'repeated infusion');
near(values[2].data['0'], 0.5, 1e-12, 'one-compartment mass balance');
near(values[3], 0.5, 1e-12, 'two-compartment mass balance');
near(values[4].data['0'], 1-4*20/30, 1e-12, 'saturable rate');
near(values[5], normal(1.1,1,0.1)+normal(1.8,2,0.1), 1e-12, 'additive observation');
near(values[6], normal(1.1,1,Math.sqrt(0.05))+normal(1.8,2,Math.sqrt(0.17)), 1e-12, 'combined observation');
near(values[7], normal(Math.log(1.1),0,0.2)-Math.log(1.1)+normal(Math.log(1.8),Math.log(2),0.2)-Math.log(1.8), 1e-12, 'lognormal observation');
console.log(JSON.stringify({checks,worstNormalizedError:worst}));
