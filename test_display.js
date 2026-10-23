// Run with Node; no packages or browser required.
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");
const html = fs.readFileSync(path.join(__dirname, "index.html"), "utf8");
const source = html.split('<script id="app">')[1].split("</script>")[0];
new vm.Script(source);
const context = vm.createContext({});
vm.runInContext(
  source.slice(0, source.indexOf("function startViewer")),
  context,
);
const rows = [
  { system: "ring5", detector: 150252, gdml_detectors: [150252] },
  {
    system: "showermax",
    plane_detector: 73003,
    gdml_detectors: [70203, 70403, 70603, 70803],
  },
  { system: "showermax", plane_detector: 73004, gdml_detectors: [70204] },
];
const data = {
  counts: { 150252: 2, 73003: 3, 70203: 9000, 70403: 8000 },
  rates: { 150252: 8, 73003: 21, 70203: 90000, 70403: 80000 },
};
for (const [mode, expected, excluded] of [
  ["counts", [2, 3, 0], 17000],
  ["rates", [8, 21, 0], 170000],
]) {
  const result = context.mapHitData(rows, data, mode);
  assert.deepEqual(Array.from(result.mapped), expected);
  assert.equal(result.excluded, excluded);
}
delete data.counts[73003];
delete data.rates[73003];
for (const mode of ["counts", "rates"]) {
  const result = context.mapHitData(rows, data, mode);
  assert.equal(result.mapped[1], null);
  assert.equal(result.mapped[2], null);
  assert.equal(result.planeSeen, false);
}
assert.equal(
  context.mapHitData(
    [{ system: "showermax", gdml_detectors: [70203] }],
    data,
    "counts",
  ).mapped[0],
  null,
);
console.log(
  "PASS: syntax; main IDs; ShowerMax planes only in Hits and Rate; layer hits excluded; missing planes unavailable.",
);
