/* ============================================================================
 * app.js — UI wiring for the sequence aligner.
 * ==========================================================================*/
"use strict";

const $ = (id) => document.getElementById(id);
const AUTO_LIMIT = 500000;      // auto-recompute only below this many DP cells
const AFFINE_LIMIT = 1000000;   // affine allocates 3 tables; cap it

const EXAMPLES = {
  "Toy DNA":                    { type: "DNA",     s1: SEQS.toy1,       s2: SEQS.toy2,          n1: "toy_1",    n2: "toy_2" },
  "HBB DNA (human/mouse)":      { type: "DNA",     s1: SEQS.hbbHumanDNA, s2: SEQS.hbbMouseDNA,  n1: "human_HBB",n2: "mouse_HBB" },
  "HBB protein (human/gorilla)":{ type: "Protein", s1: SEQS.hbbHumanAA, s2: SEQS.hbbGorillaAA,  n1: "human",    n2: "gorilla" },
  "HBB protein (human/fish)":   { type: "Protein", s1: SEQS.hbbHumanAA, s2: SEQS.hbbZebrafishAA,n1: "human",    n2: "zebrafish" },
  "SARS spike DNA":             { type: "DNA",     s1: SEQS.sars1DNA,   s2: SEQS.sars2DNA,      n1: "SARS-CoV", n2: "SARS-CoV-2" },
  "SARS spike protein":         { type: "Protein", s1: SEQS.sars1AA,    s2: SEQS.sars2AA,       n1: "SARS-CoV", n2: "SARS-CoV-2" },
};

const state = { type: "DNA", mode: "global", gapModel: "linear" };

function esc(s) { return String(s).replace(/[&<>]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;" }[c])); }

/* ---------- segmented controls ---------- */
function wireSeg(id, onChange) {
  const seg = $(id);
  seg.addEventListener("click", (e) => {
    const b = e.target.closest("button"); if (!b) return;
    [...seg.children].forEach((c) => c.classList.toggle("active", c === b));
    onChange(b.dataset.v);
  });
}

/* ---------- read current parameters ---------- */
function params() {
  return {
    type: state.type, mode: state.mode, gapModel: state.gapModel,
    match: +$("match").value, mismatch: +$("mismatch").value, gap: +$("gap").value,
    matrix: $("matrix").value, open: +$("open").value, extend: +$("extend").value,
  };
}

/* ---------- render an alignment into a container ---------- */
function renderAlignment(res, el, { width = 60, maxCols = 600 } = {}) {
  const top = res.top, bot = res.bottom, total = top.length;
  const shown = Math.min(total, maxCols);
  let html = "";
  for (let start = 0; start < shown; start += width) {
    const end = Math.min(start + width, shown);
    let t = "", mid = "", b = "";
    for (let k = start; k < end; k++) {
      const x = top[k], y = bot[k];
      let cls, m;
      if (x === "-" || y === "-") { cls = "g"; m = " "; }
      else if (x === y) { cls = "m"; m = "|"; }
      else { cls = "x"; m = "."; }
      t += `<span class="res ${cls}">${x}</span>`;
      b += `<span class="res ${cls}">${y}</span>`;
      mid += `<span class="mid">${m}</span>`;
    }
    html += `<div class="block"><div>${t}<span class="coord">  ${end}</span></div>` +
            `<div class="track">${mid}</div><div>${b}</div></div>`;
  }
  if (shown < total) html += `<div class="muted">… showing first ${shown} of ${total} columns (stats cover all ${total}).</div>`;
  el.innerHTML = html || `<span class="muted">No alignment.</span>`;
}

function renderStats(res, el) {
  const s = alignmentStats(res);
  const cells = [
    ["score", s.score.toFixed(1)],
    ["identity", s.percentIdentity.toFixed(1) + "%"],
    ["columns", s.columns],
    ["mismatches", s.mismatches],
    ["gaps", s.numGaps],
    ["longest gap", s.longestGap],
  ];
  el.innerHTML = cells.map(([k, v]) => `<div class="stat"><div class="k">${k}</div><div class="v">${v}</div></div>`).join("");
}

/* ---------- the main compute ---------- */
function buildResult(p, s1, s2, n1, n2) {
  const sub = p.type === "DNA" ? dnaMatrix(p.match, p.mismatch) : MATRICES[p.matrix];
  if (p.gapModel === "affine" && p.mode === "global") {
    return alignAffine(s1, s2, sub, p.open, p.extend, n1, n2);
  }
  return align(s1, s2, sub, p.gap, p.mode, n1, n2);
}

let pendingLarge = false;
function compute(force) {
  const p = params();
  const s1 = cleanSeq($("seq1").value, p.type);
  const s2 = cleanSeq($("seq2").value, p.type);
  if (!s1.length || !s2.length) { $("result").innerHTML = `<span class="muted">Enter two sequences.</span>`; $("statgrid").innerHTML = ""; return; }
  const cells = s1.length * s2.length;
  if (p.gapModel === "affine" && cells > AFFINE_LIMIT) {
    $("result").innerHTML = `<span class="muted">Affine gaps are capped for very large sequences — switch to Linear, or use shorter sequences.</span>`;
    return;
  }
  if (!force && cells > AUTO_LIMIT) {
    pendingLarge = true;
    $("result").innerHTML = `<span class="muted">Large sequences (${s1.length.toLocaleString()} × ${s2.length.toLocaleString()}). Click <b>Align ▶</b> to run.</span>`;
    return;
  }
  pendingLarge = false;
  const n1 = $("name1").value || "seq1", n2 = $("name2").value || "seq2";
  $("result").innerHTML = `<span class="muted">Computing…</span>`;
  setTimeout(() => {
    const res = buildResult(p, s1, s2, n1, n2);
    renderStats(res, $("statgrid"));
    renderAlignment(res, $("result"));
  }, 10);
}

const debounce = (fn, ms = 120) => { let t; return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); }; };
const autoCompute = debounce(() => compute(false));

/* ---------- sync slider value labels ---------- */
function syncLabels() {
  $("matchV").textContent = $("match").value;
  $("mismatchV").textContent = $("mismatch").value;
  $("gapV").textContent = $("gap").value;
  $("openV").textContent = $("open").value;
  $("extendV").textContent = $("extend").value;
}

/* ---------- presets & good/bad ---------- */
function loadExample(name) {
  const ex = EXAMPLES[name];
  state.type = ex.type;
  [...$("type").children].forEach((c) => c.classList.toggle("active", c.dataset.v === ex.type));
  onTypeChange(ex.type, false);
  $("seq1").value = ex.s1; $("seq2").value = ex.s2;
  $("name1").value = ex.n1; $("name2").value = ex.n2;
  compute(true);
}

function setGood() {
  if (state.type === "DNA") { $("match").value = 1; $("mismatch").value = 1; $("gap").value = 2; }
  else { $("matrix").value = "BLOSUM62"; $("gap").value = 11; }
  setGapModel("linear"); syncLabels(); compute(true);
}
function setBad() {
  if (state.type === "DNA") { $("match").value = 2; $("mismatch").value = 2; $("gap").value = 0.1; }
  else { $("matrix").value = "BLOSUM62"; $("gap").value = 0.4; }
  setGapModel("linear"); syncLabels(); compute(true);
}
function setGapModel(v) {
  state.gapModel = v;
  [...$("gapModel").children].forEach((c) => c.classList.toggle("active", c.dataset.v === v));
  $("affineControls").classList.toggle("hidden", v !== "affine");
  $("affineHelp").classList.toggle("hidden", v !== "affine");
}

function onTypeChange(v, recompute = true) {
  state.type = v;
  $("dnaControls").classList.toggle("hidden", v !== "DNA");
  $("protControls").classList.toggle("hidden", v === "DNA");
  if (recompute) compute(false);
}

/* ---------- Good vs Bad demo (compare two alignments) ---------- */
function gbCard(title, tagClass, res) {
  const s = alignmentStats(res);
  const div = document.createElement("div");
  div.className = "card";
  div.innerHTML = `<span class="tag ${tagClass}">${title}</span>
    <div class="statgrid" style="grid-template-columns:repeat(3,1fr)">
      <div class="stat"><div class="k">score</div><div class="v">${s.score.toFixed(1)}</div></div>
      <div class="stat"><div class="k">identity</div><div class="v">${s.percentIdentity.toFixed(1)}%</div></div>
      <div class="stat"><div class="k">gaps</div><div class="v">${s.numGaps}</div></div>
    </div>
    <div class="aln" style="margin-top:12px"></div>`;
  renderAlignment(res, div.querySelector(".aln"), { maxCols: 240 });
  return div;
}
function runGoodBad(which) {
  const [a, b, n1, n2] = which === "toy"
    ? [SEQS.toy1, SEQS.toy2, "toy_1", "toy_2"]
    : [SEQS.hbbHumanDNA, SEQS.hbbMouseDNA, "human_HBB", "mouse_HBB"];
  const good = align(a, b, dnaMatrix(1, 1), 2, "global", n1, n2);
  const bad = align(a, b, dnaMatrix(2, 2), 0.1, "global", n1, n2);
  const box = $("gbCompare"); box.innerHTML = "";
  box.appendChild(gbCard("GOOD · match 1, mismatch 1, gap 2", "good", good));
  box.appendChild(gbCard("BAD · gap ≈ 0 (0.1)", "bad", bad));
}

/* ---------- SARS DNA vs protein ---------- */
function bar(label, value, max, color) {
  const h = Math.max(4, Math.round(150 * value / max));
  return `<div class="bar"><div class="num">${value}</div>
    <div class="fill" style="height:${h}px;background:${color}"></div>
    <div class="lab">${label}</div></div>`;
}
function runSars() {
  $("sarsRun").textContent = "Computing… (~1s)"; $("sarsRun").disabled = true;
  setTimeout(() => {
    const dna = align(SEQS.sars1DNA, SEQS.sars2DNA, dnaMatrix(1, 1), 2, "global");
    const prot = align(SEQS.sars1AA, SEQS.sars2AA, MATRICES.BLOSUM62, 11, "global");
    const dSt = alignmentStats(dna), pSt = alignmentStats(prot);

    let c1 = 0, c2 = 0, identical = 0, synon = 0, missense = 0, indel = 0, pos = [0, 0, 0];
    for (let k = 0; k < prot.top.length; k++) {
      const x = prot.top[k], y = prot.bottom[k];
      if (x !== "-" && y !== "-") {
        const cod1 = SEQS.sars1DNA.substr(3 * c1, 3), cod2 = SEQS.sars2DNA.substr(3 * c2, 3);
        if (x === y && cod1 === cod2) identical++;
        else if (x === y) { synon++; for (let q = 0; q < 3; q++) if (cod1[q] !== cod2[q]) pos[q]++; }
        else { missense++; for (let q = 0; q < 3; q++) if (cod1[q] !== cod2[q]) pos[q]++; }
        c1++; c2++;
      } else if (x !== "-") { c1++; indel++; } else { c2++; indel++; }
    }
    const diff = synon + missense, totpos = pos[0] + pos[1] + pos[2];

    $("sarsCallout").innerHTML = [
      [`${dSt.percentIdentity.toFixed(1)}%`, "DNA identity"],
      [`${pSt.percentIdentity.toFixed(1)}%`, "protein identity — more conserved than its own gene"],
      [`${Math.round(100 * synon / diff)}%`, "of differing codons are synonymous (protein unchanged)"],
      [`${Math.round(100 * pos[2] / totpos)}%`, "of substitutions sit at the 3rd (wobble) position"],
    ].map(([n, t]) => `<div class="big"><div class="n">${n}</div><div class="t">${t}</div></div>`).join("");

    const mc = Math.max(identical, synon, missense);
    $("barsCodon").innerHTML =
      bar("identical", identical, mc, "var(--match)") +
      bar("synonymous", synon, mc, "#2b8cbe") +
      bar("missense", missense, mc, "var(--mismatch)");
    const mp = Math.max(...pos);
    $("barsPos").innerHTML =
      bar("1st", pos[0], mp, "var(--gap)") +
      bar("2nd", pos[1], mp, "#8593a1") +
      bar("3rd (wobble)", pos[2], mp, "#2b8cbe");

    $("sarsOut").classList.remove("hidden");
    $("sarsRun").textContent = "Re-run the analysis ▶"; $("sarsRun").disabled = false;
  }, 20);
}

/* ---------- quality demo: same alignment, different scores ---------- */
function scoreFlat(top, bot, match, mismatch, gap) {
  let s = 0;
  for (let k = 0; k < top.length; k++) {
    const x = top[k], y = bot[k];
    if (x === "-" || y === "-") s -= gap;
    else if (x === y) s += match;
    else s -= mismatch;
  }
  return s;
}
function runQualDemo() {
  const fixed = align(SEQS.hbbHumanDNA, SEQS.hbbMouseDNA, dnaMatrix(1, 1), 2, "global");
  const id = alignmentStats(fixed).percentIdentity.toFixed(1);
  const rows = [[1, 1, 2], [2, 1, 2], [5, 1, 2], [1, 1, 1]]
    .map(([m, mm, g]) => `<tr><td>match ${m}, mismatch ${mm}, gap ${g}</td><td>${scoreFlat(fixed.top, fixed.bottom, m, mm, g).toFixed(0)}</td></tr>`)
    .join("");
  $("qualOut").innerHTML =
    `<table class="stats"><tr><td class="muted">parameters (one fixed human/mouse HBB alignment)</td><td class="muted">score</td></tr>${rows}</table>
     <p class="muted" style="margin:8px 0 0">The alignment never changed — identity is ${id}% throughout — yet the score swings from
     ~290 to ~1770. That is why a bigger score can’t mean a better alignment.</p>`;
}

/* ---------- gap tie demo ---------- */
function runGapDemo() {
  const blockTop = "AAAAAAAAAA", blockBot = "AAAA------";
  const scatTop = "AAAAAAAAAA", scatBot = "A--A--AA--";
  const lin = (top, bot, g) => { let s = 0; for (let k = 0; k < top.length; k++) { const x = top[k], y = bot[k]; if (x === "-" || y === "-") s -= g; else if (x === y) s += 1; } return s; };
  const aff = (top, bot, o, e) => {
    let s = 0; for (let k = 0; k < top.length; k++) if (top[k] !== "-" && bot[k] !== "-" && top[k] === bot[k]) s += 1;
    for (const seq of [top, bot]) { let run = 0; for (const ch of seq + " ") { if (ch === "-") run++; else { if (run) s -= o + e * (run - 1); run = 0; } } }
    return s;
  };
  let rows = "";
  for (const g of [1, 2, 4]) rows += `<tr><td>LINEAR, gap ${g}</td><td>${lin(blockTop, blockBot, g)}</td><td>${lin(scatTop, scatBot, g)}</td><td>${lin(blockTop, blockBot, g) === lin(scatTop, scatBot, g) ? "TIE" : "differ"}</td></tr>`;
  rows += `<tr><td>AFFINE, open 6 ext 1</td><td>${aff(blockTop, blockBot, 6, 1)}</td><td>${aff(scatTop, scatBot, 6, 1)}</td><td>single gap wins</td></tr>`;
  $("gapOut").innerHTML =
    `<div class="aln" style="margin-bottom:6px">one 6-gap : ${blockTop} / ${blockBot}<br>three 2-gaps: ${scatTop} / ${scatBot}</div>
     <table class="stats"><tr><td class="muted">scoring</td><td class="muted">one 6-gap</td><td class="muted">three 2-gaps</td><td class="muted">verdict</td></tr>${rows}</table>`;
}

/* ---------- init ---------- */
function init() {
  // presets
  $("presets").innerHTML = Object.keys(EXAMPLES).map((k) => `<button class="chip" data-ex="${esc(k)}">${esc(k)}</button>`).join("");
  $("presets").addEventListener("click", (e) => { const b = e.target.closest(".chip"); if (b) loadExample(b.dataset.ex); });
  // matrix dropdown
  $("matrix").innerHTML = Object.keys(MATRICES).map((m) => `<option ${m === "BLOSUM62" ? "selected" : ""}>${m}</option>`).join("");
  // segmented controls
  wireSeg("type", (v) => onTypeChange(v));
  wireSeg("mode", (v) => { state.mode = v; compute(false); });
  wireSeg("gapModel", (v) => { setGapModel(v); compute(false); });
  // sliders + inputs
  ["match", "mismatch", "gap", "open", "extend"].forEach((id) => {
    $(id).addEventListener("input", () => { syncLabels(); autoCompute(); });
    // Don't let scrolling the page over a slider nudge its value.
    $(id).addEventListener("wheel", (e) => e.preventDefault(), { passive: false });
  });
  ["seq1", "seq2"].forEach((id) => $(id).addEventListener("input", autoCompute));
  ["name1", "name2"].forEach((id) => $(id).addEventListener("input", autoCompute));
  $("matrix").addEventListener("change", () => compute(false));
  // buttons
  $("btnGood").addEventListener("click", setGood);
  $("btnBad").addEventListener("click", setBad);
  $("btnAlign").addEventListener("click", () => compute(true));
  document.querySelectorAll("[data-demo]").forEach((b) => b.addEventListener("click", () => runGoodBad(b.dataset.demo)));
  $("sarsRun").addEventListener("click", runSars);
  $("qualDemo").addEventListener("click", runQualDemo);
  $("gapDemo").addEventListener("click", runGapDemo);

  syncLabels();
  loadExample("Toy DNA");
  runGoodBad("toy");
}
document.addEventListener("DOMContentLoaded", init);
