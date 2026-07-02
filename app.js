/* ============================================================================
 * app.js: UI wiring for the sequence aligner.
 * ==========================================================================*/
"use strict";

const $ = (id) => document.getElementById(id);
const AUTO_LIMIT = 500000;      // auto-recompute only below this many DP cells
const AFFINE_LIMIT = 1000000;   // affine allocates 3 tables; cap it

// Two made-up proteins that share one 41-letter "domain" buried inside otherwise
// different sequence: a case where local alignment shines and global flounders.
const SHARED_DOMAIN = "GHSMGGLLARYFAQERPDLVKALVLLGAVDPPHELREAFPS";
const SHARED_A = "MKQLTA" + SHARED_DOMAIN + "RR";
const SHARED_B = "AADGKLYVSAETRAWNGLMHVCQ" + SHARED_DOMAIN + "TTLHYPAPQWKDNAEFLVRESE";

const EXAMPLES = {
  "Toy DNA":                    { type: "DNA",     s1: SEQS.toy1,       s2: SEQS.toy2,          n1: "toy_1",    n2: "toy_2" },
  "Shared domain (protein)":    { type: "Protein", s1: SHARED_A,        s2: SHARED_B,           n1: "protein_A",n2: "protein_B" },
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

/* ---------- render an alignment into a container ----------
 * For LOCAL mode, if the full sequences are provided (full1/full2), the matched core
 * is shown embedded in the whole sequences: the flanks on either side are drawn as
 * unaligned overhangs (residue-vs-gap), tagged .ctx so they look dimmed, while the
 * scored core keeps the normal match/mismatch colors. Stats are always on the core. */
function renderAlignment(res, el, { width = 60, maxCols = 600, full1 = null, full2 = null } = {}) {
  const isLocal = res.mode && res.mode.indexOf("local") === 0 && res.top.length > 0;
  let top = res.top, bot = res.bottom;
  let p1 = res.start1 || 0, p2 = res.start2 || 0;   // residue position where the DISPLAY begins
  let coreStart = 0, coreEnd = res.top.length;       // core column range within the displayed rows
  let ctxNote = "";

  if (isLocal && full1 != null && full2 != null) {
    const cap = 60, s1 = full1, s2 = full2;
    const lp1 = Math.min(res.start1, cap), lp2 = Math.min(res.start2, cap);
    const rp1 = Math.min(s1.length - res.end1, cap), rp2 = Math.min(s2.length - res.end2, cap);
    const leftS1 = s1.slice(res.start1 - lp1, res.start1), leftS2 = s2.slice(res.start2 - lp2, res.start2);
    const rightS1 = s1.slice(res.end1, res.end1 + rp1), rightS2 = s2.slice(res.end2, res.end2 + rp2);
    const dash = (k) => "-".repeat(k);
    const leftTop = leftS1 + dash(leftS2.length), leftBot = dash(leftS1.length) + leftS2;
    const rightTop = rightS1 + dash(rightS2.length), rightBot = dash(rightS1.length) + rightS2;
    top = leftTop + res.top + rightTop;
    bot = leftBot + res.bottom + rightBot;
    coreStart = leftTop.length; coreEnd = coreStart + res.top.length;
    p1 = res.start1 - lp1; p2 = res.start2 - lp2;
    const hiddenL = Math.max(res.start1 - lp1, res.start2 - lp2);
    const hiddenR = Math.max((s1.length - res.end1) - rp1, (s2.length - res.end2) - rp2);
    if (hiddenL || hiddenR) ctxNote = ` (flanks capped at ${cap} residues per side${hiddenL ? `; +${hiddenL} left` : ""}${hiddenR ? `; +${hiddenR} right` : ""})`;
  }

  const total = top.length;
  const shown = Math.min(total, maxCols);

  let header = "";
  if (isLocal) {
    const l1 = res.len1 != null ? ` of ${res.len1}` : "", l2 = res.len2 != null ? ` of ${res.len2}` : "";
    const inContext = full1 != null && full2 != null;
    header = `<div style="font-family:-apple-system,'Segoe UI',Roboto,sans-serif;font-size:12.5px;`
      + `color:var(--muted);margin-bottom:8px;white-space:normal">`
      + `<b>Local alignment.</b> The scored match spans <b>${esc(res.name1)}</b> ${res.start1 + 1}–${res.end1}${l1} `
      + `and <b>${esc(res.name2)}</b> ${res.start2 + 1}–${res.end2}${l2}, and the stat tiles above describe only that region. `
      + (inContext
          ? `Below it is shown in the full context of both sequences: the <b>colored</b> columns are the match, the <b>grey</b> letters on either side are the unaligned flanks (insertions / deletions).${ctxNote} `
          : "")
      + `Switch to <b>Global</b> to make the whole sequences line up end to end.</div>`;
  }

  let html = "";
  for (let start = 0; start < shown; start += width) {
    const end = Math.min(start + width, shown);
    let t = "", mid = "", b = "";
    for (let k = start; k < end; k++) {
      const x = top[k], y = bot[k];
      const flank = k < coreStart || k >= coreEnd;   // outside the scored core
      let cls, m;
      if (x === "-" || y === "-") { cls = flank ? "ctx" : "g"; m = " "; }
      else if (x === y) { cls = "m"; m = "|"; }
      else { cls = "x"; m = "."; }
      t += `<span class="res ${cls}">${x}</span>`;
      b += `<span class="res ${cls}">${y}</span>`;
      mid += `<span class="mid">${m}</span>`;
      if (x !== "-") p1++;
      if (y !== "-") p2++;
    }
    html += `<div class="block">` +
            `<div>${t}<span class="coord" title="${esc(res.name1)} position">  ${p1}</span></div>` +
            `<div class="track">${mid}</div>` +
            `<div>${b}<span class="coord" title="${esc(res.name2)} position">  ${p2}</span></div>` +
            `</div>`;
  }
  if (shown < total) html += `<div class="muted">… showing first ${shown} of ${total} displayed columns.</div>`;
  el.innerHTML = header + (html || `<span class="muted">No alignment.</span>`);
}

function renderStats(res, el) {
  const s = alignmentStats(res);
  const cells = [
    ["score", s.score.toFixed(1), "reward for matches minus penalties for mismatches and gaps"],
    ["identity (w/ gaps)", s.percentIdentity.toFixed(1) + "%", "matches / all columns; gaps count against it (BLAST / EMBOSS convention)"],
    ["length", s.columns, "total columns in the alignment"],
    ["mismatches", s.mismatches, "columns where two different residues line up"],
    ["gaps", s.numGaps, "number of separate gaps (indel runs) across both rows"],
    ["gap columns", s.gapColumns, "total dashed columns (sum of all gap lengths)"],
    ["longest gap", s.longestGap, "length of the longest single gap"],
  ];
  el.innerHTML = cells.map(([k, v, t]) =>
    `<div class="stat" title="${t || ""}"><div class="k">${k}</div><div class="v">${v}</div></div>`).join("");
}

/* ---------- the main compute ---------- */
function buildResult(p, s1, s2, n1, n2) {
  const sub = p.type === "DNA" ? dnaMatrix(p.match, p.mismatch) : MATRICES[p.matrix];
  if (p.gapModel === "affine") {
    return alignAffine(s1, s2, sub, p.open, p.extend, p.mode, n1, n2);
  }
  return align(s1, s2, sub, p.gap, p.mode, n1, n2);
}

// Is a raw input essentially a DNA string (only A/C/G/T/U/N)?
function isMostlyDNA(raw) {
  const letters = raw.toUpperCase().match(/[A-Z*]/g);
  if (!letters || !letters.length) return false;
  const dna = letters.filter((c) => "ACGTUN".includes(c)).length;
  return dna / letters.length >= 0.9;
}

// Decide what actually gets aligned, so DNA and protein aren't silently confused.
// Protein selected + DNA input  -> translate the coding DNA to protein (frame 1).
// DNA selected + protein input  -> flag it (you cannot score amino acids as DNA bases).
function resolveSeq(raw, type) {
  if (type === "Protein") {
    if (isMostlyDNA(raw)) {
      const dna = cleanSeq(raw, "DNA");
      const aa = translate(dna).replace(/\*+$/, "");
      return { seq: aa, translated: true, hasStop: aa.includes("*") };
    }
    return { seq: cleanSeq(raw, "Protein") };
  }
  if (!isMostlyDNA(raw) && /[A-Za-z]/.test(raw)) return { seq: cleanSeq(raw, "DNA"), notDNA: true };
  return { seq: cleanSeq(raw, "DNA") };
}

let pendingLarge = false;
function compute(force) {
  const p = params();
  const r1 = resolveSeq($("seq1").value, p.type);
  const r2 = resolveSeq($("seq2").value, p.type);

  if (r1.notDNA || r2.notDNA) {
    $("statgrid").innerHTML = "";
    $("result").innerHTML = `<div class="note" style="border-left-color:var(--mismatch)">These look like
      <b>protein</b> sequences, not DNA. Switch <b>Type</b> to <b>Protein</b> to score them with a substitution
      matrix. (The four DNA bases A/C/G/T are also valid amino-acid letters, so scoring amino acids as if they
      were DNA is not biologically meaningful.)</div>`;
    return;
  }
  const s1 = r1.seq, s2 = r2.seq;
  if (!s1.length || !s2.length) { $("result").innerHTML = `<span class="muted">Enter two sequences.</span>`; $("statgrid").innerHTML = ""; return; }

  const cells = s1.length * s2.length;
  if (p.gapModel === "affine" && cells > AFFINE_LIMIT) {
    $("result").innerHTML = `<span class="muted">Affine gaps are capped for very large sequences; switch to Linear, or use shorter sequences.</span>`;
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
    renderAlignment(res, $("result"), { full1: s1, full2: s2 });
    if (r1.translated || r2.translated) {
      const stopNote = (r1.hasStop || r2.hasStop)
        ? " Heads up: reading frame 1 contains a stop codon, so this DNA is not a complete open reading frame (real genes translate cleanly)."
        : "";
      $("result").insertAdjacentHTML("afterbegin",
        `<div class="note" style="margin-bottom:10px"><b>Translated to protein.</b> Protein is selected and the
         input looks like DNA, so it was translated (reading frame 1) and the amino-acid sequences were scored
         with ${p.matrix}.${stopNote}</div>`);
    }
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
function gbCard(title, tagClass, res, full1, full2) {
  const s = alignmentStats(res);
  const div = document.createElement("div");
  div.className = "card";
  div.innerHTML = `<span class="tag ${tagClass}">${title}</span>
    <div class="statgrid" style="grid-template-columns:repeat(3,1fr)">
      <div class="stat"><div class="k">score</div><div class="v">${s.score.toFixed(1)}</div></div>
      <div class="stat" title="matches / all columns; gaps count against it"><div class="k">identity (w/ gaps)</div><div class="v">${s.percentIdentity.toFixed(1)}%</div></div>
      <div class="stat" title="number of separate gaps (indel runs)"><div class="k">gaps</div><div class="v">${s.numGaps}</div></div>
    </div>
    <div class="aln" style="margin-top:12px"></div>`;
  renderAlignment(res, div.querySelector(".aln"), { maxCols: 240, full1, full2 });
  return div;
}
function runGoodBad(which) {
  const [a, b, n1, n2] = which === "toy"
    ? [SEQS.toy1, SEQS.toy2, "toy_1", "toy_2"]
    : [SEQS.hbbHumanDNA, SEQS.hbbMouseDNA, "human_HBB", "mouse_HBB"];
  const good = align(a, b, dnaMatrix(1, 1), 2, "global", n1, n2);
  const bad = align(a, b, dnaMatrix(2, 2), 0.1, "global", n1, n2);
  const box = $("gbCompare"); box.innerHTML = "";
  box.appendChild(gbCard("GOOD · match 1, mismatch 1, gap 2", "good", good, a, b));
  box.appendChild(gbCard("BAD · gap ≈ 0 (0.1)", "bad", bad, a, b));
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
      [`${pSt.percentIdentity.toFixed(1)}%`, "protein identity (more conserved than its own gene)"],
      [`${diff ? Math.round(100 * synon / diff) : 0}%`, "of differing codons are synonymous (protein unchanged)"],
      [`${totpos ? Math.round(100 * pos[2] / totpos) : 0}%`, "of substitutions sit at the 3rd (wobble) position"],
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
     <p class="muted" style="margin:8px 0 0">The alignment never changed (identity is ${id}% throughout), yet the score swings from
     ~290 to ~1770. That is why a bigger score cannot mean a better alignment.</p>`;
}

/* ---------- quality demo 2: the quality profile (identity is not enough) ---------- */
function runProfileDemo() {
  const good = align(SEQS.hbbHumanDNA, SEQS.hbbMouseDNA, dnaMatrix(1, 1), 2, "global");
  const bad = align(SEQS.hbbHumanDNA, SEQS.hbbMouseDNA, dnaMatrix(2, 2), 0.1, "global");
  const gs = alignmentStats(good), bs = alignmentStats(bad);
  const row = (label, g, b) => `<tr><td>${label}</td><td>${g}</td><td>${b}</td></tr>`;
  $("profileOut").innerHTML =
    `<table class="stats">
       <tr><td class="muted">measure (human vs mouse HBB)</td><td class="muted">honest params</td><td class="muted">"identity-gaming" params</td></tr>
       ${row("identity, counting gaps (what this tool shows)", gs.percentIdentity.toFixed(1) + "%", bs.percentIdentity.toFixed(1) + "%")}
       ${row("identity, ignoring gaps (the naive way)", gs.identityAligned.toFixed(1) + "%", "<b>" + bs.identityAligned.toFixed(1) + "%</b>")}
       ${row("number of gaps", gs.numGaps, "<b>" + bs.numGaps + "</b>")}
       ${row("alignment length", gs.columns + " cols", bs.columns + " cols")}
     </table>
     <p class="muted" style="margin:8px 0 0">Look at the second alignment: if you <i>ignore gaps</i>, its identity is a
     perfect <b>${bs.identityAligned.toFixed(0)}%</b>, which looks better than the honest alignment. That is the trap the
     naive definition falls into. Count the gaps (the row above it) and the same alignment drops to
     <b>${bs.percentIdentity.toFixed(0)}%</b> because it is riddled with <b>${bs.numGaps}</b> of them. Always read identity
     next to the gap count.</p>`;
}

/* ---------- when local beats global: a shared domain in two different proteins ---------- */
function runLocalVsGlobal() {
  const g = align(SHARED_A, SHARED_B, MATRICES.BLOSUM62, 11, "global", "protein_A", "protein_B");
  const l = align(SHARED_A, SHARED_B, MATRICES.BLOSUM62, 11, "local", "protein_A", "protein_B");
  const box = $("lvgOut"); box.innerHTML = "";
  box.appendChild(gbCard("Global: forced to line up the whole length", "bad", g, SHARED_A, SHARED_B));
  box.appendChild(gbCard("Local: finds just the shared domain", "good", l, SHARED_A, SHARED_B));
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
  $("localVsGlobal").addEventListener("click", runLocalVsGlobal);
  $("sarsRun").addEventListener("click", runSars);
  $("qualDemo").addEventListener("click", runQualDemo);
  $("profileDemo").addEventListener("click", runProfileDemo);
  $("gapDemo").addEventListener("click", runGapDemo);

  syncLabels();
  loadExample("Toy DNA");
  runGoodBad("toy");
}
document.addEventListener("DOMContentLoaded", init);
