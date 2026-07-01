/* ============================================================================
 * aligner.js: the alignment engine, ported from my Module 2 Python code.
 *
 * Same dynamic-programming recurrence as globalAlignment.py / localAlignment.py,
 * generalized to a substitution matrix `sub` so it scores DNA (flat) and protein
 * (BLOSUM/PAM) with one function. Runs entirely in the browser.
 *
 * Scoring convention (identical to the Python): `match` rewards are ADDED; the
 * `gap` penalty is a positive magnitude that is SUBTRACTED. For proteins the
 * diagonal score comes straight from the matrix.
 * ==========================================================================*/

// Flat DNA scoring expressed as a substitution matrix.
function dnaMatrix(match, mismatch, alphabet = "ACGTN") {
  const M = {};
  for (const a of alphabet) {
    M[a] = {};
    for (const b of alphabet) M[a][b] = a === b ? match : -mismatch;
  }
  return M;
}

/* Global (Needleman–Wunsch) or local (Smith–Waterman) alignment.
 * A single Float64Array holds the score table (index i*W+j) so even the ~3800×3800
 * spike DNA table fits comfortably and fills in well under a second. */
function align(seq1, seq2, sub, gap, mode = "global", name1 = "seq1", name2 = "seq2") {
  const local = mode === "local";
  const n = seq1.length, m = seq2.length, W = m + 1;
  const T = new Float64Array((n + 1) * W);

  if (!local) {
    for (let j = 1; j <= m; j++) T[j] = -gap * j;
    for (let i = 1; i <= n; i++) T[i * W] = -gap * i;
  }
  for (let i = 1; i <= n; i++) {
    const srow = sub[seq1[i - 1]];
    const base = i * W, pbase = (i - 1) * W;
    for (let j = 1; j <= m; j++) {
      const c = seq2[j - 1];
      const diag = T[pbase + j - 1] + (srow ? (srow[c] ?? 0) : 0);
      const up = T[pbase + j] - gap;
      const left = T[base + j - 1] - gap;
      let best = diag >= up ? diag : up;
      if (left > best) best = left;
      if (local && best < 0) best = 0;
      T[base + j] = best;
    }
  }

  // Where does the traceback start?
  let i, j, score;
  if (local) {
    let best = 0, bi = 0, bj = 0;
    for (let r = 0; r <= n; r++) {
      const rb = r * W;
      for (let c = 0; c <= m; c++) {
        if (T[rb + c] > best) { best = T[rb + c]; bi = r; bj = c; }
      }
    }
    i = bi; j = bj; score = best;
  } else {
    i = n; j = m; score = T[n * W + m];
  }
  const end1 = i, end2 = j;

  const top = [], bot = [];
  while ((i > 0 && j > 0) && (!local || T[i * W + j] > 0)) {
    const cur = T[i * W + j];
    const srow = sub[seq1[i - 1]];
    const s = srow ? (srow[seq2[j - 1]] ?? 0) : 0;
    if (cur === T[(i - 1) * W + (j - 1)] + s) {
      top.push(seq1[i - 1]); bot.push(seq2[j - 1]); i--; j--;
    } else if (cur === T[(i - 1) * W + j] - gap) {
      top.push(seq1[i - 1]); bot.push("-"); i--;
    } else if (cur === T[i * W + (j - 1)] - gap) {
      top.push("-"); bot.push(seq2[j - 1]); j--;
    } else {
      break;
    }
  }
  if (!local) {
    while (i > 0) { top.push(seq1[i - 1]); bot.push("-"); i--; }
    while (j > 0) { top.push("-"); bot.push(seq2[j - 1]); j--; }
  }
  const start1 = i, start2 = j;
  top.reverse(); bot.reverse();
  return {
    name1, name2, top: top.join(""), bottom: bot.join(""),
    score, mode, gap, start1, end1, start2, end2,
  };
}

/* Alignment with AFFINE gaps (Gotoh): cost(L) = open + extend*(L-1).
 * Three matrices: M (ending aligned), X (ending in a gap in seq2, consuming
 * seq1), Y (ending in a gap in seq1, consuming seq2). Supports both global
 * (Needleman-Wunsch-Gotoh) and local (Smith-Waterman-Gotoh, M floored at 0). */
function alignAffine(seq1, seq2, sub, openPen, extendPen, mode = "global", name1 = "seq1", name2 = "seq2") {
  const local = mode === "local";
  const n = seq1.length, m = seq2.length, W = m + 1, NEG = -Infinity;
  const M = new Float64Array((n + 1) * W).fill(NEG);
  const X = new Float64Array((n + 1) * W).fill(NEG); // gap in seq2 (consume seq1)
  const Y = new Float64Array((n + 1) * W).fill(NEG); // gap in seq1 (consume seq2)

  if (local) {
    // A local alignment may start anywhere for free (and never in a gap state).
    for (let i = 0; i <= n; i++) M[i * W] = 0;
    for (let j = 0; j <= m; j++) M[j] = 0;
  } else {
    M[0] = 0;
    for (let i = 1; i <= n; i++) X[i * W] = -openPen - (i - 1) * extendPen;
    for (let j = 1; j <= m; j++) Y[j] = -openPen - (j - 1) * extendPen;
  }

  let best = local ? 0 : NEG, bi = 0, bj = 0;
  for (let i = 1; i <= n; i++) {
    const srow = sub[seq1[i - 1]];
    for (let j = 1; j <= m; j++) {
      const p = (i - 1) * W + (j - 1);
      const s = srow ? (srow[seq2[j - 1]] ?? 0) : 0;
      let mij = s + Math.max(M[p], X[p], Y[p]);
      if (local && mij < 0) mij = 0;      // the 0-floor that makes it Smith-Waterman
      M[i * W + j] = mij;
      X[i * W + j] = Math.max(M[(i - 1) * W + j] - openPen, X[(i - 1) * W + j] - extendPen);
      Y[i * W + j] = Math.max(M[i * W + (j - 1)] - openPen, Y[i * W + (j - 1)] - extendPen);
      if (local && mij > best) { best = mij; bi = i; bj = j; }
    }
  }

  // Where does the traceback start?
  let i, j, state, score;
  if (local) {
    i = bi; j = bj; state = "M"; score = best;   // best local alignment ends in a match
  } else {
    i = n; j = m;
    const mE = M[n * W + m], xE = X[n * W + m], yE = Y[n * W + m];
    score = Math.max(mE, xE, yE);
    state = score === mE ? "M" : score === xE ? "X" : "Y";
  }

  const top = [], bot = [];
  while (i > 0 || j > 0) {
    if (local && state === "M" && M[i * W + j] === 0) break;   // reached the local start
    if (state === "M") {
      top.push(seq1[i - 1]); bot.push(seq2[j - 1]);
      const p = (i - 1) * W + (j - 1);
      const mx = Math.max(M[p], X[p], Y[p]);
      i--; j--; state = mx === M[p] ? "M" : mx === X[p] ? "X" : "Y";
    } else if (state === "X") {
      top.push(seq1[i - 1]); bot.push("-");
      const opened = X[i * W + j] === M[(i - 1) * W + j] - openPen;
      i--; state = opened ? "M" : "X";
    } else {
      top.push("-"); bot.push(seq2[j - 1]);
      const opened = Y[i * W + j] === M[i * W + (j - 1)] - openPen;
      j--; state = opened ? "M" : "Y";
    }
    if (!local && i === 0 && j === 0) break;
    if (!local && i === 0 && j > 0) { while (j > 0) { top.push("-"); bot.push(seq2[j - 1]); j--; } break; }
    if (!local && j === 0 && i > 0) { while (i > 0) { top.push(seq1[i - 1]); bot.push("-"); i--; } break; }
  }
  const start1 = i, start2 = j;
  top.reverse(); bot.reverse();
  return {
    name1, name2, top: top.join(""), bottom: bot.join(""),
    score, mode: local ? "local-affine" : "global-affine", gap: openPen,
    start1: local ? start1 : 0, end1: local ? bi : n, start2: local ? start2 : 0, end2: local ? bj : m,
  };
}

// Readouts.
function alignmentStats(res) {
  const top = res.top, bottom = res.bottom;
  let matches = 0, mismatches = 0, gapCols = 0;
  for (let k = 0; k < top.length; k++) {
    const x = top[k], y = bottom[k];
    if (x === "-" || y === "-") gapCols++;
    else if (x === y) matches++;
    else mismatches++;
  }
  const countGaps = (s) => {
    let n = 0, cur = 0, longest = 0;
    for (const ch of s) {
      if (ch === "-") { cur++; if (cur === 1) n++; if (cur > longest) longest = cur; }
      else cur = 0;
    }
    return [n, longest];
  };
  const [g1n, g1l] = countGaps(top);
  const [g2n, g2l] = countGaps(bottom);
  const aligned = matches + mismatches;
  const cols = top.length;
  return {
    score: res.score,
    columns: cols,
    // Percent identity counts EVERY column, gaps included: an insertion or deletion
    // is a non-match, so indels lower identity (this is the BLAST / EMBOSS convention).
    percentIdentity: cols ? (100 * matches / cols) : 0,
    // The older "over aligned columns only" number, kept so we can show that
    // ignoring gaps lets identity be inflated toward 100%.
    identityAligned: aligned ? (100 * matches / aligned) : 0,
    matches, mismatches,
    numGaps: g1n + g2n,
    gapColumns: gapCols,
    longestGap: Math.max(g1l, g2l),
    percentGaps: cols ? (100 * gapCols / cols) : 0,
  };
}

// Clean user input: uppercase, drop unknown symbols; for DNA map U->T.
function cleanSeq(seq, type) {
  seq = seq.toUpperCase().replace(/[^A-Z*]/g, "");
  if (type === "DNA") {
    seq = seq.replace(/U/g, "T");
    return [...seq].filter((c) => "ACGTN".includes(c)).join("");
  }
  const alphabet = Object.keys(MATRICES.BLOSUM62);
  return [...seq].filter((c) => alphabet.includes(c)).join("");
}

// Translate a coding DNA string to protein (stops as '*').
function translate(dna) {
  let aa = "";
  for (let i = 0; i + 3 <= dna.length; i += 3) aa += CODON[dna.slice(i, i + 3)] ?? "X";
  return aa;
}

/* ---------- E-value by Monte-Carlo ----------
 * The E-value is a local-alignment statistic (Karlin-Altschul). Rather than
 * hard-code the lambda/K constants (which only exist for specific matrix + gap
 * combinations), we estimate it directly and honestly: shuffle one sequence many
 * times, realign, and see how the real score compares to that null distribution.
 * A fixed seed makes it deterministic, so it does not jitter as you drag sliders. */
function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
function erfc(x) { // Abramowitz-Stegun 7.1.26
  const z = Math.abs(x), t = 1 / (1 + 0.3275911 * z);
  const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-z * z);
  return x >= 0 ? 1 - y : 1 + y;
}
function evalueMonteCarlo(seq1, seq2, realScore, alignScore, mode, N, seed = 0x9e3779b9) {
  const rnd = mulberry32(seed);
  const chars = [...seq2];
  const nullScores = new Array(N);
  for (let t = 0; t < N; t++) {
    for (let k = chars.length - 1; k > 0; k--) { const j = (rnd() * (k + 1)) | 0; const tmp = chars[k]; chars[k] = chars[j]; chars[j] = tmp; }
    nullScores[t] = alignScore(seq1, chars.join(""));
  }
  const mean = nullScores.reduce((a, b) => a + b, 0) / N;
  const sd = Math.sqrt(nullScores.reduce((a, b) => a + (b - mean) ** 2, 0) / N) || 1e-9;
  const z = (realScore - mean) / sd;
  let E;
  if (mode === "local") {                 // best local scores follow an extreme-value (Gumbel) distribution
    const beta = Math.sqrt(6) * sd / Math.PI;
    const u = mean - 0.5772156649 * beta;
    const arg = (realScore - u) / beta;
    E = arg > 30 ? Math.exp(-arg) : 1 - Math.exp(-Math.exp(-arg)); // extrapolate the far tail
  } else {                                 // global scores are roughly normal
    E = 0.5 * erfc(z / Math.SQRT2);
  }
  return { E: Math.max(E, 0), z, mean, sd, N };
}

if (typeof module !== "undefined") {
  module.exports = { dnaMatrix, align, alignAffine, alignmentStats, cleanSeq, translate, evalueMonteCarlo };
}
