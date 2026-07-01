# Sequence Aligner — a browser-based, tunable pairwise aligner

An interactive **Needleman–Wunsch (global)** and **Smith–Waterman (local)** sequence aligner that runs
entirely in the browser. Tune the match reward, mismatch and gap penalties; switch between **DNA** and
**protein** (with BLOSUM / PAM matrices); and watch the alignment change. Built for the CMU Pre-College
Computational Biology **Module 2** project — a web version of the Colab notebook, running my own alignment
code ported from Python to JavaScript.

## What it does

- **The tool** — align any two sequences, live. Presets for a toy pair, HBB (hemoglobin β) DNA and protein,
  and the SARS-CoV / SARS-CoV-2 spike. Readouts: score, percent identity, gaps, longest gap, columns.
- **Good vs bad** — one click shows how a near-zero gap penalty reaches a *higher score and "100% identity"*
  by shredding the alignment into tiny gaps: technically optimal, biologically wrong.
- **DNA vs protein** — aligns the spike as DNA and as protein, then classifies every codon as
  identical / synonymous / missense and charts substitutions by codon position. Shows why the protein is
  more conserved than its own gene (synonymous mutations pile up at the wobble base).
- **The quality question** — why the raw score isn't a quality measure, and what to use instead.
- **The gap problem** — a proof that linear gaps can't tell one long indel from many short ones, and how
  affine gaps (open + extend, also built into the tool) fix it.

## Run it locally

It's plain static HTML/CSS/JS — no build step, no dependencies. Any of these work:

```bash
# option 1: just open the file
#   double-click index.html

# option 2: serve it (needed for some browsers' file:// restrictions is not an issue here,
#   but a server is cleaner)
python -m http.server 8123
#   then open http://localhost:8123
```

## Deploy to GitHub Pages

1. Create a repo and push these files (`index.html`, `styles.css`, `aligner.js`, `data.js`, `app.js`).
2. Repo **Settings → Pages → Build and deployment → Source: Deploy from a branch**, pick `main` / `/root`.
3. Your site goes live at `https://<username>.github.io/<repo>/`.

## Files

| File | What's in it |
|------|--------------|
| `index.html` | Page structure and controls |
| `styles.css` | All styling (light + dark, responsive) |
| `aligner.js` | The engine: global/local + affine alignment, readouts, translation (ported from my `globalAlignment.py` / `localAlignment.py`) |
| `data.js` | Substitution matrices (BLOSUM/PAM) and embedded example sequences |
| `app.js` | UI wiring, rendering, the good/bad and SARS analyses |

## Correctness

The JavaScript engine reproduces the Python module exactly — verified on the textbook Smith–Waterman example
(score 13), the good/bad parameter numbers, the SARS synonymous breakdown (438 synonymous / 275 missense /
609 third-position substitutions), and it agrees with Biopython's `PairwiseAligner` to the decimal on the HBB
proteins (777 / 584 / 418).

Data: HBB coding DNA (NCBI RefSeq `NM_000518.5`, `NM_008220.5`), HBB protein (UniProt), SARS-CoV &
SARS-CoV-2 spike CDS + protein. Matrices are the standard published BLOSUM / PAM tables.
