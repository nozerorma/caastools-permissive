# CAAStools 2.0: What's New (Group Discussion)

**Latest Release:** v2.0.1-caap-beta (January 17, 2026)  
**Quick Overview:** Four major features simplify CAAS analysis and enable 100-1000× faster validation

---

## Quick Version History

```
v1.0  (Original)      → Basic CAAS discovery & bootstrap
                      
2026  Pair-aware mode → Support species pairing
  ↓   Optimization    → Speed up bootstrap 100-1000×
  
v2.0 (Release)        → Stable, production-ready with all features
  ↓   
v2.0.1 (Current)      → CAAP mode added (property-based convergence)
```

**Timeline:**

- **2024-2025:** Pair-aware mode + bootstrap optimization development
- **Jan 2026:** v2.0 stable release (all features)
- **Jan 17, 2026:** CAAP mode integration complete (v2.0.1-beta)

---

## The Four Game-Changing Features

### 1️⃣ **Pair-Aware Mode** — Paired Species Comparisons

**What's the problem it solves?**

You have species pairs with contrasting traits (baboon-chimp, dog-wolf, mouse-blind mole) and want to analyze them together in a coordinated way.

**How to use it:**
```
Species    			Group(0/1)    Pair_ID
Baboon     			1             1
Chimp      			0             1
Dog        			1             2
Wolf       			0             2
Mouse       		1             3
Blind mole  		0             3
```

**What it does:**

- Validates that pairs are complete (both species present)

- Tracks convergence at pair level

- Identifies which pairs share conserved amino acids

**Command:**
```bash
ct discovery --paired_mode -t trait.cfg -a alignment.msa
```

**Key insight:**

Better for analyzing related species pairs. Especially useful if your trait naturally clusters by pairs. Useful for traceable convergence events.

**Important notes:**

- Pair-wise missing data handling: Use `--max_fg_miss` and `--max_bg_miss` flags. In paired mode, missing species from foreground and background must be paired (same pair IDs). Both flags should have same value.
- **Critical constraint:** Never set `--max_fg_miss` or `--max_bg_miss`  lower than the number of pairs minus two. This ensures at least one complete pair is present in each group for valid CAAS testing.

---

### 2️⃣ **Max Conservation Permitted** — Allow Partial Overlap

**What's the problem it solves?**

Classical CAAS requires zero overlap between foreground and background amino acids. But sometimes species converged with 1-2 amino acids in common—does that negate convergence?

**How to use it:**
```bash
ct discovery --max_conserved 1 -a alignment.msa -t trait.cfg
```

**What it does:**

- Tests CAAS allowing up to N overlapping amino acids

- Default: 0 (strict, classical CAAS)

- Typical range: 0-3 (depends on biological expectation and dataset size)

**Real example:**
```
Foreground: A, A, A      (conserved: 3 concordant AAs)
Background: A, D, D      (conserved: 2 AAs changed, 1 conserved)

Overlap: 1 (one BG species shares same AA as foreground)
Result with --max_conserved 0: NOT CAAS ✗
Result with --max_conserved 1: CAAS ✓
```

**Why this matters:**

- More realistic biology (minor exceptions happen)

- Better sensitivity without losing specificity

- Paired mode automatically tracks which pairs share conserved AAs

**Important notes:**

- The `--max_conserved` threshold should be set considering your dataset size and biological prior
- In paired mode: Conserved pairs are tracked separately. A pair is marked as conserved if all its species share overlapping amino acids
- Recommendation: Start with 0 (strict), then test with 1-2 if needed for sensitivity

---

### 3️⃣ **CAAP Mode** — Detect Property-Based Convergence

**What's the problem it solves?**

Two proteins converged to completely different amino acids—but the same *chemical properties*. Classical CAAS misses this.

**Real example:**
```
Foreground:  A (small), S (polar), T (polar)      → Small, polar
Background:  Q (polar), N (polar), D (charged)    → Polar, charged
AA overlap: NONE (classical CAAS ✓)
BUT: Both groups converged to POLAR properties
```

**How it works:**

Tests convergence at 5 different property levels:

- **GS0:** Classical identity (A→A, S→S, etc.)

- **GS1:** Basic properties (hydrophobic, polar, charged)

- **GS2:** Charge + size

- **GS3:** Functional roles

- **GS4:** Fine-grained biochemistry

**How to use it:**
```bash
ct discovery --caap_mode -a alignment.msa -t trait.cfg
```

**What the output looks like:**

One position generates 5 lines (one per scheme):
```
BRCA2  1234  GS0  p=0.001   (AA identity)
BRCA2  1234  GS1  p=0.0001  ← STRONGEST SIGNAL (property-based)
BRCA2  1234  GS2  p=0.0005
BRCA2  1234  GS3  p=0.0001
BRCA2  1234  GS4  p=0.01
```

**How to interpret:**

- Multiple schemes significant = genuine functional convergence

- Only GS0 significant = identity-driven

- GS0 fails but GS1-3 pass = property-level convergence is real

---

### 4️⃣ **Permulation Optimization** — Bootstrap 100-1000× Faster

**What's the problem it solves?**

Bootstrap on 20,000 genes × 10,000 cycles = unfeasible.

**How it works:**
```
Step 1: Run DISCOVERY (identifies positions with CAAS signal)
Step 2: Run BOOTSTRAP on ONLY those positions (skip positions without CAAS)
Result: Skip unnecessary testing → speed up proportional to CAAS % of total
```

**How to use it:**
```bash
# Step 1: Discovery (quick, tests all positions)
ct discovery -a alignment.msa -t trait.cfg -o discovery.tab

# Step 2: Bootstrap with optimization (FAST!)
ct bootstrap -s resampled/ -a alignment.msa -t trait.cfg \
    --discovery_out discovery.tab -o bootstrap.tab
```

**Technical implementation:**

- **Resample file chunking:** Load resampled traits in manageable chunks to avoid memory issues with large bootstrap datasets
- **Position filtering:** Parse discovery output, extract positions with CAAS, restrict bootstrap testing to only these positions
- **Consistency across modes:** When using CAAP mode, all 5 schemes (GS0-GS4) tested for each position to maintain discovery-bootstrap consistency
- **Throughput increase:** Architecture changes allow processing from ~1K resampled groups (v1.0) up to 1M+ in bootstrap (v2.0)
- **Speedup:** Real speedup depends on CAAS prevalence. Typical: 50-300× (not guaranteed 1000×, dataset-dependent) 

---

**Q: How does PhyloPhere use this?**  
A: PhyloPhere 2.0 (sel_algorithm_v2 branch) integrates all these features with automatic parallelization across HPC clusters.

---

## Migration from v1.0

**If you have old scripts:**
```bash
# Old command (still works)
ct discovery -a alignment.msa -t trait.cfg

# New with optimization (recommended)
ct discovery -a alignment.msa -t trait.cfg -o discovery.tab
ct bootstrap -s resampled/ -a alignment.msa -t trait.cfg --discovery_out discovery.tab
```

**What changed in outputs?**

**Discovery output new columns (paired mode):**
- `Paired_FG`: Number of foreground pairs with CAAS
- `Paired_BG`: Number of background pairs with CAAS
- `Pair_Count`: Total pairs analyzed

**Discovery output new columns (CAAP mode):**
- `Scheme`: Which GS scheme (GS0-GS4) this row tests
- `CAAP_Group`: Property group instead of individual AA

**Bootstrap output new columns (paired mode):**
- Species output sorted by pair ID for organization

- Existing columns unchanged when features not used
- Fully backward compatible with old parsing scripts

---

## Visual Workflow

```
START
  │
  ├─ Have paired species? ──→ Add --paired_mode
  │  (human-chimp, etc)      
  │
  ├─ Worried about AA overlap? ──→ Add --max_conserved N
  │  (Allow N shared AAs)       (N = number to allow)
  │
  ├─ Want property-level detection? ──→ Add --caap_mode
  │  (Similar properties, different AAs)
  │
  ├─ Testing 1000+ genes in bootstrap mode? ──→ Add --discovery_out FILE
  │  (Essential for speed)
  │
  └─ RUN ANALYSIS ──→ RESULTS
     Discovery          (positions with CAAS)
     Bootstrap          (p-values & validation)
```

---

## Repository Status

**caastools-permissive (nozerorma/caastools-permissive)**

- Current branch: main
- Latest: ffedb8c (CAAP working - Jan 17, 2026)
- All features integrated in main branch
- Status: Active development

**PhyloPhere Integration (nozerorma/PhyloPhere)**

- Branch: sel_algorithm_v2 (active development)
- Status: Synchronized with caastools main
- Adds HPC parallelization + Nextflow workflow
- Main branch: Contains separate RER analysis pipelines

---

**Questions?** Email: [miguel.ramon@upf.edu]  
**Full Technical Details:** See CAASTOOLS_2.0_TECHNICAL_SUMMARY.md  
**Repository:** https://github.com/linudz/caastools  

---

**Document Version:** 1.0 (Group Discussion Version)  
**Generated:** January 18, 2026  
**Audience:** Research group, collaborators  
