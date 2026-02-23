# CAGE analysis for zebrafish PGCs (SLIC-CAGE)

Code + notes for analyzing **SLIC-CAGE** data generated from **zebrafish Primordial Germ Cells (PGCs)** and matched **somatic cells** across two developmental stages.

---

## Dataset

Four biological conditions are used in this project:

- PGC at **high** stage  
- Somatic at **high** stage  
- PGC at **prim5** stage  
- Somatic at **prim5** stage  

Technical Replicates are merged.

---

## Pipeline overview

1. **Mapping**
2. **CAGE preprocessing + clustering (CAGEr)**
3. **Promoter architecture analyses**
   - interquantile width (sharp vs broad)
   - dinucleotide initiators + pattern heatmaps
   - genomic annotation of clusters (promoter/UTR/intron/etc.)
4. **Promoter shifting analysis** (maternal → zygotic usage changes)
5. **Enhancer discovery** (bidirectional CAGE using CAGEfightR)
6. **Integration with ATAC-seq**
   - Tn5 cutsite processing
   - promoter accessibility metaplots / heatmaps
   - footprinting (ATACseqQC)

---

## 1) Mapping

### Align reads (bowtie2)
```bash
bowtie2 --phred33-quals --threads $nCores -x <INDEX_PREFIX> \
  -1 <R1.fastq.gz> -2 <R2.fastq.gz> -S <out.sam>
```

## 2) Build CAGEr object and compute CTSS

Import BAMs into a `CAGEset` and compute CTSS signal.

## 2.1 Create `CAGEset`
- `genomeName = BSgenome.Drerio.UCSC.danRer7`
- `inputFiles = list of *.sorted.bam`
- `inputFilesType = "bam"`
- `sampleLabels = meaningful names`

## 2.2 Compute CTSS
Run `getCTSS()`

**Output**
- CTSS tracks stored inside the CAGEset
- record library sizes using `librarySizes()`

**CAGEr parameters**
### 2.3 Reverse cumulative plots
- `values = "raw"` first
- `fitInRange = c(10, 1000)`

### 2.4 Normalize tag counts (power-law)
- `method = "powerLaw"`
- `fitInRange = c(5, 1000)`

### 2.5 Cluster CTSS into tag clusters and consensus clusters
- `thresholdIsTpm = TRUE`
- `threshold = 1` TPM (you commonly use 1)
- `method = "distclu"`
- `maxDist = 20` (or alternative: 40 for sensitivity test)
- `removeSingletons = TRUE`
- `keepSingletonsAbove = 5`


## 3) Promoter architecture analyses

These steps operate on **tag clusters** 

### 3.1 Interquantile width distributions (sharp vs broad)
**Logic**
- Extract tag clusters per sample with `returnInterquantileWidth=TRUE`
- Optionally filter to standard chromosomes (`chr1..chr25`)
- Plot histograms / densities per sample

**rule used**
- sharp: `interquantile_width < 9`
- broad: `interquantile_width >= 9`

**Output**
- IQ width plots

### 3.2 Dinucleotide initiators at dominant CTSS
**Logic**
- Center on `dominant_ctss`
- Extract ±1 bp (2bp total) sequences
- Compute dinucleotide frequencies per sample

**Optional grouping**
- YR initiators: CA, TG, TA, CG
- YC initiators: CC, TC


### 3.3 Pattern heatmaps (TA/CG/SS/WW)
**Logic**
- Extract promoter-centered sequences (±400 bp)
- Order clusters by IQ width
- Compute pattern occurrence matrices
- Smooth heatmaps (kernel smoothing)
- Plot heatmaps for patterns: TA, CG, SS, WW


### 3.4 Genomic element annotation (ChIPseeker)
**Logic**
- Convert clusters to GRanges using cluster interval (e.g., `q0.1 → q0.9`)
- Annotate with `annotatePeak()` using:
  - TxDb for danRer7
  - `tssRegion=c(-3000,1000)`
  - `sameStrand=TRUE`
- Collapse categories (e.g., all exons together, all introns together)
- Plot proportions per sample


## 4) Shifting promoter analysis (maternal → zygotic)

**Goal:** quantify shifts in dominant TSS usage between early and late stages.

### 4.1 Define “early” and “late”
- Early: high stage
- Late: prim5 stage

### 4.2 Restrict to promoter-proximal clusters
**Logic**
- Annotate clusters with `annotatePeak()` around TSS
- Keep clusters within ±500 bp of TSS

### 4.3 Pair early vs late promoters
Two common strategies:
1) merge by transcript/gene identifiers
2) overlap-based pairing (resize + `findOverlaps`), then annotate

## 4.4 Compute strand-aware shift
**Definition**
- `difference = dominant_ctss_late - dominant_ctss_early`
- flip sign on `-` strand so “positive shift” is comparable across strands

**filters**
- require expression in both states (e.g. `tpm.dominant_ctss > 3` or `> 5`)
- define “strong shifting” as |difference| > 100 bp (as used)

**Outputs**
- tables of shifting promoters
- density plots / histograms of shift distributions
- subsets: PGC-only shifting, shared shifting, etc.


## 5) Enhancer analysis (CAGEfightR)

**Goal:** identify candidate enhancers using bidirectional CAGE transcription away from promoters.

### 5.1 Prepare stranded BigWigs
**Logic**
- Split BAM-derived coverage by strand (+ and -)
- Export to `.plus.bw` and `.minus.bw`
- Keep only standard chromosomes (chr1..chr25), trim out-of-bound

**Output**
- `sample.plus.bw`, `sample.minus.bw`

## 5.2 CAGEfightR bidirectional clustering
**Key parameters**
- `clusterBidirectionally(balanceThreshold = 0.8)`
- ensure TPM calculation and pooling are performed before clustering

**Post-filter**
- annotate vs TxDb and keep distal events:
  - e.g. `distanceToTSS > 100`

**Output**
- GRanges of bidirectional clusters
- filtered set of putative enhancers saved as RDS


## 6) Integrate CAGE with ATAC-seq

**Goal:** compare promoter accessibility with promoter classes/motifs/shifting.

### 6.1 Process ATAC BAMs into Tn5 cut sites
**Logic**
- adjust alignments for Tn5 offset (5 bp / -4 bp scheme in notes)
- derive cut sites by resizing to width=1
- remove blacklist overlaps
- compute fold-change coverage:
  - expected coverage = total aligned bases / genome size
  - fold-change = observed coverage / expected

**Parameters**
- genome size (you used ~1.42e9; ensure correct for danRer7)
- blacklist BED path
- width filters (exclude very long fragments)

**Output**
- strand-aware GRanges with fold-change scores
- optional export to BigWig

### 6.2 Define promoter windows from CAGE clusters
**Logic**
- use dominant CTSS
- resize to a fixed window
- optionally sort by IQ width (sharp → broad)

**Output**
- promoter window GRanges for PGC and Soma

### 6.3 ScoreMatrix and metaplots (genomation)
**Logic**
- `ScoreMatrix(target = atac_foldchange, windows = promoter_windows, weight.col="score", strand.aware=TRUE)`
- average across promoters (colMeans) for metaplot
- rescale signals 0–1 for comparability across samples

**Outputs**
- promoter accessibility metaplots
- heatmaps (optional) for subsets (e.g. WWWWW/TATA-like promoters)

## 6.4 Footprinting (ATACseqQC)
**Goal:** infer TF binding footprints (e.g. TBP/TBP2 motifs).
**Logic**
- use MotifDb PFMs
- run footprint + vPlot around motif occurrences

**Outputs**
- footprint profiles and v-plots

---

## Key parameters summary (defaults used in this project)

**CAGEr**
- reverse cumulative fit: `fitInRange = c(10, 1000)`
- normalization: `method="powerLaw"`, `fitInRange=c(5,1000)`, `alpha≈1.15`, `T=1e6`
- tag clusters: `threshold=1 TPM`, `method="distclu"`, `maxDist=20` (alt: 40)
- IQ widths: `qLow=0.1`, `qUp=0.9`
- consensus clusters: `tpmThreshold=5`, `maxDist=100`

**ChIPseeker annotation**
- promoter window: `tssRegion=c(-500,500)` (promoter-proximal filtering)
- genomic elements: `tssRegion=c(-3000,1000)` for broader categorization

**Shifting**
- expression filter: typically `tpm.dominant_ctss > 3–5` in both conditions
- strong shift cutoff: `|difference| > 100 bp`
- strand-aware sign flip for `-` strand

**CAGEfightR**
- `balanceThreshold=0.8`
- distal filter: `distanceToTSS > 100`

**ATAC**
- Tn5 shift: +5 / -4 convention
- fold-change vs expected coverage
- ScoreMatrix: `strand.aware=TRUE`, rescale metaplot to 0–1 for cross-sample plots

---
