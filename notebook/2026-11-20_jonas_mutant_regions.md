# 2026-11-20 – Jonas (lecturer), picking up the mutant group's genomic-position idea

## Question for today
The mutant group plots `diff` against genomic position for `stp_vs_anc_at_TSB` and
`stk_vs_anc_at_TSB` and finds regions where many neighbouring proteins move together.
Two questions: (a) is that real, and what are those regions? (b) do the same regions come
out of the interaction contrasts `interaction_st*_vs_anc_at_TSB_vs_SN`?

## What we did
- Wrote `analyses/look_into_mutants_jonas.qmd` (renders in ~12 s).
- Joined `read_dea("total")` to the coordinates in `read_protein_metadata()`.
- Scanned the genome with a rolling median over k = 5 and k = 25 neighbouring genes,
  with a permutation null (shuffle which gene each `diff` belongs to, 1000x, 95 % of the
  max statistic as threshold).
- Verified every region on raw per-sample values (`read_abundances`), per genotype and
  per medium — not on the fold changes.
- STRING v12 REST API via the NCTC 8325 orthologs (taxid 93061): enrichment with the
  **detected proteome** as background, plus a network with `add_nodes = 15`.
  Responses cached in `data/derived/string_cache/` (git-ignored).

## Results (key numbers)
- **stk, 2.183–2.291 Mb:** 94 proteins, median log2 FC +0.72 (= 1.65x), 85 % above +0.3
  vs 13 % in the rest of the genome. Edges sharp to the single gene and they sit exactly
  between the two rRNA operons at 2.178–2.182 and 2.294–2.299 Mb. Present in **all 9 stk
  samples in all 3 media** (block minus rest ≈ +0.9 vs ≈ +0.3 for anc and stp).
  Best explanation: a tandem genomic amplification. Contains the whole *sigB* operon
  (*rsbU-rsbV-rsbW-sigB*), the *atp* operon and *kdp*.
- **stp, 2.103–2.125 Mb:** φSa3 prophage structural/replication module (the phage is
  inserted into *hlb*, which is why *hlb* appears twice, at 2.086 and 2.129 Mb).
  Down 1.5–2.5 log2 in all 9 stp samples, but **still detected** (often 8–9 of 9), so the
  phage is not deleted — it looks less active (less spontaneous induction).
- **Interaction contrasts: both regions cancel.** stk block median goes from +0.72 in each
  medium to −0.07 / +0.03 in the interactions. The φSa3 module leaves a −0.26 residue.
  That is the arithmetic of a difference of differences: anything constitutive drops out.
- What survives the interaction contrasts (FDR < 0.1, |log2FC| > 1, ~180 hits) is enriched
  for **Secreted** (21 of 69 detected, FDR 1e-07), Extracellular region (19/94), Signal
  (24/148), Peptidoglycan-anchor (6/14), Cell wall (6/16), Potassium transport (3/3) and
  branched-chain amino acid biosynthesis (5/9).
- *kdp*: cleanly +2.2 log2 in stk in TSB and SNpulse (3/3 detected everywhere). The huge
  interaction value (kdpB +5.3) is driven by the SN side, where detection drops to 1–3 of
  3. Reported as "solid in TSB/SNpulse, not established in SN".
- σB reporters: Asp23 (B7H15_12150, at 2.320 Mb, i.e. **outside** the block) is flat in
  stk (|log2FC| ≤ 0.18) — no evidence of raised σB output despite the extra *sigB* copy.
  In stp, Asp23 and AmaP are down 0.4–0.6 in every medium with small FDRs.

## Problems / surprises
- The narrow (k = 5) scan does **not** detect the stk block at all (0 windows). The block
  is a broad, modest plateau, not a spike; only the k = 25 scan sees it. Both window sizes
  are now in the report, because the difference in *kind* between the two features is the
  clue to what they are.
- First pass of the STRING enrichment used the whole genome as background and produced a
  convincing "iron acquisition" story. With the detected proteome as background **no iron
  or haem term survives FDR < 0.05.** The Isd proteins are simply abundant and easy to
  detect. The wrong-background version is called out explicitly in the report.
- Tried to confirm the amplification independently in the phospho data: too few sites in
  the block (16–84 per sample) and the per-sample medians scatter more than the effect.
  Inconclusive, not negative.
- The PNG device drops φ, σ and the Unicode minus; all plot labels had to be plain ASCII.

## AI use
Claude Code (Opus) did the analysis and wrote the notebook. Corrections needed: the
permutation result contradicted my first draft of the text (the stk peak does not pass the
narrow scan), the enrichment background issue above, and a wrong description of the STRING
network modules. All prose claims were re-checked against the rendered tables afterwards.

## Part 2 (same day): into the phospho data

### Coverage reality check
4,245 sites on 1,149 proteins, but only **565 detected in >= 9 of 27 samples** and 1,865
seen in exactly one sample. **96 % of the phospho DEA rows are lod_imputed.** Everything
below uses only the 565 well-detected sites.

### The amplification is confirmed independently
FragPipe site intensities are not cross-run normalised, so each sample was median-centred
first. Sites with >= 2 of 3 replicates in both anc and stk: **50 sites inside the block go
up by a median +0.45 log2 (1.36x) in stk, vs +0.04 for the 728 sites elsewhere**
(Wilcoxon p = 6.0e-05). Separate enrichment, separate prep, separate normalisation, same
block. My earlier "inconclusive" test failed only because it used all sites including
one-hit wonders and no centring — that caveat in the report has been corrected.

### The interaction hits could NOT be followed
Of the 157 interaction-hit proteins, 72 have any site in the table and only **15** have a
site detected in >= 9 of 27 samples. Site-level interaction estimates for those are almost
all FDR = 1 or imputed. So the secreted/cell-surface signature cannot be confirmed as a
phosphorylation effect — a limitation of coverage, not evidence against it. Secreted and
wall-anchored proteins are simply under-represented among detected phospho-sites.

### The phospho layer has its own, orthogonal interaction signal
Restricted to the 558 well-covered sites tested: 11 (stk) and 20 (stp) pass
FDR < 0.1 & |diff| > 1. They are central carbon/energy metabolism — pyk S538, eno S196/T376,
pgi T142, ptsG S526 + **H536**, carB T916, pyrC H249, ureC T345, rplR S36/S58 — and the
**protein level is flat** for nearly all of them. Only **2 of the 18** proteins overlap the
total-proteome interaction hits. BUT almost every one is lod_imputed, i.e. presence/absence,
which is why the fold changes are +4 to +10. Candidate list for targeted MS, not a result.

### The best result in the data set: ReoM Thr7
B7H15_08845 = ReoM (NCTC 8325 ortholog Q2FXW0), "peptidoglycan synthesis regulatory
protein". pThr7 (localisation probability **1.00**, detected in 18/27 samples):
- **3/3 in every anc group and every stp group, all three media** (~22.1–22.7 log2)
- **0/9 in stk**, in all three media
- ReoM *protein* equal in all three strains: 15.81 / 15.94 / 15.10 (9/9 each)
So: missing modification, not missing protein. ReoM Thr7 is the conserved PASTA-kinase
site (PrkA phosphorylates ReoM Thr7 in *L. monocytogenes*, Wamp et al. 2020 eLife 9:e56048,
doi:10.7554/eLife.56048; the Stk–ReoM–MurA axis is also described in *S. aureus*,
PMC13692261). This is a positive control that the Fe-IMAC enrichment and the STYH search
are producing real biology.
Note: T7 and Y12 lie on the **same tryptic peptide**, so they are one observation, not two.
Note: in the *interaction* contrast ReoM T7 is FDR = 1.0, because the loss is constitutive
— the same lesson as the genomic blocks, for the third time.

**Unexpected and worth following:** stp has ~16x less Stk1 protein (12.46 vs 16.64 log2,
5/9 vs 9/9 detected) yet keeps ReoM T7 fully phosphorylated. With no Stp1 to remove the
mark, residual kinase is enough. Kinase/phosphatase balance, not kinase amount.

### A prediction that could not be tested
If ReoM is unphosphorylated in stk, ClpCP should degrade MurA → less MurA in stk. Observed:
MurA **up** in stk, +0.89 (B7H15_11665) and +0.61 (B7H15_11795) — i.e. the block median.
**Both MurA/MurZ paralogues sit inside the amplified block**, so dosage and degradation push
in opposite directions and cannot be separated here. Another argument for sequencing.

## Part 3: STRING clustering + cluster labels + fold-change halos, via the API

**Yes, the API can do the clustering and the name overlay** — but they are NOT parameters of
the `network` endpoint. They belong to the **network image** endpoint:
`network_clustering_algorithm` (kmeans / MCL / leiden),
`network_clustering_parameter_{kmeans,mcl,leiden}`, `cluster_aware_layout`,
`network_clustering_edge_visibility` (faded/dotted/solid/noshow),
`show_cluster_names`, `cluster_name_size`.

**The clusters also come back as DATA**, which is the useful part:
`/api/tsv/network_image_url` (or `json`/`xml`) takes the same parameters and returns one row
per cluster — clusterNumber, clusterColor, hexColor, proteinCount, **clusterNames**,
preferredNames, stringIds. In `json` the first record additionally carries `imageURL`
(a rendered PNG, valid 6 months); the `tsv` output drops that column.

Clusters of our 135 mapped interaction hits (Leiden, resolution 1.4, default score cut-off):
Peptidoglycan-anchor / NEAr transporter domain (11, the Isd haem cluster); valine-leucine-
isoleucine biosynthesis (8); Secreted (7); nitrogen metabolism (6); lipid degradation (6);
phosphorelay response regulator / HrtB (5); De novo IMP biosynthesis (4); potassium
transport (4); polysaccharide biosynthesis (3); citrulline/arginine (3); acetyl-CoA
carboxylase (3). These line up with the §5.1 enrichment.

**Caveat found while testing:** the clustering depends strongly on `required_score`. At the
default threshold the secreted/surface proteins form clusters; at `required_score=700` they
fall apart into singletons, because the edges holding them together are co-expression and
text-mining links below 0.7. Always report the threshold with the picture.

**Fold-change halos.** The `svg/network` output is editable text: every node is a
`<g class='nwnodecontainer' data-safe_div_label='kdpB' data-x_pos data-y_pos data-radius>`.
So we walk the node tags, look the label up in the DEA table, and inject a coloured
`<circle>` of r = radius + 7 right after the opening tag (painted *behind* the node, so it
shows as a ring). Node interior = STRING cluster colour, ring = our interaction log2FC.
All 76 nodes in the figure got a value. Added an SVG colour-bar legend and a `viewBox`.
Use `block_structure_pics_in_bubbles=1` + `flat_node_design=1`: the SVG drops from 2.6 MB
to 124 kB and the halos read much better.

**Gotcha:** pasting the SVG inline into the HTML with `output: asis` works, but pandoc
lower-cases `viewBox` to `viewbox`, which SVG then ignores. Write the SVG to a file and use
`knitr::include_graphics()` instead; Quarto base64-embeds it untouched.

## Part 4: two halos, one per interaction contrast

Jonas asked whether a node can carry **two** rings, one per interaction contrast. It can —
the injection trick generalises: draw the larger disc first, the smaller one second, and
each paints over the middle of the previous one, so what stays visible is a pair of
concentric rings. The whole control is one named vector:

```r
RINGS <- c(stk = 10, stp = 5)   # pixels added to the node radius, outer listed first
```

Outer ring = `interaction_stk_vs_anc_at_TSB_vs_SN`, inner = `interaction_stp_vs_anc_at_TSB_vs_SN`,
node fill = STRING cluster colour. Legend now sits in a 160 px strip appended below the
network (height attribute and `viewBox` both extended), not on top of it.

**Two changes of substance, not just of drawing.**

1. A ring is now coloured by the **measured** effect in that contrast, hit or not; the
   *outline* carries significance (thick dark = FDR < 0.1 and |log2FC| > 1, thin pale =
   measured but not a hit). With the old hits-only colouring, grey would have meant both
   "not measured" and "measured but not significant". All 76 nodes have a value in both
   contrasts, so no grey ring occurs in this figure.
2. The two contrasts are much less redundant than they look: of the 157 interaction hit
   proteins, **131 are a hit in exactly one contrast and only 26 in both**. One ring per
   node hid that entirely.

**A mistake the second ring exposed (corrected in the report, §5.4 callout).** The
one-ring figure took, per protein, the interaction contrast with the larger |log2FC|. I
then described `hrtA` and `hrtB` as "ringed in deep blue". Wrong: only `hrtA` was blue
(stp, -4.98, FDR 0.014); `hrtB` was orange, because its value came from the *other*
contrast (stk, +2.83, FDR 0.076). A one-number-per-node summary over two contrasts
silently mixes sources. Either say per node which contrast the colour came from, or do not
reduce.

**What the figure now shows.**
- `kdp` is almost purely a `stk` effect: `kdpB` +5.30 (FDR 0.083) and `kdpD` +5.07
  (FDR 2.2e-07) on the outer ring, +0.07 and +0.60 (both ns) on the inner one. Only
  `kdpC` is a hit in both. `kdpA` is `lod_imputed` in both contrasts - presence/absence.
- BCAA biosynthesis leans the same way but is not exclusive: `ilvA` up in both (+3.01 /
  +2.32, FDR < 0.002), `leuA` / `leuB` / `ilvD` `stk`-only.
- `hrtA` and `hrtB` are one operon pointing opposite ways in the two contrasts. Treat as
  a warning sign (thin coverage?), not a finding - `plot_protein()` both before using it.

**Verification.** Rendered SVG checked programmatically: 76 nodes, 152 injected circles
(87 thick/hit + 65 thin), legend group present, `viewBox='0 0 1392 1356'` intact after
pandoc. Rasterised with `qlmanage -t` and inspected visually: `hrtA` inner blue / outer
pale and `hrtB` outer orange / inner pale both read correctly at node size.

**Limit.** Three rings still work (add an entry with a larger offset, listed first); past
that the rings get too thin and a small pie or bar glyph per node is the better shape -
same injection trick, different element.

## Next steps
1. **Whole-genome sequencing of anc, stk and stp.** Settles the amplification directly and
   checks for other secondary mutations. Highest-value next experiment by a distance.
2. φSa3 induction assay in stp (± mitomycin C).
3. **Targeted PRM/MRM on the §6.3 candidate sites** (pyk, eno, pgi, ptsG peptides) — turn
   presence/absence calls into measurements.
4. **ReoM follow-up:** why does stp keep Thr7 phosphorylated with 16x less kinase?
   Phospho-ReoM blot across the three strains would answer it directly.
5. More σB regulon members, given the stp Asp23 result.
