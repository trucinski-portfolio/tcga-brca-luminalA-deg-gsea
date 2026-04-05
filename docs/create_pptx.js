const PptxGenJS = require("pptxgenjs");
const fs = require("fs");
const path = require("path");

const FIG = path.join(__dirname, "..", "results", "figures");
const OUT = path.join(__dirname, "pipeline_summary.pptx");

const SUBTYPE_COLORS = {
  LumA:  "2166AC",
  LumB:  "92C5DE",
  Her2:  "D6604D",
  Basal: "1A1A1A",
};

const BG    = "FFFFFF";
const DARK  = "1A1A2E";
const GRAY  = "555555";
const LGRAY = "F5F5F5";

function figPath(name) { return path.join(FIG, name); }

const pptx = new PptxGenJS();
pptx.layout = "LAYOUT_WIDE"; // 13.33 x 7.5 inches

// ─── helpers ──────────────────────────────────────────────────────────────────

function titleBar(slide, text, color = "2166AC") {
  slide.addShape(pptx.ShapeType.rect, {
    x: 0, y: 0, w: "100%", h: 0.55,
    fill: { color },
    line: { color, width: 0 },
  });
  slide.addText(text, {
    x: 0.2, y: 0, w: "95%", h: 0.55,
    fontSize: 20, bold: true, color: "FFFFFF",
    valign: "middle",
  });
}

function footerBar(slide) {
  slide.addShape(pptx.ShapeType.rect, {
    x: 0, y: 7.15, w: "100%", h: 0.35,
    fill: { color: "EEEEEE" },
    line: { color: "CCCCCC", width: 0.5 },
  });
  slide.addText("TCGA-BRCA PAM50 Pipeline  |  TOIL RNA-seq  |  limma · ElasticNet · fgsea · WGCNA · Cox PH", {
    x: 0.2, y: 7.15, w: "100%", h: 0.35,
    fontSize: 7, color: "888888", valign: "middle",
  });
}

function bodyText(slide, lines, x, y, w, h, fontSize = 11) {
  const bullets = lines.map(l => ({
    text: l.text !== undefined ? l.text : l,
    options: { bullet: l.bullet !== false, fontSize: l.size || fontSize, color: l.color || GRAY, bold: l.bold || false },
  }));
  slide.addText(bullets, { x, y, w, h, valign: "top", autoFit: true });
}

function statBox(slide, label, value, x, y, color = "2166AC") {
  slide.addShape(pptx.ShapeType.rect, {
    x, y, w: 1.9, h: 0.85,
    fill: { color: "F0F4FA" },
    line: { color, width: 1.5 },
    rectRadius: 0.05,
  });
  slide.addText(value, { x, y: y + 0.02, w: 1.9, h: 0.45, align: "center", fontSize: 22, bold: true, color });
  slide.addText(label, { x, y: y + 0.44, w: 1.9, h: 0.38, align: "center", fontSize: 8.5, color: GRAY });
}

// ─── Slide 1: Title ───────────────────────────────────────────────────────────
{
  const s = pptx.addSlide();
  // Full-width dark header
  s.addShape(pptx.ShapeType.rect, { x: 0, y: 0, w: "100%", h: 2.8, fill: { color: DARK }, line: { color: DARK } });

  s.addText("TCGA-BRCA PAM50 Multi-Subtype\nTranscriptomics Pipeline", {
    x: 0.5, y: 0.35, w: 12.3, h: 1.9,
    fontSize: 30, bold: true, color: "FFFFFF", valign: "middle",
  });
  s.addText("8-Stage All-R Analysis  ·  UCSC TOIL RNA-seq  ·  1,109 Samples  ·  18 Contrasts  ·  22 Publication Figures", {
    x: 0.5, y: 2.25, w: 12.3, h: 0.45,
    fontSize: 11, color: "AAAACC", italic: true,
  });

  // Subtype color swatches
  const subtypes = ["LumA", "LumB", "Her2", "Basal"];
  const labels = ["Luminal A  (n=420)", "Luminal B  (n=192)", "HER2-enriched  (n=66)", "Basal-like  (n=139)"];
  subtypes.forEach((st, i) => {
    const x = 0.5 + i * 3.2;
    s.addShape(pptx.ShapeType.rect, { x, y: 3.1, w: 2.7, h: 0.5, fill: { color: SUBTYPE_COLORS[st] }, line: { color: SUBTYPE_COLORS[st] }, rectRadius: 0.04 });
    s.addText(labels[i], { x, y: 3.1, w: 2.7, h: 0.5, align: "center", fontSize: 10, bold: true, color: st === "Basal" ? "DDDDDD" : "FFFFFF", valign: "middle" });
  });

  s.addText("Controls: NAT (n=113)  ·  GTEx Healthy Breast (n=179)", {
    x: 0.5, y: 3.75, w: 12.3, h: 0.4, fontSize: 11, color: GRAY, italic: true,
  });

  // Pipeline stage boxes
  const stages = ["00 Setup","01 Preprocessing","02 UMAP Clustering","03 DEG (limma)","04 ElasticNet","05 fgsea","06 WGCNA","07 Survival","08 Figures"];
  const stageColors = ["888888","2166AC","2166AC","D6604D","D6604D","4DAC26","6A4C93","E07B39","1A1A1A"];
  stages.forEach((st, i) => {
    const x = 0.3 + i * 1.42;
    s.addShape(pptx.ShapeType.rect, { x, y: 4.4, w: 1.32, h: 0.6, fill: { color: stageColors[i] }, line: { color: stageColors[i] }, rectRadius: 0.04 });
    s.addText(st, { x, y: 4.4, w: 1.32, h: 0.6, align: "center", fontSize: 7.5, bold: true, color: "FFFFFF", valign: "middle" });
    if (i < stages.length - 1) {
      s.addText("→", { x: x + 1.32, y: 4.55, w: 0.1, h: 0.3, fontSize: 10, color: "AAAAAA" });
    }
  });

  footerBar(s);
}

// ─── Slide 2: Cohort Overview ─────────────────────────────────────────────────
{
  const s = pptx.addSlide();
  titleBar(s, "Cohort Overview", "2166AC");
  footerBar(s);

  // Stat boxes
  statBox(s, "Total samples", "1,109",       0.25, 0.75, "2166AC");
  statBox(s, "Tumor (PAM50)",  "817",         2.25, 0.75, "D6604D");
  statBox(s, "NAT",            "113",         4.25, 0.75, "E07B39");
  statBox(s, "GTEx Healthy",   "179",         6.25, 0.75, "4DAC26");
  statBox(s, "Gene symbols",   "28,344",      8.25, 0.75, "6A4C93");
  statBox(s, "Contrasts",      "18",         10.25, 0.75, "888888");

  // UMAP figure
  s.addImage({ path: figPath("fig1_umap_overview.png"), x: 0.25, y: 1.75, w: 8.5, h: 5.1 });

  // Key notes
  bodyText(s, [
    { text: "3-way design — dual reference groups:", bold: true, bullet: false },
    { text: "GTEx Healthy = absolute dysregulation baseline" },
    { text: "NAT = field cancerization context" },
    { text: "Tumor = 4 PAM50 subtypes (Normal-like excluded)" },
    { text: " ", bullet: false },
    { text: "Subtype breakdown:", bold: true, bullet: false },
    { text: "LumA  420   (51.4%)", color: "2166AC" },
    { text: "LumB  192   (23.5%)", color: "557799" },
    { text: "Basal  139   (17.0%)", color: "555555" },
    { text: "Her2    66   ( 8.1%)", color: "D6604D" },
  ], 9.0, 1.75, 4.1, 5.0, 10);
}

// ─── Slide 3: Differential Expression ────────────────────────────────────────
{
  const s = pptx.addSlide();
  titleBar(s, "Differential Expression  —  18 Contrasts, limma-trend (eBayes, trend=TRUE)", "D6604D");
  footerBar(s);

  s.addImage({ path: figPath("fig2_deg_counts.png"), x: 0.25, y: 0.65, w: 6.0, h: 6.3 });

  bodyText(s, [
    { text: "Key DEG counts (FDR<0.05, |log₂FC|≥1):", bold: true, bullet: false },
    { text: "Subtype vs GTEx:  8,688–9,693 per subtype" },
    { text: "Subtype vs NAT:   5,575–7,091 per subtype" },
    { text: "NAT vs GTEx:       5,689 genes" },
    { text: "  → confirms field cancerization signal", bold: false, color: "888888" },
    { text: " ", bullet: false },
    { text: "Subtype vs subtype:", bold: true, bullet: false },
    { text: "LumA vs LumB:   1,343  (most similar)" },
    { text: "LumA vs Basal:  5,372  (most distant)" },
    { text: " ", bullet: false },
    { text: "All 28,344 genes written per file — no pre-filtering. Downstream stages apply own thresholds.", bold: false, color: "888888", size: 9 },
  ], 6.5, 0.75, 6.5, 6.2, 10.5);
}

// ─── Slide 4: Subtype Signatures ──────────────────────────────────────────────
{
  const s = pptx.addSlide();
  titleBar(s, "Subtype-Specific Gene Signatures  —  ElasticNet Bootstrap Stability Selection", "6A4C93");
  footerBar(s);

  s.addImage({ path: figPath("fig5_signature_heatmap.png"), x: 0.25, y: 0.65, w: 7.5, h: 6.3 });

  bodyText(s, [
    { text: "Method: Bootstrap ElasticNet (α=0.5)", bold: true, bullet: false },
    { text: "100 bootstrap iterations at fixed λ.1se" },
    { text: "Freq. cutoff: ≥0.80 (Meinshausen & Bühlmann 2010)" },
    { text: " ", bullet: false },
    { text: "Signature sizes:", bold: true, bullet: false },
    { text: "LumA:   30 genes", color: "2166AC" },
    { text: "LumB:   27 genes", color: "557799" },
    { text: "Her2:   36 genes", color: "D6604D" },
    { text: "Basal:  34 genes", color: "555555" },
    { text: " ", bullet: false },
    { text: "Her2/Basal slightly broader — smaller positive-class n drives wider λ.1se selection. Consistent with class-imbalanced ElasticNet.", bold: false, color: "888888", size: 9 },
    { text: " ", bullet: false },
    { text: "Row annotation: Direction (Up/Down) based on mean ElasticNet coefficient sign across bootstrap runs.", bold: false, color: "888888", size: 9 },
  ], 8.0, 0.75, 5.1, 6.2, 10.5);
}

// ─── Slide 5: Pathway Enrichment ──────────────────────────────────────────────
{
  const s = pptx.addSlide();
  titleBar(s, "Pathway Enrichment  —  fgsea, Hallmark + KEGG (236 gene sets, ranked by limma t-stat)", "4DAC26");
  footerBar(s);

  s.addImage({ path: figPath("fig6_fgsea_dotplot.png"), x: 0.25, y: 0.65, w: 8.0, h: 6.3 });

  bodyText(s, [
    { text: "Gene sets: 50 Hallmark + 186 KEGG Legacy", bold: true, bullet: false },
    { text: " ", bullet: false },
    { text: "Significant pathways (padj<0.05):", bold: true, bullet: false },
    { text: "LumA vs GTEx:   58  |  vs NAT:  66", color: "2166AC" },
    { text: "LumB vs GTEx:   73  |  vs NAT:  93", color: "557799" },
    { text: "Her2 vs GTEx:   67  |  vs NAT:  73", color: "D6604D" },
    { text: "Basal vs GTEx:  69  |  vs NAT:  83", color: "555555" },
    { text: " ", bullet: false },
    { text: "Ranking metric: limma t-statistic (not log₂FC) — integrates effect size with precision.", bold: false, color: "888888", size: 9 },
    { text: " ", bullet: false },
    { text: "Dual reference reveals field effect pathways — immune/stromal pathways enriched vs GTEx but not vs NAT suggest cancer-field-wide alterations.", bold: false, color: "888888", size: 9 },
  ], 8.5, 0.75, 4.6, 6.2, 10.5);
}

// ─── Slide 6: WGCNA ───────────────────────────────────────────────────────────
{
  const s = pptx.addSlide();
  titleBar(s, "Co-expression Network  —  WGCNA Signed Network, Tumor Samples Only (n=817)", "6A4C93");
  footerBar(s);

  s.addImage({ path: figPath("wgcna_module_trait_heatmap.png"), x: 0.25, y: 0.65, w: 6.5, h: 4.8 });
  s.addImage({ path: figPath("wgcna_soft_threshold.png"), x: 0.25, y: 5.55, w: 6.5, h: 1.4 });

  bodyText(s, [
    { text: "Network parameters:", bold: true, bullet: false },
    { text: "Soft-threshold power = 8  (R² = 0.907)" },
    { text: "Scale-free topology target: R² ≥ 0.85" },
    { text: "Gene selection: top 5,000 by MAD" },
    { text: " ", bullet: false },
    { text: "Modules identified:", bold: true, bullet: false },
    { text: "5 modules (excl. grey/unassigned)" },
    { text: "Sizes: 395–775 genes per module" },
    { text: "Grey (unassigned): 2,061 genes" },
    { text: " ", bullet: false },
    { text: "Module-trait associations:", bold: true, bullet: false },
    { text: "20 / 24 significant (p<0.05)" },
    { text: "83% of module-subtype pairs significant" },
    { text: " ", bullet: false },
    { text: "Signed network preserves co-expression direction — anti-correlated genes form separate modules.", bold: false, color: "888888", size: 9 },
  ], 7.0, 0.75, 6.1, 6.2, 10.5);
}

// ─── Slide 7: Survival ────────────────────────────────────────────────────────
{
  const s = pptx.addSlide();
  titleBar(s, "Survival Analysis  —  Cox PH, Signature Score vs Overall Survival", "E07B39");
  footerBar(s);

  s.addImage({ path: figPath("fig7_cox_forest.png"), x: 0.25, y: 0.65, w: 7.5, h: 3.2 });
  s.addImage({ path: figPath("survival_km_Basal.png"), x: 0.25, y: 3.95, w: 7.5, h: 3.0 });

  bodyText(s, [
    { text: "Cox model results:", bold: true, bullet: false },
    { text: "LumA:   HR=0.867, p=0.336  (n=370, 38 events)", color: "2166AC" },
    { text: "LumB:   HR=0.901, p=0.607  (n=167, 23 events)", color: "557799" },
    { text: "Her2:   HR=1.141, p=0.700  (n= 57, 11 events)", color: "D6604D" },
    { text: "Basal:  HR=0.696, p=0.114  (n=125, 16 events)", color: "555555" },
    { text: " ", bullet: false },
    { text: "Why null results are expected:", bold: true, bullet: false, color: "333333" },
    { text: "Signatures were optimized for subtype discrimination — a classification objective, not survival prediction." },
    { text: "A gene that identifies LumA tumors reliably ≠ a prognostic gene within LumA." },
    { text: " ", bullet: false },
    { text: "Basal trend (HR=0.70, p=0.11):", bold: true, bullet: false },
    { text: "Most plausible signal — TNBC is the highest-risk subtype. Low event count (n=16) limits power." },
    { text: "A larger Basal-enriched cohort is needed to test definitively." },
  ], 8.0, 0.75, 5.1, 6.2, 9.5);
}

// ─── Slide 8: Conclusions ─────────────────────────────────────────────────────
{
  const s = pptx.addSlide();
  titleBar(s, "Conclusions", "1A1A1A");
  footerBar(s);

  s.addImage({ path: figPath("fig4_deg_heatmap.png"), x: 7.5, y: 0.65, w: 5.6, h: 6.35 });

  const findings = [
    { text: "1.  Extensive transcriptional dysregulation across all subtypes", bold: true, bullet: false },
    { text: "8,688–9,693 significant DEGs per subtype vs healthy breast (GTEx). NAT vs GTEx shows 5,689 field cancerization genes — normal-adjacent tissue is not truly normal.", bold: false, size: 9.5, color: GRAY, bullet: false },
    { text: " ", bullet: false },
    { text: "2.  Compact, stable subtype-discriminating signatures", bold: true, bullet: false },
    { text: "27–36 genes per subtype at freq ≥ 0.80 bootstrap threshold. Suitable for external validation on matched RNA-seq data.", bold: false, size: 9.5, color: GRAY, bullet: false },
    { text: " ", bullet: false },
    { text: "3.  Biologically coherent pathway activation patterns", bold: true, bullet: false },
    { text: "58–93 Hallmark/KEGG pathways per subtype per contrast. LumA/LumB: estrogen response and metabolic programs. Basal: DNA damage, G2M, MYC targets.", bold: false, size: 9.5, color: GRAY, bullet: false },
    { text: " ", bullet: false },
    { text: "4.  Co-expression network confirms subtype identity structure", bold: true, bullet: false },
    { text: "5 WGCNA modules with 20/24 significant PAM50 trait associations. Network-based view corroborates pairwise DEG results independently.", bold: false, size: 9.5, color: GRAY, bullet: false },
    { text: " ", bullet: false },
    { text: "5.  Null survival result is methodologically expected", bold: true, bullet: false },
    { text: "Classification signatures are not prognosis signatures. Basal trend (HR=0.70) warrants follow-up in a powered TNBC cohort.", bold: false, size: 9.5, color: GRAY, bullet: false },
  ];

  s.addText(
    findings.map(f => ({ text: f.text, options: { bold: f.bold || false, fontSize: f.size || 11, color: f.color || DARK, bullet: f.bullet !== false } })),
    { x: 0.3, y: 0.65, w: 7.0, h: 6.35, valign: "top", autoFit: true }
  );
}

// ─── Save ─────────────────────────────────────────────────────────────────────
pptx.writeFile({ fileName: OUT }).then(() => {
  console.log("Saved:", OUT);
});
