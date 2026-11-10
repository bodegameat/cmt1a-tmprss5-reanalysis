# Preclinical Nerve Transcript Pharmacodynamics vs. Cutaneous Null of TMPRSS5 (Spinesin) in CMT1A

**Joshua Bryant-Gainer**  
November 2026

---

## Abstract

Charcot–Marie–Tooth disease type 1A (CMT1A) is driven by *PMP22* gene dosage, yet circulating biomarkers that report Schwann-cell pharmacodynamics remain incompletely validated across tissue compartments. Wang et al. (2020) reported plasma TMPRSS5 (spinesin) elevation in CMT1A with a **2.07-fold** increase and discriminative performance of **AUC 0.9131**, nominating TMPRSS5 as a candidate blood protein marker. Reanalysis of the Zhao/Svaren antisense oligonucleotide (ASO) sciatic-nerve RNA-seq cohort (**GSE103799**) demonstrates that *Tmprss5* transcript abundance tracks disease and treatment in peripheral nerve: untreated C22 mice exhibit a **+2.90 log2 fold-change** versus wild-type, ASO-mediated *Pmp22* suppression yields a **−1.26 log2 fold-change** relative to untreated C22, and the disease-associated TPM excess is reduced by **67.40%**. In contrast, the Svaren 2019 Dryad cutaneous (skin) microarray resource shows a cutaneous **null** for TMPRSS5 (**fold-change 1.00**, **p = 0.99**, ranking **#52 of 54** probes), indicating that skin transcript abundance does not recapitulate the nerve or plasma signal. Together, these compartments establish a coherent but incomplete biomarker cascade: nerve transcript pharmacodynamics are strong, circulating protein is elevated, and cutaneous mRNA is uninformative—motivating prospective plasma-to-nerve electrophysiology pairing studies.

---

## Multi-Compartment Summary

| Compartment | Study / Accession | Analyte | Key Quantitative Result | Interpretation |
|---|---|---|---|---|
| Sciatic nerve mRNA | GSE103799 (Zhao/Svaren ASO) | *Tmprss5* TPM | Disease **+2.90** log2FC; Treatment **−1.26** log2FC; **67.40%** rescue | Strong preclinical nerve pharmacodynamic marker |
| Cutaneous nerve / skin mRNA | Svaren 2019 (Dryad) | TMPRSS5 microarray probe | FC **1.00**, p **0.99**, rank **#52 / 54** | Cutaneous null; skin does not report disease |
| Circulating plasma protein | Wang et al., 2020 | TMPRSS5 (spinesin) protein | **2.07-fold** elevation; AUC **0.9131** | Promising blood biomarker, tissue source still inferred |

### GSE103799 group statistics (TPM)

| Group | n | Mean TPM | SD |
|---|---:|---:|---:|
| WT_PBS | 8 | 11.17 | 3.29 |
| C22_PBS | 8 | 83.60 | 39.01 |
| C22_ASO | 8 | 34.78 | 6.67 |

**Contrasts:** Disease Δ = **72.43 TPM**; ASO Δ = **48.82 TPM**; Rescue = **67.40%**.  
**Welch t-tests:** C22_PBS vs WT_PBS, **p = 0.00116**; C22_ASO vs C22_PBS, **p = 0.00926**.

---

## Mechanistic & Translational Implications

The 67.4% reduction in sciatic-nerve *Tmprss5* after *Pmp22*-directed ASO treatment places TMPRSS5 downstream of Schwann-cell gene-dosage pathology and within the dynamic range of a treatment-responsive transcript. The cutaneous null (FC 1.00, p = 0.99) argues against skin biopsy mRNA as a surrogate for this axis and cautions against assuming dermal expression mirrors endoneurial biology. Plasma protein elevation (Wang 2020) remains the most clinically accessible compartment, yet the missing link is whether circulating TMPRSS5 quantitatively tracks nerve electrophysiology and ASO engagement in the same individuals over time. Future studies should therefore prioritize **paired plasma TMPRSS5 protein assays with nerve conduction / electromyography endpoints** (and, where feasible, nerve biopsy or imaging surrogates) to establish whether blood levels report the same pharmacodynamic window observed at **+2.90 / −1.26 log2FC** in GSE103799.

---

## Reproducibility

All GSE103799 TMPRSS5 TPM values used here are distributed in:

```text
data/gse103799_tmprss5_tpm.csv
```

Recompute group means, deltas, log2 fold-changes, and Welch t-tests with:

```bash
python scripts/recompute_tmprss5.py
```

Dependencies: Python 3 standard library (`csv`, `math`) and `scipy` (`scipy.stats`).

Expected headline outputs include WT_PBS mean **11.17** (SD **3.29**), C22_PBS mean **83.60** (SD **39.01**), C22_ASO mean **34.78** (SD **6.67**), disease Δ **72.43 TPM**, ASO Δ **48.82 TPM**, rescue **67.40%**, disease log2FC **+2.90**, treatment log2FC **−1.26**, and Welch p-values **0.00116** / **0.00926**.

---

## Site & Deployment

- Static micro-site: `index.html` (Chart.js CDN)
- Netlify publish root: `.` (see `netlify.toml`)
- Repository name (intended): `cmt1a-tmprss5-reanalysis`
