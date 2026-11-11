# Tmprss5 in the Zhao 2018 C22 ASO nerve data (GSE103799)

**Reproducible Data Note & Literature Cross-Reference**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24118723.svg)](https://doi.org/10.5281/zenodo.24118723)

Joshua Bryant-Gainer (CMT1A patient and data analyst) – November 2026 – Not peer-reviewed

---

## The question

Does *Tmprss5* mRNA in mouse sciatic nerve differ between C22 and wild-type mice, and does it change after PMP22-lowering ASO treatment? And how does that sit next to the published human skin and plasma TMPRSS5 data?

## What the data show

Sciatic nerve *Tmprss5* TPM, 15-week-old mice, n = 8 per group:

| Group | n | Mean TPM | SD | Range |
|---|---:|---:|---:|---:|
| WT_PBS | 8 | 11.17 | 3.29 | 6.41–15.12 |
| C22_PBS | 8 | 83.60 | 39.01 | 30.04–153.04 |
| C22_ASO | 8 | 34.78 | 6.67 | 25.92–43.36 |
| WT_ASO | 8 | 10.86 | 3.48 | 6.76–16.01 |

- **Untreated C22 vs wild-type:** 7.48x WT (log2FC +2.90), a gap of 72.43 TPM.
- **ASO-treated vs untreated C22:** ratio 0.42 (log2FC −1.27), a shift of −48.82 TPM.
  - Reduction from untreated, (C22_PBS − C22_ASO) / C22_PBS: **58.4%**.
  - Share of the gap over WT closed, (C22_PBS − C22_ASO) / (C22_PBS − WT_PBS): **67.4%**.
  - Treated mice remain **3.11x** WT.
- **ASO control in wild-type:** wild-type + ASO 10.86 vs wild-type + PBS 11.17 TPM (ratio 0.97). Welch t-test p = 0.85907 (raw TPM) and 0.84146 (log2 TPM); Mann-Whitney U p = 0.87848. The ASO did not measurably change Tmprss5 in wild-type nerve.

Two-sided p-values:

| Contrast | Welch t (raw TPM) | Welch t (log2 TPM) | Mann-Whitney U |
|---|---:|---:|---:|
| C22_PBS vs WT_PBS | 0.00116 | 7.7e-07 | 0.00016 |
| C22_ASO vs C22_PBS | 0.00926 | 0.00222 | 0.00466 |
| C22_ASO vs WT_PBS | 3.6e-06 | 1.6e-06 | 0.00016 |
| WT_ASO vs WT_PBS | 0.85907 | 0.84146 | 0.87848 |

This is one gene with no multiple-testing correction.

Two high C22_PBS mice (125.34 and 153.04 TPM) raise the untreated mean, and one C22_PBS mouse (30.04 TPM) overlaps the ASO-treated range. The treatment difference holds in a Mann-Whitney U test (p = 0.0047) and a Welch t-test on log2 TPM (p = 0.0022).

### Methods

- Library prep and sequencing, as stated in GEO: QuantSeq 3' mRNA-Seq Library Prep Kit from Lexogen; NextSeq 500 (platform GPL19057).
- Quantification, as stated in GEO: reads were associated to transcripts using a quasi-mapping based mode of Salmon (version 0.7.2), using gene models from mouse ensembl build 81 (GRCm38.p4).
- Per-sample Salmon transcript quantification files were downloaded from GEO. Gene-level Tmprss5 TPM is the sum of TPM across its Ensembl transcripts (IDs listed in the script). Tmprss5 transcripts identified from the Ensembl release 102 (GRCm38) GTF and the current Ensembl REST lookup.
- GEO check (10 November 2026): all 32 values match the GEO files to 0.01 TPM. GEO sample titles and characteristics (genotype, treatment, 100 mg/kg, age 15 weeks) agree with the group labels.

## Study context

- **Dataset:** GSE103799, from Zhao et al., *J Clin Invest* 2018;128(1):359-368, co-authored by Dr. John Svaren. 3' poly-A mRNA sequencing of sciatic nerve.
- **Subset:** GSE103799 has 43 samples. Thirty-two are 15-week-old animals in four arms of 8: wild-type + PBS (GSM2782061–GSM2782068), C22 + PBS (GSM2782045–GSM2782052), C22 + PMP22 ASO (GSM2782053–GSM2782060), and wild-type + PMP22 ASO (GSM2782069–GSM2782076). The ASO arms received 100 mg/kg weekly by subcutaneous injection for 9 weeks, starting at 5 weeks of age. The main comparison uses the first three arms (24 mice); the wild-type + ASO arm is reported as a control. The series also contains 11 untreated 5-week animals (8 C22, 3 wild-type), not analyzed here.
- **Model:** C22 is a severe model carrying about seven copies of human PMP22; 100 mg/kg was the highest dose tested.
- Zhao et al. reported 1,056 genes differing more than twofold between 15-week C22 and wild-type nerve. A Schwann cell gene moving back toward wild-type as myelination improves is expected and is not specific to Tmprss5.

## Human skin cross-reference

Svaren et al., *Ann Neurol* 2019;85(6):887-898, used a custom 54-gene NanoString panel on skin biopsies. Only the panel's "CMT1A genes" subset was selected from the Zhao rodent ASO data. In the paper's panel table TMPRSS5 is listed with injury-related genes, and the paper's text does not discuss it.

My reading of the deposited Dryad counts (doi:10.5061/dryad.hqbzkh1st; 32 CMT1A vs 11 control biopsies) gives a flat result for TMPRSS5: FC 1.00, p = 0.99, rank 52 of 54 probes. This is my own reading of the deposited counts, not a finding stated in the paper. The code for this skin comparison is not yet included in this repository.

## Human plasma cross-reference

Wang et al., *Ann Clin Transl Neurol* 2020;7(1):69-82, was led by Sanofi with investigators from the Inherited Neuropathies Consortium, including Dr. Shy and Dr. Svaren. Plasma TMPRSS5 protein was about 2.07-fold elevated in CMT1A (47 patients, 41 controls). Only untreated patients were measured, and the authors wrote that further animal efficacy studies were needed to see whether Tmprss5 decreases with treatment.

## Limitations

- Mouse sciatic nerve mRNA does not show that human plasma TMPRSS5 protein will fall with treatment.
- Baseline human plasma TMPRSS5 did not correlate with age, nerve conduction, or clinical severity scores in Wang et al., and its stability within a person over time has not been published.
- One dataset, one model, one dose, one time point, one gene; no multiple-testing correction; no replication.
- TPM from 3' tag sequencing is a relative measure.

## What would answer the real question

Plasma TMPRSS5 measured before and after PMP22-lowering treatment, in animals or trial participants, alongside nerve conduction.

## Reproducing the numbers

```bash
# Group statistics, contrasts, and p-values from the CSV
python3 scripts/recompute_tmprss5.py

# Check every CSV value against the GEO per-sample Salmon files (needs network)
python3 scripts/recompute_tmprss5.py --fetch-geo
```

Result of the GEO check run on 10 November 2026: 32/32 samples OK. All 32 values match the GEO files to 0.01 TPM, and the GEO titles and characteristics agree with the CSV group labels.

Dependencies: Python 3.9+ and `scipy`. The `--fetch-geo` mode uses only the standard library.

Files:

- `data/gse103799_tmprss5_tpm.csv`: the 32 sample values
- `scripts/recompute_tmprss5.py`: statistics and GEO check (`TMPRSS5_TRANSCRIPTS` lists the transcript IDs and their sources)
- `index.html`: static site (Netlify, publish root `.`)
