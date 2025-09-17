# GWAS_PhenoAge-Frailty
# Phenotype and Frailty Indices

This repository contains code and workflows for calculating **PhenoAge** and the **Electronic Frailty Index (eFI-2)** from UK Biobank data.

---

## 1. PhenoAge

**PhenoAge** translates an individual's mortality risk—derived from a Gompertz proportional hazards model that combines chronological age and nine clinical biomarkers—into biological age (in years):

\[
\text{PhenoAge} = 141.50225 + \frac{\ln\left(-0.00553 \times \ln(1 - \text{mortality\_risk})\right)}{0.090165}
\]

The mortality risk is computed using the Gompertz model:

\[
\text{mortality\_risk} = 1 - \exp \Big( -\exp(xb) \times \frac{\exp(120 \cdot \gamma) - 1}{\gamma} \Big)
\]

Where:

- \(\gamma = 0.0076927\) (Gompertz hazard slope)
- \(xb\) is the linear predictor:

\[
xb = -19.9067 + \sum_{i=1}^{9} w_i \cdot x_i
\]

Here, \(w_i\) are the respective biomarker weights and \(x_i\) are the biomarker values, including chronological age.

---

## 2. Electronic Frailty Index (eFI-2)

The **electronic Frailty Index version 2 (eFI-2)** quantifies frailty as a weighted sum of specific health deficits:

\[
\text{eFI-2} = \sum_{i=1}^{n} w_i \cdot d_i
\]

Where:

- \(d_i\) = 1 if the i-th deficit is present, 0 otherwise  
- \(w_i\) = coefficient representing the relative contribution of the i-th deficit  
- \(n\) = total number of deficits included

Frailty categories:

| Category | Threshold |
|----------|-----------|
| Robust   | eFI-2 ≤ 0.0857 |
| Mild     | 0.0857 < eFI-2 ≤ 0.1624 |
| Moderate | 0.1624 < eFI-2 ≤ 0.2392 |
| Severe   | 0.2392 < eFI-2 < 1 |

### Risk prediction

The 1-year risk of the composite outcome can be calculated as:

\[
\text{1-year risk} = 0.0151 \times \exp(\text{LP})
\]

Where the **linear predictor (LP)** is:

\[
\text{LP} = \sum_{i=1}^{n} \beta_i \cdot d_i
\]

Or approximately:

\[
\text{LP} \approx 8.429 \times \text{eFI-2 score}
\]

Here, \(\beta_i\) denotes the coefficient for the i-th deficit, and \(d_i\) indicates presence (1) or absence (0).

---
