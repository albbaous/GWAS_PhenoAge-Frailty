# PhenoAge and Frailty Indices

This repository contains code and workflows for calculating **PhenoAge** and the **Electronic Frailty Index (eFI-2)** from UK Biobank data.

---

## 1. PhenoAge

**PhenoAge** translates an individual's mortality risk—derived from a Gompertz proportional hazards model that combines chronological age and nine clinical biomarkers—into biological age (in years):

```
PhenoAge = 141.50225 + ln(-0.00553 * ln(1 - mortality_risk)) / 0.090165
```

The mortality risk is computed using the Gompertz model:

```
mortality_risk = 1 - exp(-exp(xb) * ((exp(120 * gamma) - 1) / gamma))
```

Where:

- `gamma = 0.0076927` (Gompertz hazard slope)  
- `xb` = linear predictor = `-19.9067 + sum(w_i * x_i)`  
- `w_i` = biomarker weights  
- `x_i` = biomarker values (including chronological age)

Here, `w_i` are the respective biomarker weights and `x_i` are the biomarker values, including chronological age.

---

## 2. Electronic Frailty Index (eFI-2)

The **electronic Frailty Index version 2 (eFI-2)** quantifies frailty as a weighted sum of specific health deficits:

```
eFI-2 = sum(w_i * d_i)
```

Where:

- `d_i` = `1` if the i-th deficit is present, `0` otherwise  
- `w_i` = coefficient for the i-th deficit  
- `n` = total number of deficits

Frailty categories:

| Category | Threshold |
|----------|-----------|
| Robust   | eFI-2 ≤ 0.0857 |
| Mild     | 0.0857 < eFI-2 ≤ 0.1624 |
| Moderate | 0.1624 < eFI-2 ≤ 0.2392 |
| Severe   | 0.2392 < eFI-2 < 1 |

### Risk prediction

The 1-year risk of the composite outcome can be calculated as:

```
1-year risk = 0.0151 * exp(LP)
```

Where the **linear predictor (LP)** is:

```
LP = sum(beta_i * d_i)
```

Or approximately:

```
LP ≈ 8.429 * eFI-2 score
```
- `beta_i` = coefficient for the i-th deficit  
- `d_i` = presence (1) or absence (0) of deficit
