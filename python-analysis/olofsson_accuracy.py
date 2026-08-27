"""
Olofsson et al. (2014) area-adjusted accuracy assessment.
Rangamati Deforestation Monitoring — 2023 epoch, corrected spatial-block RF classifier.

Reference: Olofsson, P., Foody, G.M., Herold, M., Stehman, S.V., Woodcock, C.E.,
Wulder, M.A. (2014). Good practices for estimating area and assessing accuracy
of land change. Remote Sensing of Environment, 148, 42-57.

INPUTS:
  1. Mapped class areas (hectares) from Script 07 (corrected spatial-block RF).
  2. Confusion matrix built from the 200-point independent validation sample,
     reclassified with the corrected classifier (Script 06b).

IMPORTANT CAVEAT (read before using these numbers in the paper):
The 200 validation points were originally drawn as an equal-allocation
stratified sample (~40 points/class) using the OLD (pre-fix) classifier's
map as the stratification variable. The Olofsson estimator below groups
those same points by the NEW (corrected) classifier's predicted class to
compute stratum sample sizes (n_i.) and weights (W_i from the new mapped
areas). This is the standard practical approach when reusing an existing
reference sample against an updated map, but it is technically a departure
from a "clean" single-phase stratified design (the reference sample was not
literally drawn stratified by the map being assessed here). This should be
disclosed as a limitation in the paper's accuracy-assessment subsection.
"""

import pandas as pd
import numpy as np

# ---------------------------------------------------------------------------
# 1. Mapped class areas (hectares) — Script 07, corrected spatial-block RF
# ---------------------------------------------------------------------------
class_names = {
    1: 'Dense Forest',
    2: 'Degraded/Jhum',
    3: 'Water',
    4: 'Agriculture/Settlement',
    5: 'Bare Land',
}
classes = [1, 2, 3, 4, 5]

mapped_area_ha = {
    1: 396367.35,
    2: 99639.57,
    3: 30623.85,
    4: 35083.75,
    5: 14819.86,
}
total_area_ha = sum(mapped_area_ha.values())
W = {i: mapped_area_ha[i] / total_area_ha for i in classes}

# ---------------------------------------------------------------------------
# 2. Load validation points and build the confusion matrix
#    rows = map class (rf_class_new), columns = reference class (manual_class)
# ---------------------------------------------------------------------------
df = pd.read_csv('ValidationPoints_200_Reclassified_v2.csv')

n = pd.crosstab(df['rf_class_new'], df['manual_class'])
n = n.reindex(index=classes, columns=classes, fill_value=0)

print("=== Confusion matrix (rows=map class, cols=reference class) ===")
print(n.rename(index=class_names, columns=class_names))
print()

n_i_dot = n.sum(axis=1)          # row totals: sample size per map-class stratum
print("Sample size per stratum (n_i.):")
for i in classes:
    print(f"  {class_names[i]}: {n_i_dot[i]}")
print()

# ---------------------------------------------------------------------------
# 3. Olofsson area-weighted proportion matrix: p_ij = W_i * (n_ij / n_i.)
# ---------------------------------------------------------------------------
p = pd.DataFrame(index=classes, columns=classes, dtype=float)
for i in classes:
    for j in classes:
        p.loc[i, j] = W[i] * (n.loc[i, j] / n_i_dot[i]) if n_i_dot[i] > 0 else 0.0

p_dot_j = p.sum(axis=0)   # column sums (used for producer's accuracy & area)
p_i_dot = p.sum(axis=1)   # row sums (== W_i, used for user's accuracy)

# ---------------------------------------------------------------------------
# 4. Overall accuracy, user's accuracy, producer's accuracy
# ---------------------------------------------------------------------------
OA = sum(p.loc[i, i] for i in classes)

UA = {i: (p.loc[i, i] / p_i_dot[i]) if p_i_dot[i] > 0 else np.nan for i in classes}
PA = {j: (p.loc[j, j] / p_dot_j[j]) if p_dot_j[j] > 0 else np.nan for j in classes}

# ---------------------------------------------------------------------------
# 5. Standard errors (Olofsson et al. 2014, eqs. 5-7)
# ---------------------------------------------------------------------------
# SE of overall accuracy
se_OA = np.sqrt(sum(
    (W[i] ** 2) * (UA[i] * (1 - UA[i])) / (n_i_dot[i] - 1)
    if n_i_dot[i] > 1 else 0.0
    for i in classes
))

# SE of user's accuracy (per class i)
se_UA = {}
for i in classes:
    if n_i_dot[i] > 1:
        se_UA[i] = np.sqrt(UA[i] * (1 - UA[i]) / (n_i_dot[i] - 1))
    else:
        se_UA[i] = np.nan

# SE of producer's accuracy (per class j) — eq. 7, requires N_i (pixel counts);
# we approximate N_i with mapped area (ha), which is proportional to pixel
# count at fixed resolution and is valid for this ratio-based formula.
se_PA = {}
for j in classes:
    if p_dot_j[j] == 0:
        se_PA[j] = np.nan
        continue
    term1 = (mapped_area_ha[j] ** 2) * ((1 - PA[j]) ** 2) * UA[j] * (1 - UA[j]) / (n_i_dot[j] - 1) if n_i_dot[j] > 1 else 0.0
    term2 = sum(
        (mapped_area_ha[i] ** 2) * (n.loc[i, j] / n_i_dot[i]) * (1 - n.loc[i, j] / n_i_dot[i]) / (n_i_dot[i] - 1)
        for i in classes if i != j and n_i_dot[i] > 1
    )
    se_PA[j] = (1 / mapped_area_ha[j]) * np.sqrt(term1 + (PA[j] ** 2) * term2)

# ---------------------------------------------------------------------------
# 6. Area-adjusted class area estimates + 95% CI (eq. 9-10)
# ---------------------------------------------------------------------------
se_p_dot_j = {}
for j in classes:
    se_p_dot_j[j] = np.sqrt(sum(
        (W[i] ** 2) * (n.loc[i, j] / n_i_dot[i]) * (1 - n.loc[i, j] / n_i_dot[i]) / (n_i_dot[i] - 1)
        if n_i_dot[i] > 1 else 0.0
        for i in classes
    ))

adjusted_area_ha = {j: p_dot_j[j] * total_area_ha for j in classes}
adjusted_area_se_ha = {j: se_p_dot_j[j] * total_area_ha for j in classes}
adjusted_area_ci95 = {j: 1.96 * adjusted_area_se_ha[j] for j in classes}

# ---------------------------------------------------------------------------
# 7. Report
# ---------------------------------------------------------------------------
print("=== Area-adjusted accuracy (Olofsson et al. 2014) ===")
print(f"Overall Accuracy: {OA*100:.2f}%  (SE = {se_OA*100:.2f} pp, 95% CI +/- {1.96*se_OA*100:.2f} pp)")
print()

print(f"{'Class':<24}{'UA %':>8}{'UA SE':>8}{'PA %':>8}{'PA SE':>8}")
for i in classes:
    print(f"{class_names[i]:<24}{UA[i]*100:>8.2f}{se_UA[i]*100:>8.2f}{PA[i]*100:>8.2f}{se_PA[i]*100:>8.2f}")
print()

print(f"{'Class':<24}{'Mapped ha':>12}{'Adjusted ha':>14}{'+/- 95% CI':>12}")
for j in classes:
    print(f"{class_names[j]:<24}{mapped_area_ha[j]:>12.0f}{adjusted_area_ha[j]:>14.0f}{adjusted_area_ci95[j]:>12.0f}")

print()
print(f"Total study area: {total_area_ha:.2f} ha")
