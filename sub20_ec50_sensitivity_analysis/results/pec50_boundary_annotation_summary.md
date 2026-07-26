# pEC50 Boundary Annotation Summary

Measured concentration range: 20-640 uM.

Lower measured boundary: 20 uM -> pEC50 = 4.699.
Upper measured boundary: 640 uM -> pEC50 = 3.194.

Classification rules:

- pEC50 > 4.699: lower-boundary-limited; fitted EC50 is below the lowest measured concentration.
- 3.194 <= pEC50 <= 4.699: fitted EC50 lies within the measured 20-640 uM range.
- pEC50 < 3.194: fitted EC50 is above the highest measured concentration or weak/out-of-range.

## Counts

```text
dataset           level                       ec50_range_status  count   percent
   8171            site              in_measured_range_20_640uM  13214 97.054719
   8171            site upper_boundary_or_weak_EC50_above_640uM    401  2.945281
   8171 gene_best_pEC50              in_measured_range_20_640uM   4958 98.333994
   8171 gene_best_pEC50 upper_boundary_or_weak_EC50_above_640uM     84  1.666006
   8170            site              in_measured_range_20_640uM   5445 92.617792
   8170            site upper_boundary_or_weak_EC50_above_640uM    434  7.382208
   8170 gene_best_pEC50              in_measured_range_20_640uM   2497 96.707978
   8170 gene_best_pEC50 upper_boundary_or_weak_EC50_above_640uM     85  3.292022
```

Interpretation note: lower-boundary-limited curves should not be interpreted as precise sub-20 uM potency estimates.