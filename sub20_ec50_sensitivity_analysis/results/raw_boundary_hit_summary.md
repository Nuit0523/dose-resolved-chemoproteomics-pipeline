# Raw Data + CurveCurator Boundary Hit Analysis

20 uM lower boundary pEC50 = 4.698970.
640 uM upper boundary pEC50 = 3.193820.

A real-data site was counted as detected below 20 uM only if fitted pEC50 > 4.699.

```text
dataset                                        metric  site_count  gene_count  percent_sites
   8170                           detected_sub20_EC50           0           0       0.000000
   8170                             strict_hit_absFC2           0           0       0.000000
   8170                        strict_hit_positiveFC1           0           0       0.000000
   8170                       sub20_strict_hit_absFC2           0           0       0.000000
   8170                  sub20_strict_hit_positiveFC1           0           0       0.000000
   8170        near_lower_boundary_pEC50_4p5_to_4p699         436         352       4.056946
   8170 range_status::fitted_EC50_above_640uM_or_weak         434         376       4.038336
   8170         range_status::fitted_EC50_in_20_640uM        6405        2497      59.598027
   8170                         range_status::missing        3908        1839      36.363636
   8171                           detected_sub20_EC50           0           0       0.000000
   8171                             strict_hit_absFC2         220         189       1.435001
   8171                        strict_hit_positiveFC1          75          74       0.489205
   8171                       sub20_strict_hit_absFC2           0           0       0.000000
   8171                  sub20_strict_hit_positiveFC1           0           0       0.000000
   8171        near_lower_boundary_pEC50_4p5_to_4p699         758         681       4.944231
   8171 range_status::fitted_EC50_above_640uM_or_weak         401         379       2.615615
   8171         range_status::fitted_EC50_in_20_640uM       13214        4958      86.191377
   8171                         range_status::missing        1716        1363      11.193008
```

Important: this analysis counts detected below-range fitted EC50 values. It cannot identify true missed sub-20 uM sites from raw data alone.