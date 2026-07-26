# Sub-20 uM Detectability Simulation

Training dataset: 8170 raw TMTMosaic data.

Current design: [np.float64(20.0), np.float64(40.0), np.float64(80.0), np.float64(160.0), np.float64(320.0), np.float64(640.0)] uM.
Extended design: [np.float64(1.0), np.float64(2.5), np.float64(5.0), np.float64(10.0), np.float64(20.0), np.float64(40.0), np.float64(80.0), np.float64(160.0), np.float64(320.0), np.float64(640.0)] uM.

Simulations per EC50/log2FC/signal/design cell: 500.

Detected sub-20 uM fitted EC50 was defined as fitted pEC50 > 4.699.

## Summary by true EC50

```text
          design  true_EC50_uM    n  fit_success_rate  sub20_fit_recovery_rate  significant_detection_rate  strict_detection_rate  sub20_strict_recovery_rate  median_observed_range_log2FC
current_20_640uM             1 4500          0.995111                 0.344222                    0.000000               0.000000                    0.000000                      0.232565
current_20_640uM             5 4500          0.996444                 0.353111                    0.000000               0.000000                    0.000000                      0.234667
current_20_640uM            10 4500          0.995556                 0.374000                    0.000000               0.000000                    0.000000                      0.244634
current_20_640uM            15 4500          0.999333                 0.508222                    0.000000               0.000000                    0.000000                      0.343440
current_20_640uM            20 4500          1.000000                 0.592889                    0.000000               0.000000                    0.000000                      0.771579
current_20_640uM            40 4500          1.000000                 0.023556                    0.000000               0.000000                    0.000000                      1.947550
extended_1_640uM             1 4500          1.000000                 0.980889                    0.499333               0.000000                    0.000000                      0.806045
extended_1_640uM             5 4500          1.000000                 0.999556                    0.972444               0.178889                    0.178889                      2.125383
extended_1_640uM            10 4500          1.000000                 0.996000                    0.991111               0.307333                    0.307333                      2.161995
extended_1_640uM            15 4500          1.000000                 0.982222                    0.995111               0.379556                    0.379333                      2.171547
extended_1_640uM            20 4500          1.000000                 0.423556                    0.994000               0.402000                    0.170889                      2.168334
extended_1_640uM            40 4500          1.000000                 0.000000                    0.996000               0.504444                    0.000000                      2.177274
```

Interpretation: this simulation estimates detectability under the experimental design. It cannot identify specific real sites with missed sub-20 uM EC50 values.