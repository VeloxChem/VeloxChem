## 3pa: caffeine

Measured at `fc9b87054`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-17.

| functional | basis | nao | naux | method | 3PA (s) | speedup | SCF (s) | circular strengths (a.u.) | max rel |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HF | def2-svp | 246 | 1242 | four-centre | 432.18 | 1.00 | 12.43 | 50311.45, 327.17 | -- |
|  |  |  |  | RI-JK simd, in memory | 51.68 | 8.36 | 1.13 | 50308.72, 327.06 | 8.3e-04 |
| HF | def2-svpd | 366 | 1242 | four-centre | 1941.32 | 1.00 | 50.58 | 56742.96, 269.20 | -- |
|  |  |  |  | RI-JK simd, in memory | 111.25 | 17.45 | 2.42 | 56729.28, 269.38 | 1.8e-03 |
| HF | def2-tzvp | 494 | 1242 | four-centre | 6555.62 | 1.00 | 174.25 | 57763.86, 275.39 | -- |
|  |  |  |  | RI-JK simd, in memory | 218.09 | 30.06 | 4.44 | 57756.27, 275.32 | 2.1e-03 |
| B3LYP | def2-svp | 246 | 1242 | four-centre | 428.82 | 1.00 | 14.90 | 242960.51, 5753.01 | -- |
|  |  |  |  | RI-JK simd, in memory | 128.66 | 3.33 | 4.01 | 242946.19, 5753.03 | 1.6e-03 |
| B3LYP | def2-svpd | 366 | 1242 | four-centre | 1746.04 | 1.00 | 56.02 | 257814.68, 5646.05 | -- |
|  |  |  |  | RI-JK simd, in memory | 311.15 | 5.61 | 8.69 | 257812.15, 5645.19 | 8.0e-04 |
| B3LYP | def2-tzvp | 494 | 1242 | four-centre | 5400.94 | 1.00 | 179.75 | 239335.95, 5509.01 | -- |
|  |  |  |  | RI-JK simd, in memory | 468.60 | 11.53 | 13.35 | 239322.08, 5507.69 | 1.0e-03 |
