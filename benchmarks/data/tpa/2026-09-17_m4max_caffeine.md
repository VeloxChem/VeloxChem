## tpa: caffeine

Measured at `f899e8c78`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-17.

| functional | basis | nao | naux | method | TPA (s) | speedup | SCF (s) | circular strengths (a.u.) | max rel |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HF | def2-svp | 246 | 1242 | four-centre | 154.08 | 1.00 | 12.41 | 33.42, 0.17 | -- |
|  |  |  |  | RI-JK simd, in memory | 15.68 | 9.83 | 1.12 | 33.42, 0.17 | 1.9e-03 |
| HF | def2-svpd | 366 | 1242 | four-centre | 713.21 | 1.00 | 50.38 | 43.39, 0.20 | -- |
|  |  |  |  | RI-JK simd, in memory | 34.74 | 20.53 | 2.42 | 43.39, 0.20 | 1.3e-03 |
| HF | def2-tzvp | 494 | 1242 | four-centre | 2362.74 | 1.00 | 173.57 | 39.69, 0.21 | -- |
|  |  |  |  | RI-JK simd, in memory | 65.06 | 36.32 | 4.45 | 39.69, 0.21 | 8.9e-04 |
| B3LYP | def2-svp | 246 | 1242 | four-centre | 150.87 | 1.00 | 14.89 | 116.89, 0.92 | -- |
|  |  |  |  | RI-JK simd, in memory | 41.38 | 3.65 | 4.04 | 116.89, 0.92 | 6.7e-04 |
| B3LYP | def2-svpd | 366 | 1242 | four-centre | 626.89 | 1.00 | 56.22 | 129.62, 0.53 | -- |
|  |  |  |  | RI-JK simd, in memory | 101.27 | 6.19 | 8.66 | 129.62, 0.53 | 1.2e-03 |
| B3LYP | def2-tzvp | 494 | 1242 | four-centre | 1946.67 | 1.00 | 180.24 | 120.29, 0.58 | -- |
|  |  |  |  | RI-JK simd, in memory | 150.07 | 12.97 | 13.33 | 120.29, 0.58 | 2.1e-03 |
