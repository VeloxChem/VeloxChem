## tda: caffeine

Measured at `dfc9f00c6`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-16.

| functional | basis | nao | naux | method | TDA (s) | speedup | iter | s/iter | SCF (s) | max dE (a.u.) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HF | def2-svp | 246 | 1242 | four-centre | 58.50 | 1.00 | 15 | 3.90 | 12.46 | -- |
|  |  |  |  | RI-JK simd, in memory | 3.37 | 17.34 | 15 | 0.22 | 1.12 | 2.2e-05 |
| HF | def2-svpd | 366 | 1242 | four-centre | 289.93 | 1.00 | 17 | 17.05 | 50.89 | -- |
|  |  |  |  | RI-JK simd, in memory | 7.81 | 37.13 | 17 | 0.46 | 2.45 | 1.7e-05 |
| HF | def2-tzvp | 494 | 1242 | four-centre | 914.97 | 1.00 | 16 | 57.19 | 174.79 | -- |
|  |  |  |  | RI-JK simd, in memory | 13.16 | 69.52 | 16 | 0.82 | 4.46 | 1.3e-05 |
| HF | def2-tzvpd | 614 | 1242 | four-centre | 2554.92 | 1.00 | 18 | 141.94 | 424.33 | -- |
|  |  |  |  | RI-JK simd, in memory | 23.22 | 110.03 | 18 | 1.29 | 7.13 | 7.4e-06 |
| B3LYP | def2-svp | 246 | 1242 | four-centre | 51.64 | 1.00 | 10 | 5.16 | 15.21 | -- |
|  |  |  |  | RI-JK simd, in memory | 12.24 | 4.22 | 10 | 1.22 | 4.28 | 6.0e-06 |
| B3LYP | def2-svpd | 366 | 1242 | four-centre | 211.23 | 1.00 | 10 | 21.12 | 57.31 | -- |
|  |  |  |  | RI-JK simd, in memory | 29.41 | 7.18 | 10 | 2.94 | 8.89 | 5.9e-06 |
| B3LYP | def2-tzvp | 494 | 1242 | four-centre | 646.65 | 1.00 | 10 | 64.66 | 183.91 | -- |
|  |  |  |  | RI-JK simd, in memory | 40.26 | 16.06 | 10 | 4.03 | 13.46 | 5.6e-06 |
| B3LYP | def2-tzvpd | 614 | 1242 | four-centre | 1509.78 | 1.00 | 9 | 167.75 | 446.54 | -- |
|  |  |  |  | RI-JK simd, in memory | 64.36 | 23.46 | 9 | 7.15 | 22.64 | 4.1e-06 |
