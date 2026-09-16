## rpa: caffeine

Measured at `29419eceb`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-16.

| functional | basis | nao | naux | method | RPA (s) | speedup | iter | s/iter | SCF (s) | max dE (a.u.) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HF | def2-svp | 246 | 1242 | four-centre | 58.19 | 1.00 | 17 | 3.42 | 12.51 | -- |
|  |  |  |  | RI-JK simd, in memory | 5.32 | 10.94 | 17 | 0.31 | 1.15 | 2.3e-05 |
| HF | def2-svpd | 366 | 1242 | four-centre | 281.97 | 1.00 | 19 | 14.84 | 50.70 | -- |
|  |  |  |  | RI-JK simd, in memory | 11.75 | 24.00 | 19 | 0.62 | 2.48 | 1.9e-05 |
| HF | def2-tzvp | 494 | 1242 | four-centre | 935.20 | 1.00 | 19 | 49.22 | 174.86 | -- |
|  |  |  |  | RI-JK simd, in memory | 20.89 | 44.76 | 19 | 1.10 | 4.49 | 1.5e-05 |
| HF | def2-tzvpd | 614 | 1242 | four-centre | 2666.00 | 1.00 | 22 | 121.18 | 424.60 | -- |
|  |  |  |  | RI-JK simd, in memory | 39.17 | 68.06 | 22 | 1.78 | 7.23 | 1.3e-05 |
| B3LYP | def2-svp | 246 | 1242 | four-centre | 49.25 | 1.00 | 11 | 4.48 | 15.15 | -- |
|  |  |  |  | RI-JK simd, in memory | 12.94 | 3.81 | 11 | 1.18 | 4.13 | 6.5e-06 |
| B3LYP | def2-svpd | 366 | 1242 | four-centre | 205.82 | 1.00 | 11 | 18.71 | 56.24 | -- |
|  |  |  |  | RI-JK simd, in memory | 31.90 | 6.45 | 11 | 2.90 | 8.78 | 5.9e-06 |
| B3LYP | def2-tzvp | 494 | 1242 | four-centre | 656.47 | 1.00 | 12 | 54.71 | 181.17 | -- |
|  |  |  |  | RI-JK simd, in memory | 46.16 | 14.22 | 12 | 3.85 | 13.52 | 5.7e-06 |
| B3LYP | def2-tzvpd | 614 | 1242 | four-centre | 1639.05 | 1.00 | 11 | 149.00 | 433.46 | -- |
|  |  |  |  | RI-JK simd, in memory | 80.01 | 20.49 | 11 | 7.27 | 22.70 | 4.5e-06 |
