## gradient: caffeine

Measured at `f8c24e0ee`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-17.

| functional | basis | nao | naux | method | gradient (s) | speedup | SCF (s) | vs four-centre |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HF | def2-svp | 246 | 1242 | four-centre | 21.17 | 1.00 | 73.02 | -- |
|  |  |  |  | RI-JK simd, in memory | 0.84 | 25.32 | 3.91 | 1.0e-04 |
| HF | def2-svpd | 366 | 1242 | four-centre | 94.75 | 1.00 | 299.57 | -- |
|  |  |  |  | RI-JK simd, in memory | 1.59 | 59.55 | 8.47 | 9.8e-05 |
| HF | def2-tzvp | 494 | 1242 | four-centre | 333.08 | 1.00 | 1054.83 | -- |
|  |  |  |  | RI-JK simd, in memory | 2.89 | 115.29 | 15.48 | 8.0e-05 |
| HF | def2-tzvpd | 614 | 1242 | four-centre | 834.48 | 1.00 | 2546.53 | -- |
|  |  |  |  | RI-JK simd, in memory | 4.41 | 189.23 | 25.91 | 8.2e-05 |
| B3LYP | def2-svp | 246 | 1242 | four-centre | 22.19 | 1.00 | 56.10 | -- |
|  |  |  |  | RI-JK simd, in memory | 1.88 | 11.79 | 9.28 | 1.5e-05 |
| B3LYP | def2-svpd | 366 | 1242 | four-centre | 96.79 | 1.00 | 218.52 | -- |
|  |  |  |  | RI-JK simd, in memory | 3.88 | 24.97 | 23.17 | 1.8e-05 |
| B3LYP | def2-tzvp | 494 | 1242 | four-centre | 335.82 | 1.00 | 794.56 | -- |
|  |  |  |  | RI-JK simd, in memory | 5.84 | 57.47 | 33.14 | 6.4e-05 |
| B3LYP | def2-tzvpd | 614 | 1242 | four-centre | 842.87 | 1.00 | 1912.56 | -- |
|  |  |  |  | RI-JK simd, in memory | 9.37 | 89.93 | 56.75 | 6.5e-05 |

### Scaling, cost proportional to nao to the p

| functional | method | p |
| --- | --- | ---: |
| HF | four-centre | 4.02 |
| HF | RI-JK simd | 1.82 |
| B3LYP | four-centre | 3.98 |
| B3LYP | RI-JK simd | 1.71 |
