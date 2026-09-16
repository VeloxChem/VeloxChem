## gradient: caffeine

Measured at `15b33b62e`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank, veloxchem 1.0rc4, 2026-09-16.

| functional | basis | nao | naux | method | gradient (s) | speedup | SCF (s) | vs four-centre |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HF | def2-svp | 246 | 1242 | four-centre | 8.61 | 1.00 | 12.38 | -- |
|  |  |  |  | RI-JK simd, in memory | 2.22 | 3.87 | 1.11 | 7.7e-05 |
| HF | def2-svpd | 366 | 1242 | four-centre | 41.68 | 1.00 | 50.26 | -- |
|  |  |  |  | RI-JK simd, in memory | 2.91 | 14.33 | 2.45 | 7.7e-05 |
| HF | def2-tzvp | 494 | 1242 | four-centre | 136.34 | 1.00 | 173.64 | -- |
|  |  |  |  | RI-JK simd, in memory | 4.13 | 33.00 | 4.43 | 7.1e-05 |
| HF | def2-tzvpd | 614 | 1242 | four-centre | 352.64 | 1.00 | 423.08 | -- |
|  |  |  |  | RI-JK simd, in memory | 5.46 | 64.61 | 7.16 | 7.3e-05 |
| B3LYP | def2-svp | 246 | 1242 | four-centre | 9.26 | 1.00 | 14.97 | -- |
|  |  |  |  | RI-JK simd, in memory | 2.78 | 3.33 | 4.04 | 1.6e-05 |
| B3LYP | def2-svpd | 366 | 1242 | four-centre | 42.80 | 1.00 | 55.86 | -- |
|  |  |  |  | RI-JK simd, in memory | 4.09 | 10.47 | 8.72 | 1.7e-05 |
| B3LYP | def2-tzvp | 494 | 1242 | four-centre | 138.23 | 1.00 | 179.76 | -- |
|  |  |  |  | RI-JK simd, in memory | 5.75 | 24.06 | 13.45 | 6.2e-05 |
| B3LYP | def2-tzvpd | 614 | 1242 | four-centre | 356.52 | 1.00 | 428.84 | -- |
|  |  |  |  | RI-JK simd, in memory | 7.92 | 45.04 | 22.49 | 6.3e-05 |

### Scaling, cost proportional to nao to the p

| functional | method | p |
| --- | --- | ---: |
| HF | four-centre | 4.04 |
| HF | RI-JK simd | 0.98 |
| B3LYP | four-centre | 3.97 |
| B3LYP | RI-JK simd | 1.13 |
