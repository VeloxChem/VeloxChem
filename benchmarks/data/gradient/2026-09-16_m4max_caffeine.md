## gradient: caffeine

Measured at `68e136edd`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-16.

| functional | basis | nao | naux | method | gradient (s) | speedup | SCF (s) | vs four-centre |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HF | def2-svp | 246 | 1242 | four-centre | 8.84 | 1.00 | 12.54 | -- |
|  |  |  |  | RI-JK simd, in memory | 0.74 | 12.01 | 1.14 | 7.7e-05 |
| HF | def2-svpd | 366 | 1242 | four-centre | 42.11 | 1.00 | 50.36 | -- |
|  |  |  |  | RI-JK simd, in memory | 1.36 | 30.98 | 2.45 | 7.7e-05 |
| HF | def2-tzvp | 494 | 1242 | four-centre | 136.48 | 1.00 | 174.08 | -- |
|  |  |  |  | RI-JK simd, in memory | 2.52 | 54.18 | 4.43 | 7.1e-05 |
| HF | def2-tzvpd | 614 | 1242 | four-centre | 352.46 | 1.00 | 421.68 | -- |
|  |  |  |  | RI-JK simd, in memory | 3.84 | 91.69 | 7.14 | 7.3e-05 |
| B3LYP | def2-svp | 246 | 1242 | four-centre | 9.17 | 1.00 | 15.02 | -- |
|  |  |  |  | RI-JK simd, in memory | 1.25 | 7.31 | 4.06 | 1.6e-05 |
| B3LYP | def2-svpd | 366 | 1242 | four-centre | 42.79 | 1.00 | 55.86 | -- |
|  |  |  |  | RI-JK simd, in memory | 2.53 | 16.89 | 8.70 | 1.7e-05 |
| B3LYP | def2-tzvp | 494 | 1242 | four-centre | 138.15 | 1.00 | 180.23 | -- |
|  |  |  |  | RI-JK simd, in memory | 4.01 | 34.43 | 13.43 | 6.2e-05 |
| B3LYP | def2-tzvpd | 614 | 1242 | four-centre | 353.66 | 1.00 | 428.91 | -- |
|  |  |  |  | RI-JK simd, in memory | 6.64 | 53.29 | 22.46 | 6.3e-05 |

### Scaling, cost proportional to nao to the p

| functional | method | p |
| --- | --- | ---: |
| HF | four-centre | 4.01 |
| HF | RI-JK simd | 1.82 |
| B3LYP | four-centre | 3.97 |
| B3LYP | RI-JK simd | 1.78 |
