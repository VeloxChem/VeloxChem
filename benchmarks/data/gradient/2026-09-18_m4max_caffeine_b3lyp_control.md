## gradient: caffeine

Measured at `255140804`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-18.

| functional | basis | nao | naux | method | gradient (s) | speedup | SCF (s) | vs four-centre |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| B3LYP | def2-svp | 246 | 1242 | four-centre | 9.29 | 1.00 | 15.60 | -- |
|  |  |  |  | RI-JK simd, in memory | 1.29 | 7.21 | 4.08 | 1.6e-05 |
| B3LYP | def2-svpd | 366 | 1242 | four-centre | 43.67 | 1.00 | 57.67 | -- |
|  |  |  |  | RI-JK simd, in memory | 2.68 | 16.30 | 8.83 | 1.7e-05 |
| B3LYP | def2-tzvp | 494 | 1242 | four-centre | 139.76 | 1.00 | 181.81 | -- |
|  |  |  |  | RI-JK simd, in memory | 4.29 | 32.62 | 13.65 | 6.2e-05 |
| B3LYP | def2-tzvpd | 614 | 1242 | four-centre | 359.18 | 1.00 | 434.90 | -- |
|  |  |  |  | RI-JK simd, in memory | 6.99 | 51.40 | 22.82 | 6.3e-05 |

### Scaling, cost proportional to nao to the p

| functional | method | p |
| --- | --- | ---: |
| B3LYP | four-centre | 3.97 |
| B3LYP | RI-JK simd | 1.81 |
