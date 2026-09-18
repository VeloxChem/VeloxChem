## gradient: caffeine

Measured at `255140804`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-18.

| functional | basis | nao | naux | method | gradient (s) | speedup | SCF (s) | vs four-centre |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CAM-B3LYP | def2-svp | 246 | 1242 | four-centre | 16.72 | 1.00 | 27.19 | -- |
|  |  |  |  | RI-JK simd, in memory | 1.97 | 8.49 | 5.23 | 1.5e-05 |
| CAM-B3LYP | def2-svpd | 366 | 1242 | four-centre | 77.01 | 1.00 | 107.00 | -- |
|  |  |  |  | RI-JK simd, in memory | 4.07 | 18.90 | 11.65 | 1.7e-05 |
| CAM-B3LYP | def2-tzvp | 494 | 1242 | four-centre | 257.77 | 1.00 | 357.61 | -- |
|  |  |  |  | RI-JK simd, in memory | 6.90 | 37.36 | 17.09 | 6.1e-05 |
| CAM-B3LYP | def2-tzvpd | 614 | 1242 | four-centre | 655.69 | 1.00 | 868.74 | -- |
|  |  |  |  | RI-JK simd, in memory | 11.01 | 59.55 | 28.56 | 6.2e-05 |

### Scaling, cost proportional to nao to the p

| functional | method | p |
| --- | --- | ---: |
| CAM-B3LYP | four-centre | 4.00 |
| CAM-B3LYP | RI-JK simd | 1.86 |
