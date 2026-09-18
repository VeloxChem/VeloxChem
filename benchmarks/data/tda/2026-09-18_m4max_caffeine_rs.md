## tda: caffeine

Measured at `cf796efd6`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-18.

| functional | basis | nao | naux | method | TDA (s) | speedup | iter | s/iter | SCF (s) | max dE (a.u.) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| B3LYP | def2-svp | 246 | 1242 | four-centre | 51.12 | 1.00 | 10 | 5.11 | 15.32 | -- |
|  |  |  |  | RI-JK simd, in memory | 12.18 | 4.20 | 10 | 1.22 | 4.09 | 6.0e-06 |
| CAM-B3LYP | def2-svp | 246 | 1242 | four-centre | 102.32 | 1.00 | 12 | 8.53 | 27.64 | -- |
|  |  |  |  | RI-JK simd, in memory | 17.07 | 5.99 | 12 | 1.42 | 5.28 | 9.1e-06 |
| WB97X-D4 | def2-svp | 246 | 1242 | four-centre | 108.94 | 1.00 | 13 | 8.38 | 27.31 | -- |
|  |  |  |  | RI-JK simd, in memory | 18.12 | 6.01 | 13 | 1.39 | 5.18 | 8.1e-06 |
