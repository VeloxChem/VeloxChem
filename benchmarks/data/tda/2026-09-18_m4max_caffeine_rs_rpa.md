## rpa: caffeine

Measured at `cf796efd6`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-18.

| functional | basis | nao | naux | method | RPA (s) | speedup | iter | s/iter | SCF (s) | max dE (a.u.) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| B3LYP | def2-svp | 246 | 1242 | four-centre | 49.65 | 1.00 | 11 | 4.51 | 15.17 | -- |
|  |  |  |  | RI-JK simd, in memory | 13.63 | 3.64 | 11 | 1.24 | 4.10 | 6.5e-06 |
| CAM-B3LYP | def2-svp | 246 | 1242 | four-centre | 99.30 | 1.00 | 14 | 7.09 | 27.32 | -- |
|  |  |  |  | RI-JK simd, in memory | 20.55 | 4.83 | 14 | 1.47 | 5.24 | 9.2e-06 |
| WB97X-D4 | def2-svp | 246 | 1242 | four-centre | 105.51 | 1.00 | 14 | 7.54 | 27.01 | -- |
|  |  |  |  | RI-JK simd, in memory | 21.66 | 4.87 | 14 | 1.55 | 5.24 | 8.4e-06 |
