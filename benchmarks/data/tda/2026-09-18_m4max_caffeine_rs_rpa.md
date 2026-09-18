## rpa: caffeine

Measured at `df4b99a27`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-18.

| functional | basis | nao | naux | method | RPA (s) | speedup | iter | s/iter | SCF (s) | max dE (a.u.) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CAM-B3LYP | def2-svp | 246 | 1242 | four-centre | 110.60 | 1.00 | 14 | 7.90 | 29.91 | -- |
|  |  |  |  | RI-JK simd, in memory | 21.41 | 5.17 | 14 | 1.53 | 5.55 | 9.2e-06 |
| WB97X-D4 | def2-svp | 246 | 1242 | four-centre | 118.56 | 1.00 | 14 | 8.47 | 29.70 | -- |
|  |  |  |  | RI-JK simd, in memory | 22.73 | 5.22 | 14 | 1.62 | 5.54 | 8.4e-06 |
