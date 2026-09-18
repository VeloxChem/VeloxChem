## tda: caffeine

Measured at `df4b99a27`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-18.

| functional | basis | nao | naux | method | TDA (s) | speedup | iter | s/iter | SCF (s) | max dE (a.u.) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CAM-B3LYP | def2-svp | 246 | 1242 | four-centre | 113.72 | 1.00 | 12 | 9.48 | 29.86 | -- |
|  |  |  |  | RI-JK simd, in memory | 17.73 | 6.41 | 12 | 1.48 | 5.64 | 9.1e-06 |
| WB97X-D4 | def2-svp | 246 | 1242 | four-centre | 122.72 | 1.00 | 13 | 9.44 | 29.76 | -- |
|  |  |  |  | RI-JK simd, in memory | 19.06 | 6.44 | 13 | 1.47 | 5.57 | 8.1e-06 |
