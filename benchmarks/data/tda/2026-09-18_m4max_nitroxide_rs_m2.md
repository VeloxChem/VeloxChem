## tda: nitroxide, multiplicity 2

Measured at `298840912`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-18.

| functional | basis | nao | naux | method | TDA (s) | speedup | iter | s/iter | SCF (s) | max dE (a.u.) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| B3LYP | def2-svp | 219 | 1067 | four-centre | 151.40 | 1.00 | 13 | 11.65 | 35.00 | -- |
|  |  |  |  | RI-JK simd, in memory | 30.43 | 4.98 | 13 | 2.34 | 6.65 | 7.1e-06 |
| CAM-B3LYP | def2-svp | 219 | 1067 | four-centre | 213.84 | 1.00 | 11 | 19.44 | 55.08 | -- |
|  |  |  |  | RI-JK simd, in memory | 31.04 | 6.89 | 11 | 2.82 | 7.77 | 7.7e-06 |
| WB97X-D4 | def2-svp | 219 | 1067 | four-centre | 199.21 | 1.00 | 10 | 19.92 | 54.34 | -- |
|  |  |  |  | RI-JK simd, in memory | 28.82 | 6.91 | 10 | 2.88 | 7.72 | 7.9e-06 |
