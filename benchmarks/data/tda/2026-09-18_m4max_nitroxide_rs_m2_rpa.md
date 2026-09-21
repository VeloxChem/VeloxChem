## rpa: nitroxide, multiplicity 2

Measured at `298840912`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-18.

| functional | basis | nao | naux | method | RPA (s) | speedup | iter | s/iter | SCF (s) | max dE (a.u.) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| B3LYP | def2-svp | 219 | 1067 | four-centre | 127.08 | 1.00 | 12 | 10.59 | 35.06 | -- |
|  |  |  |  | RI-JK simd, in memory | 28.39 | 4.48 | 12 | 2.37 | 6.65 | 7.1e-06 |
| CAM-B3LYP | def2-svp | 219 | 1067 | four-centre | 206.47 | 1.00 | 12 | 17.21 | 54.56 | -- |
|  |  |  |  | RI-JK simd, in memory | 35.36 | 5.84 | 12 | 2.95 | 7.87 | 7.7e-06 |
| WB97X-D4 | def2-svp | 219 | 1067 | four-centre | 205.92 | 1.00 | 12 | 17.16 | 54.84 | -- |
|  |  |  |  | RI-JK simd, in memory | 35.27 | 5.84 | 12 | 2.94 | 7.78 | 7.9e-06 |
