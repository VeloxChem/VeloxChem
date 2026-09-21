## tda: caffeine

Measured at `dfc9f00c6`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-16.

| functional | basis | nao | naux | method | TDA (s) | speedup | iter | s/iter | SCF (s) | max dE (a.u.) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HF | def2-svp | 246 | 1242 | four-centre | 58.50 | 1.00 | 15 | 3.90 | 12.46 | -- |
|  |  |  |  | RI-JK simd, in memory | 3.37 | 17.34 | 15 | 0.22 | 1.12 | 2.2e-05 |
| B3LYP | def2-svp | 246 | 1242 | four-centre | 51.64 | 1.00 | 10 | 5.16 | 15.21 | -- |
|  |  |  |  | RI-JK simd, in memory | 12.24 | 4.22 | 10 | 1.22 | 4.28 | 6.0e-06 |
