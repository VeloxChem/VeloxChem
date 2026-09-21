## tda: tagrisso

Measured at `633e3a8ef`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-16.

| functional | basis | nao | naux | method | TDA (s) | speedup | iter | s/iter | SCF (s) | max dE (a.u.) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HF | def2-svp | 683 | 3387 | four-centre | 866.14 | 1.00 | 19 | 45.59 | 151.18 | -- |
|  |  |  |  | RI-JK simd, in memory | 138.38 | 6.26 | 19 | 7.28 | 29.89 | 7.4e-06 |
| B3LYP | def2-svp | 683 | 3387 | four-centre | 693.42 | 1.00 | 13 | 53.34 | 151.20 | -- |
|  |  |  |  | RI-JK simd, in memory | 167.62 | 4.14 | 13 | 12.89 | 46.74 | 5.0e-06 |
