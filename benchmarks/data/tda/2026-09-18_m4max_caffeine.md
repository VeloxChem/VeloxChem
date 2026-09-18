## tda: caffeine

Measured at `df4b99a27`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-18.

| functional | basis | nao | naux | method | TDA (s) | speedup | iter | s/iter | SCF (s) | max dE (a.u.) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| B3LYP | def2-svp | 246 | 1242 | four-centre | 57.88 | 1.00 | 10 | 5.79 | 17.39 | -- |
|  |  |  |  | RI-JK simd, in memory | 12.80 | 4.52 | 10 | 1.28 | 4.34 | 6.0e-06 |
