## rpa: caffeine

Measured at `df4b99a27`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-18.

| functional | basis | nao | naux | method | RPA (s) | speedup | iter | s/iter | SCF (s) | max dE (a.u.) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| B3LYP | def2-svp | 246 | 1242 | four-centre | 56.36 | 1.00 | 11 | 5.12 | 16.76 | -- |
|  |  |  |  | RI-JK simd, in memory | 14.25 | 3.95 | 11 | 1.30 | 4.42 | 6.5e-06 |
