## optimization: nitroxide, multiplicity 2

Measured at `65d178052`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-17.

| functional | basis | nao | naux | method | total (s) | speedup | steps | s/step | energy (a.u.) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HF | def2-svp | 219 | 1067 | four-centre | 347.0 | 1.00 | 8 | 43.37 | -530.72492668 |
|  |  |  |  | RI-JK simd, in memory | 15.2 | 22.82 | 8 | 1.90 | -530.72453903 |
| B3LYP | def2-svp | 219 | 1067 | four-centre | 411.2 | 1.00 | 9 | 45.68 | -534.00249701 |
|  |  |  |  | RI-JK simd, in memory | 60.2 | 6.83 | 9 | 6.68 | -534.00258958 |
