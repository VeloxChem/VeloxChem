## optimization: nitroxide, multiplicity 2

Measured at `06617da2a`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-17.

| functional | basis | nao | naux | method | total (s) | speedup | steps | s/step | energy (a.u.) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HF | def2-svpd | 330 | 1067 | four-centre | 1198.6 | 1.00 | 7 | 171.24 | -530.75639276 |
|  |  |  |  | RI-JK simd, in memory | 26.9 | 44.48 | 7 | 3.85 | -530.75600318 |
| B3LYP | def2-svpd | 330 | 1067 | four-centre | 1521.9 | 1.00 | 9 | 169.10 | -534.04656378 |
|  |  |  |  | RI-JK simd, in memory | 135.7 | 11.22 | 9 | 15.08 | -534.04664616 |
