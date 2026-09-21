## optimization: caffeine

Measured at `ac2b22a5f`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank, veloxchem 1.0rc4, 2026-09-16.

| functional | basis | nao | naux | method | total (s) | speedup | steps | s/step | energy (a.u.) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HF | def2-svp | 246 | 1242 | four-centre | 619.8 | 1.00 | 33 | 18.78 | -675.83732476 |
|  |  |  |  | RI-JK simd, in memory | 90.2 | 6.87 | 27 | 3.34 | -675.83673457 |
| HF | def2-svpd | 366 | 1242 | four-centre | 2961.8 | 1.00 | 35 | 84.62 | -675.86719050 |
|  |  |  |  | RI-JK simd, in memory | 174.2 | 17.01 | 35 | 4.98 | -675.86660602 |
| B3LYP | def2-svp | 246 | 1242 | four-centre | 628.5 | 1.00 | 29 | 21.67 | -679.88509552 |
|  |  |  |  | RI-JK simd, in memory | 185.6 | 3.39 | 29 | 6.40 | -679.88514659 |
| B3LYP | def2-svpd | 366 | 1242 | four-centre | 2471.5 | 1.00 | 27 | 91.53 | -679.92921001 |
|  |  |  |  | RI-JK simd, in memory | 322.4 | 7.67 | 27 | 11.94 | -679.92926228 |
