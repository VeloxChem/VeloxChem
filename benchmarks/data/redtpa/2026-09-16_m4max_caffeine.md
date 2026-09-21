## redtpa: caffeine

Measured at `ad3b36cec`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-16.

| functional | basis | nao | naux | method | TPA (s) | speedup | SCF (s) | Re gamma | Im gamma | max rel |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HF | def2-svp | 246 | 1242 | four-centre | 521.42 | 1.00 | 12.49 | 13978.0 | 2441.8 | -- |
|  |  |  |  | RI-JK simd, in memory | 54.72 | 9.53 | 1.16 | 13981.5 | 2445.5 | 8.9e-04 |
| HF | def2-svpd | 366 | 1242 | four-centre | 2465.55 | 1.00 | 51.13 | 19061.8 | 7474.3 | -- |
|  |  |  |  | RI-JK simd, in memory | 118.25 | 20.85 | 2.45 | 19054.7 | 7477.3 | 6.5e-04 |
| B3LYP | def2-svp | 246 | 1242 | four-centre | 515.02 | 1.00 | 15.04 | 32374.7 | 20396.0 | -- |
|  |  |  |  | RI-JK simd, in memory | 136.57 | 3.77 | 4.13 | 32370.3 | 20392.8 | 5.3e-04 |
| B3LYP | def2-svpd | 366 | 1242 | four-centre | 2281.42 | 1.00 | 56.42 | 56171.7 | 26439.5 | -- |
|  |  |  |  | RI-JK simd, in memory | 358.93 | 6.36 | 8.89 | 56187.4 | 26445.3 | 6.4e-04 |
