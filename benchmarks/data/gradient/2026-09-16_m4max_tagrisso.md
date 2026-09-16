## gradient: tagrisso

Measured at `59ee1892c`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-16.

| functional | basis | nao | naux | method | gradient (s) | speedup | SCF (s) | vs four-centre |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HF | def2-svp | 683 | 3387 | four-centre | 114.17 | 1.00 | 151.21 | -- |
|  |  |  |  | RI-JK simd, in memory | 92.92 | 1.23 | 29.86 | 6.0e-05 |
| B3LYP | def2-svp | 683 | 3387 | four-centre | 117.31 | 1.00 | 151.02 | -- |
|  |  |  |  | RI-JK simd, in memory | 96.05 | 1.22 | 46.55 | 1.7e-05 |
