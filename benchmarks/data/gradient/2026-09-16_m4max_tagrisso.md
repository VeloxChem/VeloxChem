## gradient: tagrisso

Measured at `68e136edd`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-16.

| functional | basis | nao | naux | method | gradient (s) | speedup | SCF (s) | vs four-centre |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HF | def2-svp | 683 | 3387 | four-centre | 114.50 | 1.00 | 151.18 | -- |
|  |  |  |  | RI-JK simd, in memory | 9.95 | 11.51 | 29.79 | 6.0e-05 |
| B3LYP | def2-svp | 683 | 3387 | four-centre | 117.10 | 1.00 | 151.51 | -- |
|  |  |  |  | RI-JK simd, in memory | 13.01 | 9.00 | 46.49 | 1.7e-05 |
