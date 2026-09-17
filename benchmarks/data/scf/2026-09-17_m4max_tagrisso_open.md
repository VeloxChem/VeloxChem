## SCF: tagrisso, multiplicity 3

Measured at `71bbf63da`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-17.

#### HF

| basis | nao | fitting set | naux | method | wall | B vectors | 2e build | XC | rest | iters | energy | build x | whole x |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| def2-svp | 683 | def2-universal-jkfit | 3387 | full four-centre | 1737.89 | -- | 1704.45 | 0.00 | 33.45 | 100 | did not converge | -- | -- |
|  |  |  |  | RI-JK veloxchem | 456.85 | 43.20 | 405.64 | 0.00 | 8.01 | 100 | did not converge | 4.20 | 3.80 |
|  |  |  |  | RI-JK simd, in memory | 187.20 | 7.86 | 171.29 | 0.00 | 8.04 | 100 | did not converge | 9.95 | 9.28 |

#### B3LYP

| basis | nao | fitting set | naux | method | wall | B vectors | 2e build | XC | rest | iters | energy | build x | whole x |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| def2-svp | 683 | def2-universal-jkfit | 3387 | full four-centre | 836.20 | -- | 749.58 | 56.81 | 29.82 | 45 | -1619.19746283 | -- | -- |
|  |  |  |  | RI-JK veloxchem | 320.75 | 42.97 | 208.79 | 63.95 | 5.05 | 50 | -1619.19877019 | 3.98 | 2.61 |
|  |  |  |  | RI-JK simd, in memory | 164.02 | 7.91 | 87.62 | 63.30 | 5.19 | 50 | -1619.19877019 | 9.48 | 5.10 |

