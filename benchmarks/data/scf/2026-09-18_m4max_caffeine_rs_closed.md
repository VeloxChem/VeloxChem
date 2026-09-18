## SCF: caffeine

Measured at `ba5fdd191`, working tree dirty on m4max (Apple M4 Max, 14 cores), 1 rank of 14 threads, veloxchem 1.0rc4, 2026-09-18.

#### CAM-B3LYP

| basis | nao | fitting set | naux | method | wall | B vectors | 2e build | XC | rest | iters | energy | build x | whole x |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| def2-svp | 246 | def2-universal-jkfit | 1242 | full four-centre | 27.54 | -- | 23.45 | 2.96 | 1.14 | 20 | -679.55269065 | -- | -- |
|  |  |  |  | RI-JK simd, in memory | 5.20 | 0.60 | 1.30 | 3.06 | 0.23 | 21 | -679.55275724 | 18.86 | 5.30 |
| def2-svpd | 366 | def2-universal-jkfit | 1242 | full four-centre | 107.72 | -- | 99.21 | 6.86 | 1.66 | 20 | -679.59776076 | -- | -- |
|  |  |  |  | RI-JK simd, in memory | 11.63 | 1.22 | 2.83 | 7.24 | 0.34 | 21 | -679.59782715 | 36.73 | 9.26 |
| def2-tzvp | 494 | def2-universal-jkfit | 1242 | full four-centre | 361.56 | -- | 348.51 | 8.87 | 4.18 | 20 | -680.30697360 | -- | -- |
|  |  |  |  | RI-JK simd, in memory | 16.95 | 2.11 | 5.10 | 9.12 | 0.62 | 21 | -680.30712927 | 71.62 | 21.34 |
| def2-tzvpd | 614 | def2-universal-jkfit | 1242 | full four-centre | 868.74 | -- | 847.79 | 15.28 | 5.67 | 20 | -680.30989647 | -- | -- |
|  |  |  |  | RI-JK simd, in memory | 28.31 | 3.21 | 8.57 | 15.64 | 0.89 | 21 | -680.31005482 | 103.65 | 30.69 |

#### WB97X-D4

| basis | nao | fitting set | naux | method | wall | B vectors | 2e build | XC | rest | iters | energy | build x | whole x |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| def2-svp | 246 | def2-universal-jkfit | 1242 | full four-centre | 27.18 | -- | 23.20 | 2.86 | 1.12 | 20 | -680.11304033 | -- | -- |
|  |  |  |  | RI-JK simd, in memory | 5.12 | 0.60 | 1.31 | 2.98 | 0.23 | 21 | -680.11312348 | 18.54 | 5.31 |
| def2-svpd | 366 | def2-universal-jkfit | 1242 | full four-centre | 106.92 | -- | 98.49 | 6.75 | 1.68 | 20 | -680.15182844 | -- | -- |
|  |  |  |  | RI-JK simd, in memory | 11.46 | 1.21 | 2.82 | 7.08 | 0.35 | 21 | -680.15191055 | 36.56 | 9.33 |
| def2-tzvp | 494 | def2-universal-jkfit | 1242 | full four-centre | 360.31 | -- | 347.09 | 8.83 | 4.39 | 20 | -680.85206707 | -- | -- |
|  |  |  |  | RI-JK simd, in memory | 17.54 | 2.11 | 5.35 | 9.45 | 0.63 | 22 | -680.85222816 | 71.03 | 20.54 |
| def2-tzvpd | 614 | def2-universal-jkfit | 1242 | full four-centre | 868.04 | -- | 846.68 | 15.49 | 5.87 | 20 | -680.85474498 | -- | -- |
|  |  |  |  | RI-JK simd, in memory | 29.23 | 3.21 | 8.91 | 16.24 | 0.87 | 22 | -680.85490763 | 104.05 | 29.70 |

