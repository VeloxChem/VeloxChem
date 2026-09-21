## two-center Coulomb: the plain driver against the range separated one

Measured at `befe15723`, working tree dirty on m4max (Apple M4 Max, 14 cores), veloxchem 1.0rc4, 2026-09-18.

omega 0.3, best of three after a cold call, each case in its own process.

| molecule | fitting set | nao | lmax | packed GB | threads | plain (s) | both (s) | ratio | settled |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| tagrisso | jfit | 2176 | 4 | 0.02 | 14 | 0.0022 | 0.0038 | 1.74 | **no, 6% spread** |
| tagrisso | jfit | 2176 | 4 | 0.02 | 1 | 0.0131 | 0.0254 | 1.94 | yes |
| tagrisso | jkfit | 3387 | 4 | 0.04 | 14 | 0.0037 | 0.0067 | 1.82 | **no, 10% spread** |
| tagrisso | jkfit | 3387 | 4 | 0.04 | 1 | 0.0255 | 0.0487 | 1.91 | yes |
| taxol | jfit | 3528 | 4 | 0.05 | 14 | 0.0046 | 0.0086 | 1.85 | yes |
| taxol | jfit | 3528 | 4 | 0.05 | 1 | 0.0351 | 0.0698 | 1.99 | yes |
| taxol | jkfit | 5489 | 4 | 0.11 | 14 | 0.0085 | 0.0160 | 1.87 | yes |
| taxol | jkfit | 5489 | 4 | 0.11 | 1 | 0.0678 | 0.1322 | 1.95 | yes |
| crambin | jfit | 19500 | 4 | 1.42 | 14 | 0.1362 | 0.2787 | 2.05 | yes |
| crambin | jfit | 19500 | 4 | 1.42 | 1 | 1.1712 | 2.2915 | 1.96 | yes |
| crambin | jkfit | 30751 | 4 | 3.52 | 14 | 0.2670 | 0.5395 | 2.02 | yes |
| crambin | jkfit | 30751 | 4 | 3.52 | 1 | 2.3826 | 4.8323 | 2.03 | yes |
| ubiquitin | jfit | 36419 | 4 | 4.94 | 14 | 0.5770 | 1.2082 | 2.09 | yes |
| ubiquitin | jfit | 36419 | 4 | 4.94 | 1 | 4.5558 | 8.9174 | 1.96 | yes |
| ubiquitin | jkfit | 56971 | 4 | 12.09 | 14 | 1.1865 | 2.4359 | 2.05 | **no, 354% spread** |
| ubiquitin | jkfit | 56971 | 4 | 12.09 | 1 | 9.3864 | 34.6513 | 3.69 | **no, 15% spread** |

The range separated driver holds two packed matrices where the plain one holds a
single one, so the last case is 24.18 GB of matrix alone. Its three timed calls
were 11.06, 2.47 and 2.44 seconds at fourteen threads against plain calls tight to
under one per cent in the same run, which is memory pressure and not arithmetic.
Those rows carry a ratio in the table and are marked as what they are.

