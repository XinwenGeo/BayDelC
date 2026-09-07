# Contributing to BayDelC

Bug reports and focused pull requests are welcome. Before proposing a change:

1. open an issue describing the scientific or software problem;
2. keep changes small and do not replace calibrated model assets without a
   documented model-version decision;
3. add or update a regression test for numerical behavior; and
4. run `R CMD build .` and `R CMD check BayDelC_*.tar.gz`.

Please do not commit private data, raw MCMC working files, or machine-specific
paths. By contributing code, you agree that it can be distributed under the
package's MIT license.
