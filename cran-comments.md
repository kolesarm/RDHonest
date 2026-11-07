## Submission note

- This update fixes the test failures reported on CRAN's MKL additional-issues
  check. The failing tests compared separately computed results at or near
  machine precision. Some of these comparisons passed through numerical
  bandwidth optimization, which amplifies rounding differences (such as those
  between BLAS implementations) to about sqrt(.Machine$double.eps).
- The tests now fix the bandwidth and smoothness constant in cases where we are
  not testing their selection, so they isolate the computation being tested and
  agree to rounding error. Comparisons that do involve optimization use a
  tolerance of 1e-6. I checked the margins by perturbing the computations in
  ways that leave the results mathematically unchanged: differences are about
  1e-15 with tuning parameters fixed and below 5e-8 with them re-selected.
- The tests also pass with MKL 2025.3 on R-hub and with MKL 2026.2 on GitHub
  Actions, although neither setup reproduced the failures in the previous
  version, presumably due to hardware differences.

## Test environments
* local Ubuntu 24.04.5 LTS, R 4.6.1
* Github actions
  - Ubuntu 24.04.5 LTS, R-devel with BLIS
  - macOS 26.7.2, 4.6.2
  - Windows Server 2025, R 4.6.2
  - Ubuntu 24.04.5 LTS, R-devel
  - Ubuntu 24.04.5 LTS, R 4.6.2
  - Ubuntu 24.04.5 LTS, R-oldrel 4.5.3
* win-builder, R-devel and R-release
* Rhub
  - mkl (Intel MKL 2025.3), R-devel (2026-11-06 r90845)
  - ubuntu-gcc12 24.04.5 LTS, R-devel (2026-04-13 r89874)
  - windows server 2025, R-devel (2026-11-06 r90845 ucrt)
  - ubuntu-release 24.04.5 LTS, R 4.6.2
  - macOS-arm64 26.7.2, R-devel (2026-11-06 r90845)
  - macOS 15.7.9 R-devel (2026-11-06 r90845)




## R CMD check results
There were no ERRORs, WARNINGs or NOTEs


## Downstream dependencies
There are currently no downstream dependencies for this package
