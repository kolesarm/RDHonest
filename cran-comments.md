## Submission note

- This update fixes a broken check on BLIS. The unit tests have been adjusted to
  pass. I checked this with a github actions workflow that uses blis:

  === R Session Info ===
  R Under development (unstable) (2026-09-13 r90534)
  Platform: x86_64-pc-linux-gnu
  Running under: Ubuntu 24.04.5 LTS

  Matrix products: default
  BLAS:   /usr/lib/x86_64-linux-gnu/blis-openmp/libblis.so.4
  LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0


## Test environments
* local Debian GNU/Linux 13 (trixie) install, R 4.5.0
* Github actions
  - macOS 26.6.2, 4.6.1
  - Windows Server 2025, R 4.6.1
  - Ubuntu 24.04.5 LTS, R 4.6.1
  - Ubuntu 24.04.5 LTS, R-devel
  - Ubuntu 24.04.5 LTS, R-oldrel 4.5.3
  - Ubuntu 24.04.5 LTS, R-devel with BLIS
* Rhub
  - macOS 13.7.1 R-devel (2024-12-15 r87442)
  - macOS-arm64 14.7.1, R-devel (2024-12-15 r87442)
  - ubuntu-gcc12 22.04.5 LTS R-devel (2024-12-15 r87442)
  - ubuntu-nold 22.04.5 LTS, R-devel (2024-12-15 r87442)
  - ubuntu-release 22.04.5 LTS, R 4.4.2
* macbuilder macOS 13.3.1 (22E261) R4.4.0 (2024-04-24)
* win-builder, R-devel and R-release


## R CMD check results
There were no ERRORs, WARNINGs or NOTEs


## Downstream dependencies
There are currently no downstream dependencies for this package
