## Submission note

- This update fixes a broken check on BLIS. The unit tests have been updated,
  and they now pass on BLIS as well as other BLAS implementations. I checked
  this with a GitHub actions workflow that uses BLIS, here is the session info:

  R Under development (unstable) (2026-09-13 r90534)
  Platform: x86_64-pc-linux-gnu
  Running under: Ubuntu 24.04.5 LTS

  Matrix products: default
  BLAS:   /usr/lib/x86_64-linux-gnu/blis-openmp/libblis.so.4
  LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0


## Test environments
* local Ubuntu 24.04.5 LTS, R 4.6.1
* Github actions
  - macOS 26.6.2, 4.6.1
  - Windows Server 2025, R 4.6.1
  - Ubuntu 24.04.5 LTS, R 4.6.1
  - Ubuntu 24.04.5 LTS, R-devel
  - Ubuntu 24.04.5 LTS, R-oldrel 4.5.3
  - Ubuntu 24.04.5 LTS, R-devel with BLIS
* win-builder, R-devel and R-release
* Rhub
  - macOS 15.7.9 R-devel (2026-09-13 r90534)
  - windows server 2025, R-devel (2026-09-13 r90534)
  - macOS-arm64 26.6.2, R-devel (2026-09-13 r90534)
  - ubuntu-gcc12 24.04.5 LTS, R-devel (2026-04-13 r89874)
  - ubuntu-release 24.04.5 LTS, R 4.6.1


## R CMD check results
There were no ERRORs, WARNINGs or NOTEs


## Downstream dependencies
There are currently no downstream dependencies for this package
