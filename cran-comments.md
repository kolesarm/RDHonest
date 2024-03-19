## Test environments
* local Debian 12 (bookworm) install, R 4.3.1
* Github actions
  - macOS 12.7.3, R 4.3.3
  - Windows Server 2022, R 4.3.3
  - Ubuntu 22.04.4 LTS, R 4.3.3
  - Ubuntu 20.04.4 LTS, R-devel


* win-builder, R-devel and R-release
* macbuilder macOS 13.3.1 (22E261) R4.3.0 Patched (2023-05-18 r84451)
* Rhub
  - Ubuntu Linux 20.04.1 LTS, R-release, GCC
  - Debian Linux, R-devel, GCC, no long double

## R CMD check results
There were no ERRORs or WARNINGs.
There was 1 NOTE:

  Maintainer: ‘Michal Kolesár <kolesarmi@googlemail.com>’

  New submission

Possibly misspelled words in DESCRIPTION:
  covariates (16:33)

THe misspelled word is correct, it refers to independent variables (covariates)
in a regression.

## Downstream dependencies
There are currently no downstream dependencies for this package