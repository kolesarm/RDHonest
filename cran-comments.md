## Test environments
* macbuilder macOS 13.3.1 (22E261) R4.4.0 (2024-04-24)
* Github actions
  - macOS 14.7.1, R 4.4.2
  - Windows Server 2022, R 4.4.2
  - Ubuntu 22.04.5 LTS, R 4.4.2
  - Ubuntu 22.04.5 LTS, R-devel
* local Ubuntu 24.04.1 LTS install, R 4.3.3


* win-builder, R-devel and R-release

* Rhub
  - Ubuntu Linux 20.04.1 LTS, R-release, GCC
  - Debian Linux, R-devel, GCC, no long double

## R CMD check results
There were no ERRORs or WARNINGs.
There was 1 NOTE:

Maintainer: ‘Michal Kolesár <kolesarmi@googlemail.com>’

New submission

Possibly misspelled words in DESCRIPTION:
  Rothe (17:47)
  covariates (18:42)

The misspelled words are correct, "covariates" refers to independent variables
in a regression. "Rothe" is a surname.

## Downstream dependencies
There are currently no downstream dependencies for this package
