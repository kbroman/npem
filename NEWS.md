## Revision history for the R/npem package

### Version 0.53-1, 2023-02-22

- Changelog -> NEWS.md

- In CITATION file, citEntry() -> bibentry()


### Version 0.52, 2021-07-14

- Remove Numerical Recipes in C code for sort and rand


### Version 0.51-6, 2019-06-15

- Add "LazyData: true" in DESCRIPTION file to avoid some warnings.


### Version 0.51-5, 2019-03-20

- Convert documentation to roxygen2.


### Version 0.51-4, 2017-05-21

- Small changes to avoid note about "R_registerRoutines"


### Version 0.51-3, 2015-07-18

- Changed license to GPL-3.


### Verison 0.51-2, 2011-11-07

- Added NAMESPACE.


### Version 0.50-10, 2009-04-22

- Added a keyword to the help file for npem.sem


### Version 0.50-9, 2007-10-09

- Minuscule changes to the help files, to conform to a change in R.


### Version 0.50-7, 2006-10-16

- Revised argument defaults in documentation to match the code.

- Data sets were in old .RData format; replaced with new format.


### Version 0.50-6, 2003-06-22

- In calls to C code (using the .C() function), now make use of the
  PACKAGE argument.

- Replaced calls to printf and fprintf with calls to Rprintf in
  the C code.


### Version 0.50-5, 2002-11-01

- Fixed a few minor errors in some of the help pages.


### Version 0.50-4, 2001-11-19

- Fixed a few minor errors in some of the help pages.


### Version 0.50-3, 2001-05-28

- I've revised some of the help pages so that they work with LaTeX.
  I also added examples to the help pages for npem.sem, npem.start,
  npem.ll, and npsim.


### Version 0.50, 2000-08-20

- This is a set of functions for analyzing cell proliferation assays
  with a normal-Poisson mixture model.
