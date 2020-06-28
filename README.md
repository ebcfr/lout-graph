# A slightly modified graph module

This project provides a document for the Lout typesetting system (https://savannah.nongnu.org/projects/lout/).

It contains experiences on a slightly modified "graph" module.

This module 
* is compliant with the standard graph module found in the Lout distribution
* adds an explicit grid option tied to ticks usable for all graph styles
* adds a subgrid automatically used when using log scales, and optionnal when using linear scales.

The comands used to generate the example document are

```
lout ex_doc >ex_doc.ps
ps2pdf -sPAPERSIZE=a4 -dPDFSETTINGS=/prepress -dAutoRotatePages=/All ex_doc.ps
```
