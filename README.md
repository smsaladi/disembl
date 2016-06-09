## disembl
fork of [DisEMBL v1.4](http://dis.embl.de/)

This fork of DisEMBL is an almost-complete rewrite of the code. The regions
identified as disordered do not change.

Much of the calculation is done with NumPy/Scipy instead of the previously written,
custom C code since it ends up being faster and significantly simplifies
maintainance.

The original calculation will be available via '--old' but will be removed in
further versions.

### Authors:

Originally -

Copyright (C) 2004 Rune Linding & Lars Juhl Jensen - EMBL

Rewrite (v2.0+) -

Copyright (C) 2016 Shyam Saladi - Caltech

The DisEMBL is licensed under the [GNU General Public License]
(http://www.opensource.org/licenses/gpl-license.php)
