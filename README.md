## disembl
fork of [DisEMBL v1.4](http://dis.embl.de/)

This fork of DisEMBL is an almost-complete rewrite of the Python code. The regions
identified as disordered do not change.

[TISEAN](www.mpipks-dresden.mpg.de/~tisean/) Savitzky–Golay filter is replaced by 
Scipy's since it ends up being faster (no stdin/stdout necessary).

The original calculation will be available via '--old' but will be removed in
further versions.

### Authors:

Originally -

Copyright (C) 2004 Rune Linding & Lars Juhl Jensen - EMBL

Rewrite (v2.0+) -

Copyright (C) 2016 Shyam Saladi & Nadine Bradbury - Caltech

The DisEMBL is licensed under the [GNU General Public License]
(http://www.opensource.org/licenses/gpl-license.php)
