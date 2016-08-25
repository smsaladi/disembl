 ____  _     _____ __  __ ____  _        ____      ______
|  _ \(_)___| ____|  \/  | __ )| |      /___  )   (  __  )
| | | | / __|  _| | |\/| |  _ \| |         / /    | |  | |
| |_| | \__ \ |___| |  | | |_) | |___     / /_    | |__| |
|____/|_|___/_____|_|  |_|____/|_____|  /_____|(_)(______)

This fork of [DisEMBL v1.4](http://dis.embl.de/) represents a significant
rewriting of the code without changing the scores/regions identified as
disordered.

Much of the calculation is done with NumPy/Scipy instead of previously written,
custom Python code. C code is wrapped directly instead of through stdin/out and
subprocess. These changes result in significant speed increase and will
simplify maintenance.

## Authors:

Originally -
Copyright (C) 2004 Rune Linding & Lars Juhl Jensen - EMBL

Rewrite (2.0+) -
Copyright (C) 2016 Nadine Bradbury & Shyam Saladi - Caltech

## License
The DisEMBL is licensed under the [GNU General Public License]
(http://www.opensource.org/licenses/gpl-license.php)
