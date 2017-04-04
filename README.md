[![Build Status](https://travis-ci.org/smsaladi/disembl.svg?branch=master)](https://travis-ci.org/smsaladi/disembl)

 ____  _     _____ __  __ ____  _        ____      ______
|  _ \(_)___| ____|  \/  | __ )| |      /___  )   (  __  )
| | | | / __|  _| | |\/| |  _ \| |         / /    | |  | |
| |_| | \__ \ |___| |  | | |_) | |___     / /_    | |__| |
|____/|_|___/_____|_|  |_|____/|_____|  /_____|(_)(______)

This fork of [DisEMBL v1.4](http://dis.embl.de/) represents a significant
rewriting of the code without changing the scores/regions identified as
disordered.

Much of the Python code is revised.
[TISEAN](www.mpipks-dresden.mpg.de/~tisean/)'s Savitzky–Golay filter is replaced
by [Scipy](http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.signal.savgol_filter.html)'s
since it ends up being faster (no stdin/stdout necessary). C code is wrapped
directly instead of through stdin/out and subprocess. These changes result in
significant speed increase and will future simplify maintenance.

## Authors:

Originally -
Copyright (C) 2004 Rune Linding & Lars Juhl Jensen - EMBL

Rewrite (v2.0+) -
Copyright (C) 2016 Shyam Saladi & Nadine Bradbury - Caltech

## Installation

* Clone repo

```shell
git clone git@github.com:smsaladi/disembl.git
```

* Install with `pip`

```shell
cd disembl
pip install .
```

* If you've already installed the package and want to reinstall, try

```shell
pip install . -I
```

All dependencies should be checked for and, if necessary, installed
automatically by `pip`.

## License
The DisEMBL is licensed under the [GNU General Public License Version 2]
(https://opensource.org/licenses/GPL-2.0).



## Making releases

1. Reconcile package requirements in `requirements.txt` with those listed at
`install_requires`, `setup_requires`, and `test_requires` in `setup.py`.

2. Confirm build passes all tests

2. Update the version number and `download_url` in `setup.py`

3. Tag a release on GitHub with the new version number
