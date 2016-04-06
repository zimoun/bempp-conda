

[Bem++](http://www.bempp.org)
is really sensitive about the version of the dependencies.

Therefore, it is possible that these recipes do not work as they are.

How to download and install Miniconda
=====================================

Please, give a look at:
 - http://conda.pydata.org/miniconda.html
 - http://conda.pydata.org/docs/install/quick.html

In short, just run

```
 $ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
  $ bash Miniconda3-latest-Linux-x86_64.sh
```

then ENTER, then q, then yes then ENTER then yes
and that's all.


How to build
============

Once `conda-build` installed by,

```
    $ conda install conda-build
```

you only have to run,

```
    $ cd VERSION/
    $ conda build .
```

where `VERSION` is ... the version that you want to build.
And then pray.

All the logging file are in you `conda` install,
by defautl `path/to/{mini}conda/conda-bld/work/{VERSION/}build/`
especially the `config.log`, `cmake.log` and `make.log`.


Note
----

The files `build.sh`
and all the python test files are created by this script.

```
#!/bin/sh

cat eg/laplace.py eg/scattering.py eg/maxwell.py > eg/run_test.py

ln eg/run_test.py VERSION/run_test.py
ln eg/build.sh VERSION/build.sh
```
