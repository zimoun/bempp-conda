

[Bem++](http://www.bempp.org)
is really sensitive about the version of the dependencies.

Therefore, it is possible that these recipes do not work as they are.
(tested with Debian stable and testing, and Ubuntun LTS 14)


Works on Linux: Ubuntu LTS and Debian stable/testing.


How to download and install Miniconda
=====================================

Please, give a look at:
 - http://conda.pydata.org/miniconda.html
 - http://conda.pydata.org/docs/install/quick.html

In short, just run,

```
    $ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    $ bash Miniconda3-latest-Linux-x86_64.sh
```

then ENTER, then q, then yes then ENTER then yes
and that's all.

Then you need the special build package from conda,

```
    $ conda install conda-build --yes
```

About version of Python
-----------------------

At the building step, you can specify the version of Python 2 or 3.
However, you might then remember when you install it.

By default, everything is Python 3.
And you should move to Pyhton 3 (see
[PEP 373](https://www.python.org/dev/peps/pep-0373)
and [PEP 404](https://www.python.org/dev/peps/pep-0404))


Note
----

You can also `update` the conda environment before doing anything,

```
    $ conda update --all
```


How to build
============

Once `conda-build` installed, you only have to run,

```
    $ cd VERSION/
    $ conda build .
```

where `VERSION` is ... the version that you want to build.
And then pray.

All the logging files are in you `conda` install,
by default `path/to/{mini}conda/conda-bld/work/{VERSION/}build/`
especially the `config.log`, `cmake.log` and `make.log`.

It is really sensitive and should not work, see below.

Note
----

By default, `conda` fetches `numpy` with the `MKL` support
(see [there](https://docs.continuum.io/mkl-optimizations/index)).
However, you need to have Intel processor, otherwise it will be
slower.

To desactivate the `MKL` support, you need to edit the `meta.yaml`
file of the considered version, then you need to add `nomkl` at both
place, i.e., write `- nomkl` after `- numpy`.

How to check your processor, type in a terminal:
```
   $ cat /proc/cpuinfo | grep -e 'model name'
```
otherwise, look for around the flags.

Note (reminder to me)
----

The files `build.sh`
and all the python test files are created by this script,

```
#!/bin/sh

cat eg/bempp_tests.py eg/laplace.py eg/scattering.py eg/maxwell.py > eg/run_test.py

cd VERSION/
ln -s ../eg/run_test.py .
ln -s ../eg/build.sh .
```


How to install the new built `bempp` package
============================================

Once bempp built, and if `TEST` passed,
then you can locally use it.

Recommended
-----------

Create `conda environment`, e.g.,

```
    $ conda create -n my-name bempp --use-local --yes
```

Then, you switch to this environment to use it,

```
    $ source activate my-name
```

And just try to launch `ipython` and import `bempp.api`.

Further readings about conda
----------------------------

Give a look at http://conda.pydata.org/docs/index.html


Less recommended
----------------

Install in the root conda environment,

```
    $ conda install bempp --use-local --yes
```

Issues (versions)
======

- Using these versions, we are able to compile and pass the tests
[![Build Status](https://travis-ci.org/zimoun/bempp.svg?branch=eg-travis)](https://travis-ci.org/zimoun/bempp)

 + CMake needs to be < 3.5
 + Boost: ?
 + TBB: release tbb43_20150209oss_lin.tgz
 + Dune: 2.3.1

If the error message is cryptic when you are trying to compile,
then give a look to these dependencies.

- Special files are created but not deployed from the building tree,
e.g., `build/bin/ibempp`
Running these files allows the access to special import, such that
some low-level quadrature rules. I am working on it...
