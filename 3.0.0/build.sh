#!/bin/sh

mycmake=/usr/bin/cmake

mkdir -p build
cd build

echo '# configuration output' | tee config.log

which $mycmake | tee config.log
$mycmake --version | tee config.log

which make | tee config.log
make --version | tee config.log

which cc | tee config.log
cc --version | tee config.log

which c++ | tee config.log
c++ --version | tee config.log

which ld | tee config.log
ld --version | tee config.log

echo 'check config done.'


$mycmake \
    -DCMAKE_INSTALL_PREFIX=$PREFIX \
    .. | tee cmake.log 2>&1

make -j \
    | tee make.log 2>&1

make install \
    | tee install.log 2>&1

make test \
    | tee test.log 2>&1
