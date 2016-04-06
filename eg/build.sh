#!/bin/sh

mycmake=/usr/bin/cmake

mkdir -p build
cd build

echo '# configuration output' >> config.log

which $mycmake >> config.log
$mycmake --version >> config.log

which make >> config.log
make --version >> config.log

which cc >> config.log
cc --version >> config.log

which c++ >> config.log
c++ --version >> config.log

which ld >> config.log
ld --version >> config.log

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
