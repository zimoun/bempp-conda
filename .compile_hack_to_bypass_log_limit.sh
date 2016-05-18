#!/bin/bash

# hack
error_handler() {
    echo "##############################################"
    echo "ERROR !! tail -n 1000 out.log"
    tail -n 1000 out.log
    kill $PING_LOOP_PID
    echo "#while true killed."
    echo "##############################################"
}
trap 'error_handler' ERR

bash -c "while true; do echo $(date); sleep 9m; done" &
PING_LOOP_PID=$!
# end

# ###### Compile bempp from scratch
#
#
conda install python=2.7 numpy scipy -q -y
python -c 'import numpy'
python -c 'import scipy'
## same test about cython ?
conda install cython=0.23.4 -y -q

wget https://github.com/bempp/bempp/archive/master.tar.gz -q
tar -xzf master.tar.gz && cd bempp-master

# ###### step 1: build
#    ## travis_wait avoids time out when downloading
#    ## and building external dependencies

mkdir -p build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/miniconda

# Building all the dependencies
echo "###################TBB..."
make TBB 1> out.log
echo "###################TBB done."
echo "###################Boost..."
make Boost 1> out.log
echo "###################Boost done."
echo "###################Dune..."
make Dune 1> out.log
echo "###################Dune done."

# the effective build
make 1> out.log

# ###### step 2: install

sudo make install > out.log

# ###### step 3: test
#    ## travis_wait has troubles with ""
#    ## but travis_wait should be useful for intensive tests

# ## Disable this test, because  it does not necessary pass.
make test ARGS="-E python_tests"

# ## Run test suite from python
python -c "import bempp.api ; bempp.api.test()"

# ## Re-Run with Gmsh installed
sudo apt-get install gmsh -q -y
python -c "import bempp.api ; bempp.api.test()"

pwd
cd ../..
#
#
# ##### Done from scratch.

kill $PING_LOOP_PID
