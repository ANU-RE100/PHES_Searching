#!/bin/bash

mkdir -p build
pushd build || return
cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 ..
make install -j
cp compile_commands.json ..
popd || return
echo "Finished building PHES Searching :-)"


