#!/usr/bin/bash

INITIAL_WD=$(pwd)
cd "$(dirname "$0")"

python ./setup.py build_ext --build-lib ./

if [ $? -eq 0 ]
then
  echo -e "\e[32mCompilation successful!"
else
  echo -e "\e[31mCompilation exited with error!" >&2
  cd $INITIAL_WD
  exit 1
fi

mv ./cy_stringscan.*.so ../scripts/cy_stringscan.so
mv ./cy_bamindexer.*.so ../scripts/cy_bamindexer.so 
rm -r ./build 

cd $INITIAL_WD
exit 0
