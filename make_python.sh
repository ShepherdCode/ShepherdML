#!/bin/sh

echo "Convert each python notebook to python"

for D in Geron_*; do
    cd $D 
    pwd
    jupyter nbconvert --to script *.ipynb
    cd ..
done
