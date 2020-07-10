#!/bin/sh

export SRC="/Users/jasonmiller/Source/Python/ShepherdML/Strings"

date
echo "K=2 ML"
python3 ${SRC}/features_to_train_set.py ncRNA.2mer_co.features.csv pcRNA.2mer_co.features.csv

python3 ${SRC}/machine_learning.py
date



