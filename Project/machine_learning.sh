#!/bin/sh

export SRC="/Users/jasonmiller/Source/Python/ShepherdML/Strings"

date
echo "K=2 ML"
python3 ${SRC}/features_to_train_set.py ncRNA.2mer.features.csv pcRNA.2mer.features.csv
python3 ${SRC}/machine_learning.py
date

date
echo "K=3 ML"
python3 ${SRC}/features_to_train_set.py ncRNA.3mer.features.csv pcRNA.3mer.features.csv
python3 ${SRC}/machine_learning.py
date

date
echo "K=4 ML"
python3 ${SRC}/features_to_train_set.py ncRNA.4mer.features.csv pcRNA.4mer.features.csv
python3 ${SRC}/machine_learning.py
date

date
echo "K=5 ML"
python3 ${SRC}/features_to_train_set.py ncRNA.5mer.features.csv pcRNA.5mer.features.csv
python3 ${SRC}/machine_learning.py
date

date
echo "K=6 ML"
python3 ${SRC}/features_to_train_set.py ncRNA.6mer.features.csv pcRNA.6mer.features.csv
python3 ${SRC}/machine_learning.py
date




