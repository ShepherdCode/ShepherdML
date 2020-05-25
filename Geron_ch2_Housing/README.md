# Machine Learning - Housing
Hands-On Machine Learning by Aurelien Geron.\
Chapter 2 - End-to-End Project with California housing prices.

Housing 1-6 was first try. Followed book. Didn't understand it all.

Step 1-5 was second try. Experimented more. Changed things.

## Installation

Update (use -U for update):\
`pip install -U jupyter ipython`\
`pip install -U sklearn scikit-learn scipy numpy pandas matplotlib`\
`pip install -U openssl`

Launch:\
`jupyter notebook`\
Starts a Jupyter kernel.
Starts a web server on localhost:8888.
Opens new tab in default browser showing directory.

Convert jupyter notebook (\*.ipynb) to python (\*.py):\
`jupyter nbconvert --to script Housing5.ipynb`\
Note that `ipython nbconvert` is deprecated.

## Components

### IPython
Enhanced python interpreter.
See ipython.readthedocs.io.
Some ipython features have since been added to the standard python interpreter.

IPython supports parallel computing.
IPython runs its own python kernel.
IPython adds meta commands

* Get help with `?thing` and `??thing`
* Store data with `%store`
* Change directory with `%cd dir`
* Execute system commands with `!!cmd` or `var=!cmd`

IPython magic refers to a set of convenience functions.
Line magic looks like %cmd and operate on this line only.
Cell magic looks like %%cmd and involve subsequent lines.

### matplotlib
This library draws plots.
Its frontend is the user or the python code.
Its backend is settable: new window, postscript, Jupyter.

The python command `get_ipython().run_line_magic('matplotlib', 'inline')` is used in case the code runs in standard python.
The ipython equivalent `%matplotlib inline` means
the output of plotting commands will be displayed inline
(directly below the code cell that produced it);
with Jupyter, the image will get encoded in the notebook too.

### Jupyter
See jupyter.readthedocs.io.
Extension to ipython.
Jupyter is a metapackage of packages.

* Jupyter notebook package is a frontend to ipython kernel.
Saves sessions to files called notebooks.
Notebook native format is JSON.
The nbconvert tool converts notebook to python, HTML, LaTeX.
* jupyter_console terminal application
* jupyterhub for collaboration
* nbgrader for creating and grading homework

### Pandas
Python library.

* [DataFrame](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html)
Can hold dict, list, list of dict, etc.
Can usually incorporate column nmes.
Support for I/O from dict, csv file, numpy ndarray.
Hundreds of methods like describe(), groupby(), kurtosis().
We use cut() with bins as parameter.

### SKLearn
Machine Learning library.
We use OneHotEncoder, StratifiedShuffleSplit,
LinearRegression, mean_squared_error, DecisionTreeRegressor,
SimpleImputer.

## Pipeline
Partition learning and validation data.

Partition learning data int train and test data.
* dataframe = pandas.read_csv()
* histogram = pandas.cut(dataframe, bins=, labels=)
* splitter = StratifiedShuffleSplit(test+size=, random_state=)\
* splitter.split() iterates (train,test) pairs
* Evaluate stratification, find undifferentiated features
* for each (train,test), drop() useless columns

Data preparation.
* plot(scatter,longitude,lattitude) to get map!
* plot(colorbar) to color by housing price
* Measure correlation with dataframe.corr().sort_values()
* 4x4 array of scatter plots for 4x4 different features
* Plot again to zoom in on selected features.
* Temporarily separate categoric and numeric data.
* Use SimpleImputer.fit() to fix numeric missing values.
* Apply OneHotEncoder to one feature: ocean proximity.
* Encoder.fit_transform(categoric) builds sparse array.
* Recombine numeric and categoric data into one dataframe.

Learning
* LinearRegression.fit(prepared_data,labels)
* lin_reg.predict() on a subset of data
* lin_reg.coef_ to view coefficient per feature (want positive)

Validation
* pred.mean_squared_error(labels,predictions)
* cross_val_score(lin_reg,prepared_data,lables,cost_func)

Compare to other models
* DecisionTreeRegressor
* cross_val_score()
