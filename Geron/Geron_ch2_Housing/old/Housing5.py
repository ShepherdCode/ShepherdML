#!/usr/bin/env python
# coding: utf-8

# In[5]:


import pandas as pd
import os
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import numpy
from sklearn.preprocessing import OrdinalEncoder
from sklearn.preprocessing import OneHotEncoder
from sklearn.model_selection import StratifiedShuffleSplit
import numpy as np

HOUSING_PATH=os.path.join("../../datasets","housing")
HOUSING_FILE=os.path.join(HOUSING_PATH,"housing.csv")

# returns a Pandas DataFrame
def load_data (data_file=HOUSING_FILE):
    data = pd.read_csv(data_file)
    return data

housing=load_data()

# Create random test and training sets with stratification by income
# Pandas cut() retuns a temporary column containing bin labels 
housing["temp_income_bin"] = pd.cut(housing["median_income"],bins=[0.0,1.5,3.0,4.5,6.0,np.inf],labels=[1,2,3,4,5])
# Splitter contains indexes to a stratified-random subset
splitter = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
for train_index, test_index in splitter.split(housing, housing["temp_income_bin"]):
    # Pandas DataFrame.loc() retrieves from array by location 
    train_set = housing.loc[train_index]
    test_set = housing.loc[test_index]
# Remove the temporary column containing bin labels.
for s in (train_set,test_set,housing):
    # Pandas DataFrame.drop() removes a column 
    s.drop("temp_income_bin",axis=1,inplace=True)
    
# From now on, look at the training set only.
# From now on, separate training data from training labels. Look at data only.
train_set.sort_index(inplace=True)  # in case they got shuffled, sort by original index
train_set.reset_index(drop=True,inplace=True)  # in case indices are no longer continuous, reindex
housing_data = train_set.drop("median_house_value",axis=1) # Pandas drop() creates a copy minus one column.
housing_labels = train_set["median_house_value"].copy()

# Establish all possible labels in our one categoric field.
temp_categoric = housing_data[["ocean_proximity"]]   # this is a Pandas DataFrame object
encoder = OneHotEncoder()
encoded = encoder.fit_transform(temp_categoric)   # 1 bit per row, stored in sparse matrix
OCEAN_VALUES = encoder.categories_


# In[6]:


def full_pipeline (housing_data, ocean_values):

    # Fix missing features. Use a repeatable process using scikit. 
    # Apply the Imputer. Unfortunately, it only works on numeric data.
    temp_numeric = housing_data.drop("ocean_proximity",axis=1)
    from sklearn.impute import SimpleImputer
    imputer = SimpleImputer(strategy="median")
    imputer.fit(temp_numeric)
    array_with_imputes = imputer.transform(temp_numeric)
    housing_with_imputes = pd.DataFrame(array_with_imputes, columns=temp_numeric.columns, index=temp_numeric.index)

    # Apply Encoder to categoric data.
    # Apply one-hot encoding rather than sequential numbers.
    temp_categoric = housing_data[["ocean_proximity"]]   # this is a Pandas DataFrame object
    encoder = OneHotEncoder(ocean_values)
    encoded = encoder.fit_transform(temp_categoric)   # 1 bit per row, stored in sparse matrix
    temp_df = pd.DataFrame(encoded.toarray())
    temp_df.columns = encoder.categories_

    # Combine numeric and categoric data
    housing_prepared = housing_with_imputes.join(temp_df)
    return (housing_prepared)


# In[7]:


housing_prepared = full_pipeline(housing_data,OCEAN_VALUES) 
housing_prepared


# In[8]:


# Linear Regression on first five data rows
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error

lin_reg = LinearRegression()
lin_reg.fit(housing_prepared,housing_labels)

some_data = housing_data.iloc[:5]
some_prepared = full_pipeline(some_data,OCEAN_VALUES)
some_labels = housing_labels.iloc[:5]

predictions = lin_reg.predict(some_prepared)
some_labels.tolist(), predictions.round()
lin_reg.coef_   # largest positive coefficients on house_age, area_income, island

# Linear Regression on test set
#housing_prepared = full_pipeline(housing_data,OCEAN_VALUES) 
predictions = lin_reg.predict(housing_prepared)
lin_mse = mean_squared_error(housing_labels,predictions)
lin_rmse = np.sqrt(lin_mse)
lin_rmse.round()     # typical house price prediction off by $69,000 on the same test data!

# Cross Validation
from sklearn.model_selection import cross_val_score
lin_scores = cross_val_score(lin_reg,housing_prepared,housing_labels,
                         scoring="neg_mean_squared_error",cv=10)
lin_rscores = np.sqrt(-lin_scores)
lin_rscores, lin_rscores.mean().round(), lin_rscores.std().round()


# In[9]:


# Decision Tree Regressor
from sklearn.tree import DecisionTreeRegressor
tree_reg = DecisionTreeRegressor()
tree_reg.fit(housing_prepared,housing_labels)
predictions = tree_reg.predict(housing_prepared)
tree_mse = mean_squared_error(housing_labels,predictions)
tree_rmse = np.sqrt(tree_mse)
tree_rmse.round()     # typical house price prediction is perfect i.e. overfit

# Cross Validation
tree_scores = cross_val_score(tree_reg,housing_prepared,housing_labels,
                         scoring="neg_mean_squared_error",cv=10)
tree_rscores = np.sqrt(-tree_scores)
tree_rscores, tree_rscores.mean().round(), tree_rscores.std().round()


# In[ ]:




