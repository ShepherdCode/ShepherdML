#!/usr/bin/env python
# coding: utf-8

# In[5]:


import pandas as pd
import os

HOUSING_PATH=os.path.join("../../datasets","housing")
HOUSING_FILE=os.path.join(HOUSING_PATH,"housing.csv")

# returns a Pandas DataFrame
def load_data (data_file=HOUSING_FILE):
    data = pd.read_csv(data_file)
    return data

housing=load_data()

# Data exploration
# housing.head()
# housing.info()
# housing["ocean_proximity"].value_counts()
# housing.describe()

# Jupyter directive to draw in Jupyter
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
# housing.hist(bins=50,figsize=(20,10))
# plt.show()

# Create random test and training sets with stratification by income
from sklearn.model_selection import StratifiedShuffleSplit
import numpy as np
# Pandas cut() retuns a column containing a histogram
housing["temp_income_cat"] = pd.cut(housing["median_income"],bins=[0.0,1.5,3.0,4.5,6.0,np.inf],labels=[1,2,3,4,5])
splitter = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
for train_index, test_index in splitter.split(housing, housing["temp_income_cat"]):
    # Pandas DataFrame.loc() retrieves from array by location 
    train_set = housing.loc[train_index]
    test_set = housing.loc[test_index]

# Evaluate statification
# housing["temp_income_cat"].value_counts()/len(housing["temp_income_cat"])
# train_set["temp_income_cat"].value_counts()/len(strat_train_set)

# Remove the temporary column income_category.
for s in (train_set,test_set,housing):
    # Pandas DataFrame.drop() removes a column 
    s.drop("temp_income_cat",axis=1,inplace=True)
  
# Plot a distribution for each column of the training set.
# %matplotlib inline
# import matplotlib.pyplot as plt
# train_set.hist(bins=50,figsize=(20,10))
# plt.show()   # optional in Jupyter

# From now on, look at the training set only
housing = train_set.copy()

# Visualize geographical distribution
housing.plot(kind="scatter",x="longitude",y="latitude",alpha=0.1)
housing.plot(kind="scatter",x="longitude",y="latitude",alpha=0.4,
    s=housing["population"]/100, label="population", figsize=(10,7),
    c="median_house_value", cmap=plt.get_cmap("jet"), colorbar=True)
plt.show()


# In[6]:


# Pearson's
corr_matrix = housing.corr()
corr_matrix["median_house_value"].sort_values()


# In[57]:


# Scatter plots of Pearson's
# Skew (long tails to right) suggests log transformation before training.
from pandas.plotting import scatter_matrix
attributes=["median_house_value","median_income","total_rooms","housing_median_age"]
scatter_matrix(housing[attributes], figsize=(12,8))
plt.show()


# In[60]:


housing.plot(kind="scatter",x="median_income",y="median_house_value",alpha=0.1)
plt.show()


# In[ ]:




