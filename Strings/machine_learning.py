#!/usr/bin/env python
# coding: utf-8

# # Attempt nc/pc classification with K-mers

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# These files could countain K-mers for any K or any combination of K.
X_train=pd.read_pickle("ncRNA.pcRNA.X_train.pkl")
y_train=pd.read_pickle("ncRNA.pcRNA.y_train.pkl")

# Normalize by row sum.
# Effectively convert k-mer counts to k-mer frequencies.
X_norm=X_train.div(X_train.sum(axis=1), axis=0)

# Feature Scaling by column.
from sklearn.preprocessing import StandardScaler
scaler=StandardScaler()
scaled = scaler.fit_transform(X_norm)
X_scaled = pd.DataFrame(scaled,columns=X_train.columns)
X_train = X_scaled
X_scaled = None
X_norm = None

from sklearn.model_selection import cross_val_score
from sklearn.metrics import accuracy_score
from sklearn.linear_model import LogisticRegression
split=20000

from sklearn.svm import SVC
svc = SVC(gamma='auto')
svc.fit( X_train[:split], np.ravel(y_train[:split]) )
y_pred=svc.predict(X_train[split:])
svc_score=svc.score(X_train[split:],y_train[split:]) # train on 2/3, test on 1/3

from sklearn.ensemble import RandomForestClassifier
rfc = RandomForestClassifier()
rfc.fit( X_train[:split], np.ravel(y_train[:split]) )
y_pred=rfc.predict(X_train[split:])
rfc_score=rfc.score(X_train[split:],y_train[split:]) # train on 2/3, test on 1/3

print("%10f %s"%(svc_score,'SupportVectorMachine'))
print("%10f %s"%(rfc_score,'RandomForestClassifier'))
