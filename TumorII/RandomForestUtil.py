import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.metrics import confusion_matrix
class RF_Util:
    def __init__(self):
        self.model=RandomForestClassifier()
    def get_model(self):
        return self.model
    def set_train(self,X,y):
        self.Xtr = X
        self.ytr = y
    def set_validation(self,X,y):
        self.Xval = X
        self.yval = y
    def fit(self):
        self.model=RandomForestClassifier()
        self.model.fit(self.Xtr,self.ytr)
    def cross_validation(self,fold=5):   # shuffle required!
        self.model=RandomForestClassifier()
        # Crashes on large data when n_jobs is set for parallel operation.
        cv_scores = cross_val_score(self.model, self.Xtr, self.ytr, cv=fold, verbose=1)
        return cv_scores
    def validation_accuracy(self):
        ypred = self.model.predict(self.Xval)
        matches = np.count_nonzero(self.yval==ypred)
        accuracy = 100.0 * matches / len(ypred)  
        return accuracy
    def validation_confusion(self):
        ypred = self.model.predict(self.Xval)
        cm = confusion_matrix(self.yval, ypred)
        return cm
    def important_features(self):
        names = self.model.feature_names_in_
        importances = self.model.feature_importances_
        pairs = np.column_stack( (names,importances) )
        top_array = sorted(pairs, key = lambda e:e[1], reverse=True)
        # There must be a way to do this witout a loop!
        top_list = []
        for i in top_array:
             top_list.append((i[1],i[0]))  # 0=feature_name, 1=importance
        top_df = pd.DataFrame(top_list)
        return top_df