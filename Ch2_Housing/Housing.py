import pandas as pd
import os

HOUSING_PATH=os.path.join("datasets","housing")
HOUSING_FILE=os.path.join(HOUSING_PATH,"housing.csv")

def load_data (data_file=HOUSING_FILE):
    data = pd.read_csv(data_file)
    return data
