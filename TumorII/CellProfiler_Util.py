import random as rnd
import numpy as np
import pandas as pd
class CP_Util():
    def __init__(self,filepath='./'):
        self.FILEPATH  =filepath
        self.NUCLEI_FN ='Process100_Nucleus.csv'
        self.PATCH_FN  ='Process100_Image.csv'
        self.TEST_SET_ASIDE =0.20
        self.num_patches   =0
        self.test_patches  =None
        self.train_patches =None
    def _load_patches(self):
        """Load dataframe from CellProfiler Image.csv file.
        Change CP's word "Image" to "Patch".
        Parse strings like FILE-NAME_XCOORD_YCOORD.PNG
        example: TCGA-DB-A4XF-01Z-00-DX1_10200_20700.png"""
        filename = self.FILEPATH+self.PATCH_FN
        cols=['ImageNumber','FileName_Tumor']
        image_info = pd.read_csv(filename,usecols=cols)
        self.num_patches = len(image_info)
        column_rename = {'ImageNumber':'PatchNumber','FileName_Tumor':'FileName'}
        patch_info = image_info.rename(column_rename,axis=1)
        patch_info.set_index('PatchNumber',inplace=True)
        tumor_prefix = []
        tumor_x=[]
        tumor_y=[]
        for index,row in patch_info.iterrows():
            tumor = row['FileName']
            delim = tumor.index('_')
            prefix = tumor[:delim]
            suffix = tumor[delim+1:]
            x = suffix[:suffix.index('_')]
            y = suffix[suffix.index('_')+1:suffix.index('.')]
            tumor_prefix.append(prefix)
            tumor_x.append(x)
            tumor_y.append(y)
        patch_info.insert(0,'TumorName',tumor_prefix)
        patch_info['PatchX']=tumor_x
        patch_info['PatchY']=tumor_y
        return patch_info
    def _subset(self,df,tumors):
        return df.loc[df['TumorName'].isin(tumors)]
    def train_test_split(self):
        patch_info = self._load_patches()  # dataframe
        tumor_names = patch_info['TumorName'].unique()  # ndarray
        num_tumors = len(tumor_names)
        population = range(num_tumors)
        test_size = int(num_tumors*self.TEST_SET_ASIDE+0.5)
        test_indices = rnd.sample(population,test_size)
        train_indices = np.setdiff1d(population,test_indices)
        test_tumor_names = tumor_names[test_indices]
        train_tumor_names = tumor_names[train_indices]
        self.test_patches = self._subset(patch_info,test_tumor_names)
        self.train_patches = self._subset(patch_info,train_tumor_names)
    def get_train_patches(self):
        return self.train_patches  # dataframe
    def get_test_patches(self):
        return self.test_patches  # dataframe
    def validate_split(self):
        num_train_patches = len(self.train_patches)
        num_test_patches = len(self.test_patches)
        if not num_train_patches + num_test_patches == self.num_patches:
            raise Exception('Apparent data loss.')
        df1=pd.DataFrame(self.train_patches['TumorName'])
        df2=pd.DataFrame(self.test_patches['TumorName'])
        in_common = df1.merge(df2,how='inner',indicator=False)
        if len(in_common)>0:
            raise Exception('Sets not mutually exclusive.')
