import random as rnd
import numpy as np
import pandas as pd
class CP_Util():
    def __init__(self,filepath='./'):
        self.FILEPATH  =filepath
        self.NUCLEI_FN      ='Process100_Nucleus.csv'
        self.RBC_FN         ='Process100_MergeRBC.csv'
        self.PATCH_FN       ='Process100_Image.csv'
        self.RBC_ROLLUP_FN  ='RBC_Rollup.csv'
        self.NUC_ROLLUP_FN  ='Nucleus_Rollup.csv'
        self.TEST_SET_ASIDE =0.20
        self.num_patches    =0
        self.test_patches   =None
        self.train_patches  =None
        self.reproducible   = rnd.Random(123456)
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
        participants = []
        tumor_prefix = []
        tumor_x      = []
        tumor_y      = []
        for index,row in patch_info.iterrows():
            tumor = row['FileName']
            delim = tumor.index('_')
            prefix = tumor[:delim]
            suffix = tumor[delim+1:]
            x = suffix[:suffix.index('_')]
            y = suffix[suffix.index('_')+1:suffix.index('.')]
            rest = prefix
            delim = rest.index('-')
            rest = rest[delim+1:] # remove project
            delim = rest.index('-')
            rest = rest[delim+1:] # remove tissue
            delim = rest.index('-')
            participant = rest[:delim] # remove everything after participant
            #
            participants.append(participant)
            tumor_prefix.append(prefix)
            tumor_x.append(x)
            tumor_y.append(y)
        patch_info.insert(0,'Participant',participants)
        patch_info['TumorName']=tumor_prefix
        patch_info['PatchX']=tumor_x
        patch_info['PatchY']=tumor_y
        return patch_info
    def subset_(self,df,col,tumors):
        return df.loc[df[col].isin(tumors)]
    def train_test_split(self):
        patch_info = self._load_patches()  # type dataframe
        patient_names = patch_info['Participant'].unique()  # type ndarray
        num_patients = len(patient_names)
        population = range(num_patients)
        test_size = int(num_patients*self.TEST_SET_ASIDE+0.5)
        test_indices = self.reproducible.sample(population,test_size)
        train_indices = np.setdiff1d(population,test_indices)
        test_names = patient_names[test_indices]
        train_names = patient_names[train_indices]
        self.test_patches = self.subset_(patch_info,'Participant',test_names)
        self.train_patches = self.subset_(patch_info,'Participant',train_names)
        print('Num patients in test/train sets:',len(test_indices),len(train_indices))
        print('Num patches in test/train sets:',len(self.test_patches),len(self.train_patches))
    def get_train_patches(self) -> pd.DataFrame:
        return self.train_patches  
    def get_test_patches(self) -> pd.DataFrame:
        return self.test_patches  
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
    def get_nuclei(self,test_set=False):
        return self.get_objects_(self.NUCLEI_FN,test_set)
    def get_RBC(self,test_set=False):
        return self.get_objects_(self.RBC_FN,test_set)
    def get_RBC_rollup(self,test_set=False):
        return self.get_objects_(self.RBC_ROLLUP_FN,test_set)
    def get_nucleus_rollup(self,test_set=False):
        return self.get_objects_(self.NUC_ROLLUP_FN,test_set)
    def get_objects_(self,FN,test_set):
        good_patches = self.get_train_patches()
        if test_set:
            good_patches = self.get_test_patches()
        filename = self.FILEPATH+FN
        object_info = pd.read_csv(filename)
        column_rename = {'ImageNumber':'PatchNumber'}
        object_info = object_info.rename(column_rename,axis=1)
        object_info = object_info.set_index('PatchNumber')
        good_objects = object_info[object_info.index.isin(good_patches.index)]
        return good_objects
