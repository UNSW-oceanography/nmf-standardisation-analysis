"""
#################################################################
File: PickleStuff.py

Author: Michael Hemming 
Date: 19/03/2024
Description: Functions to load and save data as a pickle
#################################################################
"""

# %% -------------------------------------------------------------------------
# Import necessary libraries

import pickle

# %% -------------------------------------------------------------------------
# Functions

def PickleSave(file_path,data2save):
    
    print('saving data as a pickle')
    
    with open(file_path, 'wb') as file:
        pickle.dump(data2save, file)
    
def PickleLoad(file_path):
    
    print('Loading pickled data')
    
    with open(file_path, 'rb') as file:
        data = pickle.load(file)
    
    return data