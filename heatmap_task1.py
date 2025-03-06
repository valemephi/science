import pandas as pd 
import numpy as np
import seaborn as sns
dataframe1 = pd.read_csv("heatmap_task1_data.txt", delimiter = ' ') 
dataframe1.to_csv('heatmap_task1.csv',  
                  index = None) 