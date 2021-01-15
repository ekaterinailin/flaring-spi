# utf-8
# python3
# basics for jupyter notebooks

import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

import matplotlib 
matplotlib.rc('xtick', labelsize=14) 
matplotlib.rc('ytick', labelsize=14) 

font = {'weight' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)