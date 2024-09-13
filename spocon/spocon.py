import numpy as np
import pandas as pd
from scipy.interpolate import interp2d
import os

keys = np.arange(1700, 5100, 100)
df_contrasts = {}

module_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(module_dir, 'data')

for key in keys:
    csv_file = os.path.join(data_dir, f'contrast_{key}_new.csv')
    cs = pd.read_csv(csv_file).drop('Unnamed: 0', axis=1)
    df_contrasts[key] = cs

def contrast(Ts, Tp, b):
    keys = np.array(sorted(df_contrasts.keys()))
    ts_values = np.array(df_contrasts[keys[0]]['Tspot'])
    Ts_mesh, Tp_mesh = np.meshgrid(ts_values, keys)

    contrast_values = np.array([df_contrasts[key][b] for key in keys])

    interp_func = interp2d(ts_values, keys, contrast_values, kind='cubic')
    interpolated_values = interp_func(Ts, Tp)

    return interpolated_values

def estimate_Tspot_emp(Tp, model='Herbst2'):
    
    # empirical relations from Herbst+21
    # model 'Herbst1' : Eq.(3) in Herbst+21, derived from the sample in Berdyugina05
    # model 'Herbst2' : Eq.(5) in Herbst+21, derived from the sample in Berdyugina05, Biazzo+06, Valio17
    
    w1 = -3.58e-5
    
    if model == 'Herbst1':
        w2 = 0.801
        e_w2 = 0.065
        w3 = 666.5
        e_w3 = 280.3
    
    if model == 'Herbst2':
        w2 = 1.0188
        e_w2 = 0.068
        w3 = 239.3
        e_w3 = 317.8
    
    Ts = w1 * Tp**2 + w2 * Tp + w3
    Ts_low = w1 * Tp**2 + (w2-e_w2) * Tp + (w3-e_w3)
    Ts_up = w1 * Tp**2 + (w2+e_w2) * Tp + (w3+e_w3)
    
    return Ts_low, Ts, Ts_up

    
