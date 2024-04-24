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
