import numpy as np
import pandas as pd
import os

keys = np.arange(1700, 5100, 100)
df_contrasts = {}

module_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(module_dir, 'data')

for key in keys:
    csv_file = os.path.join(data_dir, f'contrast_{key}_new.csv')
    cs = pd.read_csv(csv_file).drop('Unnamed: 0', axis=1)
    df_contrasts[key] = cs


def _bilinear_interp_2d(x, y, xgrid, ygrid, values):
    """
    Bilinear interpolation on a regular 2D grid.

    Parameters
    ----------
    x, y : float or array-like
        Query points.
    xgrid : 1D ndarray
        Grid along x-axis (must be strictly increasing).
    ygrid : 1D ndarray
        Grid along y-axis (must be strictly increasing).
    values : 2D ndarray
        Shape must be (len(ygrid), len(xgrid)).

    Returns
    -------
    ndarray
        Interpolated values with broadcasted shape of x and y.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    xgrid = np.asarray(xgrid, dtype=float)
    ygrid = np.asarray(ygrid, dtype=float)
    values = np.asarray(values, dtype=float)

    if values.shape != (len(ygrid), len(xgrid)):
        raise ValueError(
            f"values shape must be ({len(ygrid)}, {len(xgrid)}), "
            f"but got {values.shape}"
        )

    if np.any(np.diff(xgrid) <= 0):
        raise ValueError("xgrid must be strictly increasing.")
    if np.any(np.diff(ygrid) <= 0):
        raise ValueError("ygrid must be strictly increasing.")

    x_b, y_b = np.broadcast_arrays(x, y)
    x_flat = x_b.ravel()
    y_flat = y_b.ravel()

    if np.any((x_flat < xgrid[0]) | (x_flat > xgrid[-1])):
        raise ValueError("Some x values are outside the interpolation grid.")
    if np.any((y_flat < ygrid[0]) | (y_flat > ygrid[-1])):
        raise ValueError("Some y values are outside the interpolation grid.")

    ix = np.searchsorted(xgrid, x_flat, side='right') - 1
    iy = np.searchsorted(ygrid, y_flat, side='right') - 1

    ix = np.clip(ix, 0, len(xgrid) - 2)
    iy = np.clip(iy, 0, len(ygrid) - 2)

    x0 = xgrid[ix]
    x1 = xgrid[ix + 1]
    y0 = ygrid[iy]
    y1 = ygrid[iy + 1]

    z00 = values[iy, ix]
    z10 = values[iy, ix + 1]
    z01 = values[iy + 1, ix]
    z11 = values[iy + 1, ix + 1]

    tx = (x_flat - x0) / (x1 - x0)
    ty = (y_flat - y0) / (y1 - y0)

    out = (
        (1 - tx) * (1 - ty) * z00
        + tx * (1 - ty) * z10
        + (1 - tx) * ty * z01
        + tx * ty * z11
    )

    return out.reshape(x_b.shape)


def contrast(Ts, Tp, b):
    min_temp, max_temp = 1700, 5000

    Ts = np.asarray(Ts, dtype=float)
    Tp = np.asarray(Tp, dtype=float)

    if np.any(Ts < min_temp) or np.any(Ts > max_temp):
        bad = Ts[(Ts < min_temp) | (Ts > max_temp)]
        raise ValueError(
            f"Tspot value(s) out of bound. Values: {bad} "
            f"are out of bound ({min_temp} to {max_temp} K)!"
        )

    if np.any(Tp < min_temp) or np.any(Tp > max_temp):
        bad = Tp[(Tp < min_temp) | (Tp > max_temp)]
        raise ValueError(
            f"Tphot value(s) out of bound. Values: {bad} "
            f"are out of bound ({min_temp} to {max_temp} K)!"
        )

    grid_tp = np.array(sorted(df_contrasts.keys()), dtype=float)
    grid_ts = np.array(df_contrasts[grid_tp[0]]['Tspot'], dtype=float)

    if b not in df_contrasts[grid_tp[0]].columns:
        raise ValueError(f"Unknown filter '{b}'.")

    # shape = (len(grid_tp), len(grid_ts))
    contrast_values = np.array(
        [df_contrasts[key][b].to_numpy(dtype=float) for key in grid_tp]
    )

    result = _bilinear_interp_2d(
        x=Ts,
        y=Tp,
        xgrid=grid_ts,
        ygrid=grid_tp,
        values=contrast_values,
    )

    return result.item() if result.shape == () else result


def estimate_Tspot_emp(Tp, model='Herbst2'):
    """
    Empirical relations from Herbst+21.

    model 'Herbst1' : Eq.(4) in Herbst+21
    model 'Herbst2' : Eq.(5) in Herbst+21
    """

    Tp = np.asarray(Tp, dtype=float)

    w1 = -3.58e-5

    if model == 'Herbst1':
        w2 = 0.801
        e_w2 = 0.065
        w3 = 666.5
        e_w3 = 280.3
    elif model == 'Herbst2':
        w2 = 1.0188
        e_w2 = 0.068
        w3 = -239.3
        e_w3 = 317.8
    else:
        raise ValueError("model must be 'Herbst1' or 'Herbst2'.")

    Ts = w1 * Tp**2 + w2 * Tp + w3
    Ts_low = w1 * Tp**2 + (w2 - e_w2) * Tp + (w3 - e_w3)
    Ts_up = w1 * Tp**2 + (w2 + e_w2) * Tp + (w3 + e_w3)

    return np.round([Ts_low, Ts, Ts_up]).astype(int)
