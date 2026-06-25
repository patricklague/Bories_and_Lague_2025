"""Compute order-parameter ABC deviation (%) for each analog vs POPC."""

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

NAMES = [
    'sca', 'scv', 'scl', 'sci', 'scc', 'scm', 'scs', 'sct', 'scn', 'scq',
    'scf', 'scy', 'scw', 'scp', 'glyd', 'sche', 'schd', 'scdn', 'scen', 'sckn', 'scrn',
    'schp', 'scd', 'sce', 'sccm', 'scym', 'sck', 'scr', 'scrn-1', 'scw-1',
    'scf-1', 'scy-1', 'scv-1', 'sci-1', 'scl-1', 'scm-1'
]
LABELS = [
    'ALA', 'VAL', 'LEU', 'ILE', 'CYS', 'MET', 'SER', 'THR', 'ASN', 'GLN',
    'PHE', 'TYR', 'TRP', 'PRO', 'GLYD', 'HSE', 'HSD', 'ASP0', 'GLU0', 'LYS0',
    'ARG0', 'HSP+', 'ASP-', 'GLU-', 'CYS-', 'TYR-', 'LYS+', 'ARG+', 'ARG0 (0.1M)', 'TRP (0.1M)',
    'PHE(0.1M)', 'TYR(0.1M)', 'VAL(0.1M)', 'ILE(0.1M)', 'LEU(0.1M)', 'MET(0.1M)'
]

CHAINS = [2, 3]
CHAIN_LABELS = ['chain_18_1', 'chain_16_0']

ORDER_POPC = '../../../../figures/data/SUPP_membrane_parm/order_parameter/popc-chain{}.dat'    # format with chain number
ORDER_ANALOG = '../../../../figures/data/SUPP_membrane_parm/order_parameter/{}-chain{}.dat'     # format with (analog, chain)
OUTPUT = '../../../../figures/data/SUPP_membrane_parm/computed_order_deviation.csv'

BLOCK_COLS = [f'SCD_traj{t}_{b}' for t in range(1, 4) for b in ['401-600', '601-800', '801-1000']]


def abc_deviation(x_ref, y_ref, x_test, y_test):
    """Area Between Curves normalized by reference area (%)."""
    x_lo = max(x_ref.min(), x_test.min())
    x_hi = min(x_ref.max(), x_test.max())
    n_pts = max(len(x_ref), len(x_test))
    x_common = np.linspace(x_lo, x_hi, n_pts)
    f_ref = interp1d(x_ref, y_ref, kind='linear', fill_value=0, bounds_error=False)
    f_test = interp1d(x_test, y_test, kind='linear', fill_value=0, bounds_error=False)
    P_i = f_ref(x_common)
    A_i = f_test(x_common)
    denom = np.sum(np.abs(P_i))
    return (np.sum(np.abs(P_i - A_i)) / denom) * 100 if denom != 0 else np.nan


def main():
    rows = []

    for name, label in zip(NAMES, LABELS):
        row = {'name': name, 'label': label}

        for chain, clabel in zip(CHAINS, CHAIN_LABELS):
            df_ref = pd.read_csv(ORDER_POPC.format(chain), sep='\t')
            df_t = pd.read_csv(ORDER_ANALOG.format(name, chain), sep='\t')
            xr = df_ref['Carbon'].values
            xt = df_t['Carbon'].values

            devs = []
            for col in BLOCK_COLS:
                yr = df_ref[col].values
                yt = df_t[col].values
                devs.append(abc_deviation(xr, yr, xt, yt))

            row[f'{clabel}_mean'] = np.mean(devs)
            row[f'{clabel}_se'] = np.std(devs, ddof=1) / np.sqrt(len(devs))

        rows.append(row)

    result = pd.DataFrame(rows)
    result.to_csv(OUTPUT, index=False)
    print(f"Saved {OUTPUT}  ({len(result)} rows)")


if __name__ == '__main__':
    main()
