"""Compute area compressibility modulus (K_A) for POPC and all analogs."""

import numpy as np
import pandas as pd

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

BOX_POPC = '../../../../figures/data/SUPP_membrane_parm/area_per_lipid/popc-apl.dat'
BOX_ANALOG = '../../../../figures/data/SUPP_membrane_parm/area_per_lipid/{}-apl.dat'
OUTPUT = '../../../../figures/data/SUPP_membrane_parm/computed_acm.csv'
BLOCKS = [(401, 600), (601, 800), (801, 1000)]  # frame ranges per replica


def acm(x, y, z):
    """Area compressibility modulus K_A (mN/m ≡ dyn/cm)."""
    k_B = 1.380649   # Å²·g·s⁻²·K⁻¹
    T = 303.15        # K
    x, y, z = (np.asarray(v, dtype=float) for v in (x, y, z))
    A = 2.0 * (x * y + x * z + y * z)
    var_A = np.var(A, ddof=0)
    return k_B * T * np.mean(A) / var_A if var_A != 0 else np.nan


def compute_for_file(filepath):
    """Return (mean_KA, std_error) over 9 blocks (3 replicas × 3 time blocks)."""
    df = pd.read_csv(filepath, sep='\t')
    vals = []
    for i in range(1, 4):
        for start, end in BLOCKS:
            block = df[(df['#section'] >= start) & (df['#section'] <= end)]
            vals.append(acm(block[f'x-{i}'], block[f'y-{i}'], block[f'z-{i}']))
    return np.mean(vals), np.std(vals, ddof=1) / np.sqrt(len(vals))


def main():
    rows = []

    # POPC reference
    mean_ka, se_ka = compute_for_file(BOX_POPC)
    rows.append({'name': 'popc', 'label': 'POPC', 'acm': mean_ka, 'std_error': se_ka})

    # Analogs
    for name, label in zip(NAMES, LABELS):
        mean_ka, se_ka = compute_for_file(BOX_ANALOG.format(name))
        rows.append({'name': name, 'label': label, 'acm': mean_ka, 'std_error': se_ka})

    result = pd.DataFrame(rows)
    result.to_csv(OUTPUT, index=False)
    print(f"Saved {OUTPUT}  ({len(result)} rows)")


if __name__ == '__main__':
    main()
