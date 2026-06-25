"""Compute and export bilayer thickness data sorted by analog order."""

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

INPUT = '../../../../figures/data/SUPP_membrane_parm/thickness/all_thicknesses.dat'
OUTPUT = '../../../../figures/data/SUPP_membrane_parm/computed_thickness.csv'


def main():
    df = pd.read_csv(INPUT, sep='\t', usecols=['name', 'thickness', 'std_error'])

    df_popc = df[df['name'].str.lower() == 'popc'].copy()
    df_analogs = df[df['name'].str.lower() != 'popc'].copy()

    order_map = {name.lower(): i for i, name in enumerate(NAMES)}
    df_analogs = (
        df_analogs
        .assign(order=df_analogs['name'].str.lower().map(order_map))
        .dropna(subset=['order'])
        .sort_values('order')
        .drop(columns='order')
        .reset_index(drop=True)
    )
    df_analogs['label'] = LABELS
    
    df_popc['label'] = 'POPC'
    result = pd.concat([df_popc, df_analogs], ignore_index=True)

    result.to_csv(OUTPUT, index=False)
    print(f"Saved {OUTPUT}  ({len(result)} rows)")


if __name__ == '__main__':
    main()
