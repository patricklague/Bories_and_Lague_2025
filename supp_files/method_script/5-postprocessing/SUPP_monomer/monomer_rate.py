##############################################################################
# Monomer / Multimer rate from raw_data, 4.5A cutoff, 9 batches (3 traj x 3 blocs of 20,000 frames)
##############################################################################
import numpy as np
import pandas as pd
import os

RAW_DIR = '../../../../figures/data/distribution_data/raw_data'
OUTPUT = '../../../../figures/data/SUPP_monomer/monomer_rates_45A_9batches.dat'

aa_all = ['GLYD', 'SCP', 'SCA', 'SCV', 'SCL', 'SCI', 'SCC', 'SCM', 'SCS', 'SCT', 'SCN', 'SCQ',
          'SCF', 'SCY', 'SCW', 'SCCM', 'SCYM', 'SCE', 'SCEN', 'SCD',
          'SCDN', 'SCK', 'SCKN', 'SCR', 'SCRN', 'SCHD', 'SCHE', 'SCHP',
          'SCW-1', 'SCRN-1', 'SCV-1', 'SCI-1', 'SCL-1', 'SCF-1', 'SCY-1', 'SCM-1']
batches = [(0, 19999), (20000, 39999), (40000, 59999)]

results = []

for acid in aa_all:
    acid_lower = acid.lower()
    batch_mono_rates = []
    batch_multi_rates = []

    for traj in [1, 2, 3]:
        filepath = os.path.join(RAW_DIR, acid_lower, f'{acid_lower}_contacts_{traj}.dat')
        if not os.path.isfile(filepath):
            print(f"Warning: {filepath} not found, skipping")
            continue

        df = pd.read_csv(filepath, sep=r'\s+')

        for f_start, f_end in batches:
            df_batch = df[(df['#frame'] >= f_start) & (df['#frame'] <= f_end)]
            n_total = len(df_batch)
            if n_total == 0:
                continue
            n_mono = (df_batch['4.5A_cutoff'] == 0).sum()
            mono_rate = n_mono / n_total * 100.0
            multi_rate = 100.0 - mono_rate
            batch_mono_rates.append(mono_rate)
            batch_multi_rates.append(multi_rate)

    if len(batch_mono_rates) > 0:
        mono_mean = np.mean(batch_mono_rates)
        mono_std = np.std(batch_mono_rates, ddof=1) #/ np.sqrt(len(batch_mono_rates))
        multi_mean = np.mean(batch_multi_rates)
        multi_std = np.std(batch_multi_rates, ddof=1) #/ np.sqrt(len(batch_multi_rates))
    else:
        mono_mean = mono_std = multi_mean = multi_std = np.nan

    results.append({
        'SC': acid,
        'mono_mean': mono_mean,
        'mono_std': mono_std,
        'multi_mean': multi_mean,
        'multi_std': multi_std,
        'n_batches': len(batch_mono_rates)
    })

monomer_df = pd.DataFrame(results)
monomer_rates = dict(zip(monomer_df['SC'], monomer_df['mono_mean']))

os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
with open(OUTPUT, 'w', newline='') as f:
    f.write('# Monomer rates (4.5A cutoff, 9 batches)\n')
    monomer_df.to_csv(f, sep='\t', index=False, float_format='%.2f')
print(f"Saved: {OUTPUT}")
