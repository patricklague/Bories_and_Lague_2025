########SCRIPT TO EXTRACT DATAS################
import pandas as pd
import numpy as np
import glob
from functools import reduce

# Lire les fichiers scd*-t*-*.dat
file_list = sorted(glob.glob("scd*-t*-*.dat"))
print("Fichiers trouvés :", file_list)

# Liste pour stocker les DataFrames transformés
dfs = []
col_names = []

for filename in file_list:
    # Extraire traj et bloc du nom de fichier: scd{chain}-t{traj}-{bloc}.dat
    parts = filename.replace('.dat', '').split('-')
    traj = parts[1].replace('t', '')
    bloc = parts[2]+'-'+parts[3]

    col_name = f"SCD_traj{traj}_{bloc}"
    col_names.append(col_name)

    # Lecture du fichier
    df = pd.read_csv(filename, sep=r'\s+')
    df_tmp = df[['Carbon', '-SCD']].copy()
    df_tmp = df_tmp.rename(columns={'-SCD': col_name})

    dfs.append(df_tmp)

# Fusionner tous les DataFrames sur la colonne 'Carbon'
merged = reduce(lambda left, right: pd.merge(left, right, on='Carbon'), dfs)

# Calcul de la moyenne et de l'erreur-type (SE) sur les 9 blocs
scd_cols = [col for col in merged.columns if col.startswith("SCD_")]
n = len(scd_cols)
merged['scd_mean'] = merged[scd_cols].mean(axis=1)
merged['se']   = merged[scd_cols].std(axis=1, ddof=1) / np.sqrt(n)

# Sauvegarde du résultat
output_filename = "trajectory_scd.dat"
merged.to_csv(output_filename, sep="\t", index=False, float_format="%.6f")
print(f"Fichier '{output_filename}' généré avec succès.")



