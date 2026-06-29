########SCRIPT TO EXTRACT DATAS################
import pandas as pd
import numpy as np
import glob
from functools import reduce

# Lire les fichiers trajectory{traj}-{bloc}.dat
file_list = sorted(glob.glob("trajectory*-*.dat"))
print("Fichiers trouvés :", file_list)

# Liste pour stocker les DataFrames transformés
dfs = []
col_names = []

for filename in file_list:
    # Extraire traj et bloc du nom: trajectory{traj}-{bloc}.dat
    base = filename.replace('.dat', '').replace('trajectory', '')
    parts = base.split('-')
    traj = parts[0]
    bloc = parts[1]+'-'+parts[2]

    col_name = f"dens_traj{traj}_{bloc}"
    col_names.append(col_name)

    # Lecture du fichier (pas d'en‑tête, séparation par espaces)
    df = pd.read_csv(filename, sep=r'\s+', header=None)

    # Ne garder que la colonne z (colonne 0) et la densité (dernière colonne)
    df_tmp = df[[0, df.shape[1] - 1]].copy()
    df_tmp.columns = ['z', col_name]

    dfs.append(df_tmp)

# Fusionner tous les DataFrames sur la colonne 'z'
merged = reduce(lambda left, right: pd.merge(left, right, on='z'), dfs)

# Calcul de la densité moyenne et de l'erreur-type (SE) sur les 9 blocs
dens_cols = [col for col in merged.columns if col.startswith("dens_")]
n = len(dens_cols)
merged['dens_mean'] = merged[dens_cols].mean(axis=1)
merged['se']   = merged[dens_cols].std(axis=1, ddof=1) / np.sqrt(n)

# Sauvegarde du résultat
output_filename = "trajectory-dens.dat"
merged.to_csv(output_filename, sep="\t", index=False, float_format="%.6f")
print(f"Fichier '{output_filename}' généré avec succès.")



