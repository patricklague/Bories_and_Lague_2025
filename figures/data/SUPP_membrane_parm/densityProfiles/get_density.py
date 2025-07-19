########SCRIPT TO EXTRACT DATAS################
import pandas as pd
import numpy as np
import glob
from functools import reduce

# Lire les fichiers trajectory*.dat
file_list = sorted(glob.glob("trajectory*.dat"))
print("Fichiers trouvés :", file_list)

# Liste pour stocker les DataFrames transformés
dfs = []

for i, filename in enumerate(file_list, start=1):
    # Lecture du fichier (pas d'en‑tête, séparation par espaces)
    df = pd.read_csv(filename, sep=r'\s+', header=None)
    
    # Ne garder que la colonne z (colonne 0) et la densité (dernière colonne)
    df_tmp = df[[0, df.shape[1] - 1]].copy()
    df_tmp.columns = ['z', f'dens_{i}']
    
    # Stocker dans la liste
    dfs.append(df_tmp)

# Fusionner tous les DataFrames sur la colonne 'z'
merged = reduce(lambda left, right: pd.merge(left, right, on='z'), dfs)

# Calcul de la densité moyenne et de l'erreur-type (SE)
n = len(dfs)
dens_cols = [f'dens_{i}' for i in range(1, n+1)]
merged['dens_moy'] = merged[dens_cols].mean(axis=1)
merged['se_moy']   = merged[dens_cols].std(axis=1) / np.sqrt(n)

# Sauvegarde du résultat
output_filename = "trajectory-dens.dat"
merged.to_csv(output_filename, sep="\t", index=False, float_format="%.6f")
print(f"Fichier '{output_filename}' généré avec succès.")



