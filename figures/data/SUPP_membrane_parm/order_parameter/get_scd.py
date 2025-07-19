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
    df = pd.read_csv(filename, sep=r'\s+')
    
    # Ne garder que la colonne z (colonne 0) et la densité (dernière colonne)
    df_tmp = df[['Carbon', '-SCD']].copy()
    
    # Stocker dans la liste
    dfs.append(df_tmp)

# Fusionner tous les DataFrames sur la colonne 'z'
merged = reduce(lambda left, right: pd.merge(left, right, on='Carbon'), dfs)
merged.rename(columns={'-SCD': '-SCD_3', '-SCD_x': '-SCD_1', '-SCD_y': '-SCD_2'}, inplace=True)

# Calcul de la densité moyenne et de l'erreur-type (SE)
n = len(dfs)
scd_cols = [col for col in merged.columns if col.startswith("-SCD")]
merged['scd_moy'] = merged[scd_cols].mean(axis=1)
merged['se_moy']   = merged[scd_cols].std(axis=1) / np.sqrt(n)

# Sauvegarde du résultat
output_filename = "trajectory_scd.dat"
merged.to_csv(output_filename, sep="\t", index=False, float_format="%.6f")
print(f"Fichier '{output_filename}' généré avec succès.")



