import pandas as pd
import numpy as np
import glob
from functools import reduce

# Liste des fichiers *‑PROFILE.dat
file_list = sorted(glob.glob("*-PROFILE.dat"))
print("Fichiers trouvés :", file_list)

# Liste pour stocker les DataFrames
dfs = []

for filename in file_list:
    # Nom du système avant le tiret
    systeme_name = filename.split("-")[0]
    
    # Lecture du profil : deux colonnes (z et densité moyenne)
    # Si vos fichiers ont un en‑tête, remplacez header=None et names=… par header=0
    df = pd.read_csv(filename, sep=r'\s+')
    
    # Renommer la colonne dens_moy avec le nom du système
    df = df[['z', 'dens_moy']].rename(columns={'dens_moy': systeme_name})
    
    # Ajouter à la liste
    dfs.append(df)

# Fusionner tous les DataFrames sur la colonne 'z' (outer pour inclure tous les z)
summary_df = reduce(lambda left, right: pd.merge(left, right, on='z', how='outer'), dfs)

# Remplacer tous les NaN par 0
summary_df = summary_df.fillna(0)

# Sauvegarde du tableau résumé
summary_df.to_csv("all_PROFILE.dat", sep="\t", index=False, float_format="%.4f")
print("Fichier 'all_PROFILE.dat' généré avec succès.")


