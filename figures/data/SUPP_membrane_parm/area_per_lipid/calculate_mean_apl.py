import pandas as pd
import numpy as np
import glob

# Liste des fichiers système
file_list = sorted(glob.glob("*-apl.dat"))

# Liste pour stocker les résultats
results = []

for filename in file_list:
    # Nom du système sans l'extension
    systeme_name = filename.split("-")[0]

    # Lecture du fichier
    df = pd.read_csv(filename, sep='\s+')

    # Filtrage des sections entre 400 et 1000
    df_filtered = df[(df["#section"] >= 400) & (df["#section"] <= 1000)]

    # Récupérer toutes les colonnes "thickness-*"
    apl_cols = [col for col in df_filtered.columns if col.startswith("apl")]

    # Calcul de la moyenne et de l'erreur standard (std / sqrt(n))
    all_values = df_filtered[apl_cols].values.flatten()
    all_values = all_values[~np.isnan(all_values)]  # retirer les éventuels NaNs

    mean_apl = np.mean(all_values)
    std_error = np.std(all_values, ddof=1) / np.sqrt(len(all_values))

    # Ajouter à la liste
    results.append([systeme_name, mean_apl, std_error])

# Création du DataFrame final
summary_df = pd.DataFrame(results, columns=["name", "apl", "std_error"])

# Sauvegarde
summary_df.to_csv("all_apl.dat", sep="\t", index=False, float_format="%.4f")

