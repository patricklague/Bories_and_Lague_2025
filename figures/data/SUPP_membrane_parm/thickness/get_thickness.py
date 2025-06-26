########SCRIPT TO EXTRACT DATAS################
import pandas as pd
import glob

# Liste triée des fichiers d'entrée
file_list = sorted(glob.glob("thickness*.dat"))
print("Fichiers détectés :", file_list)

# Liste pour stocker les DataFrames
dfs = []

# Lecture et renommage des colonnes
for i, filename in enumerate(file_list, start=1):
    df = pd.read_csv(filename, sep='\s+')
    df = df.rename(columns={"#section": "section", "thickness": f"thickness-{i}"})
    dfs.append(df[["section", f"thickness-{i}"]])

# Fusion progressive sur 'section'
result = dfs[0]
for df in dfs[1:]:
    result = pd.merge(result, df, on="section", how="outer")

# Réordonner les colonnes si nécessaire
result = result.sort_values(by="section")

# Ajouter le '#' à 'section' dans l'en-tête
result = result.rename(columns={"section": "#section"})

# Sauvegarde
result.to_csv("FILENAME-thickness.dat", sep="\t", index=False, float_format="%.2f")

