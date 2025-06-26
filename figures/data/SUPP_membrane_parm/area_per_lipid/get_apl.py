########SCRIPT TO EXTRACT DATAS################
import pandas as pd
import glob

# Lire les fichiers cell*.dat
file_list = sorted(glob.glob("cell*.dat"))
print("Fichiers trouvés :", file_list)

# Liste pour stocker les dataframes transformés
dfs = []

for i, filename in enumerate(file_list, start=1):
    df = pd.read_csv(filename, sep='\s+')

    # Renommer colonnes
    df = df.rename(columns={"#section": "section"})
    df[f"x{i}"] = df["x"]
    df[f"y{i}"] = df["y"]
    df[f"z{i}"] = df["z"]
    df[f"apl{i}"] = (df[f"x{i}"] * df[f"y{i}"]) / 32

    dfs.append(df[["section", f"x{i}", f"y{i}", f"z{i}", f"apl{i}"]])

# Fusion progressive sur 'section'
merged = dfs[0]
for df in dfs[1:]:
    merged = pd.merge(merged, df, on="section", how="outer")

# Trier par section et renommer pour sortie
merged = merged.sort_values(by="section")
merged = merged.rename(columns={"section": "#section"})

# Sauvegarde
merged.to_csv("FILENAME-apl.dat", sep="\t", index=False, float_format="%.6f")


