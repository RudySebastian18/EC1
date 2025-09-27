import itertools

def generar_estereoisomeros(smiles: str):
    posiciones = []
    i = 0
    while i < len(smiles):
        if smiles[i] == "@":
            if i + 1 < len(smiles) and smiles[i+1] == "@":
                posiciones.append((i, True))  # ya es @@
                i += 2
            else:
                posiciones.append((i, False)) # es @ simple
                i += 1
        else:
            i += 1

    n = len(posiciones)
    print(f"🔎 Se encontraron {n} centros quirales (@).")

    # Verificación: aceptar 1, 2 o 3; rechazar > 3
    if n == 0:
        print("⚠️ El SMILES no tiene centros quirales. No se generarán isómeros.")
        return []
    elif n > 3:
        print("❌ El SMILES tiene más de 3 centros quirales. No se generarán isómeros.")
        return []

    # Generar todas las combinaciones posibles
    combinaciones = list(itertools.product(["@", "@@"], repeat=n))

    resultados = []
    for comb in combinaciones:
        chars = list(smiles)
        offset = 0
        for (pos, era_doble), val in zip(posiciones, comb):
            real_pos = pos + offset
            if era_doble:
                chars[real_pos:real_pos+2] = list(val)
                offset += len(val) - 2
            else:
                chars[real_pos:real_pos+1] = list(val)
                offset += len(val) - 1
        resultados.append("".join(chars))

    return resultados


# -------------------
# INTERACTIVO
# -------------------
smiles = input("👉 Ingresa el código SMILES: ")

isomeros = generar_estereoisomeros(smiles)

if isomeros:
    print(f"\n✅ Total estereoisómeros generados: {len(isomeros)}")
    print("Ejemplos:")
    for s in isomeros[:5]:
        print(s)
import itertools
from google.colab import files

def generar_estereoisomeros(smiles: str):
    posiciones = []
    i = 0
    while i < len(smiles):
        if smiles[i] == "@":
            if i + 1 < len(smiles) and smiles[i+1] == "@":
                posiciones.append((i, True))  # ya es @@
                i += 2
            else:
                posiciones.append((i, False)) # es @ simple
                i += 1
        else:
            i += 1

    n = len(posiciones)
    print(f"🔎 Se encontraron {n} centros quirales (@).")

    combinaciones = list(itertools.product(["@", "@@"], repeat=n))

    resultados = []
    for comb in combinaciones:
        chars = list(smiles)
        offset = 0
        for (pos, era_doble), val in zip(posiciones, comb):
            real_pos = pos + offset
            if era_doble:
                chars[real_pos:real_pos+2] = list(val)
                offset += len(val) - 2
            else:
                chars[real_pos:real_pos+1] = list(val)
                offset += len(val) - 1
        resultados.append("".join(chars))

    return resultados


# -------------------
# INTERACTIVO
# -------------------
smiles = input("👉 Ingresa el código SMILES: ")

isomeros = generar_estereoisomeros(smiles)

print(f"\n✅ Total estereoisómeros generados: {len(isomeros)}")
print("Ejemplos:")
for s in isomeros[:5]:
    print(s)

# Guardar archivo
filename = "archivo.smi"
with open(filename, "w") as f:
    for s in isomeros:
        f.write(s + "\n")

print(f"\n📁 Archivo '{filename}' guardado con éxito.")

# Descargar archivo en Colab
files.download(filename)
# ========================
# 1. Instalar RDKit en Colab
# ========================
#!pip install rdkit

# ========================
# 2. Importar librerías
# ========================
from rdkit import Chem
from rdkit.Chem import AllChem
from google.colab import files
import os, zipfile

# ========================
# 3. Subir archivo.smi
# ========================
uploaded = files.upload()  # selecciona tu archivo.smi
input_file = "archivo.smi"

output_folder = "xyz_files_rdkit"
zip_name = "xyz_results_rdkit.zip"
os.makedirs(output_folder, exist_ok=True)

# ========================
# 4. Función para convertir SMILES → XYZ
# ========================
def smiles_to_xyz(smiles, filename):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"❌ Error: SMILES inválido {smiles}")
        return

    # Agregar hidrógenos
    mol = Chem.AddHs(mol)

    # Generar geometría inicial con ETKDG
    params = AllChem.ETKDGv3()
    params.randomSeed = 42  # reproducible
    if AllChem.EmbedMolecule(mol, params) != 0:
        print(f"⚠️ No se pudo generar conformación 3D para {smiles}")
        return

    # Optimizar con MMFF94 (si está disponible)
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol)
    else:
        AllChem.UFFOptimizeMolecule(mol)

    # Guardar en formato XYZ
    conf = mol.GetConformer()
    with open(filename, "w") as f:
        f.write(f"{mol.GetNumAtoms()}\n{smiles}\n")
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            f.write(f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n")

# ========================
# 5. Procesar cada línea del archivo.smi
# ========================
with open(input_file, "r") as f:
    smiles_list = [line.strip() for line in f if line.strip()]

print(f"🔎 Se encontraron {len(smiles_list)} moléculas")

for i, smi in enumerate(smiles_list, start=1):
    out_file = os.path.join(output_folder, f"mol_{i}.xyz")
    smiles_to_xyz(smi, out_file)

# ========================
# 6. Comprimir en ZIP y descargar
# ========================
with zipfile.ZipFile(zip_name, "w") as zipf:
    for file in os.listdir(output_folder):
        zipf.write(os.path.join(output_folder, file), file)

files.download(zip_name)

print(f"\n📦 Archivo '{zip_name}' generado con éxito 🎉")
