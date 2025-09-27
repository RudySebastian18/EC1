import streamlit as st
import itertools
import os
import zipfile
import tempfile
import io
import sys
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import py3Dmol
from PIL import Image
import warnings

# Evitar warnings de RDKit
warnings.filterwarnings('ignore')

# ---------------- FUNCIONES ---------------- #

def mostrar_imagen_2d(smiles: str):
    """Genera imagen 2D de un SMILES"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    img = Draw.MolToImage(mol, size=(300, 300))
    return img

def mostrar_imagen_3d(smiles: str):
    """Genera visualizador 3D con py3Dmol"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.UFFOptimizeMolecule(mol)

    conf = mol.GetConformer()
    xyz = f"{mol.GetNumAtoms()}\n\n"
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        xyz += f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n"

    viewer = py3Dmol.view(width=400, height=400)
    viewer.addModel(xyz, "xyz")
    viewer.setStyle({"stick": {}})
    viewer.setBackgroundColor("white")
    viewer.zoomTo()
    viewer.render()  # <- Importante
    return viewer

def detectar_quiralidad(smiles: str):
    """Detecta si la molécula tiene centros quirales"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "SMILES inválido", []
        
        centros = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        
        if len(centros) == 0:
            return False, "Su molécula no es quiral", []
        else:
            return True, f"Su molécula es quiral. Se detectaron {len(centros)} posibles centros", centros
            
    except Exception as e:
        return False, f"Error al analizar la molécula: {str(e)}", []

def analizar_centros_existentes(smiles: str):
    """Cuenta cuántos centros quirales están definidos con @ o @@ en el SMILES"""
    centros_especificados = 0
    posiciones_at = []
    i = 0
    
    while i < len(smiles):
        if smiles[i] == "@":
            if i + 1 < len(smiles) and smiles[i+1] == "@":
                centros_especificados += 1
                posiciones_at.append(i)
                i += 2
            else:
                centros_especificados += 1
                posiciones_at.append(i)
                i += 1
        else:
            i += 1
    
    return centros_especificados, posiciones_at

def generar_estereoisomeros(smiles: str):
    """Genera combinaciones de estereoisómeros cambiando @ y @@"""
    posiciones = []
    i = 0
    while i < len(smiles):
        if smiles[i] == "@":
            if i + 1 < len(smiles) and smiles[i+1] == "@":
                posiciones.append((i, True))  # ya es @@
                i += 2
            else:
                posiciones.append((i, False))  # es @ simple
                i += 1
        else:
            i += 1
    
    n = len(posiciones)
    
    if n == 0:
        st.warning("⚠️ El SMILES no tiene centros quirales especificados con @ o @@. No se generarán isómeros.")
        return [], n
    elif n > 3:
        st.error("❌ El SMILES tiene más de 3 centros quirales. No se generarán isómeros.")
        return [], n
    
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
    
    return resultados, n

def smiles_to_xyz(smiles, mol_id):
    """Convierte un SMILES a coordenadas XYZ"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, f"❌ Error: SMILES inválido {smiles}"
        
        mol = Chem.AddHs(mol)
        
        params = AllChem.ETKDGv3()
        params.randomSeed = 42  
        
        embed_result = AllChem.EmbedMolecule(mol, params)
        if embed_result != 0:
            params.useRandomCoords = True
            embed_result = AllChem.EmbedMolecule(mol, params)
            if embed_result != 0:
                return None, f"⚠️ No se pudo generar conformación 3D para {smiles}"
        
        try:
            if AllChem.MMFFHasAllMoleculeParams(mol):
                AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            else:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
        except:
            pass
        
        conf = mol.GetConformer()
        xyz_content = f"{mol.GetNumAtoms()}\n{smiles}\n"
        
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            xyz_content += f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n"
        
        return xyz_content, f"✅ Molécula {mol_id} procesada correctamente"
        
    except Exception as e:
        return None, f"❌ Error procesando {smiles}: {str(e)}"

def crear_archivo_zip(archivos_xyz):
    """Crea un ZIP con los archivos XYZ"""
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
        for filename, content in archivos_xyz.items():
            zip_file.writestr(filename, content)
    return zip_buffer.getvalue()

# ---------------- APP PRINCIPAL ---------------- #

def main():
    st.set_page_config(
        page_title="Inchiral - Generador de Estereoisómeros",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    
    st.title("🧬 Generador de Estereoisómeros")
    st.markdown("**Genera todos los estereoisómeros posibles y convierte a formato XYZ**")
    
    with st.sidebar:
        try:
            st.image("imagenes1/inchiral final.png", width=200)
        except:
            st.markdown("**🧬 Inchiral**")
        
        st.markdown("---")
        st.title("ℹ️ Información")
        st.markdown("""
        **Instrucciones:**
        1. Ingresa un código SMILES (con o sin quiralidad especificada)
        2. El sistema detecta automáticamente si la molécula es quiral
        3. Si tiene centros quirales especificados (@ o @@), genera todos los estereoisómeros
        4. Máximo 3 centros quirales para evitar demasiados isómeros
        5. Opcionalmente convierte a formato XYZ para visualización 3D
        """)

    st.subheader("📝 Entrada de Datos")
    smiles_input = st.text_input(
        "👉 Ingresa el código SMILES:",
        placeholder="Ejemplo: C[C@H](O)[C@@H](N)C"
    )
    
    if smiles_input:
        st.subheader("🔍 Análisis de Quiralidad")
        
        es_quiral, mensaje_quiralidad, centros_detectados = detectar_quiralidad(smiles_input)
        centros_especificados, posiciones_at = analizar_centros_existentes(smiles_input)

        col1, col2 = st.columns(2)
        
        with col1:
            if es_quiral:
                st.success(f"✅ {mensaje_quiralidad}")
                if centros_detectados:
                    st.write("**Centros detectados:**")
                    for i, (idx, chirality) in enumerate(centros_detectados):
                        tipo_quiralidad = str(chirality) if chirality else "Sin asignar"
                        st.write(f"• Átomo {idx}: {tipo_quiralidad}")
            else:
                st.warning(mensaje_quiralidad)
        
        with col2:
            if centros_especificados > 0:
                st.success(f"✅ {centros_especificados} centros con @ o @@ especificados")
            else:
                st.warning("⚠️ No hay centros especificados con @ o @@")
        
        # Tabs
        tab1, tab2, tab3, tab4 = st.tabs(["📋 Lista", "💾 Descargar SMI", "🧪 Convertir a XYZ", "🖼️ Visualizar Molécula"])

        # Generación de estereoisómeros
        isomeros, n_centros = [], 0
        if centros_especificados > 0:
            with st.spinner("🔄 Generando estereoisómeros..."):
                isomeros, n_centros = generar_estereoisomeros(smiles_input)

        if isomeros:
            with tab1:
                for i, isomero in enumerate(isomeros, 1):
                    st.code(f"{i}. {isomero}")

            with tab2:
                smi_content = "\n".join(isomeros)
                st.download_button(
                    label="📥 Descargar archivo.smi",
                    data=smi_content,
                    file_name="estereoisomeros.smi",
                    mime="text/plain"
                )

            with tab3:
                if st.button("🚀 Convertir todos a XYZ", type="primary"):
                    archivos_xyz = {}
                    for i, smiles in enumerate(isomeros, 1):
                        xyz_content, mensaje = smiles_to_xyz(smiles, i)
                        if xyz_content:
                            archivos_xyz[f"mol_{i}.xyz"] = xyz_content
                    if archivos_xyz:
                        zip_data = crear_archivo_zip(archivos_xyz)
                        st.download_button(
                            label="📦 Descargar archivos XYZ (ZIP)",
                            data=zip_data,
                            file_name="estereoisomeros_xyz.zip",
                            mime="application/zip"
                        )

        with tab4:
            st.header("🖼️ Visualizador de Moléculas")
            smiles_vis = st.text_input("Introduce un código SMILES para visualizar")
            if smiles_vis:
                # 2D
                img2d = mostrar_imagen_2d(smiles_vis)
                if img2d:
                    st.image(img2d, caption="Estructura 2D", use_column_width=False)

                # 3D
                st.subheader("Vista 3D Interactiva")
                viewer = mostrar_imagen_3d(smiles_vis)
                if viewer:
                    st.components.v1.html(viewer._make_html(), height=500, width=500)
                else:
                    st.warning("⚠️ No se pudo generar la vista 3D para este SMILES")

    st.markdown("---")
    st.markdown(
        """
        <div style='text-align: center'>
            <small>🧬 <strong>Inchiral</strong><br>
            Generador de Estereoisómeros | Desarrollado con Streamlit y RDKit</small>
        </div>
        """,
        unsafe_allow_html=True
    )

if __name__ == "__main__":
    main()
