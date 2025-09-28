# 🧬 Inchiral - Generador de Estereoisómeros

Este proyecto permite **analizar códigos SMILES** para detectar quiralidad, generar todos los estereoisómeros posibles (máximo 3 centros quirales) y exportarlos en formato **SMI** o **XYZ**.  
Además, incluye la opción de visualizar las moléculas en **2D y 3D** de forma interactiva.

---

## 🚀 Funcionalidades principales
- Detección de centros quirales en moléculas.
- Generación automática de estereoisómeros a partir de SMILES.
- Exportación de resultados en:
  - 📄 `.smi` (lista de isómeros)
  - 📦 `.xyz` (conformaciones 3D en ZIP)
- Visualización molecular en:
  - 🖼️ 2D (estructura plana)
  - 🔬 3D (vista interactiva con py3Dmol)

---

## ⚙️ Tecnologías usadas
- [Streamlit](https://streamlit.io/) → Interfaz web interactiva.
- [RDKit](https://www.rdkit.org/) → Procesamiento químico y generación de conformaciones 3D.
- [py3Dmol](https://pypi.org/project/py3Dmol/) → Visualización molecular en 3D.
- [Pillow (PIL)](https://pillow.readthedocs.io/en/stable/) → Manejo de imágenes.

---

## 📦 Instalación
Clona este repositorio y asegúrate de tener instaladas las dependencias:

```bash
pip install streamlit rdkit-pypi py3Dmol pillow

---
## 📄 Licencia del Proyecto

-Este repositorio utiliza la **GNU General Public License v3.0 (GPL v3)**.  
A continuación, se compara con la **MIT License**, que también podría ser pertinente para un proyecto de este tipo:

| Aspecto                        | GNU General Public License v3.0 (GPL v3)            | MIT License                                                                 |
|--------------------------------|-----------------------------------------------------|------------------------------------------------------------------------------|
| **Tipo de licencia**           | Copyleft fuerte (derivados también deben ser GPL).  | Permisiva (permite reutilización incluso en software privativo).             |
| **Uso en proyectos**           | Mantiene el proyecto siempre libre y abierto.       | Facilita máxima difusión, incluso en software cerrado.                       |
| **Distribución de código**     | Obliga a publicar el código fuente de cualquier modificación o derivado. | No obliga a publicar cambios ni derivados, solo mantener atribución.         |
| **Protección de libertad**     | Alta: defiende que siempre se mantenga open source. | Moderada: prioriza la flexibilidad sobre la obligación de compartir.         |
| **Impacto en este repositorio**| Garantiza que mejoras al software sigan siendo libres. | Permitiría que universidades, laboratorios o empresas usen/adapten sin liberar sus cambios. |

---

### ✅ Justificación de la elección
- Se optó por la **GPL v3** porque este es un trabajo académico y científico donde es importante que el conocimiento y las mejoras se compartan de forma abierta.  
- Sin embargo, la **MIT License** también es pertinente si se quisiera priorizar la adopción del código en la comunidad, sin importar si los derivados son cerrados o no.

