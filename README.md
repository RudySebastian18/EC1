# 🧬 Estereoisomería en Moléculas - Análisis con Streamlit y RDKit  

Este proyecto es una aplicación web construida con **Streamlit** que permite:  

1. ✅ Detectar si una molécula (dada en formato **SMILES**) presenta **estereoisomería**.  
2. 🔄 Generar **todos los estereoisómeros posibles** en código **SMILES**.  
3. 📂 Exportar las moléculas en formato **XYZ** para obtener las coordenadas cartesianas de sus átomos.  

El objetivo es proporcionar una herramienta sencilla, interactiva y didáctica para **química computacional** y **química orgánica**.  

---

## 🚀 Características  

- 📌 **Entrada de moléculas** mediante cadenas SMILES.  
- 🔍 **Detección de estereocentros** y evaluación de estereoisomería.  
- 🧬 **Enumeración automática de estereoisómeros** usando RDKit.  
- 📦 **Exportación a formato `.xyz`**, listo para usarse en programas de química computacional.  
- 🌐 Interfaz web desarrollada con **Streamlit** para fácil interacción.  

---

## 📂 Estructura del Repositorio  

```bash
├── app.py              # Aplicación principal de Streamlit
├── requirements.txt    # Dependencias del proyecto
├── images/             # Recursos estáticos (imágenes, íconos)
└── README.md           # Documentación del proyecto

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

