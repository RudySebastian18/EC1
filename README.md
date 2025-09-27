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
├── utils/              # Funciones auxiliares (detección, conversión, exportación)
├── outputs/            # Carpeta donde se guardan los archivos .xyz generados
├── assets/             # Recursos estáticos (imágenes, íconos)
└── README.md           # Documentación del proyecto
