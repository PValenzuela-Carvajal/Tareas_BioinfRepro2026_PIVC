#!/bin/bash 

# Ruta donde estan tus archivos
src=$HOME/research/project 

# Lista de muestras con su nombre
files="sample_01 sample_02 sample_03"

# Ruta de la carpeta donde estan los subscripts
SCRIPT_DIR="./"

#Iniciar pipeline
echo "INICIANDO PIPELINE DE STACKS"


# Paso 1.1: Alineamiento
bash "$SCRIPT_DIR/1.1Alignment_Script.sh" "$src" "$files"
echo "Alineamiento terminado"
# Paso 1.2: pstacks
bash "$SCRIPT_DIR/1.2pstacks.sh" "$src" "$files"
echo "Pstacks terminado"
# Paso 1.3: Catálogo
bash "$SCRIPT_DIR/1.3build_catalogs.sh" "$src" "$files"
echo "Catalogo listo"
# Paso 1.4: Mapeo contra catálogo
bash "$SCRIPT_DIR/1.4sstacks.sh" "$src" "$files"

# Paso 1.5: Estadísticas poblacionales
echo "Calculando population stats"
bash "$SCRIPT_DIR/1.5populations_statistics.sh" "$src"


echo "PIPELINE FINALIZADA CON ÉXITO

