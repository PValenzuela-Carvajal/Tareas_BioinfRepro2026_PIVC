#!/bin/bash

archivo="$1"

# Número de reads
echo "NUMERO DE READS"
echo $(($(zcat "$archivo" | wc -l) / 4))

# Primeras 40 líneas 
echo "PRIMERAS 40 LÍNEAS DEL ARCHIVO"
zcat "$archivo" | head -40

echo "TERCER READ"
# Tercer read
zcat "$archivo" | sed -n '9,12p'

# Primeras 10 calidades
echo "PRIMERAS 10 CALIDADES DEL TERCER READ"
zcat "$archivo" | sed -n '12p' | cut -c1-10
