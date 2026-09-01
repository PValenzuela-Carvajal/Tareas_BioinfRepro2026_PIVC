# Respuestas Ejercicio 2

## 1. ¿Cuantos pasos tiene este script?

EL script tiene 5 pasos de analsisi y un preambulo donde se le da la direccion donde trabajar, los pasos de analisis estan delimitados por comentarios con #

1. Alinea las lecturas con gsnap y convierte el resultado a BAM con samtools y borra el archivo SAM intermedio
2. Corre Stacks con los datos generados en el paso anterior 
3. Se utiliza un bucle para crear una lista de archivos que posteriormente se suministraran a cstacks 
4. Construye el catálogo de loci con cstacks y luego busca las muestras contra ese catálogo usando sstacks y se imprimen a un log 
5. Calcula estadisticos poblacionales y exporta en formato phylip y genepop

## 2. ¿Si quisieras correr este script y que funcionara en tu propio equipo, qué línea deberías cambiar y a qué?

SE debe cambiar la linea 3, que es la ruta a la carpeta de trabajo y es donde estan los archivos que utilizara el script: 

```
3 - src=$HOME/research/project
```

Tambien se deben cambiar las comillas cursivas que van mas adelante

## 3. ¿A qué equivale $HOME?

$Home equivale a la ruta absoluta del usuario, es el equivalente al signo `~`

## 4. ¿Qué paso del análisis hace el programa gsnap?

GSnap se encarga del alineamiento a un genoma de referencia, recibe los archivos .fq (FASTQ) y los mapea contra gac_gen_broads1_e64 con un minimo de cobertura del read de 0.9 

## 5. ¿Qué paso del análisis hace el programa gsnap?

1. Toma cada una de las muestras en cada iteracion y las alinea contra el genoma de referencia 
2. Le asigna un ID a cada muestra que recorre en la variable files
3. Construye una cadena de texto, junta los nombres de todas las muestras en una sola linea (guardada en la varible $samps) 
4. 4. Para cada file, guarda un log file
