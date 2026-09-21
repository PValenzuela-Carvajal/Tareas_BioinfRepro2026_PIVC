## Buena práctica

# que su WD sea DONDE VIVE EL SCRIPT
# hacer esto con 
# a) abrir R dando doble click en el script
# b) En el menú de RStudio: Session > Set Working directory > To surce file location

# cargar el archivo
# Prac_Uni3/maices/meta/maizteocintle_SNP50k_meta_extended.txt

tabla_maices<- read.delim(file="../meta/maizteocintle_SNP50k_meta_extended.txt",
               header=TRUE, sep = "\t", fileEncoding = "latin1")

#¿Qué tipo de objeto creamos al cargar la base?
class(tabla_maices)
#¿Cómo se ven las primeras 6 líneas del archivo?
head(tabla_maices)
#¿Cuántas muestras hay?
nrow(tabla_maices)
# ¿De cuántos estados se tienen muestras?
unique(tabla_maices$Estado)
length(unique(tabla_maices$Estado))
# ¿Cuántas muestras fueron colectadas antes de 1980?
sum(tabla_maices$AÌ.o._de_colecta[!is.na(tabla_maices$AÌ.o._de_colecta)]<1980)
# ¿Cuántas muestras hay de cada raza?
table(tabla_maices$Raza)
# En promedio ¿a qué altitud fueron colectadas las muestras?
mean(tabla_maices$Altitud)
# ¿Y a qué altitud máxima y mínima fueron colectadas?
c(min=min(tabla_maices$Altitud),
  max=max(tabla_maices$Altitud))
## Crea una nueva df de datos sólo con las muestras de la raza Olotillo

Olotillo_df<-tabla_maices[tabla_maices$Raza=="Olotillo",]
class(Olotillo_df)

## Crea una nueva df de datos sólo con las muestras de la raza Reventador, Jala y Ancho

Razas_df <- tabla_maices[tabla_maices$Raza %in% c("Ancho", "Jala", "Reventador"), ]
Razas_df
## Escribe la matriz anterior a un archivo llamado "submat.cvs" en /meta.
write.csv(Razas_df, file = "../meta/submat.csv", row.names = FALSE)
