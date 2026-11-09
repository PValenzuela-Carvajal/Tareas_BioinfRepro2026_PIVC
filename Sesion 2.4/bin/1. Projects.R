setwd("")

#Cargamos librerias 

library(admixtools)
library(tidyverse)

rm(list = ls())

##Cargamos la data
f2_blocks <- f2_from_precomp("../data/aadr_1000G_f2_all")
metadf <- read.table(
  "../data/v62.0_1240k_public_metadata2.csv",
  header = TRUE,
  sep = ",")


##Proyecto 1
## Definimos target y demas
target<- c("GBR.DG")

sources <- c("Turkey_Marmara_Barcin_N.AG","Russia_Samara_EBA_Yamnaya.AG","Luxembourg_Mesolithic.DG")

outgroups <- c("Mbuti.DG","CHB.DG","Papuan.DG","Russia_UstIshim_IUP.DG","Denisova.DG")

#f3 statistics
res_f3 <- f3(
  f2_blocks,
  pop1 = target,
  pop2 = sources,
  pop3 = sources
)

#Qpwave
res_qpwave <- qpwave(
  f2_blocks,
  left  = c(target, sources),
  right = outgroups
)

#Qpadm
res_qpadm <- qpadm(
  f2_blocks,
  left   = c(target, sources),
  right  = outgroups,
  target = target
)


res_f3
res_qpwave
res_qpadm



## Proyecto 2. Steppe formation model.

#Can the Steppe Early/Middle Bronze Age (Steppe_EMBA) population be modeled as a mixture of an Iranian Neolithic–related group (Iran_N) and an Eastern Hunter-Gatherer (EHG: Sidelkino) population?
  

target2 <- c("Russia_MLBA_Sintashta.AG", "Kazakhstan_Maitan_MLBA_Alakul.AG", "Russia_LBA_Srubnaya_Alakul.SG") # Steppe_EMBA groups

source2 <- c("Iran_GanjDareh_N.AG", "Russia_Sidelkino_HG.SG")

outgroup2 <- c("Mbuti.DG", "CHB.DG", "Papuan.DG", "Russia_UstIshim_IUP.DG", "Denisova.DG")


f3_results <- f3(f2_blocks, pop1 = target2, pop2 = source2[1], pop3 = source2[2]) # Iran_N and EHG
f3_results 

## Test number of ancestry streams explaining Steppe_EMBA + sources

## Model Steppe_EMBA as mixture of Iran_N and EHG
resultados_qpadm <- list()

# El bucle para procesar target
for(poblacion in target2) {
  wave <- qpwave(f2_blocks, left = c(source2,target2), right = outgroup2)
  
  
  resultados_qpadm[[poblacion]] <- qpadm(
    data = f2_blocks,
    left = source2,
    right = outgroup2,
    target = poblacion
  )
}

#Mirar resultados
for (poblacion in target2) {
  print(resultados_qpadm[[poblacion]])}





##Proyecto 3 Población de las americas

# Bucle para todos los targets 

targets3 <- c("Pima.DG","Karitiana.DG","Surui.DG","PEL.DG")
sources3 <- c("USA_Anzick_realigned.SG",
              "USA_Ancient_Beringian.SG",
              "USA_Nevada_SpiritCave_11000BP.SG")
outgroup3 <- c("Mbuti.DG","CHB.DG","Papuan.DG",
               "Russia_UstIshim_IUP.DG","Denisova.DG")  

#Streams
qpwave(f2_blocks, left = c(targets3[1], sources3), right = outgroup3)
qpwave(f2_blocks, left = c(targets3[2], sources3), right = outgroup3)
qpwave(f2_blocks, left = c(targets3[3], sources3), right = outgroup3)
qpwave(f2_blocks, left = c(targets3[4], sources3), right = outgroup3)

#qpadm y qpwave
for (t in targets3) {
  
  # 2-way: Anzick + Beringian
  m2a <- qpadm(f2_blocks, left = sources3[1:2],
               right = outgroup3, target = t)
  print(m2a)
  print(m2a$weights)
  
  # 2-way: Beringian + SpiritCave
  m2b <- qpadm(f2_blocks, left = sources3[2:3],
               right = outgroup3, target = t)
  print(m2b)
  print(m2b$weights)
  
  # 3-way: Anzick + Beringian + SpiritCave
  m3  <- qpadm(f2_blocks, left = sources3,
               right = outgroup3, target = t)
  print(m3)
  print(m3$weights)
}



## Proyecto 4 

target4 <- "IBS.DG" 

sources4 <- c("Spain_SE_Iberia_CA.AG",
              "Algeria_NumidoRoman_Berber.SG", "YRI.DG")#"Algeria_NumidoRoman_Berber.SG") 

outgroups4 <- c("Mbuti.DG", "Han.DG", "Papuan.DG", "Russia_UstIshim_IUP.DG")

f3_iberia <- map_dfr(sources4, function(src) {
  f3(
    f2_blocks,
    pop1 = "Mbuti.DG",
    pop2 = target4,
    pop3 = src
  ) %>%
    mutate(source = src)
})

print(f3_iberia)

# f4 tests

f4_iberia1 <- f4(
  f2_blocks,
  pop1 = target4,
  pop2 = "Spain_SE_Iberia_CA.AG",
  pop3 = "Algeria_NumidoRoman_Berber.SG",
  pop4 = "Mbuti.DG"
)

f4_iberia2 <- f4(
  f2_blocks,
  pop1 = target4,
  pop2 = "Spain_SE_Iberia_CA.AG",
  pop3 = "YRI.DG",
  pop4 = "Mbuti.DG"
)

print(f4_iberia1)
print(f4_iberia2)

# qpWave

wave4 <- qpwave(
  f2_blocks,
  left = c(target4, sources4),
  right = outgroups4
)

print(wave4$rankdrop)

# qpAdm

# Modelo 2-way: European + North African

admix4_2way <- qpadm(
  f2_blocks,
  left = c(target4, sources4[c(1, 2)]),
  right = outgroups4,
  target = target4
)

# Modelo 3-way: European + North African + Sub-Saharan

admix4_3way <- qpadm(
  f2_blocks,
  left = c(target4, sources4[c(1, 2, 3)]),
  right = outgroups4,
  target = target4
)

admix4_2way
admix4_3way