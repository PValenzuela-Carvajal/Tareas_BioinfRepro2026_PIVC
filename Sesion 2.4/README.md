# Informe: f-statistics y modelado de mestizaje con admixtools2

**Estudiante:** Pablo Valenzuela-Carvajal

**Profesores:** Dr. Ricardo Verdugo y Dra. Constanza de la Fuente

**Curso:** Bioinformática e investigación reproducible para análisis genómicos

**Sesión:** 2.4

**Fecha:** 2026-09-11


## Proyecto 1: ¿Muestran todas las poblaciones europeas actuales la misma mezcla tripartita?

### Objetivo

Evaluar si la población europea actual analizada (`GBR.DG`) puede explicarse mediante un modelo de mezcla con tres fuentes candidatas:

- `Turkey_Marmara_Barcin_N.AG` 
- `Russia_Samara_EBA_Yamnaya.AG` 
- `Luxembourg_Mesolithic.DG` 

y comparar con lo esperado para poblaciones europeas actuales (3-way admixture).

### Materiales y métodos

Se utilizaron bloques `f2` precomputados y los siguientes grupos:

- **Target:** `GBR.DG`
- **Fuentes:** `Turkey_Marmara_Barcin_N.AG`, `Russia_Samara_EBA_Yamnaya.AG`, `Luxembourg_Mesolithic.DG`
- **Outgroups:** `Mbuti.DG`, `CHB.DG`, `Papuan.DG`, `Russia_UstIshim_IUP.DG`, `Denisova.DG`

Se ejecutaron tres tipos de análisis:

1. Estadísticos `f3`
2. `qpWave` para estimar el número de olas de mezcla
3. `qpAdm` para ajustar modelos de mezcla de dos y tres fuentes 

```r
target   <- c("GBR.DG")
sources  <- c("Turkey_Marmara_Barcin_N.AG",
              "Russia_Samara_EBA_Yamnaya.AG",
              "Luxembourg_Mesolithic.DG")
outgroups <- c("Mbuti.DG", "CHB.DG", "Papuan.DG",
               "Russia_UstIshim_IUP.DG", "Denisova.DG")

res_f3 <- f3(f2_blocks, pop1 = target, pop2 = sources, pop3 = sources)
res_qpwave <- qpwave(f2_blocks, left = c(target, sources), right = outgroups)
res_qpadm  <- qpadm(f2_blocks, left = c(target, sources),
                    right = outgroups, target = target)
```

### Resultados

#### 1. f3

Se ejecutó lo siguiente:

```r
f3(f2_blocks, pop1 = target, pop2 = sources, pop3 = sources)
```

Lo que corresponde a un f3-admixture, donde un estadístico f3 negativo aporta evidencia de que la población objetivo posee una ascendencia relacionada con ambas poblaciones de origen. 

**Tabla 1.  :** f3-admixture de `GBR.DG` para cada par de fuentes.

| pop2 (fuente i) | pop3 (fuente j) | est | se | z | p | 
|---|---|---:|---:|---:|---:|
| Luxembourg_Mesolithic.DG | Russia_Samara_EBA_Yamnaya.AG | -0.000753 | 0.000196 | -3.84 | 1.25e-4 | 
| Luxembourg_Mesolithic.DG | Turkey_Marmara_Barcin_N.AG | -0.00282 | 0.000164 | -17.2 | 3.32e-66 | 
| Russia_Samara_EBA_Yamnaya.AG | Turkey_Marmara_Barcin_N.AG | -0.00151 | 0.0000806 | -18.7 | 4.69e-78 |

Los tres pares únicos de fuentes (Luxembourg-Yamnaya, Luxembourg-Barcin, Yamnaya-Barcin) dan `f3` negativo y significativo. En un admixture f3, esto indica de forma directa que `GBR.DG` no se comporta como un simple outgroup respecto a cada par, sino que es compatible con ser una mezcla de las dos fuentes del par.

#### 2. qpWave

Resultado de `rankdrop`:

**Tabla 2.  :** qpWave (`rankdrop`) para `GBR.DG` y tres fuentes.

| f4rank | dof | chisq | p | Interpretación |
|---:|---:|---:|---:|---|
| 2 | 2 | 0.301 | 8.60e-1 | No se rechaza |
| 1 | 6 | 39.8 | 4.94e-7 | Se rechaza |
| 0 | 12 | 316 | 2.25e-60 | Se rechaza |

El modelo de rango 2 es el primer modelo no rechazado (`p = 0.86`). Se rechazan los modelos de rango 1 y 0, lo que indica que, para el conjunto `GBR.DG + Barcin_N + Yamnaya + Luxembourg_Mesolithic`, se requieren al menos dos olas de mezcla para explicar los estadísticos `f4` frente a los outgroups.

#### 3. qpAdm

Modelo completo con tres fuentes:

**Tabla 3.  :** qpAdm de tres fuentes para `GBR.DG`.

| Fuente | Peso | SE | z |
|---|---:|---:|---:|
| Turkey_Marmara_Barcin_N.AG | 0.293 | 0.0565 | 5.18 |
| Russia_Samara_EBA_Yamnaya.AG | 0.834 | 0.102 | 8.16 |
| Luxembourg_Mesolithic.DG | -0.127 | 0.0846 | -1.50 |

La suma de pesos es 1, pero el peso de `Luxembourg_Mesolithic.DG` es negativo, por lo que el modelo de tres fuentes es no factible. Un peso negativo no tiene interpretación biológica directa y suele indicar que la fuente no es necesaria o que el modelo está mal especificado.

El análisis `popdrop` muestra que el mejor modelo factible es el que elimina `Luxembourg_Mesolithic.DG`:

**Tabla 4.  :** `popdrop` con modelos de dos fuentes para `GBR.DG`.

| Modelo | Fuentes incluidas | Pesos | p | 
|---|---|---|---:|
| 001 | Barcin_N + Yamnaya | Barcin_N = 0.284; Yamnaya = 0.716 | 0.243 | 
| 010 | Yamnaya + Luxembourg_Mesolithic | Yamnaya = 0.246; Luxembourg = 0.754 | 1.23e-7 | 
| 100 | Barcin_N + Luxembourg_Mesolithic | Barcin_N = NA; Yamnaya = 1.60; Luxembourg = -0.603 | — | 

El modelo de dos fuentes Barcin_N + Yamnaya es factible y presenta un ajuste aceptable (`p = 0.243`). Estima aproximadamente:

- 28.4% de ascendencia relacionada con `Turkey_Marmara_Barcin_N.AG`
- 71.6% de ascendencia relacionada con `Russia_Samara_EBA_Yamnaya.AG`

Para `GBR.DG`, los datos no apoyan un modelo de tres vías con `Barcin_N + Yamnaya + Luxembourg_Mesolithic`. La inclusión de `Luxembourg_Mesolithic.DG` produce un peso negativo y el modelo completo es no factible. El mejor modelo factible es de dos vías: `Barcin_N` + `Yamnaya`.

Esto no significa necesariamente que `GBR.DG` no tenga ascendencia correspondiente al grupo que `Luxembourg_Mesolithic.DG` represente, sino que, con estas fuentes y estos outgroups, la señal no se estima como un componente positivo, dentro de las posibles causas: que `Luxembourg_Mesolithic.DG` no es la fuente óptima para este modelo o que los outgroups utilizados no permiten resolver bien ese componente.

Con los resultados actuales, no se puede concluir que todas las poblaciones europeas actuales muestren la misma mezcla tripartita. Para `GBR.DG`, el modelo preferido es de dos fuentes (`Barcin_N` + `Yamnaya`), no de tres bajo las fuentes y outgroups especificados. 


## Proyecto 2: Modelo de formación de la estepa (Steppe_EMBA)

### Objetivo

Evaluar si las poblaciones de la estepa de la Edad del Bronce Medio/Tardío (Steppe_EMBA) pueden modelarse como una mezcla de dos fuentes: un grupo relacionado con el Neolítico iraní (`Iran_GanjDareh_N.AG`) y un cazador-recolector oriental (`Russia_Sidelkino_HG.SG`, EHG).

### Materiales y métodos

Se utilizaron bloques `f2` precomputados y los siguientes grupos:

- **Targets:** `Russia_MLBA_Sintashta.AG`, `Kazakhstan_Maitan_MLBA_Alakul.AG`, `Russia_LBA_Srubnaya_Alakul.SG`
- **Fuentes:** `Iran_GanjDareh_N.AG`, `Russia_Sidelkino_HG.SG`
- **Outgroups:** `Mbuti.DG`, `CHB.DG`, `Papuan.DG`, `Russia_UstIshim_IUP.DG`, `Denisova.DG`

Se ejecutaron:

1. `f3` de mezcla con la orientación `f3(target; Iran_N, EHG)`
2. `qpWave` 
3. `qpAdm` para cada target con ambas fuentes y los outgroups

```r
target2   <- c("Russia_MLBA_Sintashta.AG",
               "Kazakhstan_Maitan_MLBA_Alakul.AG",
               "Russia_LBA_Srubnaya_Alakul.SG")
source2   <- c("Iran_GanjDareh_N.AG", "Russia_Sidelkino_HG.SG")
outgroup2 <- c("Mbuti.DG", "CHB.DG", "Papuan.DG",
               "Russia_UstIshim_IUP.DG", "Denisova.DG")

f3_results <- f3(f2_blocks, pop1 = target2,
                 pop2 = source2[1], pop3 = source2[2])

wave2 <- qpwave(f2_blocks, left = c(source2, target2), right = outgroup2)

resultados_qpadm <- list()
for (poblacion in target2) {
  resultados_qpadm[[poblacion]] <- qpadm(
    data    = f2_blocks,
    left    = source2,
    right   = outgroup2,
    target  = poblacion
  )
}
```

### Resultados

#### 1. f3 de mezcla: `f3(target; Iran_N, EHG)`

**Tabla 5.  :** f3(target; Iran_N, EHG) en tres grupos esteparios.

| Target | est | se | z | p |
|---|---:|---:|---:|---:|
| Kazakhstan_Maitan_MLBA_Alakul.AG | 0.00378 | 0.000337 | 11.2 | 2.65e-29 |
| Russia_LBA_Srubnaya_Alakul.SG | -0.00176 | 0.000316 | -5.57 | 2.60e-8 | 
| Russia_MLBA_Sintashta.AG | 0.0000900 | 0.000286 | 0.314 | 0.753 | 

El resultado es heterogéneo, ya que solo `Russia_LBA_Srubnaya_Alakul.SG` muestra una señal significativa de mezcla bajo `f3`, ya que arroja un valor negativo significativo, luego está `Kazakhstan_Maitan_MLBA_Alakul.AG`, que tiene un valor positivo y significativo, lo que indica que, con estas dos fuentes concretas, el modelo de mezcla simple no queda respaldado por `f3` y `Russia_MLBA_Sintashta.AG` queda en una zona intermedia sin señal detectable.

#### 2. qpWave: número de olas de mezcla

**Tabla 6.  :** qpWave conjunto para Iran_N/EHG y tres targets esteparios.

| Rango (`f4rank`) | gl (`dof`) | χ² | p-valor | 
|---:|---:|---:|---:|
| 3 | 1 | 0.460 | 0.497 |
| 2 | 4 | 1.27 | 0.866 | 
| 1 | 9 | 4.85 | 0.847 | 
| 0 | 16 | 146 | 4.99 × 10⁻²³ | 

El análisis conjunto de qpWave no rechaza los modelos de rango 3 (p = 0.497), rango 2 (p = 0.866) ni rango 1 (p = 0.847). En cambio, el modelo de rango 0 es rechazado fuertemente (p = 4.99 × 10⁻²³). Por lo tanto, el rango 1 es el modelo de menor complejidad que no se rechaza, indicando que las relaciones de afinidad genética entre las dos fuentes (Iran_N y EHG) y los tres grupos esteparios pueden explicarse mediante un único componente ancestral diferenciado respecto del conjunto de outgroups.

#### 3. qpAdm: modelo de mezcla de dos fuentes

##### Russia_MLBA_Sintashta.AG

**Tabla 7.  :** qpAdm de dos fuentes para `Russia_MLBA_Sintashta.AG`.

| Fuente | Peso | SE | z |
|---|---:|---:|---:|
| Iran_GanjDareh_N.AG | 0.319 | 0.0683 | 4.66 |
| Russia_Sidelkino_HG.SG | 0.681 | 0.0683 | 9.97 |

`rankdrop`: f4rank = 1, p = 0.461 (modelo no rechazado).

##### Kazakhstan_Maitan_MLBA_Alakul.AG

**Tabla 8.  :** qpAdm de dos fuentes para `Kazakhstan_Maitan_MLBA_Alakul.AG`.

| Fuente | Peso | SE | z |
|---|---:|---:|---:|
| Iran_GanjDareh_N.AG | 0.354 | 0.0727 | 4.87 |
| Russia_Sidelkino_HG.SG | 0.646 | 0.0727 | 8.89 |

`rankdrop`: f4rank = 1, p = 0.319 (modelo no rechazado).

##### Russia_LBA_Srubnaya_Alakul.SG

**Tabla 9.  :** qpAdm de dos fuentes para `Russia_LBA_Srubnaya_Alakul.SG`.

| Fuente | Peso | SE | z |
|---|---:|---:|---:|
| Iran_GanjDareh_N.AG | 0.215 | 0.0793 | 2.71 |
| Russia_Sidelkino_HG.SG | 0.785 | 0.0793 | 9.89 |

`rankdrop`: f4rank = 1, p = 0.482 (modelo no rechazado).

Todos los pesos son positivos, suman 1 y son significativos (z > 2.5), por lo que los tres modelos son biológicamente interpretables.

`qpWave` y `qpAdm` coinciden en respaldar un modelo de dos componentes (Iran_N + EHG) para los tres targets de la estepa, con un único evento de mezcla, en todos los casos el modelo de rango 1 no se rechaza, con pesos positivos y significativos.

Las proporciones estimadas son:

- **Russia_MLBA_Sintashta.AG:** ~32% Iran_N, ~68% EHG
- **Kazakhstan_Maitan_MLBA_Alakul.AG:** ~35% Iran_N, ~65% EHG
- **Russia_LBA_Srubnaya_Alakul.SG:** ~22% Iran_N, ~78% EHG

En los tres casos predomina el componente EHG (tipo `Russia_Sidelkino_HG.SG`), con una contribución iraní-neolítica menor (~20–35%). Los tres grupos comparten la misma estructura cualitativa (mezcla Iran_N + EHG), pero difieren en proporciones.

Los tres grupos de la estepa considerados (`Russia_MLBA_Sintashta.AG`, `Kazakhstan_Maitan_MLBA_Alakul.AG`, `Russia_LBA_Srubnaya_Alakul.SG`) son compatibles con un modelo de mezcla de dos fuentes entre un grupo relacionado con el Neolítico iraní (`Iran_GanjDareh_N.AG`) y un cazador-recolector (`Russia_Sidelkino_HG.SG`), con mayoría EHG (~65–80%) y minoría Iran_N (~20–35%). 

### Conclusión

Los tres grupos esteparios analizados son compatibles con un modelo de mezcla de dos fuentes (Iran_N + EHG) con un único evento de mezcla. El componente EHG domina (~65–78%), con un aporte Iran_N relevante (~22–35%). La hipótesis de formación de la estepa como mezcla entre un grupo relacionado con el Neolítico iraní y un cazador-recolector queda respaldada por el análisis `qpWave` + `qpAdm`, aunque el `f3` aislado produce resultados heterogéneos que sugieren que la elección de fuentes y outgroups puede refinarse.


## Proyecto 3: Poblamiento de América

### Objetivo

Evaluar modelos de poblamiento de las Américas considerando el número de olas o linajes independientes necesarios para explicar la estructura genética de poblaciones americanas (`Pima`, `Karitiana`, `Surui`, `PEL`).

### Materiales y métodos

Se utilizaron bloques `f2` precomputados y los siguientes grupos:

- **Targets:** `Pima.DG`, `Karitiana.DG`, `Surui.DG`, `PEL.DG`
- **Fuentes:** `USA_Anzick_realigned.SG`, `USA_Ancient_Beringian.SG`, `USA_Nevada_SpiritCave_11000BP.SG`
- **Outgroups:** `Mbuti.DG`, `CHB.DG`, `Papuan.DG`, `Russia_UstIshim_IUP.DG`, `Denisova.DG`

Se ejecutaron:

1. `qpWave` por target con `left = c(target, sources)` y `right = outgroups`
2. `qpAdm` de dos vías con las tres combinaciones posibles de pares de fuentes
3. `qpAdm` de tres vías con las tres fuentes simultáneas

```r
targets3   <- c("Pima.DG","Karitiana.DG","Surui.DG","PEL.DG")
sources3   <- c("USA_Anzick_realigned.SG",
                "USA_Ancient_Beringian.SG",
                "USA_Nevada_SpiritCave_11000BP.SG")
outgroup3  <- c("Mbuti.DG","CHB.DG","Papuan.DG",
                "Russia_UstIshim_IUP.DG","Denisova.DG")

for (t in targets3) {
  qpwave(f2_blocks, left = c(t, sources3), right = outgroup3)

  # 2-way: Anzick + Beringian
  m2a <- qpadm(f2_blocks, left = sources3[1:2], right = outgroup3, target = t)
  # 2-way: Beringian + SpiritCave
  m2b <- qpadm(f2_blocks, left = sources3[2:3], right = outgroup3, target = t)
  # 3-way: Anzick + Beringian + SpiritCave
  m3  <- qpadm(f2_blocks, left = sources3, right = outgroup3, target = t)
}
```

### Resultados

#### 1. qpWave: número de olas o linajes independientes

**Tabla 10.  :** qpWave por target americano (rangos 2 y 1).

| Target | rank 2 (p) | rank 1 (p) | Conclusión |
|---|---:|---:|---|
| Pima.DG | 0.771 | 0.359 | Rango 1 suficiente |
| Karitiana.DG | 0.0496 | 0.0436 | Rango 2 necesario |
| Surui.DG | 0.0572 | 0.0424 | Rango 2 necesario |
| PEL.DG | 0.124 | 0.0236 | Rango 2 necesario |

Pima es la única población para la que una única ola de mezcla con dos fuentes es suficiente, mientras que con Karitiana, Surui y PEL se requieren al menos dos olas (`f4rank = 2`), en Karitiana el rango 2 queda al borde del rechazo (`p = 0.0496`), lo que sugiere que incluso el modelo de dos olas podría ser insuficiente con este conjunto de fuentes.

#### 2. qpAdm: modelos de dos y tres vías

##### Pima.DG

**Tabla 11.  :** qpAdm de dos/tres vías para `Pima.DG`.

| Modelo | Pesos | p | 
|---|---|---:|
| Anzick + Beringian | Anzick = -3.09; Beringian = 4.09 | 0.131 | 
| Beringian + SpiritCave | Beringian = 0.103; SpiritCave = 0.897 | 0.797 |
| Anzick + Beringian + SpiritCave | Anzick = -0.195; Beringian = 0.274; SpiritCave = 0.921 | 0.771 |

El mejor modelo factible para Pima es Beringian + SpiritCave, con ~10% Beringian y ~90% SpiritCave. El modelo de tres vías es no factible porque el peso de Anzick es negativo. El modelo Anzick + Beringian produce pesos no biológicos (`-3.09` y `4.09`).

##### Karitiana.DG

**Tabla 12.  :** qpAdm de dos/tres vías para `Karitiana.DG`.

| Modelo | Pesos | p | 
|---|---|---:|
| Anzick + Beringian | Anzick = -7.00; Beringian = 8.00 | 0.109 | 
| Beringian + SpiritCave | Beringian = -0.412; SpiritCave = 1.41 | 0.087 | 
| Anzick + Beringian + SpiritCave | Anzick = 16.8; Beringian = -20.9; SpiritCave = 5.02 | 0.050 | 

Ningún modelo es factible, los pesos del modelo de tres vías son enormes, lo que podría indicar una mala especificación de las fuentes. En `popdrop`, el único modelo factible es de rango 0 (una sola fuente), lo que significa que ninguno de los pares probados explica a Karitiana.

##### Surui.DG

**Tabla 13.  :** qpAdm de dos/tres vías para `Surui.DG`.

| Modelo | Pesos | p | 
|---|---|---:|
| Anzick + Beringian | Anzick = -6.19; Beringian = 7.19 | 0.121 | 
| Beringian + SpiritCave | Beringian = -0.981; SpiritCave = 1.98 | 0.055 | 
| Anzick + Beringian + SpiritCave | Anzick = -3.34; Beringian = 3.45; SpiritCave = 0.898 | 0.057 | 

Surui.DG sigue el mismo patrón que Karitiana.

##### PEL.DG

**Tabla 14.  :** qpAdm de dos/tres vías para `PEL.DG`.

| Modelo | Pesos | p | 
|---|---|---:|
| Anzick + Beringian | Anzick = -5.81; Beringian = 6.81 | 0.123 | 
| Beringian + SpiritCave | Beringian = 2.60; SpiritCave = -1.60 | 0.198 | 
| Anzick + Beringian + SpiritCave | Anzick = -1.09; Beringian = 3.55; SpiritCave = -1.46 | 0.124 | 

Sigue el mismo patrón que Karitiana y Surui. En `popdrop`, el único modelo factible es de rango 0 (una sola fuente), es decir, ninguno de los pares o tríos probados explica a PEL (población contemporánea de Perú).

Los resultados de las dos aproximaciones utilizadas no son convergentes: mientras que qpWave sugiere que Karitiana, Surui y PEL requieren dos olas de mezcla (rango 2) y que Pima se explica con una sola (rango 1), qpAdm solo encuentra un modelo factible para Pima (Beringian + SpiritCave). Para Karitiana, Surui y PEL, ni los pares ni el trío de fuentes producen pesos positivos, lo que podría explicarse por la especificación no correcta de las fuentes o los outgroups. 

Los datos sugieren al menos dos olas o linajes independientes para Karitiana, Surui y PEL (qpWave, rango 2). Para Pima, una sola ola con dos fuentes es suficiente. Esto es compatible con escenarios de múltiples pulsos, pero deben refinarse las fuentes y los outgroups. 

### Conclusión

Los datos analizados son compatibles con un escenario de múltiples olas o linajes en el poblamiento de América: Karitiana, Surui y PEL requieren al menos dos olas (`qpWave` rango 2), mientras que Pima se ajusta a una sola. Sin embargo, solo Pima encuentra un modelo `qpAdm` factible con las tres fuentes americanas consideradas (`~90% SpiritCave + ~10% Beringian`). Para Karitiana, Surui y PEL, ninguno de los pares o del trío de fuentes produce pesos positivos, probablemente por la inadecuación de las fuentes o los outgroups para representar la estructura de estas poblaciones. 


## Proyecto 4: Mestizaje medieval / ibérico

### Objetivo

Evaluar si la población ibérica moderna (`IBS.DG`) puede modelarse como una mezcla de un componente ibérico antiguo (`Spain_SE_Iberia_CA.AG`), un componente norteafricano (`Algeria_NumidoRoman_Berber.SG`) y, opcionalmente, como un componente subsahariano (`YRI.DG`).

### Materiales y métodos

Se utilizaron bloques `f2` precomputados y los siguientes grupos:

- **Target:** `IBS.DG`
- **Fuentes:** `Spain_SE_Iberia_CA.AG`, `Algeria_NumidoRoman_Berber.SG`, `YRI.DG`
- **Outgroups:** `Mbuti.DG`, `Han.DG`, `Papuan.DG`, `Russia_UstIshim_IUP.DG`

> Nota: el diseño usado difiere ligeramente del propuesto en el guion del proyecto. Se usó `Spain_SE_Iberia_CA.AG` (en lugar de un `IBS` medieval como fuente, ya que el que seleccionamos está datado cerca de 7000 años antes del presente) y no se incluyó `Onge.DG` entre los outgroups ya que no se encontraba en la base de datos de Allen Ancient DNA Resource (AADR) utilizada (v62). 

Se realizaron los siguientes análisis: 

1. Outgroup-f3 con `pop1 = "Mbuti.DG"`, `pop2 = target`, `pop3 = fuente`
2. Dos tests `f4` para detectar asimetrías
3. `qpWave` con `left = c(target, sources)` y `right = outgroups`
4. `qpAdm` de 2 vías (`Spain_CA + NA_Berber`) y de 3 vías (`Spain_CA + NA_Berber + YRI`)

```r
target4   <- "IBS.DG"
sources4  <- c("Spain_SE_Iberia_CA.AG",
               "Algeria_NumidoRoman_Berber.SG",
               "YRI.DG")
outgroups4 <- c("Mbuti.DG", "Han.DG", "Papuan.DG",
                "Russia_UstIshim_IUP.DG")

f3_iberia <- map_dfr(sources4, function(src) {
  f3(f2_blocks, pop1 = "Mbuti.DG", pop2 = target4, pop3 = src) %>%
    mutate(source = src)
})

f4_iberia1 <- f4(f2_blocks, pop1 = target4, pop2 = "Spain_SE_Iberia_CA.AG",
                 pop3 = "Algeria_NumidoRoman_Berber.SG", pop4 = "Mbuti.DG")
f4_iberia2 <- f4(f2_blocks, pop1 = target4, pop2 = "Spain_SE_Iberia_CA.AG",
                 pop3 = "YRI.DG", pop4 = "Mbuti.DG")

wave4 <- qpwave(f2_blocks, left = c(target4, sources4), right = outgroups4)

admix4_2way <- qpadm(f2_blocks, left = c(target4, sources4[c(1,2)]),
                     right = outgroups4, target = target4)
admix4_3way <- qpadm(f2_blocks, left = c(target4, sources4),
                     right = outgroups4, target = target4)
```

### Resultados

#### 1. Outgroup-f3 (`f3(Mbuti; IBS, fuente)`)

**Tabla 15.  :** outgroup-f3 de `IBS.DG` con cada fuente.

| Fuente | est | se | z | p |
|---|---:|---:|---:|---:|
| Spain_SE_Iberia_CA.AG | 0.0706 | 0.000476 | 148 | <0.0001 |
| Algeria_NumidoRoman_Berber.SG | 0.0656 | 0.000475 | 138 | <0.0001 |
| YRI.DG | 0.0234 | 0.000268 | 87.4 | <0.0001 |

Las tres fuentes muestran deriva compartida significativa con `IBS.DG`. La magnitud es alta para `Spain_SE_Iberia_CA.AG` (0.071), muy similar para `Algeria_NumidoRoman_Berber.SG` (0.066), y mucho menor para `YRI.DG` (0.023). Este patrón es el esperado si `IBS` tiene ascendencia ibérica antigua, una contribución norteafricana menor y una señal subsahariana aún menor.

#### 2. Tests f4 (asimetrías)

**Tabla 16.  :** tests f4 de asimetría para `IBS.DG`.

| Test | est | se | z | p | 
|---|---:|---:|---:|---:|
| `f4(IBS; Spain_CA, NA_Berber, Mbuti)` | -0.00148 | 0.000117 | -12.6 | 2.46e-36 | 
| `f4(IBS; Spain_CA, YRI, Mbuti)` | 0.0000668 | 0.0000505 | 1.32 | 0.186 | 

El primer `f4` es altamente significativo: `IBS.DG` no se comporta como un simple intermedio entre `Spain_CA` y `NA_Berber` respecto a `Mbuti`. El signo negativo indica que la afinidad `IBS ↔ NA_Berber` no es mayor que la de `Spain_CA ↔ NA_Berber`, es decir, el genoma ibérico moderno no muestra más afinidad norteafricana que el ibérico en este eje. El segundo `f4` no es significativo (`p = 0.186`): no se detecta asimetría entre `Spain_CA` e `YRI` respecto a `Mbuti`, lo que sugiere que la señal subsahariana, si existe, es demasiado débil para resolverse con este diseño.

#### 3. qpWave

**Tabla 17.  :** qpWave para `IBS.DG` con tres fuentes.

| f4rank | dof | chisq | p | 
|---:|---:|---:|---:|
| 2 | 1 | 6.67 | 0.00983 |
| 1 | 4 | 46.7 | 1.74e-9 | 
| 0 | 9 | 11332 | <0.0001 | 

Todos los rangos testables son rechazados, incluido el rango 2 (`p = 0.0098`). Esto implica que, con estas tres fuentes y estos cuatro outgroups, se requieren más de dos olas o linajes independientes para explicar la estructura de `IBS.DG`. Como el rango 3 no es testable con tres fuentes (0 grados de libertad), el resultado significa que el conjunto `Spain_CA + NA_Berber + YRI` no basta para modelar `IBS`, sin embargo, esto puede ocurrir debido a la no inclusión del outgroup `ONGE.DG`.

#### 4. qpAdm de 2 vías: `Spain_CA + NA_Berber`

**Tabla 18.  :** qpAdm de dos vías (`Spain_CA + NA_Berber`) para `IBS.DG`.

| Fuente | Peso | SE | z |
|---|---:|---:|---:|
| Spain_SE_Iberia_CA.AG | 0.919 | 0.0392 | 23.5 |
| Algeria_NumidoRoman_Berber.SG | 0.0808 | 0.0392 | 2.06 |

`rankdrop`: rango 1 rechazado (`p = 1.73e-8`). Los pesos son positivos y significativos, pero el modelo se rechaza globalmente. Es decir: la estimación puntual es ~92% ibérico antiguo y ~8% norteafricano.

#### 5. qpAdm de 3 vías: `Spain_CA + NA_Berber + YRI`

**Tabla 19.  :** qpAdm de tres vías (`Spain_CA + NA_Berber + YRI`) para `IBS.DG`.

| Fuente | Peso | SE | z |
|---|---:|---:|---:|
| Spain_SE_Iberia_CA.AG | 0.157 | 0.313 | 0.501 |
| Algeria_NumidoRoman_Berber.SG | 0.913 | 0.341 | 2.68 |
| YRI.DG | -0.0703 | 0.0278 | -2.52 |

El peso de `YRI.DG` es negativo, por lo que el modelo de tres fuentes es no factible. El análisis `popdrop` demuestra que el mejor modelo de 2 vías factible (sin `YRI`) sigue teniendo `p = 1.66e-8`, es decir, rechazado.

El outgroup-f3 muestra que `IBS.DG` comparte deriva significativa con `Algeria_NumidoRoman_Berber.SG` (`f3 = 0.066`, `z = 138`). El `qpAdm` de 2 vías asigna ~8% de peso a esta fuente con `z = 2.06`. Sin embargo, el ajuste global del modelo es rechazado (`p = 1.7e-8`), lo que sugiere que se debe especificar mejor el modelo, es decir, las fuentes, ya que los tiempos en los que están datados los datos son mucho más allá que la época medieval, lo que sugiere que las fuentes no estarían representando bien lo ocurrido durante ese periodo.

### Conclusión

`IBS.DG` muestra deriva compartida significativa con `Spain_SE_Iberia_CA.AG` (`f3 = 0.071`), `Algeria_NumidoRoman_Berber.SG` (`f3 = 0.066`) y, en menor medida, `YRI.DG` (`f3 = 0.023`). El modelo `qpAdm` de 2 vías estima ~92% ibérico antiguo + ~8% norteafricano con pesos positivos y significativos, pero es rechazado globalmente (`p = 1.7e-8`). El modelo de 3 vías añadiendo `YRI` produce un peso negativo y no es factible. `qpWave` rechaza incluso el rango 2, indicando que se necesitan más de dos olas o componentes para explicar `IBS`. En conjunto, los resultados son compatibles con un aporte norteafricano en el ibérico moderno, y no permiten confirmar ni cuantificar un aporte subsahariano independiente. La inclusión de un componente como fuente adicional, así como proxies ibéricos y norteafricanos más adecuados, es necesaria para obtener un modelo factible y biológicamente interpretable.