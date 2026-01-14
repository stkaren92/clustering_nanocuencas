# Análisis de agrupamiento nanocuenca Bosque de Agua

En este repositorio se encuentran los _script_ necesarios para realizar el
análisis de agrupamiento de la regiones (nanocuencas) utilizando sus 
características (suelos) de composición.

El repositorio corresponde a un proyectos de R el cual maneja sus dependencias
usando [`renv`](https://rstudio.github.io/renv/). 

## Preparación del entorno

Lo más sencillo es usar [RStudio](https://posit.co/download/rstudio-desktop/). 
En RStudio creamos un nuevo proyecto de control de versiones `git`. En la 
dirección del repositorio usamos la de este:

```
https://github.com/stkaren92/clustering_nanocuencas.git
```

y se escoge la carpeta donde se guardará el proyecto. 

Después de esto en la termina de R ejecutamos:

```
> renv::restore()
```

para instalar las dependencias que usa el proyecto.

## Organización del proyecto

Existen tres carpertas principales en el proyecto. Una para guardar los datos 
originales (`00-raw-data`), la segunda (`01-scripts`) donde se encuentran los 
_scripts_ de procesamiento  y la tercera (`02-processed_data`) para salvar los 
resultados de intermedios y finales de los análisis. 

```
00-raw_data
├── DEM
├── geologia
├── matriz_v0.xlsx
├── nanocuencas
├── suelos
└── uso_suelo_vegetacion
01-scripts
├── 01-preprocessing.R # Extrae covariable y crea el dataset
└── 02-clustering.R # Hace el análisis de agrupamiento
02-processed_data
└── 2026-01-12_dataset.csv
```
