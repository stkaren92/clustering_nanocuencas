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

Después de esto en la terminal de R ejecutamos:

```
> renv::restore()
```

para instalar las dependencias que usa el proyecto.

## Organización del proyecto

Existen cuatro carpertas principales en el proyecto. Una para guardar los datos 
originales (`00-raw_data`), la segunda (`01-scripts`) donde se encuentran los 
_scripts_ de procesamiento, la tercera (`02-processed_data`) para salvar los 
datos procesados y la cuarta (`03-clustering_output`) para los resultados finales 
de los análisis. 

```
00-raw_data
├── DEM
├── geologia
├── nanocuencas
├── suelos
└── uso_suelo_vegetacion
01-scripts
├── 01-preprocessing.R # Extrae covariable y crea el dataset
└── 02-clustering.R # Hace el análisis de agrupamiento
02-processed_data
└── YYYY-MM-DD_dataset.csv
03-clustering_output
└── YYYY-MM-DD_dendrogram.jpg
└── YYYY-MM-DD_cluster_description.jpg
└── YYYY-MM-DD_nanocuecas_cluster.shp
```

## Metodología

El análisis se realiza en dos etapas: primero se construye una matriz de
covariables por nanocuenca y después se realiza el clustering.

### Preprocesamiento

El _script_ `01-scripts/01-preprocessing.R` lee las nanocuencas y las
covariables espaciales desde `00-raw_data`. Para las covariables vectoriales se calcula el área de intersección en hectáreas entre cada nanocuenca y cada
categoría de la covariable. Además, a partir del _raster_ de elevación, se
estima la covariable _pendiente_ y se extrae la media y desviación estándar por nanocuenca. La tabla resultante se guarda como `02-processed_data/YYYY-MM-DD_dataset.csv`.

### Matriz de distancia

El _script_ `01-scripts/02-clustering.R` lee los datos procesados. Las categorías se normalizan para que sumen 1 por cada covariable dentro de cada nanocuenca. De esta forma, cada covariable se interpreta como una _distribución_ que representa la composición relativa de la nanocuenca.

La distancia entre dos nanocuencas se calcula por tipo de variable:

- Las variables de _distribución_ se comparan con la distancia
  _Jensen-Shannon_.
- Las variables continuas se comparan con la distancia _Hellinger_ a partir de su media y desviación estándar.

La distancia total entre dos nanocuencas es la suma de las distancias calculadas para cada covariable. El resultado es una matriz de distancia donde la entrada _i,j_ representa la distancia entre la nanocuenca _i_ y la _j_.

### Clustering

El _clustering_ se hace a partir de la matriz de distancia con _hierarchical agglomerative clustering_.

El número de grupos se evalúa usando los diagnósticos
WSS _(within cluster sum of squares)_ y _silhouette_. Además de la interpretación del dendrograma.

Los resultados del agrupamiento se guardan en `03-clustering_output` con la
fecha de ejecución y el número de grupos en el nombre del archivo.
