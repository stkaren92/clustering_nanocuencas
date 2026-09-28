# Análisis de agrupamiento de nanocuencas del Bosque de Agua

Este repositorio contiene el flujo de trabajo en R para agrupar nanocuencas del
Bosque de Agua a partir de sus características, en especial del tipo de suelo.

El proyecto usa [`renv`](https://rstudio.github.io/renv/) para registrar y
restaurar las versiones de sus dependencias.

## Preparación del entorno

Se recomienda trabajar desde RStudio mediante el archivo
`clustering_nanocuencas.Rproj`, aunque los scripts también se pueden ejecutar
con `Rscript` desde la raíz del repositorio.

Para clonar el proyecto:

```bash
git clone https://github.com/stkaren92/clustering_nanocuencas.git
cd clustering_nanocuencas
```

Después, restaure las dependencias desde una sesión de R:

```r
renv::restore()
```

Los paquetes geoespaciales, en particular `sf` y `stars`, pueden requerir
bibliotecas del sistema como GDAL, GEOS y PROJ.

## Organización del proyecto

```text
00-raw_data/
├── DEM/
│   └── wc2.1_30s_elev.tif
└── 03_Entrega_21092026/
    └── Shapefiles/
        ├── Estadisticas/
        ├── Nanocuencas/
        ├── Variables/
        ├── categorical_variables_catalog.csv
        ├── stats_variables_catalog.csv
        ├── usv_rename.csv
        └── variables_ordering_catalog.csv
01-scripts/
├── 00-rename_usv.R
├── 01-preprocessing.R
└── 02-clustering.R
02-processed_data/
├── YYYY-MM-DD_dataset.csv
└── aoi_shp/
    └── YYYY-MM-DD_dataset.gpkg
03-clustering_output/
├── YYYY-MM-DD_wws_plot.png
├── YYYY-MM-DD_silhouette_plot.png
├── YYYY-MM-DD_dendrogram_<k>_k.jpg
├── YYYY-MM-DD_cluster_description_<k>_k.jpg
├── YYYY-MM-DD_cluster_mean_boxplots_<k>_k.jpg
└── YYYY-MM-DD_nanocuecas_cluster_<k>_k.{shp,shx,dbf,prj}
```

La ruta de entrada utilizada actualmente por los scripts es
`00-raw_data/03_Entrega_21092026/Shapefiles`.

## Variables

Las variables se dividen en composiciones categóricas y variables continuas.
Los prefijos de las columnas se definen mediante catálogos, por lo que pueden
actualizarse sin modificar la lógica principal de procesamiento.

### Variables de distribución

Para cada nanocuenca se calcula el área en hectáreas ocupada por cada categoría.
Posteriormente, estas áreas se normalizan para formar una distribución que suma
1 dentro de cada variable.

| Prefijo | Variable | Archivo y campo de origen |
|---|---|---|
| `glg` | Geología | `r250k_ccl_dissolveTipo.shp`, campo `TIPO` |
| `suelo1` | Unidad de suelo principal | `contedafo_cw.shp`, campo `NOM_SUE1` |
| `suelo2` | Unidad de suelo secundaria | `contedafo_cw.shp`, campo `NOM_SUE2` |
| `usv` | Uso del suelo y vegetación reclasificado | `usv250s7cw_renamed.shp`, campo `ETIQUETA` |
| `sueloinifap` | Descripción y calificador del suelo | `eda251mcw.shp`, campo `DESCRIPCIO` |

### Variables representadas por parámetros

Cada una de estas variables se representa por un par de columnas con sufijos
`_mean` y `_sd`.

| Prefijo | Variable | Origen |
|---|---|---|
| `alt` | Altitud | Estadísticas precomputadas en `shtv4_stats_altitud.shp` |
| `densff` | Densidad de fallas y fracturas | Estadísticas precomputadas en `shtv4_stats_densidadff.shp` |
| `prom` | Precipitación | Estadísticas precomputadas en `shtv4_stats_promedio.shp` |
| `slope` | Pendiente | Calculada a partir de `wc2.1_30s_elev.tif` |

## Catálogos

Los archivos de catálogo ubicados en la carpeta `Shapefiles` controlan la carga
y presentación de las variables:

- `categorical_variables_catalog.csv`: relaciona cada archivo vectorial con el
  campo categórico que se debe procesar y el prefijo de sus columnas.
- `stats_variables_catalog.csv`: identifica los archivos con estadísticas
  precomputadas, la llave de unión y el prefijo de cada variable continua.
- `variables_ordering_catalog.csv`: define el orden y las etiquetas utilizadas
  en la gráfica de composiciones por cluster.
- `usv_rename.csv`: relaciona las descripciones originales de uso del suelo y
  vegetación con las etiquetas agrupadas utilizadas en el análisis.

## Ejecución

Los comandos siguientes deben ejecutarse desde la raíz del proyecto.

### 1. Reclasificar uso del suelo y vegetación

Ejecute este paso cuando sea necesario regenerar
`usv250s7cw_renamed.shp` a partir de los datos originales y de
`usv_rename.csv`:

```bash
Rscript 01-scripts/00-rename_usv.R
```

Las descripciones sin correspondencia reciben la etiqueta
`usv_sin_clasificar` y se reportan en la consola.

### 2. Construir el conjunto de datos

```bash
Rscript 01-scripts/01-preprocessing.R
```

Este script genera:

- `02-processed_data/YYYY-MM-DD_dataset.csv`, sin geometría;
- `02-processed_data/aoi_shp/YYYY-MM-DD_dataset.gpkg`, con geometría.

### 3. Ejecutar el agrupamiento

```bash
Rscript 01-scripts/02-clustering.R
```

El script genera los diagnósticos WSS y silhouette, construye el agrupamiento
jerárquico con el valor de `k` configurado dentro del archivo y guarda el
dendrograma, las descripciones por cluster y el shapefile final. Para seleccionar
el número de grupos, revise los diagnósticos y el dendrograma, ajuste `k` en
`01-scripts/02-clustering.R` y vuelva a ejecutar el script para producir la
salida definitiva.

> **Importante:** `01-preprocessing.R` escribe el dataset con la fecha de
> ejecución y `02-clustering.R` busca un dataset con esa misma fecha. Ambos
> scripts deben ejecutarse el mismo día, a menos que se modifique
> explícitamente `current_date`.

## Metodología

### Preprocesamiento espacial

`01-preprocessing.R` carga las nanocuencas y procesa las covariables declaradas
en los catálogos. Los datos vectoriales se transforman al sistema de referencia
de las nanocuencas y se filtran al área de interés. Para cada categoría se
calcula el área de intersección con cada nanocuenca en hectáreas.

Las estadísticas precomputadas se incorporan mediante `idNano`, conservando el
número y el orden de las nanocuencas. La pendiente se calcula a partir del DEM y
se resume mediante su media y desviación estándar dentro de cada polígono.

### Matriz de distancia

Antes de calcular distancias, las columnas categóricas de cada variable se
normalizan por nanocuenca para representar composiciones relativas.

La distancia entre dos nanocuencas se calcula según el tipo de variable:

- las composiciones categóricas se comparan con la distancia de
  Jensen–Shannon;
- las variables continuas se comparan con la distancia de Hellinger a partir de
  su media y varianza, donde la varianza se obtiene como `_sd^2`.

La distancia total es la suma de las contribuciones de todas las variables. El
resultado se convierte en un objeto de distancias para el agrupamiento.

### Agrupamiento jerárquico

El agrupamiento se realiza mediante clustering jerárquico aglomerativo con
enlace completo (`complete`). El número de grupos se evalúa mediante WSS
(_within-cluster sum of squares_), silhouette y la interpretación del
dendrograma.

## Productos

`02-clustering.R` guarda los siguientes productos en
`03-clustering_output`:

- `YYYY-MM-DD_wws_plot.png`: diagnóstico WSS;
- `YYYY-MM-DD_silhouette_plot.png`: diagnóstico silhouette;
- `YYYY-MM-DD_dendrogram_<k>_k.jpg`: dendrograma con los grupos señalados;
- `YYYY-MM-DD_cluster_description_<k>_k.jpg`: gráfica de barras con la
  composición promedio de cada variable categórica por cluster;
- `YYYY-MM-DD_cluster_mean_boxplots_<k>_k.jpg`: boxplots de las medias de las variables continuas por cluster;
- `YYYY-MM-DD_nanocuecas_cluster_<k>_k.*`: shapefile de nanocuencas con la
  asignación de cluster.

En las gráficas de composición, cada nanocuenca tiene el mismo peso, se promedian las proporciones dentro de cada
cluster. Los boxplots también utilizan una observación por nanocuenca.
