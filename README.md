# alphaearthnemeton

Embeddings satellitaires annuels **AlphaEarth Foundations** (Google DeepMind,
64 dimensions, 10 m, 2017-2025) pour le paquet [`nemeton`](https://github.com/pobsteta/nemeton) :
indicateur **N5 — distance à la référence** (spec 064 de `nemeton`).

R pur : ni Python, ni compte Earth Engine. Les COG publics du bucket
`gs://alphaearth_foundations` sont lus par fenêtre (requêtes HTTP partielles,
GDAL `/vsicurl/`), après une requête spatiale sur l'index GeoPackage publié.

```r
# remotes::install_github("pobsteta/alphaearthnemeton")
library(alphaearthnemeton)

emb <- aef_fetch(units, years = 2023:2024)   # liste nommée par année
res <- nemeton::indicateur_n5_distance_reference(units, emb,
                                                 reference = reserves)
```

## Ce que fait `aef_fetch()`

1. interroge l'index `aef_index.gpkg` (une ligne par COG) sur l'emprise ;
2. lit la seule fenêtre utile de chaque COG (int8, NoData -128) ;
3. dé-quantifie : `sign(x) * (x / 127.5)^2` (chaque pixel redevient un vecteur
   de norme 1) ;
4. reprojette en Lambert-93 (EPSG:2154 ; 3035 hors France) sur une grille de
   10 m calée sur les multiples de 10 m, au **plus proche voisin** (une
   interpolation casserait la norme des vecteurs) ;
5. mosaïque les tuiles et écrit un COG float32 dans un cache adressé par
   contenu (version, année, CRS, emprise, tuiles sources).

Mesure du 2026-11-05 : 2 × 2 km près de Dabo, années 2023 et 2024, 52 s ;
même appel ensuite : 4 s (cache).

## Accès

Google déclare le bucket en « requester pays » ; les requêtes anonymes étaient
servies le 2026-11-05. Si cela change :
`options(alphaearthnemeton.base_url = "https://data.source.coop/tge-labs/aef/v1/annual/")`
(miroir Source Cooperative, même arborescence) ou identifiants GCS de GDAL.

## Licences

- Paquet : GPL (>= 3).
- Données : **CC-BY 4.0**. Attribution obligatoire : *« The AlphaEarth
  Foundations Satellite Embedding dataset is produced by Google and Google
  DeepMind. »* (`aef_info()$attribution`).

## Tests

`testthat` sur des tuiles simulées (aucun réseau). Le test sur le bucket
public est optionnel : `AEF_TEST_NETWORK=true`.
