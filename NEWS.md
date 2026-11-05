# alphaearthnemeton 2026.11.1 (2026-11-05)

Première version (spec 064 de `nemeton`, chantier A).

- `aef_fetch()` : embeddings AlphaEarth sur une emprise, par année, lus par
  fenêtre dans les COG publics, dé-quantifiés, reprojetés au plus proche
  voisin sur une grille de 10 m, mosaïqués, mis en cache (COG adressé par
  contenu).
- `aef_index()` : tuiles couvrant une emprise (index GeoPackage distant ou
  index en mémoire).
- `aef_dequantize()`, `aef_info()` (licence CC-BY 4.0 et attribution),
  `aef_cache_dir()`.
