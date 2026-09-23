# Guide d'exécution -- Mémoire Boeny (fouille de données spatio-temporelles)

Ce guide couvre toute la chaîne technique, de l'installation à l'export des
résultats. Suivez les étapes dans l'ordre : chacune dépend des sorties de la
précédente.

---

## Étape 1 -- Préparer l'environnement

### Python
```bash
pip install earthengine-api geopandas rasterio matplotlib cartopy folium plotly
earthengine authenticate
```
`earthengine authenticate` ouvre une page web pour lier votre compte Google.
Vous devez aussi disposer d'un **projet Google Cloud** associé à Earth
Engine (obligatoire depuis 2022) -- son ID va dans `EE_PROJECT` du script
d'extraction.

### R
```r
install.packages(c("terra", "arules", "geodata"))
```

### Contour de la région Boeny
Le nom "Boeny" (région créée en 2004) n'existe pas dans les jeux de limites
administratives intégrés à Earth Engine. Le package R `geodata` permet de le
récupérer proprement :

```r
library(geodata)
library(terra)

mdg1 <- gadm(country = "MDG", level = 1, path = "donnees/gadm")

# Vérifiez d'abord la liste exacte des noms disponibles : selon la version
# de GADM, la région peut être orthographiée différemment ou absente à ce
# niveau (auquel cas il faudra chercher au niveau 2, ou une source officielle
# malgache comme le BNGRC/INSTAT).
print(unique(mdg1$NAME_1))

boeny <- mdg1[mdg1$NAME_1 == "Boeny", ]   # ajustez le nom selon ce qui s'affiche ci-dessus
writeVector(boeny, "donnees/boeny_boundary.shp", overwrite = TRUE)
```

Puis importez ce contour comme asset Earth Engine :
```bash
earthengine upload table --asset_id=users/VOTRE_NOM_UTILISATEUR/boeny_boundary donnees/boeny_boundary.shp
```
Suivez l'avancement avec `earthengine task list`. Une fois terminé, notez
l'`asset_id` complet -- c'est lui qui va dans `AOI_ASSET_ID`.

---

## Étape 2 -- Extraire les composites annuels (Python / Earth Engine)

Dans `dynamicworldBoeny.py`, renseignez :
```python
EE_PROJECT = "projet-sigd-01"
AOI_ASSET_ID = "projects/projet-sigd-01/assets/boeny-boundary"
```
(Ces valeurs sont déjà renseignées sur ce poste.) Puis exécutez :
```bash
cd /home/bass57/Downloads/script_thesis
./venv/bin/python dynamicworldBoeny.py
```
Cela lance 10 tâches d'export (2016-2025) vers votre Google Drive, dans le
dossier `boeny_dynamicworld`. Suivez leur avancement dans l'onglet **Tasks**
du [Code Editor Earth Engine](https://code.earthengine.google.com/), ou :
```python
import ee
ee.Initialize(project="votre-projet-gcp")
print(ee.data.getTaskList())
```
Ces exports peuvent prendre de quelques minutes à plusieurs heures selon la
charge des serveurs Earth Engine.

---

## Étape 3 -- Télécharger et agréger à 250 m (Python)

1. Téléchargez les 10 GeoTIFF depuis Google Drive vers
   `donnees/boeny_dynamicworld_10m/` (chemin attendu par le script).
2. Exécutez :
```bash
./venv/bin/python agregation_trajectoires_boeny.py
```
Sorties dans `donnees/boeny_250m/` : le raster multi-bandes 250 m, et les
tables `trajectoires_large_250m.csv` / `trajectoires_longue_250m.csv`.

> ⚠️ IMPORTANT : GEE exporte la boîte englobante de Boeny et écrit les pixels
> masqués (mer, hors-polygone) en `0` = "eau". Le script exclut donc
> automatiquement les cellules dont le centre tombe hors du polygone GADM
> (`donnees/boeny_boundary/boeny_boundary.shp`) : comptez ~480 000 cellules
> (et non ~957 000) après la passe correcte.

> Les versions R (`agregation_trajectoires_boeny.R`, `inspection_donnees.R`,
> `clustering_regles_association_boeny.R`) restent disponibles comme
> référence, mais les paquets R requis (`terra`, `arules`) ne sont pas
> installés sur ce poste : utilisez les versions Python.

---

## Étape 4 -- Diagnostiquer la qualité des trajectoires (Python)

```bash
./venv/bin/python inspection_donnees.py
```
Regardez le pourcentage de cellules signalées comme bruit
(`diagnostic_bruit_cellules.csv`). Repère :
- **< 5 %** : négligeable, continuez directement à l'étape 5.
- **5-15 %** : excluez ces cellules avant le clustering (filtrez sur
  `bruit_suspect == FALSE`).
- **> 15 %** : il faut lisser les séries (étape 4 bis ci-dessous).

---

## Étape 4 bis -- Lissage temporel glissant 3 ans (si > 15 % de bruit)

```bash
./venv/bin/python lissage_trajectoires_boeny.py
```
Remplace chaque année intérieure par la classe majoritaire sur 3 ans
(t-1, t, t+1), uniquement en cas de consensus (≥ 2 années sur 3) ; les bords
2016 et 2025 sont inchangés. Sorties dans `donnees/boeny_250m/` :
`trajectoires_large_250m_lisse.csv` (nouvelle table d'entrée de l'étape 5) et
`diagnostic_bruit_apres_lissage.csv` (comparaison par cellule).

---

## Étape 5 -- Clustering des trajectoires + règles d'association (Python)

```bash
./venv/bin/python clustering_regles_association_boeny.py
```
Ce script :
1. calcule des indicateurs de trajectoire par cellule (entropie, nombre de
   changements, classes de début/fin/dominante) ;
2. affiche un graphique du coude (`coude_choix_k.png`) -- **ouvrez-le et
   ajustez `K_FINAL` dans le script** avant de considérer le clustering
   définitif (5 est une valeur de départ, pas un choix validé) ;
3. exporte `trajectoires_avec_clusters.csv` (nécessaire à l'étape 6) et
   `profil_clusters.csv` (pour interpréter chaque cluster) ;
4. extrait les règles d'association entre transitions et exporte
   `regles_association_transitions.csv` et `regles_association_chaines.csv`
   (les enchaînements de conversion en plusieurs étapes).

Si le nombre de règles obtenu est très élevé (peu informatif) ou très faible
(pas assez de signal), ajustez `supp` et `conf` dans l'appel à `apriori()`.

---

## Étape 6 -- Visualiser les résultats (Python)

Dans `visualisation_spatiale_boeny.py`, vérifiez que le chemin du CSV pointe
vers `trajectoires_avec_clusters.csv` produit à l'étape précédente, puis :
```bash
./venv/bin/python visualisation_spatiale_boeny.py
```
Sorties : cartes statiques (pour le document), carte interactive HTML, et
graphique d'évolution interactif.

---

## Étape 7 -- Rédiger le chapitre 3 et assembler le mémoire

Une fois les étapes 5 et 6 terminées, vous aurez tout le matériel nécessaire
(cartes, profils de clusters, règles d'association) pour rédiger le chapitre
Résultats et discussion. On pourra alors assembler introduction, état de
l'art, méthodologie et résultats en un seul document.
