# Trajectoires spatio-temporelles du changement d'occupation du sol — Boeny (2016–2025)

Analyse des **trajectoires de changement** de l'occupation du sol sur la région
Boeny (nord-ouest de Madagascar, estuaire de la Betsiboka) à partir des données
**Dynamic World à 10 m (Google Earth Engine)**, agrégées à 250 m, lissées puis
regroupées en **clusters de trajectoires** et analysées par **règles
d'association**.

**Mémoire de master 2 — Analyse spatiale.** Pipeline reproductible
(scripts Python + R, graine fixe, coude + silhouette, LaTeX).

## 📊 Résultats clés

| Résultat | Valeur |
|---|---|
| Cellules analysées (250 m, lissé) | 482 359 |
| Clusters de trajectoires (k-means) | **K = 6** (silhouette 0,765) |
| Clusters stables (entropie 0) | Savane stable (269 063), Forêt stable (100 773) |
| Cluster le plus dynamique | Mosaïque forêt–savane (73 390) |
| Règle la plus forte | eau ⇄ sol nu (lift **34,9** — vasières estuaire Betsiboka) |
| Règles extraites (Apriori) | 10 (support ≥ 0,01, confiance ≥ 0,5) |

## 🗂️ Structure

```
script_thesis/
├── clustering_regles_association_boeny.py   # clustering K=6 + règles Apriori
├── agregation_trajectoires_boeny.py         # 10 m → 250 m + masque polygone
├── lissage_trajectoires_boeny.py            # lissage (bruit 16,4 % → 2,0 %)
├── visualisation_spatiale_boeny.py          # cartes statique/interactive/contexte
├── dynamicworldBoeny.py                     # exports GEE 2016–2025
├── R/…                                      # équivalents R
├── donnees/boeny_250m/*.csv                 # profils, matrice, règles (résultats)
├── carte_clusters_trajectoires.png          # figure mémoire (chapitre 3)
├── carte_contexte_boeny.png                 # contexte géographique
├── *.html                                   # cartes/figures interactives
├── resume_methodologique_boeny.md           # annexe méthodo K=6 (1 page)
├── interpretation_regles_association.md     # interprétation des 10 règles
└── chapitre3_resultats_clustering.tex       # chapitre 3 LaTeX (compilable)
```

## 🚀 Reproduire

```bash
python agregation_trajectoires_boeny.py              # donnees brutes requises (GEE)
python lissage_trajectoires_boeny.py                 # trajectoires_large_250m_lisse.csv
python clustering_regles_association_boeny.py        # -> trajectoires_avec_clusters.csv (K=6)
python visualisation_spatiale_boeny.py               # cartes + html
pdflatex chapitre3_resultats_clustering.tex          # chapitre 3
```

> **Note** : les données brutes (Dynamic World 10 m ~160 Mo, GADM) ne sont **pas
> incluses** dans ce dépôt — à télécharger via `dynamicworldBoeny.py` (GEE) selon
> le `GUIDE_EXECUTION.md`.

## 📄 Licence

CC‑BY‑NC‑4.0 — attribution + usage non commercial (mémoire universitaire).

## 🔬 Reproductibilité

- Graine fixe `SEED = 1` (`random_state = 1`, `nstart = 25`)
- Coude : `methode_du_coude` + silhouette de Rousseeuw (échantillon 5 000)
- Lissage : fenêtre 3×3 + médian temporel, seuil de stabilité 5 %
- Fichiers : voir `donnees/boeny_250m/diagnostic_bruit_*.csv`
