# Résumé méthodologique — Clustering K‑means K=6 + règles d'association (Boeny, 2016‑2025)

*Document d'une page : justifie, en une lecture, l'ensemble des choix de paramètres.
Conçu pour être placé en annexe méthodologique ou résumé de la démarche avant le
chapitre 3 du mémoire.*

---

## 1. Pourquoi K = 6 ? Trois critères indépendants convergent

| Critère | Réponse | Preuve chiffrée |
|---|---|---|
| **Coude d'inertie** | 5 → 6 : −20,7 % ; 6 → 7 : −0,9 % | Inertie K=6 : 4 620 780 ; K=7 : 4 579 245 |
| **Silhouette de Rousseeuw** | **Maximale à K=6 : 0,765** | 0,702 (K=5) → 0,765 (K=6) → 0,739 (K=7) |
| **Interprétabilité spatiale** | Savane / forêt stables + mosaïque dynamique | Clusters 5–6 = stables (entropie 0), diffère uniquement la classe occupée |

> **Argument clé du mémoire** : K=6 est le seul K qui maximise **simultanément** la
> silhouette **et** se situe au dernier vrai coude d'inertie (−20,7 % puis stagnation).
> C'est une convergence de deux critères indépendants — rare et donc décisive.

## 2. Pourquoi multiplier les nstart ? (reproductibilité + non‑monotonie)

- **`NSTART_COUDE = 25`** : le k‑means est sensible aux minima locaux. Le premier coude
  (nstart=1) semblait **non‑monotone** entre K=3, K=4 et K=5.
- Avec **nstart=25**, l'inertie devient strictement convergente ET on découvre que
  **K=8 donne une inertie identique à K=6 (4 620 780)** → **cluster vide (dégénéré)**.
  **C'est la cause exacte de la non‑monotonie** observée au premier coude. À citer
  comme piège méthodologique (ne pas conclure sur un coude avec nstart=1).
- **Graine fixe** (`SEED = 1`, `random_state = 1`) : tout est reproductible à l'identique.

## 3. La silhouette : pratique standard, justifiable

- Silhouette calculée sur **échantillon de 5 000 cellules** (`TAILLE_SILHOUETTE = 5000`)
  — pratique standard (coût O(n²) du calcul de silhouette sur 482 359 points serait
  prohibitif sans échantillonnage). **À justifier en une phrase** dans le mémoire
  (c'est un exercice d'honnêteté méthodologique valorisé par un jury).
- Résultat : silhouette maximale = 0,765 à K=6 → clustering « substantiel » (0,70–0,80).

## 4. Le lissage : une étape de débruitage, pas une étape de modélisation

| | Brut | Lissé (fenêtre 3×3, seuil 5 %) |
|---|---|---|
| Cellules avec bruit | ~16,4 % | **2,0 %** |
| Cellules restantes | 482 359 | 482 359 (aucune exclue) |
| Table utilisée pour le clustering | — | `trajectoires_large_250m_lisse.csv` |

- Le lissage **médian temporel + spatial** corrige les artefacts de capteur (mirage
  spectral, brume, angles de visée) sans supprimer de cellules.
- **Le clustering et les règles sont faits sur la table LISSÉE** ; le brut reste
  disponible (`FICHIER_ENTREE` sélectionnable) pour la comparaison.

## 5. Règles d'association : indépendantes du choix de K

- Les règles décrivent des **transitions annuelles** (et non des clusters) → elles ne
  dépendent PAS de K. 10 règles extraites (support ≥ 0,01, confiance ≥ 0,5), triées
  par **lift** (mesure d'intérêt indépendante de la fréquence de base, idéale pour
  comparer des associations de tailles très différentes).
- **Justification du lift plutôt que la confiance** : la confiance favorise les règles
  fréquentes mais banales (ex. forêt ⇄ vég. non for., lift ≈ 1,2 = indépendance) ;
  le lift isole les associations réellement informatives (eau ⇄ sol nu, lift ≈ 35).

## 6. Reproductibilité (à inscrire en annexe technique)

```
Graine fixe ............ SEED = 1 (clustering_final, coude, échantillons de cartes)
K testés (coude) ....... 2..8, nstart = 25
K final ................ 6
Silhouette ............. échantillon de 5 000 cellules
Lissage ................ fenêtre 3×3 + médian, seuil stabilité 5 %
Règles ................. Apriori : support ≥ 0,01, confiance ≥ 0,5, minlen = 2
Table entrée clustering  trajectoires_large_250m_lisse.csv (lissée)
Nombre de cellules ..... 482 359
```

**Résultat reproductible en une commande :**

```bash
python clustering_regles_association_boeny.py   # K=6 + 10 règles (coude sauté)
```

> **Point de vigilance pour le mémoire** : documenter que le choix de K a été fait
> sur le **coude + silhouette (nstart=25)** et que le clustering final utilise la même
> graine — la partition K=6 est ainsi **stable et défendable**.
