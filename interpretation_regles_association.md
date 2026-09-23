# Interprétation des règles d'association — Boeny (2016‑2025)

> **Source** : `regles_association_transitions.csv` (extraction Apriori, support ≥ 0,01,
> confiance ≥ 0,5) sur les 86 988 cellules ayant au moins une transition, parmi les
> 482 359 trajectoires. **10 lignes = 5 associations bidirectionnelles** (chaque paire
> est comptée deux fois, une fois par sens d'item). Les règles de transitions sont
> **structurellement indépendantes du choix de K** : elles ne dépendent que des
> transitions elles‑mêmes, pas de la partition.

---

## Synthèse en un tableau (prêt à coller au chapitre 3)

| # | Association (anti‑symétrique) | Support | Confiance | **Lift** | Effectif | Lecture géographique |
|---|---|---|---|---|---|---|
| 1 | eau ⇄ sol nu | 1,03 % | 0,56–0,65 | **34,9** | 899 | Vasières / estuaire : alternance eau ⇄ sol nu sous l'effet de la marée et de la décrue (Betsiboka) |
| 2 | vég. non forestière ⇄ vég. inondée | 1,05 % | 0,53 | 26,7 | 914 | Frange littorale : bascule saisonnière entre végétation inondée et non inondée |
| 3 | vég. inondée ⇄ forêt | 1,16 % | 0,53 | 22,3 | 1007 | Bordure forêt / zone inondée (secteurs de transition humide‑forestier) |
| 4 | cultures ⇄ vég. non forestière | 6,18 % | 0,54–0,62 | 5,4 | 5 380 | **Front agricole** : défriche de savane et retour en végétation non forestière |
| 5 | sol nu ⇄ cultures | ~1,0 % | 0,52 | 5,2 | ~894 | Alternance parcellaire sol nu ⇄ cultures (cycle agricole annuel) |
| 6 | forêt ⇄ vég. non forestière | 26,8 % | 0,54–0,60 | **1,2** | 23 307 | Co‑occurrence **massive mais lift ≈ 1 → peu informative** (les deux classes coexistent partout à 250 m) |

> ⚠️ **Ne pas lire « forêt ⇄ vég. non forestière » comme une transition écologique** :
> son lift (1,2) est proche de l'indépendance. C'est une **co‑occurrence spatiale**
> (mosaïque forêt‑savane à l'échelle 250 m), pas une dynamique systématique.

---

## Paragraphes d'interprétation (prêts à coller, format LaTeX)

### 1. La dynamique estuarienne : eau ⇄ sol nu (lift 34,9)

L'association de loin la plus forte (lift ≈ 35, confiance 0,56–0,65, 899 cellules) associe
la classe *eau* et la classe *sol nu*. Sur le littoral Boeny, cette alternance correspond
aux **vasières et aux bancs de l'estuaire de la Betsiboka**, où le marnage et la décrue
font basculer une cellule d'une classe à l'autre d'une année sur l'autre. Le profil du
cluster 4 (*végétation inondée*, côte) et les cartes montrent que ces cellules sont
**concentrées en une frange côtière étroite** : c'est une dynamique **spatialement
localisée et intense**, représentative des zones intertidales. Cette paire illustre le cas
d'école d'une règle d'association parfaitement **symétrique** (eau → sol nu et
sol nu → eau avec les mêmes lift), signe que les deux états s'alternent sans état
dominant privilégié.

### 2. Les franges de végétation inondée (lifts 26,7 et 22,3)

Deux associations de niveau intermédiaire (lift ≈ 22–27, ~900–1 000 cellules) relient la
*végétation inondée* d'une part à la *végétation non forestière*, d'autre part à la
*forêt* : **vég. non forestière ⇄ vég. inondée** (lift 26,7) et **vég. inondée ⇄ forêt**
(lift 22,3). Elles décrivent la **bordure humide** du paysage : des zones où la présence
d'eau ou de sols inondés alterne avec de la végétation herbacée/arborée. Sur la carte de
contexte, ces franges bordent la côte et les grandes vallées inondables. Géographiquement,
elles prolongent le même système que la règle n°1 mais en domaine **végétalisé**, formant
une transition douce entre l'estuaire nu et les terres hautes.

### 3. Le front agricole : cultures ⇄ vég. non forestière (lift 5,4)

Avec un support nettement plus élevé (6,18 %, soit ~5 400 cellules) et un lift de 5,4,
l'association **cultures ⇄ vég. non forestière** est la plus **étendue** des règles
significatives. Elle matérialise le **front pionnier agricole** : des cellules qui passent
de savane / végétation non forestière à cultures, et inversement. Combiné aux clusters 3
(mosaïque dynamique) et 5/6 (savane et forêt stables), ce résultat chiffre la pression
agricole en bordure des espaces naturels : la conversion savane → cultures n'est pas
partout unidirectionnelle, une partie des parcelles « retourne » en végétation non
forestière (rotation, jachère, déprise). C'est l'une des associations **les plus
importantes à discuter** dans le mémoire, car elle relie directement dynamique
observationnelle et hypothèse socio‑économique (agriculture extensive).

### 4. Le cycle parcellaire : sol nu ⇄ cultures (lift 5,2)

L'association sol nu ⇄ cultures (lift ≈ 5,2, ~890 cellules) complète la précédente : il
s'agit du **cycle agricole intra‑annuel**, où une parcelle passe par une phase de sol
dénudé (préparation du sol, fin de cycle) puis est remis en culture. Plus localisée que
le front agricole, elle signale un **système culturaux intensif** dans les zones cultivées
stables (cluster 2). Dans le mémoire, on peut la présenter comme l'échelle « intra‑
parcellaire » du même phénomène que la règle n°3 (échelle « paysage »).

### 5. La grande co‑occurrence forêt / vég. non forestière (lift 1,2) — à interpréter avec prudence

La règle numériquement la plus volumineuse (support 26,8 %, **23 307 cellules**) est aussi,
structurellement, la **moins informative** : son lift (≈ 1,2) est très proche de 1, c'est‑
à‑dire de l'indépendance statistique entre les deux items. Autrement dit, la co‑occurrence
forêt–végétation non forestière à 250 m est surtout due au **fait que ces deux classes sont
omniprésentes** dans le paysage Boeny (mosaïque savane‑forêt), et non à une relation de
transition privilégiée. **Sur le plan méthodologique, ce résultat est précieux à citer** :
le lift permet précisément de distinguer une association *fréquente mais banale* d'une
association *fréquente ET rare* (comme eau ⇄ sol nu). C'est un argument pédagogique fort
pour justifier le tri par lift.

---

## Règles « en chaîne » (`regles_association_chaines.csv`)

Les règles en chaîne concatènent plusieurs transitions consécutives (ex. une cellule passant
par sol nu puis eau). Elles confortent les interprétations ci‑dessus en montrant que les
associations ne sont pas isolées année par année mais s'enchaînent : l'alternance
eau ⇄ sol nu se répète sur plusieurs années consécutives dans l'estuaire, et la conversion
savane → cultures est souvent suivie d'une phase de culture stable. **À utiliser comme
annexe ou comme illustration des trajectoires typiques** (complément des clusters).
