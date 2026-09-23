## ============================================================================
## Clustering des trajectoires + règles d'association - Boeny
## ============================================================================
## Entrées :
##   - trajectoires_large_250m.csv  (une ligne par cellule, une colonne par année)
##   - trajectoires_longue_250m.csv (une ligne par cellule x année)
## Sorties :
##   - trajectoires_avec_clusters.csv  (attendu par visualisation_spatiale_boeny.py)
##   - regles_association_transitions.csv
##
## IMPORTANT - choix méthodologique :
## Une distance d'appariement optimal (Optimal Matching, ex. package TraMineR),
## comme utilisée par Mas et al. (2019) qu'on cite au chapitre 1, suppose de
## calculer une matrice de dissimilarité entre TOUTES les paires de cellules.
## Sur ~480 000 cellules (grille 250 m de Boeny), cela représente environ
## 1,15.10^11 paires : totalement infaisable en temps et en mémoire sur un
## poste de travail standard.
## On utilise donc une alternative qui reste fidèle à l'esprit de Mas et al.
## (2019) : extraire, PAR CELLULE, des indicateurs résumant sa trajectoire
## (proches de leurs mesures d'entropie et de "turbulence"), puis appliquer un
## K-Means sur ces indicateurs. Le coût de calcul devient linéaire en nombre
## de cellules, donc gérable à cette échelle.
## ============================================================================

library(arules)

dossier <- "donnees/boeny_250m"
annees <- 2016:2025
colonnes_annees <- as.character(annees)

# ==============================================================================
# PARTIE A : CLUSTERING DES TRAJECTOIRES (K-Means sur indicateurs par cellule)
# ==============================================================================

table_large <- read.csv(file.path(dossier, "trajectoires_large_250m.csv"))
matrice_codes <- as.matrix(table_large[, colonnes_annees])

## ---------------------------------------------------------------------------
## A.1 - Extraction des indicateurs de trajectoire, par cellule
## ---------------------------------------------------------------------------

entropie_shannon <- function(x) {
  p <- table(x) / length(x)
  p <- p[p > 0]
  -sum(p * log2(p))
}

classe_dominante <- function(x) {
  tab <- table(x)
  names(tab)[which.max(tab)]
}

indicateurs <- data.frame(
  cell_id          = table_large$cell_id,
  entropie         = apply(matrice_codes, 1, entropie_shannon),
  nb_changements   = apply(matrice_codes, 1, function(x) sum(diff(x) != 0)),
  nb_classes_visitees = apply(matrice_codes, 1, function(x) length(unique(x))),
  classe_debut     = factor(matrice_codes[, 1]),
  classe_fin       = factor(matrice_codes[, ncol(matrice_codes)]),
  classe_dominante = factor(apply(matrice_codes, 1, classe_dominante))
)

## ---------------------------------------------------------------------------
## A.2 - Mise en forme : indicateurs numériques + indicatrices (one-hot) pour
## les variables catégorielles, puis standardisation
## ---------------------------------------------------------------------------

caracteristiques <- model.matrix(
  ~ entropie + nb_changements + nb_classes_visitees +
    classe_debut + classe_fin + classe_dominante - 1,
  data = indicateurs
)
caracteristiques_std <- scale(caracteristiques)

## ---------------------------------------------------------------------------
## A.3 - Choix du nombre de clusters (méthode du coude)
## ---------------------------------------------------------------------------
## Le K-Means complet est peu coûteux ici (indicateurs agrégés, pas de matrice
## de distance entre cellules) : on peut le lancer directement sur l'ensemble
## des ~480 000 cellules pour plusieurs valeurs de K.

set.seed(1)
k_teste <- 2:8
inertie <- sapply(k_teste, function(k) {
  kmeans(caracteristiques_std, centers = k, nstart = 5, iter.max = 50)$tot.withinss
})

png(file.path(dossier, "coude_choix_k.png"), width = 700, height = 500)
plot(k_teste, inertie, type = "b", pch = 19,
     xlab = "Nombre de clusters (K)", ylab = "Inertie intra-cluster",
     main = "Choix du nombre de clusters (méthode du coude)")
dev.off()
cat("Graphique du coude exporté : coude_choix_k.png\n")
cat("Inspectez-le pour choisir K (le 'coude' de la courbe), puis ajustez K_FINAL ci-dessous.\n")

## Valeur par défaut à ajuster après lecture du graphique du coude :
K_FINAL <- 5

## ---------------------------------------------------------------------------
## A.4 - Clustering final et export
## ---------------------------------------------------------------------------

kmeans_final <- kmeans(caracteristiques_std, centers = K_FINAL, nstart = 10, iter.max = 100)
table_large$cluster <- kmeans_final$cluster

cat("\nRépartition des cellules par cluster :\n")
print(table(table_large$cluster))

write.csv(table_large, file.path(dossier, "trajectoires_avec_clusters.csv"), row.names = FALSE)
cat("Table exportée : trajectoires_avec_clusters.csv\n")

## Profil moyen de chaque cluster, pour aider à l'interprétation (ex. cluster
## à forte entropie + beaucoup de changements = zone instable/en conversion
## active ; cluster à entropie nulle = zone stable sur toute la période)
profil_clusters <- aggregate(
  indicateurs[, c("entropie", "nb_changements", "nb_classes_visitees")],
  by = list(cluster = table_large$cluster),
  FUN = mean
)
cat("\nProfil moyen des clusters (pour interprétation) :\n")
print(profil_clusters)
write.csv(profil_clusters, file.path(dossier, "profil_clusters.csv"), row.names = FALSE)

# ==============================================================================
# PARTIE B : RÈGLES D'ASSOCIATION ENTRE TRANSITIONS (Apriori)
# ==============================================================================
## Principe (dans l'esprit de Mennis & Liu, 2003, cité au chapitre 1) :
## chaque cellule est une "transaction" contenant l'ensemble des transitions
## de classe (ex. "foret->sol_nu") qu'elle a connues sur 2016-2025. Apriori
## cherche alors quelles transitions co-apparaissent fréquemment dans les
## mêmes cellules - révélant des enchaînements de conversion.

noms_classes <- c(
  "0" = "eau", "1" = "foret", "2" = "vegetation_non_forestiere",
  "3" = "vegetation_inondee", "4" = "cultures", "5" = "bati", "6" = "sol_nu"
)

## B.1 - Construction des transactions (une liste de transitions par cellule)

construire_transitions_cellule <- function(sequence_codes) {
  transitions <- character(0)
  for (t in seq_len(length(sequence_codes) - 1)) {
    de <- sequence_codes[t]
    vers <- sequence_codes[t + 1]
    if (de != vers) {
      transitions <- c(transitions, paste0(noms_classes[as.character(de)], "->",
                                            noms_classes[as.character(vers)]))
    }
  }
  unique(transitions)
}

liste_transactions <- apply(matrice_codes, 1, construire_transitions_cellule)

# Retirer les cellules sans aucune transition (rien à apprendre pour Apriori)
liste_transactions <- liste_transactions[lengths(liste_transactions) > 0]
cat("\n", length(liste_transactions), "cellules avec au moins une transition, sur",
    nrow(matrice_codes), "au total.\n", sep = "")

transactions <- as(liste_transactions, "transactions")

## B.2 - Extraction des règles (seuils de départ, à ajuster selon le volume
## de règles obtenu - remonter les seuils si trop de règles peu informatives)

regles <- apriori(
  transactions,
  parameter = list(supp = 0.01, conf = 0.5, minlen = 2)
)

cat("\n", length(regles), "règles extraites (support >= 0.01, confiance >= 0.5).\n", sep = "")

regles_triees <- sort(regles, by = "lift")
cat("\nTop 10 règles par lift :\n")
inspect(head(regles_triees, 10))

## B.3 - Filtrage des règles "en chaîne" : celles où la classe d'arrivée de
## l'antécédent correspond à la classe de départ du conséquent (ex.
## {foret->sol_nu} => {sol_nu->cultures}), qui révèlent une conversion en
## plusieurs étapes plutôt qu'une simple co-occurrence.

df_regles <- as(regles, "data.frame")
df_regles$lhs <- as.character(labels(lhs(regles)))
df_regles$rhs <- as.character(labels(rhs(regles)))

extraire_vers <- function(item) sub(".*->", "", gsub("[{}]", "", item))
extraire_de   <- function(item) sub("->.*", "", gsub("[{}]", "", item))

est_une_chaine <- mapply(function(l, r) {
  vers_lhs <- extraire_vers(l)
  de_rhs   <- extraire_de(r)
  vers_lhs == de_rhs
}, df_regles$lhs, df_regles$rhs)

regles_chaine <- df_regles[est_une_chaine, ]
cat("\n", nrow(regles_chaine), "règles \"en chaîne\" identifiées (sur", nrow(df_regles), ").\n", sep = "")

write.csv(df_regles, file.path(dossier, "regles_association_transitions.csv"), row.names = FALSE)
write.csv(regles_chaine, file.path(dossier, "regles_association_chaines.csv"), row.names = FALSE)
cat("Fichiers exportés : regles_association_transitions.csv, regles_association_chaines.csv\n")
