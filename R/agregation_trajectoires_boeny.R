## ============================================================================
## Construction de la table de trajectoires multi-temporelles - Boeny
## ============================================================================
## Entrée  : les 10 GeoTIFF annuels (2016-2025) exportés depuis Google Earth
##           Engine par extraction_dynamicworld_boeny.py (résolution 10 m,
##           nomenclature à 7 classes déjà fusionnée)
## Sortie  : - un raster multi-bandes agrégé à 250 m (une bande par année)
##           - une table "large" (1 ligne = 1 cellule de 250 m, 1 colonne par
##             année)
##           - une table "longue" (1 ligne = 1 cellule x 1 année), prête pour
##             la construction des règles d'association et du clustering de
##             trajectoires à l'étape suivante
## ============================================================================

library(terra)

## ---------------------------------------------------------------------------
## 0. CONFIGURATION - à adapter
## ---------------------------------------------------------------------------

dossier_entree   <- "donnees/boeny_dynamicworld_10m"   # dossier des GeoTIFF 10 m
dossier_sortie   <- "donnees/boeny_250m"
annee_debut      <- 2016
annee_fin        <- 2025
facteur_agrega   <- 25   # 10 m x 25 = 250 m

noms_classes <- c(
  "0" = "eau",
  "1" = "foret",
  "2" = "vegetation_non_forestiere",
  "3" = "vegetation_inondee",
  "4" = "cultures",
  "5" = "bati",
  "6" = "sol_nu"
)

dir.create(dossier_sortie, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 1. AGREGATION 10 m -> 250 m (résampling majoritaire, pas sous-échantillonnage)
## ---------------------------------------------------------------------------

annees <- annee_debut:annee_fin

agregation_modale <- function(annee) {
  chemin <- file.path(dossier_entree, sprintf("boeny_occupation_sol_%d.tif", annee))
  if (!file.exists(chemin)) {
    stop(sprintf("Fichier manquant pour %d : %s", annee, chemin))
  }
  r <- rast(chemin)

  # fun = "modal" : conserve la classe majoritaire dans chaque bloc de 25x25
  # pixels (10 m) -> cellule de 250 m. C'est l'équivalent d'un vote majoritaire,
  # pas une simple décimation qui perdrait de l'information.
  r_250 <- aggregate(r, fact = facteur_agrega, fun = "modal", na.rm = TRUE)
  names(r_250) <- as.character(annee)
  r_250
}

cat("Agrégation des", length(annees), "composites annuels (10 m -> 250 m)...\n")
liste_rasters <- lapply(annees, agregation_modale)

# Vérification : toutes les couches doivent partager la même grille avant
# l'empilement (elles le devraient, puisqu'elles proviennent de la même AOI)
pile <- rast(liste_rasters)
names(pile) <- as.character(annees)

writeRaster(
  pile,
  file.path(dossier_sortie, "boeny_occupation_sol_250m_2016_2025.tif"),
  overwrite = TRUE
)
cat("Raster multi-bandes 250 m exporté.\n")

## ---------------------------------------------------------------------------
## 2. TABLE "LARGE" : 1 ligne = 1 cellule de 250 m, 1 colonne par année
## ---------------------------------------------------------------------------

# xy = TRUE conserve les coordonnées du centre de chaque cellule, utiles pour
# une cartographie ultérieure des motifs de trajectoire
table_large <- as.data.frame(pile, xy = TRUE, na.rm = TRUE)
table_large$cell_id <- seq_len(nrow(table_large))

write.csv(
  table_large,
  file.path(dossier_sortie, "trajectoires_large_250m.csv"),
  row.names = FALSE
)
cat("Table large exportée :", nrow(table_large), "cellules.\n")

## ---------------------------------------------------------------------------
## 3. TABLE "LONGUE" : 1 ligne = 1 cellule x 1 année (format pour la fouille
##    de données : séquences par cellule, règles d'association entre années)
## ---------------------------------------------------------------------------

table_longue <- do.call(rbind, lapply(as.character(annees), function(an) {
  data.frame(
    cell_id = table_large$cell_id,
    x       = table_large$x,
    y       = table_large$y,
    annee   = as.integer(an),
    classe_code = table_large[[an]],
    classe  = noms_classes[as.character(table_large[[an]])]
  )
}))

table_longue <- table_longue[order(table_longue$cell_id, table_longue$annee), ]

write.csv(
  table_longue,
  file.path(dossier_sortie, "trajectoires_longue_250m.csv"),
  row.names = FALSE
)
cat("Table longue exportée :", nrow(table_longue), "lignes.\n")

## ---------------------------------------------------------------------------
## 4. APERCU : matrice de transition entre deux années (exemple)
## ---------------------------------------------------------------------------

matrice_transition <- table(
  debut = noms_classes[as.character(table_large[[as.character(annee_debut)]])],
  fin   = noms_classes[as.character(table_large[[as.character(annee_fin)]])]
)

cat("\nMatrice de transition", annee_debut, "->", annee_fin, ":\n")
print(matrice_transition)

write.csv(
  as.data.frame.matrix(matrice_transition),
  file.path(dossier_sortie, sprintf("matrice_transition_%d_%d.csv", annee_debut, annee_fin))
)
