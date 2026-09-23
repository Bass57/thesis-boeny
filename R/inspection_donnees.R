## ============================================================================
## Inspection de la qualité des trajectoires - Boeny (250 m, 2016-2025)
## ============================================================================
## Entrée  : trajectoires_large_250m.csv (produit par agregation_trajectoires_boeny.R)
## Sortie  : - évolution du pourcentage de chaque classe par année (CSV + PNG)
##           - diagnostic de bruit par cellule (nb de changements, "aller-retour")
##           - liste des cellules à traiter avec prudence avant la fouille de données
## ============================================================================

dossier <- "donnees/boeny_250m"
table_large <- read.csv(file.path(dossier, "trajectoires_large_250m.csv"))

annees <- 2016:2025
colonnes_annees <- as.character(annees)

noms_classes <- c(
  "0" = "eau", "1" = "foret", "2" = "vegetation_non_forestiere",
  "3" = "vegetation_inondee", "4" = "cultures", "5" = "bati", "6" = "sol_nu"
)

## ---------------------------------------------------------------------------
## 1. EVOLUTION DES CLASSES DANS LE TEMPS
## ---------------------------------------------------------------------------
## Objectif : vérifier qu'aucune classe n'est quasi-vide (trop peu de cellules
## pour être exploitable en fouille de données) et que les tendances observées
## sont plausibles (pas de saut brutal isolé qui trahirait une erreur).

proportions <- sapply(colonnes_annees, function(an) {
  codes <- table_large[[an]]
  tab <- table(factor(codes, levels = names(noms_classes)))
  round(100 * tab / sum(tab), 2)
})
rownames(proportions) <- noms_classes[rownames(proportions)]
colnames(proportions) <- annees

cat("Pourcentage de chaque classe par année :\n")
print(proportions)

write.csv(proportions, file.path(dossier, "evolution_classes_pourcentage.csv"))

png(file.path(dossier, "evolution_classes.png"), width = 900, height = 600)
matplot(annees, t(proportions), type = "l", lty = 1, lwd = 2,
        col = rainbow(nrow(proportions)),
        xlab = "Année", ylab = "% de la superficie",
        main = "Évolution de l'occupation du sol - Boeny (250 m)")
legend("topright", legend = rownames(proportions),
       col = rainbow(nrow(proportions)), lty = 1, lwd = 2, cex = 0.8)
dev.off()

classes_quasi_vides <- rownames(proportions)[apply(proportions, 1, max) < 1]
if (length(classes_quasi_vides) > 0) {
  cat("\nATTENTION - classes représentant moins de 1% sur toute la période :\n")
  print(classes_quasi_vides)
  cat("A envisager : les fusionner avec une classe voisine avant la fouille de données.\n")
}

## ---------------------------------------------------------------------------
## 2. VOLATILITE PAR CELLULE (nombre de changements sur la série)
## ---------------------------------------------------------------------------
## Une cellule qui change de classe presque chaque année est suspecte : les
## conversions réelles d'occupation du sol (forêt -> culture, culture -> bâti,
## etc.) sont rarement aussi fréquentes sur des composites annuels.

matrice_codes <- as.matrix(table_large[, colonnes_annees])

nb_changements <- apply(matrice_codes, 1, function(x) sum(diff(x) != 0))

cat("\nDistribution du nombre de changements par cellule (sur 9 transitions) :\n")
print(summary(nb_changements))
print(table(nb_changements))

## ---------------------------------------------------------------------------
## 3. DETECTION DES "ALLER-RETOUR" (bruit typique : A -> B -> A en 1 an)
## ---------------------------------------------------------------------------
## Un aller-retour sur une seule année (ex : forêt en 2019, sol nu en 2020,
## forêt en 2021) est presque toujours un artefact (nuage résiduel, confusion
## du classifieur) plutôt qu'une vraie dynamique territoriale.

detecte_flipflop <- function(x) {
  n <- length(x)
  if (n < 3) return(0)
  sum(x[1:(n - 2)] == x[3:n] & x[1:(n - 2)] != x[2:(n - 1)])
}

nb_flipflops <- apply(matrice_codes, 1, detecte_flipflop)

diagnostic <- data.frame(
  cell_id       = table_large$cell_id,
  x             = table_large$x,
  y             = table_large$y,
  nb_changements = nb_changements,
  nb_flipflops   = nb_flipflops,
  bruit_suspect  = nb_flipflops > 0 | nb_changements >= 5
)

write.csv(diagnostic, file.path(dossier, "diagnostic_bruit_cellules.csv"), row.names = FALSE)

pct_bruit <- round(100 * mean(diagnostic$bruit_suspect), 1)
cat("\n", pct_bruit, "% des cellules sont signalées comme bruit potentiel ",
    "(", sum(diagnostic$bruit_suspect), " cellules sur ", nrow(diagnostic), ").\n", sep = "")
cat("Fichier détaillé : diagnostic_bruit_cellules.csv\n")
