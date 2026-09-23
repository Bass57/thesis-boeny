"""
Agrégation 10 m -> 250 m et construction des tables de trajectoires - Boeny
===========================================================================
Équivalent Python de `agregation_trajectoires_boeny.R` (mêmes entrées/sorties,
même nomenclature à 7 classes).

Entrée  : les 10 GeoTIFF annuels (2016-2025) exportés par `dynamicworldBoeny.py`
          dans `donnees/boeny_dynamicworld_10m/boeny_occupation_sol_<annee>.tif`
Sorties (dans donnees/boeny_250m/) :
  - boeny_occupation_sol_250m_2016_2025.tif  (raster multi-bandes à 250 m)
  - trajectoires_large_250m.csv              (1 ligne = 1 cellule de 250 m)
  - trajectoires_longue_250m.csv             (1 ligne = 1 cellule x 1 année)
  - matrice_transition_<deb>_<fin>.csv       (exemple 2016 -> 2025)

Dépendance : rasterio (pip install rasterio)

Principe de l'agrégation : chaque bloc de 25x25 pixels de 10 m devient une
cellule de 250 m par VOTE MAJORITAIRE (classe la plus fréquente du bloc),
pas par simple décimation -- on conserve l'information intra-bloc.

MASQUE POLYGONE (important) : le raster exporté par Earth Engine couvre la
boîte englobante de Boeny, et GEE écrit les pixels masqués (hors polygone,
mer...) en 0 = "eau" car la bande uint8 n'a pas de NoData. Sans précaution,
~la moitié de la grille se retrouve en "eau" fantôme. On retranche donc les
cellules dont le CENTRE est en dehors du polygone GADM (CHEMIN_LIMITE).
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_origin

# ============================ CONFIGURATION ============================
dossier_entree = "donnees/boeny_dynamicworld_10m"  # GeoTIFF 10 m téléchargés de Drive
dossier_sortie = "donnees/boeny_250m"
annee_debut, annee_fin = 2016, 2025
facteur_agrega = 25            # 10 m x 25 = 250 m
VALEUR_INVALIDE = 255          # classe "hors emprise / non renseignée" (exclue des tables)

# Polygone GADM de Boeny : les cellules dont le centre est en DEHORS de ce
# polygone sont exclues (voir l'explication en tête de fichier).
CHEMIN_LIMITE = "donnees/boeny_boundary/boeny_boundary.shp"

# Nomenclature des 7 classes (les GeoTIFF exportés sont déjà reclassés 0-6)
NOMS_CLASSES = {
    0: "eau", 1: "foret", 2: "vegetation_non_forestiere",
    3: "vegetation_inondee", 4: "cultures", 5: "bati", 6: "sol_nu",
}
# Variante à clés textuelles : pratique pour les opérations pandas (map, crosstab)
NOMS_CLASSES_STR = {str(k): v for k, v in NOMS_CLASSES.items()}
PAS_BLOCS = 100_000  # nombre de blocs traités par morceau (limite la mémoire)
# =======================================================================


def agregation_modale(chemin: str, facteur: int):
    """
    Lit un GeoTIFF 10 m et renvoie la grille 250 m (classe majoritaire par bloc)
    + les métadonnées géoréférencées nécessaires.

    Règle du jeu (équivalent de `terra::aggregate(fun="modal", na.rm=TRUE)`) :
      - on ne compte que les pixels dont la valeur est dans 0..6 ;
      - tout autre valeur (hors emprise, no-data, neige...) est ignorée ;
      - si un bloc ne contient AUCUN pixel valide -> VALEUR_INVALIDE.
    """
    with rasterio.open(chemin) as src:
        bande = src.read(1)
        crs = src.crs
        res_x = src.res[0]  # résolution horizontale du raster d'entrée (10 m ici)
        # Coin supérieur gauche (origine) du raster d'entrée
        x0, y0 = src.transform * (0, 0)

    ny, nx = bande.shape
    n_lignes = int(np.ceil(ny / facteur))     # nb de cellules 250 m en hauteur
    n_cols = int(np.ceil(nx / facteur))       # ... en largeur

    # On ajoute un liseré de VALEUR_INVALIDE pour que tout soit divisible par
    # `facteur` (les blocs du bord droit/bas sont donc comptés partiellement).
    hauteur, largeur = n_lignes * facteur, n_cols * facteur
    pad = np.full((hauteur, largeur), VALEUR_INVALIDE, dtype=bande.dtype)
    pad[:ny, :nx] = bande

    # (n_lignes, facteur, n_cols, facteur) -> une ligne par bloc de 25x25
    blocs = (
        pad.reshape(n_lignes, facteur, n_cols, facteur)
        .transpose(0, 2, 1, 3)
        .reshape(n_lignes * n_cols, facteur * facteur)
    )

    grille = np.full(n_lignes * n_cols, VALEUR_INVALIDE, dtype=np.uint8)
    classes = np.arange(7, dtype=np.uint8)
    for debut in range(0, blocs.shape[0], PAS_BLOCS):
        b = blocs[debut : debut + PAS_BLOCS]
        # Pour chaque bloc : nb d'occurrences de chaque classe 0..6
        totaux = np.stack([(b == c).sum(axis=1) for c in classes], axis=1)
        effectif = totaux.sum(axis=1)                 # pixels valides du bloc
        classe_majoritaire = totaux.argmax(axis=1)    # en cas d'égalité : classe la plus petite
        grille[debut : debut + PAS_BLOCS] = np.where(
            effectif > 0, classe_majoritaire, VALEUR_INVALIDE
        )

    grille = grille.reshape(n_lignes, n_cols)
    # Transformée de sortie : même origine, résolution multipliée par `facteur`
    out_res = res_x * facteur
    transformee = from_origin(x0, y0, out_res, out_res)
    return grille.astype(np.uint8), transformee, crs


def coordonnees_centres(transformee, n_lignes, n_cols):
    """Coordonnées UTM (x, y) du CENTRE de chaque cellule, déduites de la
    transformée affine (résolution lue depuis les coefficients a et e)."""
    x0, y0 = transformee * (0, 0)
    a = transformee.a  # pas en x (positif, north-up)
    e = transformee.e  # pas en y (négatif si north-up)
    xs = x0 + (np.arange(n_cols) + 0.5) * a
    ys = y0 + (np.arange(n_lignes) + 0.5) * e
    grille_x, grille_y = np.meshgrid(xs, ys)
    return grille_x, grille_y


def masque_boeny(crs, transformee, n_lignes, n_cols):
    """
    Grille booléenne 2D : True si le centre de la cellule est DANS le polygone
    Boeny. Utilise shapely.contains_xy (test géométrique vectorisé, très rapide
    sur les ~10^6 cellules). Renvoie None si aucun fichier limite n'est fourni.
    """
    if CHEMIN_LIMITE is None or not Path(CHEMIN_LIMITE).exists():
        print("AVERTISSEMENT : pas de fichier limite -> aucun masque appliqué.")
        return None

    import geopandas as gpd
    from shapely import contains_xy

    poly = gpd.read_file(CHEMIN_LIMITE).to_crs(crs).geometry.union_all()
    grille_x, grille_y = coordonnees_centres(transformee, n_lignes, n_cols)
    return contains_xy(poly, grille_x, grille_y)


def main():
    os.makedirs(dossier_sortie, exist_ok=True)
    annees = list(range(annee_debut, annee_fin + 1))

    print(f"Agrégation des {len(annees)} composites annuels (10 m -> 250 m)...")
    empile = []           # grilles (lignes, cols) par année
    metadonnees = None
    for annee in annees:
        chemin = Path(dossier_entree) / f"boeny_occupation_sol_{annee}.tif"
        if not chemin.exists():
            raise FileNotFoundError(f"Fichier manquant pour {annee} : {chemin}")
        grille, transformee, crs = agregation_modale(str(chemin), facteur_agrega)
        if metadonnees is None:
            metadonnees = (transformee, crs, grille.shape)
        elif grille.shape != metadonnees[2]:
            raise ValueError(f"Grille incohérente pour {annee} ({grille.shape} != {metadonnees[2]})")
        empile.append(grille)

    n_lignes, n_cols = metadonnees[2]
    transformee, crs = metadonnees[:2]

    # ---------------------------------------------------------------
    # 0 bis. Masque polygone Boeny (avant écriture du raster et des tables)
    # ---------------------------------------------------------------
    masque = masque_boeny(crs, transformee, n_lignes, n_cols)
    if masque is not None:
        n_tot = masque.size
        n_dans = int(masque.sum())
        print(f"Mascage polygone : {n_dans} cellules dans Boeny sur {n_tot} "
              f"({100 * n_dans / n_tot:.1f} %) ; {n_tot - n_dans} exclues.")
        # Hors-polygone -> VALEUR_INVALIDE dans TOUTES les bandes
        for i in range(len(annees)):
            empile[i][~masque] = VALEUR_INVALIDE

    # ---------------------------------------------------------------
    # 1. Raster multi-bandes 250 m (une bande par année)
    # ---------------------------------------------------------------
    chemin_raster = Path(dossier_sortie) / f"boeny_occupation_sol_250m_{annee_debut}_{annee_fin}.tif"
    empile_arr = np.stack(empile, axis=0)  # (annees, lignes, cols)
    with rasterio.open(
        chemin_raster, "w", driver="GTiff", height=n_lignes, width=n_cols,
        count=len(annees), dtype="uint8", crs=crs, transform=transformee,
        nodata=VALEUR_INVALIDE,
    ) as dst:
        for i in range(len(annees)):
            dst.write(empile_arr[i], i + 1)
    print(f"Raster multi-bandes 250 m exporté : {chemin_raster}")

    # ---------------------------------------------------------------
    # 2. Table "LARGE" : 1 ligne = 1 cellule, 1 colonne par année
    # ---------------------------------------------------------------
    grille_x, grille_y = coordonnees_centres(transformee, n_lignes, n_cols)
    # Matrice (n_cellules, 10) des codes
    codes = empile_arr.reshape(len(annees), -1).T.astype(np.int16)  # colonnes = années
    # On ne garde que les cellules 100 % renseignées (équivalent na.rm=TRUE)
    valide = np.all(np.isin(codes, list(NOMS_CLASSES.keys())), axis=1)
    if masque is not None:
        valide &= masque.ravel()  # en plus, le centre de cellule doit être dans Boeny
    codes_valides = codes[valide]
    xs = grille_x.ravel()[valide]
    ys = grille_y.ravel()[valide]

    table_large = pd.DataFrame(codes_valides, columns=[str(a) for a in annees])
    table_large.insert(0, "cell_id", np.arange(1, len(codes_valides) + 1))
    table_large.insert(1, "x", np.round(xs, 2))
    table_large.insert(2, "y", np.round(ys, 2))
    chemin_large = Path(dossier_sortie) / "trajectoires_large_250m.csv"
    table_large.to_csv(chemin_large, index=False)
    n_cellules = len(table_large)
    print(f"Table large exportée : {n_cellules} cellules ({chemin_large})")

    # ---------------------------------------------------------------
    # 3. Table "LONGUE" : 1 ligne = 1 cellule x 1 année
    # ---------------------------------------------------------------
    lignes = []
    for i, annee in enumerate(annees):
        lignes.append(pd.DataFrame({
            "cell_id": table_large["cell_id"],
            "x": table_large["x"],
            "y": table_large["y"],
            "annee": annee,
            "classe_code": codes_valides[:, i],
            "classe": pd.Series(codes_valides[:, i]).astype(str).map(NOMS_CLASSES_STR),
        }))
    table_longue = pd.concat(lignes, ignore_index=True).sort_values(["cell_id", "annee"])
    table_longue.to_csv(Path(dossier_sortie) / "trajectoires_longue_250m.csv", index=False)
    print(f"Table longue exportée : {len(table_longue)} lignes.")

    # ---------------------------------------------------------------
    # 4. Matrice de transition exemple : annee_debut -> annee_fin
    # ---------------------------------------------------------------
    md = pd.crosstab(
        table_large[str(annee_debut)].astype(str).map(NOMS_CLASSES_STR),
        table_large[str(annee_fin)].astype(str).map(NOMS_CLASSES_STR),
    )
    print(f"\nMatrice de transition {annee_debut} -> {annee_fin} :")
    print(md)
    md.to_csv(Path(dossier_sortie) / f"matrice_transition_{annee_debut}_{annee_fin}.csv")


if __name__ == "__main__":
    main()