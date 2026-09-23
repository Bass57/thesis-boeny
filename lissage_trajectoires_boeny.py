"""
Lissage temporel glissant 3 ans - Boeny
========================================
Étape 4 bis : préparation des trajectoires avant la fouille de données.

Le diagnostic (inspection_donnees.py) a signalé > 15 % de cellules bruitées,
majoritairement des ALLER-RETOURS en 1 an (ex. foret 2019 -> veg_non_for 2020
-> foret 2021) : archétype du bruit de composite annuel (nuage résiduel, feu
de savane, confusion du classifieur pendant une seule année).

Principe du filtre (filtre majoritaire glissant sur 3 ans) :
    pour chaque année intérieure t, on regarde (t-1, t, t+1) :
      - si DEUX années sur trois partagent la même classe -> t prend cette
        classe (consensus majoritaire) ;
      - sinon (trois classes différentes, pas de consensus) -> on GARDE la
        valeur observée en t (on n'invente pas de classe).
    Les années de bord 2016 et 2025 sont laissées inchangées (fenêtre
    incomplète).

Entrée  : donnees/boeny_250m/trajectoires_large_250m.csv
Sorties (dans donnees/boeny_250m/) :
  - trajectoires_large_250m_lisse.csv        (table lissée, mêmes colonnes)
  - diagnostic_bruit_apres_lissage.csv       (statistiques bruit par cellule)
  + comparaison avant/après dans la console
"""

from pathlib import Path

import numpy as np
import pandas as pd

# ============================ CONFIGURATION ============================
dossier = "donnees/boeny_250m"
FICHIER_ENTREE = "trajectoires_large_250m.csv"
FICHIER_SORTIE = "trajectoires_large_250m_lisse.csv"
FICHIER_DIAGNOSTIC = "diagnostic_bruit_apres_lissage.csv"
annees = list(range(2016, 2026))
SEUIL_CHANGEMENTS = 5  # même définition de "bruit" que l'inspection
# =======================================================================


def mode_glissant(tab):
    """
    (n, 10) -> (n, 10). Les années intérieures (indices 1..8) sont remplacées
    par la classe majoritaire de la fenêtre (t-1, t, t+1) SI le consensus
    existe (>= 2 années sur 3), sinon inchangées. Les bords sont inchangés.
    """
    M = np.asarray(tab, dtype=np.int16)
    lisse = M.copy()
    # Fenêtres : colonnes (t-1, t, t+1) pour t = 1..8
    fenetres = np.stack([M[:, 0:8], M[:, 1:9], M[:, 2:10]], axis=2)  # (n, 8, 3)
    # Comptage des 7 classes dans chaque fenêtre
    comptes = np.stack([(fenetres == c).sum(axis=2) for c in range(7)], axis=2)  # (n, 8, 7)
    majorite = comptes.max(axis=2) >= 2   # vraie majorité (>= 2 sur 3)
    mode = comptes.argmax(axis=2)         # classe majoritaire (plus petite si égalité)
    lisse[:, 1:9] = np.where(majorite, mode, M[:, 1:9])
    return lisse


def stats_bruit(mat):
    """Statistiques de bruit identiques à inspection_donnees.py."""
    mat = np.asarray(mat)
    nb_changements = (np.diff(mat, axis=1) != 0).sum(axis=1)
    if mat.shape[1] >= 3:
        nb_flipflops = ((mat[:, :-2] == mat[:, 2:]) & (mat[:, :-2] != mat[:, 1:-1])).sum(axis=1)
    else:
        nb_flipflops = np.zeros(mat.shape[0], dtype=int)
    bruit = (nb_flipflops > 0) | (nb_changements >= SEUIL_CHANGEMENTS)
    return nb_changements, nb_flipflops, bruit


def main():
    chemin = Path(dossier) / FICHIER_ENTREE
    if not chemin.exists():
        raise FileNotFoundError(
            f"{chemin} introuvable -- lancez d'abord agregation_trajectoires_boeny.py"
        )
    table = pd.read_csv(chemin)
    colonnes = [str(a) for a in annees]
    mat = table[colonnes].to_numpy()

    nb_ch_avant, nb_ff_avant, bruit_avant = stats_bruit(mat)
    print(f"AVANT lissage : {100 * bruit_avant.mean():.1f} % de cellules bruitées "
          f"({int(bruit_avant.sum())} / {len(mat)})")

    mat_lisse = mode_glissant(mat)
    nb_ch_apres, nb_ff_apres, bruit_apres = stats_bruit(mat_lisse)
    print(f"APRÈS lissage : {100 * bruit_apres.mean():.1f} % de cellules bruitées "
          f"({int(bruit_apres.sum())} / {len(mat_lisse)})")

    # Détail avant / après
    table_comp = pd.DataFrame({
        "avant_nb_changements": nb_ch_avant,
        "avant_nb_flipflops": nb_ff_avant,
        "avant_bruit": bruit_avant,
        "apres_nb_changements": nb_ch_apres,
        "apres_nb_flipflops": nb_ff_apres,
        "apres_bruit": bruit_apres,
    })
    print("\nÉvolution moyenne par cellule :")
    print("- changements :", round(nb_ch_avant.mean(), 3), "->", round(nb_ch_apres.mean(), 3))
    print("- aller-retours :", round(nb_ff_avant.mean(), 3), "->", round(nb_ff_apres.mean(), 3))
    print(table_comp.describe().loc[["mean", "50%", "75%", "max"]].round(2).to_string())

    # Table lissée : mêmes colonnes, mêmes identifiants
    table_lisse = table.copy()
    table_lisse[colonnes] = mat_lisse
    table_lisse.to_csv(Path(dossier) / FICHIER_SORTIE, index=False)
    print(f"\nTable lissée exportée : {Path(dossier) / FICHIER_SORTIE}")

    # Diagnostic comparatif par cellule
    comparison = pd.DataFrame({
        "cell_id": table["cell_id"],
        "x": table["x"],
        "y": table["y"],
        **{col: table_comp[col] for col in table_comp.columns},
    })
    comparison.to_csv(Path(dossier) / FICHIER_DIAGNOSTIC, index=False)
    print(f"Diagnostic comparatif exporté : {Path(dossier) / FICHIER_DIAGNOSTIC}")


if __name__ == "__main__":
    main()