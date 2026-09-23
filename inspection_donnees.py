"""
Inspection de la qualité des trajectoires - Boeny (250 m, 2016-2025)
====================================================================
Équivalent Python de `inspection_donnees.R` (mêmes entrées/sorties).

Entrée  : donnees/boeny_250m/trajectoires_large_250m.csv
          (produit par agregation_trajectoires_boeny.py)
Sorties :
  - evolution_classes_pourcentage.csv  (% de chaque classe par année)
  - evolution_classes.png              (courbes)
  - diagnostic_bruit_cellules.csv      (changements, aller-retours par cellule)
  - messages console : distribution de la volatilité, % de cellules "bruit"

Objectifs pédagogiques du diagnostic :
  1. Une classe quasi vide (< 1 %) sur toute la période n'est pas exploitable
     en fouille -> envisager une fusion avec une voisine.
  2. Une cellule qui change de classe presque chaque année est SUSPECTE : les
     vraies conversions (forêt -> culture, culture -> bâti...) sont rarement
     annuelles sur des composites. C'est souvent un artefact (nuage résiduel,
     confusion du classifieur).
  3. Un "aller-retour" en 1 an (A -> B -> A) est l'archétype du bruit.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")  # mode non-interactif (aucune fenêtre)
import matplotlib.pyplot as plt

# ============================ CONFIGURATION ============================
dossier = "donnees/boeny_250m"
annees = list(range(2016, 2026))
NOMS_CLASSES = {
    0: "eau", 1: "foret", 2: "vegetation_non_forestiere",
    3: "vegetation_inondee", 4: "cultures", 5: "bati", 6: "sol_nu",
}
SEUIL_CLASSE_QUASI_VIDE = 1.0   # % minimal (1 % d'après le guide)
SEUIL_CHANGEMENTS = 5           # au-delà : cellule jugée trop volatile
# =======================================================================


def evolution_classes(table_large):
    """% de chaque classe par année (CSV + graphique PNG)."""
    colonnes = [str(a) for a in annees]
    proportions = pd.DataFrame(index=sorted(NOMS_CLASSES.values()))
    for annee in annees:
        decompte = table_large[str(annee)].value_counts(normalize=True) * 100
        for code, nom in NOMS_CLASSES.items():
            proportions.loc[nom, annee] = decompte.get(code, 0.0)
    proportions = proportions.round(2)

    print("Pourcentage de chaque classe par année :")
    print(proportions)
    (Path(dossier) / "evolution_classes_pourcentage.csv").write_text(
        proportions.to_csv(index_label="classe"), encoding="utf-8"
    )

    fig, ax = plt.subplots(figsize=(9, 6))
    couleurs = plt.cm.rainbow(np.linspace(0, 1, len(proportions)))
    for i, (nom, serie) in enumerate(proportions.iterrows()):
        ax.plot(annees, serie.values, label=nom, lw=2, color=couleurs[i])
    ax.set_xlabel("Année")
    ax.set_ylabel("% de la superficie")
    ax.set_title("Évolution de l'occupation du sol - Boeny (250 m)")
    ax.legend(loc="best", fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(Path(dossier) / "evolution_classes.png", dpi=150)
    plt.close(fig)

    quasi_vides = [nom for nom, s in proportions.iterrows() if s.max() < SEUIL_CLASSE_QUASI_VIDE]
    if quasi_vides:
        print(f"\nATTENTION - classes représentant moins de {SEUIL_CLASSE_QUASI_VIDE:.0f}% sur toute la période :")
        print("  ", quasi_vides)
        print("A envisager : les fusionner avec une classe voisine avant la fouille.")
    return proportions


def diagnostique_bruit(table_large):
    """Volatilité (nb de changements) et aller-retours par cellule."""
    matrice = table_large[[str(a) for a in annees]].to_numpy()
    nb_changements = (np.diff(matrice, axis=1) != 0).sum(axis=1)

    print("\nDistribution du nombre de changements par cellule (sur 9 transitions) :")
    print(pd.Series(nb_changements).describe().round(2).to_string())
    print(pd.Series(nb_changements).value_counts().sort_index().to_string())

    # Aller-retour A -> B -> A en deux transitions consécutives
    if matrice.shape[1] >= 3:
        nb_flipflops = (
            (matrice[:, :-2] == matrice[:, 2:]) & (matrice[:, :-2] != matrice[:, 1:-1])
        ).sum(axis=1)
    else:
        nb_flipflops = np.zeros(len(matrice), dtype=int)

    diagnostic = pd.DataFrame({
        "cell_id": table_large["cell_id"],
        "x": table_large["x"],
        "y": table_large["y"],
        "nb_changements": nb_changements,
        "nb_flipflops": nb_flipflops,
        "bruit_suspect": (nb_flipflops > 0) | (nb_changements >= SEUIL_CHANGEMENTS),
    })
    diagnostic.to_csv(Path(dossier) / "diagnostic_bruit_cellules.csv", index=False)

    pct = round(100 * diagnostic["bruit_suspect"].mean(), 1)
    n_bruit = int(diagnostic["bruit_suspect"].sum())
    print(f"\n{pct}% des cellules signalées comme bruit potentiel "
          f"({n_bruit} cellules sur {len(diagnostic)}).")
    print("Fichier détaillé : diagnostic_bruit_cellules.csv")
    return diagnostic


def main():
    chemin = Path(dossier) / "trajectoires_large_250m.csv"
    if not chemin.exists():
        raise FileNotFoundError(f"{chemin} introuvable -- lancez d'abord agregation_trajectoires_boeny.py")
    table_large = pd.read_csv(chemin)
    evolution_classes(table_large)
    diagnostique_bruit(table_large)


if __name__ == "__main__":
    main()