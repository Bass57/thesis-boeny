"""
Clustering des trajectoires + règles d'association - Boeny
===========================================================
Équivalent Python de `clustering_regles_association_boeny.R` (mêmes
entrées/sorties, mêmes hypothèses méthodologiques).

Entrée  : donnees/boeny_250m/trajectoires_large_250m_lisse.csv (table LISSÉE,
          produite par lissage_trajectoires_boeny.py ; passer par
          FICHIER_ENTREE pour utiliser la table brute)
Sorties :
  - coude_choix_k.png                     (graphique du coude -> choisir K_FINAL)
  - trajectoires_avec_clusters.csv        (attendu par visualisation_spatiale_boeny.py)
  - profil_clusters.csv                   (moyennes des indicateurs par cluster)
  - regles_association_transitions.csv    (règles Apriori, triées par lift)
  - regles_association_chaines.csv        (règles "en chaîne" seulement)

Hypothèse méthodologique (identique au script R, cf. commentaire d'origine) :
Une distance d'appariement optimal (TraMineR, cf. Mas et al. 2019) exigerait
une matrice de dissimilarité entre TOUTES les paires de cellules, ce qui est
infeasible (~4,8e5 cellules -> ~1,15e11 paires). On résume donc chaque
trajectoire par des indicateurs (entropie, volatilité, classes de début/fin/
dominante), puis on applique un K-Means sur ces indicateurs standardisés.

Les deux algorithmes (K-Means et Apriori) sont implémentés "à la main" en
NumPy/NumPy pur, sans scikit-learn : pédagogique et suffisant ici
(respectivement <= 7 classes et <= 42 transitions distinctes).
"""

import itertools
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")  # mode non-interactif
import matplotlib.pyplot as plt

# ============================ CONFIGURATION ============================
dossier = "donnees/boeny_250m"
# Table d'entrée : la table LISSÉE par défaut (elle remplace la brute pour la
# fouille). Remettre "trajectoires_large_250m.csv" pour travailler sur la brute.
FICHIER_ENTREE = "trajectoires_large_250m_lisse.csv"
annees = list(range(2016, 2026))
NOMS_CLASSES = {
    0: "eau", 1: "foret", 2: "vegetation_non_forestiere",
    3: "vegetation_inondee", 4: "cultures", 5: "bati", 6: "sol_nu",
}
K_FINAL = 6  # coude + silhouette : inertie 5->6 = -20,7 % (dernière grosse
             # chute), silhouette MAX à K=6 (0,765) ; K=7 ne gagne que -0,9 % ;
             # K=8 est DÉGÉNÉRÉ (inertie identique à K=6 -> cluster vide, à
             # écarter -- d'où la non-monotonie vue au premier coude).
K_TESTES = range(2, 9)      # valeurs de K pour la méthode du coude
NSTART_COUDE = 25           # nb d'initialisations du coude (5 -> minima locaux possibles)
TAILLE_SILHOUETTE = 5000    # taille de l'échantillon pour la silhouette
SEED = 1                    # reproductibilité (équivalent de set.seed(1))
FAIRE_COUDE = False         # True : re-trace le coude (4-5 min, nstart=25).
                            # False : NE refait PAS le coude (K=6 déjà validé) et
                            # va directement au clustering final + règles (rapide).
                            # -> Mettre True une dernière fois pour la figure
                            #    "production" du mémoire si besoin.
SEUIL_SUPPORT = 0.01        # support minimal des règles (paramètre 'supp')
SEUIL_CONFIANCE = 0.5       # confiance minimale (paramètre 'conf')
# =======================================================================


# ---------------------------------------------------------------------------
# OUTILS GÉNÉRIQUES
# ---------------------------------------------------------------------------

def entropie_shannon(x):
    """H = -sum(p * log2(p)) sur les classes présentes (0 si une seule classe)."""
    _, comptes = np.unique(x, return_counts=True)
    p = comptes / comptes.sum()
    return float(-(p * np.log2(p)).sum())


def classe_dominante(x):
    """Classe la plus fréquente ; en cas d'égalité, la plus petite (comme R)."""
    valeurs, comptes = np.unique(x, return_counts=True)
    return int(valeurs[np.argmax(comptes)])


def zscore(mat):
    """Standardisation colonne par colonne (équivalent de R `scale`, ddof=1)."""
    mat = mat.astype(np.float64)
    moyennes = mat.mean(axis=0)
    ecarts = mat.std(axis=0, ddof=1)
    if np.any(ecarts == 0):
        nuls = int(np.any(ecarts == 0))
        print(f"AVERTISSEMENT : {nuls} colonne(s) à variance nulle supprimée(s) avant scaling.")
        garder = ecarts != 0
        return mat[:, garder] - moyennes[garder], garder
    return (mat - moyennes) / ecarts, np.ones(mat.shape[1], dtype=bool)


def kmeans_numpy(X, k, nstart=1, iter_max=50, seed=SEED):
    """
    K-Means classique (style kmeans() de R) : nstart initialisations aléatoires,
    on garde celle de plus faible inertie intra-cluster (tot.withinss).
    """
    rng = np.random.default_rng(seed)
    n = X.shape[0]
    normes_x = (X ** 2).sum(axis=1)  # précalcul pour les distances

    meilleur = None
    for _ in range(nstart):
        indices = rng.choice(n, size=k, replace=False)
        centres = X[indices].copy()
        labels_prev = None
        for _ in range(iter_max):
            # ||x - c||^2 = ||x||^2 + ||c||^2 - 2 x.c  (décomposition évitant 3D)
            distances = (
                normes_x[:, None]
                + (centres ** 2).sum(axis=1)[None, :]
                - 2.0 * (X @ centres.T)
            )
            labels = np.argmin(distances, axis=1)
            if labels_prev is not None and np.array_equal(labels, labels_prev):
                break
            labels_prev = labels
            nouveaux = np.zeros_like(centres)
            for c in range(k):
                masque = labels == c
                if masque.any():
                    nouveaux[c] = X[masque].mean(axis=0)
                else:  # cluster vide : on conserve l'ancien centre
                    nouveaux[c] = centres[c]
            centres = nouveaux
        tot_withinss = float(np.take_along_axis(distances, labels[:, None], axis=1).sum())
        if meilleur is None or tot_withinss < meilleur[0]:
            meilleur = (tot_withinss, labels.copy(), centres.copy())
    return meilleur  # (tot_withinss, labels, centres)


# ---------------------------------------------------------------------------
# PARTIE A : CLUSTERING DES TRAJECTOIRES
# ---------------------------------------------------------------------------

def construire_indicateurs(table_large):
    """Indicateurs de trajectoire par cellule (entropie, changements, classes...)."""
    matrice_codes = table_large[[str(a) for a in annees]].to_numpy()
    indicateurs = pd.DataFrame({
        "cell_id": table_large["cell_id"],
        "entropie": np.apply_along_axis(entropie_shannon, 1, matrice_codes),
        "nb_changements": (np.diff(matrice_codes, axis=1) != 0).sum(axis=1),
        "nb_classes_visitees": np.apply_along_axis(lambda x: len(np.unique(x)), 1, matrice_codes),
        "classe_debut": matrice_codes[:, 0],
        "classe_fin": matrice_codes[:, -1],
        "classe_dominante": np.apply_along_axis(classe_dominante, 1, matrice_codes),
    })
    return matrice_codes, indicateurs


def caracteristiques_model_matrix(indicateurs):
    """
    Équivalent de model.matrix(~ entropie + nb_changements + nb_classes_visitees
    + classe_debut + classe_fin + classe_dominante - 1) : les variables
    numériques + une colonne indicatrice par modalité des variables catégorielles,
    PUIS standardisation globale (comme scale() dans le script R).
    """
    base = indicateurs[["entropie", "nb_changements", "nb_classes_visitees"]].reset_index(drop=True)
    dummies = pd.concat([
        pd.get_dummies(indicateurs[c], prefix=c, dtype=float)
        for c in ("classe_debut", "classe_fin", "classe_dominante")
    ], axis=1)
    # Reconstruire le nominage R (classe_debut_0, ... classe_dominante_6)
    dummies_noms = {}
    for c in ("classe_debut", "classe_fin", "classe_dominante"):
        for code in range(7):
            dummies_noms[f"{c}_{code}"] = 1.0 * (indicateurs[c].to_numpy() == code)
    colonnes = pd.DataFrame(dummies_noms)
    return zscore(pd.concat([base, colonnes], axis=1).to_numpy())


def silhouette_echantillon(X, labels, taille=TAILLE_SILHOUETTE, seed=SEED):
    """
    Silhouette moyenne (Rousseeuw, 1987) estimée sur un échantillon aléatoire.
    s(i) = (b - a) / max(a, b) avec :
      a(i) = distance moyenne de i aux autres points de SON cluster ;
      b(i) = distance moyenne de i aux points du cluster VOISIN le plus proche.
    s proche de 1 : clusters compacts et bien séparés ; ~0 : chevauchement ;
    négatif : points mal classés. Le calcul exact sur 482 000 points serait en
    O(n^2) (infaisable), d'où l'échantillon (pratique standard).
    """
    rng = np.random.default_rng(seed)
    n = X.shape[0]
    m_ech = min(taille, n)
    indices = rng.choice(n, size=m_ech, replace=False)
    Xe = X[indices]
    lab = labels[indices]
    k = int(lab.max()) + 1
    ens_clusters = {c: np.where(lab == c)[0] for c in range(k)}
    # Matrice de distances (échantillon x échantillon) via la décomposition
    # ||x-y||^2 = ||x||^2 + ||y||^2 - 2 x.y
    normes = (Xe ** 2).sum(axis=1)
    D = normes[:, None] + normes[None, :] - 2.0 * (Xe @ Xe.T)
    np.maximum(D, 0.0, out=D)

    silhouette = np.zeros(m_ech)
    for c, idx_c in ens_clusters.items():
        nc = len(idx_c)
        if nc < 2:
            continue  # convention : silhouette = 0 pour un cluster d'un point
        # a(i) : moyenne des distances vers les autres points du cluster
        # (la diagonale D[i,i] = 0, donc diviser par (nc - 1) exclut i lui-même)
        a = D[:, idx_c].sum(axis=1) / (nc - 1)
        # b(i) : min des distances moyennes vers chaque autre cluster
        b = np.full(m_ech, np.inf)
        for c2, idx2 in ens_clusters.items():
            if c2 == c or len(idx2) == 0:
                continue
            np.minimum(b, D[:, idx2].mean(axis=1), out=b)
        denom = np.maximum(a, b)
        s = np.where(denom > 0, (b - a) / denom, 0.0)
        silhouette[idx_c] = s[idx_c]
    return float(silhouette.mean())


def methode_du_coude(X, k_testes):
    """
    Inertie intra-cluster + silhouette moyenne pour différentes valeurs de K
    (double axe, sauvegardé en PNG). nstart élevé pour éviter les minima locaux.
    """
    inerties, silhouettes = [], []
    for k in k_testes:
        debut = time.time()
        tot_withinss, labels, _ = kmeans_numpy(X, k, nstart=NSTART_COUDE, iter_max=50)
        inerties.append(tot_withinss)
        silhouettes.append(silhouette_echantillon(X, labels))
        print(f"K={k} : inertie = {tot_withinss:,.0f} | silhouette = "
              f"{silhouettes[-1]:.3f} | {time.time() - debut:.0f}s")

    ks = list(k_testes)
    fig, ax1 = plt.subplots(figsize=(8, 5))
    ax1.plot(ks, inerties, "o-", color="tab:blue", label="Inertie intra-cluster")
    ax1.set_xlabel("Nombre de clusters (K)")
    ax1.set_ylabel("Inertie intra-cluster", color="tab:blue")
    ax1.tick_params(axis="y", labelcolor="tab:blue")
    ax2 = ax1.twinx()
    ax2.plot(ks, silhouettes, "s--", color="tab:red",
             label=f"Silhouette (échantillon {TAILLE_SILHOUETTE})")
    ax2.set_ylabel("Silhouette moyenne", color="tab:red")
    ax2.tick_params(axis="y", labelcolor="tab:red")
    lignes1, labels1 = ax1.get_legend_handles_labels()
    lignes2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lignes1 + lignes2, labels1 + labels2, loc="best", fontsize=8)
    ax1.set_title("Choix du nombre de clusters (coude + silhouette)")
    ax1.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(Path(dossier) / "coude_choix_k.png", dpi=150)
    plt.close(fig)
    print(f"\nGraphique exporté : coude_choix_k.png")
    print("Critère : K au coude de l'inertie ET/OU silhouette maximale (avec un "
          "compromis d'interprétabilité).")


def clustering_final(X, table_large, indicateurs):
    """K-Means final sur les indicateurs, export des clusters et profils."""
    tot_withinss, labels, _ = kmeans_numpy(X, K_FINAL, nstart=25, iter_max=100)
    table_large["cluster"] = labels + 1  # R numérote les clusters de 1 à K
    print(f"\nRépartition des cellules par cluster :")
    print(table_large["cluster"].value_counts().sort_index().to_string())

    table_large.to_csv(Path(dossier) / "trajectoires_avec_clusters.csv", index=False)
    print("Table exportée : trajectoires_avec_clusters.csv")

    profil = (
        indicateurs[["entropie", "nb_changements", "nb_classes_visitees"]]
        .assign(cluster=labels + 1)
        .groupby("cluster")
        .mean()
    )
    print("\nProfil moyen des clusters (pour interprétation) :")
    print(profil.round(3).to_string())
    profil.to_csv(Path(dossier) / "profil_clusters.csv")
    print("Exporté : profil_clusters.csv")


# ---------------------------------------------------------------------------
# PARTIE B : RÈGLES D'ASSOCIATION ENTRE TRANSITIONS (Apriori)
# ---------------------------------------------------------------------------
# Principe (Mennis & Liu, 2003) : chaque cellule est une "transaction" contenant
# l'ensemble des transitions de classe qu'elle a connues (ex. "foret->sol_nu").
# Apriori cherche les transitions qui co-apparaissent fréquemment dans les
# mêmes cellules -> enchaînements de conversion en plusieurs étapes.

def construire_transactions(matrice_codes):
    """Une liste d'items (transitions uniques) par cellule."""
    def transitions_cellule(seq):
        transitions = set()
        for t in range(len(seq) - 1):
            de, vers = int(seq[t]), int(seq[t + 1])
            if de != vers:
                transitions.add(f"{NOMS_CLASSES[de]}->{NOMS_CLASSES[vers]}")
        return transitions

    lignes = [transitions_cellule(matrice_codes[i]) for i in range(matrice_codes.shape[0])]
    actifs = [t for t in lignes if t]  # cellules sans transition = rien à apprendre
    print(f"\n{len(actifs)} cellules avec au moins une transition, "
          f"sur {matrice_codes.shape[0]} au total.")
    return actifs


def frequent_itemsets(transactions, min_supp):
    """Retourne (itemsets_freq: liste de frozenset, support: dict, n)."""
    n = len(transactions)
    support = {}

    def est_frequent(itemset):
        if itemset not in support:
            support[itemset] = sum(1 for t in transactions if itemset.issubset(t))
        return support[itemset] / n >= min_supp

    # 1-itemsets fréquents
    unaires = {it for t in transactions for it in t}
    L = {1: [frozenset({it}) for it in unaires if est_frequent(frozenset({it}))]}

    k = 2
    while L[k - 1]:
        lus = set(L[k - 1])  # itemsets fréquents de taille k-1 (pour l'élagage)
        candidats = []
        for i, a in enumerate(L[k - 1]):
            for b in L[k - 1][i + 1:]:
                if sorted(a)[:-1] == sorted(b)[:-1]:  # fusion par préfixe commun
                    cand = a | b
                    # Élagage Apriori : tout sous-ensemble de taille k-1 d'un
                    # candidat doit lui-même être fréquent.
                    if all((cand - {x}) in lus for x in cand):
                        candidats.append(cand)
        L[k] = [c for c in dict.fromkeys(candidats) if est_frequent(c)]
        k += 1

    itemsets_freq = [itemset for niveau in L.values() for itemset in niveau]
    return itemsets_freq, support, n


def regles_association(transactions, min_supp=SEUIL_SUPPORT, min_conf=SEUIL_CONFIANCE):
    """
    Génère toutes les règles antécédent -> conséquent (itemsets fréquents de
    taille >= 2) avec support, confiance, lift -- équivalent de arules::apriori
    puis as(rules, "data.frame") avec minlen=2.
    """
    itemsets_freq, support, n = frequent_itemsets(transactions, min_supp)

    regles = []
    for itemset in itemsets_freq:
        if len(itemset) < 2:
            continue
        supp_zs = support[itemset] / n
        for taille_a in range(1, len(itemset)):
            for A in itertools.combinations(itemset, taille_a):
                A, B = frozenset(A), itemset - frozenset(A)
                supp_a, supp_b = support[A] / n, support[B] / n
                confiance = supp_zs / supp_a
                if confiance >= min_conf:
                    regles.append({
                        "lhs": ",".join(sorted(A)),
                        "rhs": ",".join(sorted(B)),
                        "support": supp_zs,
                        "confidence": confiance,
                        "coverage": supp_a,
                        "lift": supp_zs / (supp_a * supp_b),
                        "count": support[itemset],
                    })
    df = pd.DataFrame(regles).sort_values("lift", ascending=False, ignore_index=True)
    return df


def extraire_vers(item):
    """Classe d'arrivée de (la dernière) transition d'un item -- équivalent
    de `sub(".*->", "", ...)` en R (greedy)."""
    return re.sub(r".*->", "", item, count=1)


def extraire_de(item):
    """Classe de départ de la première transition -- équivalent de
    `sub("->.*", "", ...)` en R."""
    return re.sub(r"->.*", "", item, count=1)


def identifier_chaines(df_regles):
    """Règles "en chaîne" : classe d'arrivée du précurseur = classe de départ
    du conséquent (ex. {foret->sol_nu} => {sol_nu->cultures})."""
    vers_lhs = df_regles["lhs"].map(extraire_vers)
    de_rhs = df_regles["rhs"].map(extraire_de)
    return df_regles[vers_lhs == de_rhs].copy()


# ---------------------------------------------------------------------------
# SCRIPT PRINCIPAL
# ---------------------------------------------------------------------------

def main():
    chemin = Path(dossier) / FICHIER_ENTREE
    if not chemin.exists():
        raise FileNotFoundError(
            f"{chemin} introuvable -- lancez d'abord lissage_trajectoires_boeny.py "
            f"(ou modifiez FICHIER_ENTREE)"
        )
    table_large = pd.read_csv(chemin)

    # ---- A. Clustering ------------------------------------------------
    print("Construction des indicateurs de trajectoire...")
    matrice_codes, indicateurs = construire_indicateurs(table_large)
    X, _ = caracteristiques_model_matrix(indicateurs)
    print(f"Matrice de caractéristiques : {X.shape} (cellules x variables standardisées)")

    if FAIRE_COUDE:  # coude déjà validé (K=6) : le sauter pour gagner 4-5 min
        methode_du_coude(X, K_TESTES)
    clustering_final(X, table_large, indicateurs)

    # ---- B. Règles d'association --------------------------------------
    transactions = construire_transactions(matrice_codes)
    df_regles = regles_association(transactions)
    print(f"\n{len(df_regles)} règles extraites "
          f"(support >= {SEUIL_SUPPORT}, confiance >= {SEUIL_CONFIANCE}).")
    print("\nTop 10 règles par lift :")
    print(df_regles.head(10).to_string(index=False))

    df_regles.to_csv(Path(dossier) / "regles_association_transitions.csv", index=False)
    chaines = identifier_chaines(df_regles)
    print(f"{len(chaines)} règles \"en chaîne\" identifiées (sur {len(df_regles)}).")
    chaines.to_csv(Path(dossier) / "regles_association_chaines.csv", index=False)
    print("Fichiers exportés : regles_association_transitions.csv, regles_association_chaines.csv")


if __name__ == "__main__":
    main()