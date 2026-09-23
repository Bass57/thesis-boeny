"""
Visualisation spatiale des résultats - Boeny
=============================================

Ce script sert de REFERENCE : il suppose que vous disposez déjà, à l'issue
de la fouille de données (R), d'un fichier CSV avec au minimum :
    cell_id, x, y, cluster        (résultat du clustering de trajectoires)
et éventuellement les colonnes de classes par année (2016..2025) pour les
graphiques d'évolution.

Quatre outils, quatre usages différents dans votre mémoire :
  1. GeoPandas + Matplotlib -> cartes STATIQUES pour le document LaTeX
  2. Cartopy                -> cartes statiques avec contexte géographique
                               (littoral, graticule) quand une simple carte
                               de points ne suffit pas
  3. Folium                 -> carte INTERACTIVE HTML, utile pour explorer
                               vos résultats ou pour la démonstration en
                               soutenance (pas pour le document imprimé)
  4. Plotly                 -> graphiques interactifs (évolution des
                               classes, support/confiance des règles
                               d'association)

Installation :
    pip install geopandas matplotlib cartopy folium plotly
"""

import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd

# Pour une ébauche rapide, limitez le nombre de points des cartes statiques
# (None = tous les points, ~482 000 -> figure plus longue)
ECHANTILLON_CARTE = None  # ex. 150_000

# -----------------------------------------------------------------------
# 0. CHARGEMENT DES DONNEES (à adapter à vos fichiers réels)
# -----------------------------------------------------------------------

# Exemple : table_large (trajectoires_large_250m.csv) enrichie d'une
# colonne "cluster" après le clustering fait en R (à exporter depuis R
# avec write.csv, ou à recalculer ici en Python selon votre préférence).
df = pd.read_csv("donnees/boeny_250m/trajectoires_avec_clusters.csv")

# -----------------------------------------------------------------------
# NOMS DES CLUSTERS (à ajuster selon votre interprétation)
# -----------------------------------------------------------------------
# Proposition ISSUE DE L'ANALYSE des profils (voir profil_clusters.csv +
# composition par classe dominante, K=6 final) :
NOMS_CLUSTERS = {
    1: "Eau stable (vasières / estuaire)",
    2: "Cultures / sol nu (agriculture)",
    3: "Transition forêt-savane (mosaïque, dynamique)",
    4: "Végétation inondée (côte)",
    5: "Savane stable",
    6: "Forêt stable",
}
df["type"] = df["cluster"].map(NOMS_CLUSTERS).fillna("Autre")

# Les coordonnées x, y issues du raster R sont en UTM 38S (EPSG:32738),
# cohérent avec l'export du script d'agrégation.
gdf = gpd.GeoDataFrame(
    df,
    geometry=gpd.points_from_xy(df.x, df.y),
    crs="EPSG:32738",
)

# -----------------------------------------------------------------------
# 1. GEOPANDAS + MATPLOTLIB - carte statique pour le document
# -----------------------------------------------------------------------
# Usage : figure à insérer telle quelle dans le chapitre 3 (LaTeX).
# NB : 482 359 points s'affichent bien mais la figure peut prendre ~1 min ;
# réglez ECHANTILLON_CARTE (ex. 150_000) pour une première ébauche rapide.

fig, ax = plt.subplots(figsize=(8, 8))
ech_carte = gdf.sample(ECHANTILLON_CARTE, random_state=1) if ECHANTILLON_CARTE else gdf
ech_carte.plot(
    column="type",
    categorical=True,
    legend=True,
    markersize=3,
    cmap="tab10",
    ax=ax,
)
ax.set_title("Typologie des trajectoires de changement - Boeny (2016-2025)")
ax.set_axis_off()
plt.tight_layout()
plt.savefig("carte_clusters_trajectoires.png", dpi=300)
plt.close()
print("Carte statique exportée : carte_clusters_trajectoires.png")

# -----------------------------------------------------------------------
# 2. CARTOPY - carte avec contexte géographique (côte, graticule)
# -----------------------------------------------------------------------
# Usage : quand une carte de points seule manque de repères géographiques
# (utile en figure d'introduction pour situer Boeny, moins nécessaire pour
# les cartes de résultats détaillées où GeoPandas seul suffit).

import cartopy.crs as ccrs
import cartopy.feature as cfeature

gdf_wgs84 = gdf.to_crs(epsg=4326)  # Cartopy travaille en lat/lon

fig = plt.figure(figsize=(8, 8))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([45.5, 49.5, -18.3, -15.0], crs=ccrs.PlateCarree())  # bbox Boeny
ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax.add_feature(cfeature.BORDERS, linestyle=":")
ax.gridlines(draw_labels=True, linewidth=0.3)

ax.scatter(
    gdf_wgs84.geometry.x,
    gdf_wgs84.geometry.y,
    c=gdf_wgs84["cluster"],
    cmap="tab10",
    s=2,
    transform=ccrs.PlateCarree(),
)
ax.set_title("Localisation des trajectoires de changement - Boeny")
plt.savefig("carte_contexte_boeny.png", dpi=300, bbox_inches="tight")
plt.close()
print("Carte avec contexte géographique exportée : carte_contexte_boeny.png")

# -----------------------------------------------------------------------
# 3. FOLIUM - carte interactive HTML (exploration / soutenance)
# -----------------------------------------------------------------------
# Usage : à ouvrir dans un navigateur pour explorer vos résultats, ou à
# montrer en direct pendant la soutenance. NE PAS l'inclure telle quelle
# dans le document imprimé (ce n'est pas un format statique).

import folium

centre = [gdf_wgs84.geometry.y.mean(), gdf_wgs84.geometry.x.mean()]
carte = folium.Map(location=centre, zoom_start=9, tiles="OpenStreetMap")

couleurs = ["darkred", "darkblue", "darkgreen", "purple", "orange",
            "brown", "cyan", "magenta", "olive", "teal"]
couleurs_cluster = {
    c: couleurs[i % len(couleurs)]
    for i, c in enumerate(sorted(gdf_wgs84["cluster"].dropna().unique()))
}

# Pour de gros volumes de points, échantillonner avant d'afficher
# (Folium devient lent au-delà de quelques milliers de marqueurs)
echantillon = gdf_wgs84.sample(min(3000, len(gdf_wgs84)), random_state=1)

for _, row in echantillon.iterrows():
    folium.CircleMarker(
        location=[row.geometry.y, row.geometry.x],
        radius=2,
        color=couleurs_cluster.get(row["cluster"], "gray"),
        fill=True,
        fill_opacity=0.7,
    ).add_to(carte)

carte.save("carte_interactive_clusters.html")
print("Carte interactive exportée : carte_interactive_clusters.html")

# -----------------------------------------------------------------------
# 4. PLOTLY - graphiques interactifs (évolution des classes, règles)
# -----------------------------------------------------------------------
# Usage : exploration personnelle ou annexe interactive ; pour le document
# imprimé, préférez une version statique (matplotlib) de ces graphiques.

import plotly.express as px

annees = [str(a) for a in range(2016, 2026)]
noms_classes = {
    "0": "eau",
    "1": "foret",
    "2": "vegetation_non_forestiere",
    "3": "vegetation_inondee",
    "4": "cultures",
    "5": "bati",
    "6": "sol_nu",
}

# Reconstruction d'un tableau long (année, classe, % de la superficie)
# à partir des colonnes annuelles de df, si elles sont présentes
if all(a in df.columns for a in annees):
    lignes = []
    for annee in annees:
        proportions = df[annee].value_counts(normalize=True) * 100
        for code, pct in proportions.items():
            lignes.append(
                {
                    "annee": int(annee),
                    "classe": noms_classes.get(str(code), str(code)),
                    "pourcentage": pct,
                }
            )
    evolution = pd.DataFrame(lignes)

    fig = px.line(
        evolution,
        x="annee",
        y="pourcentage",
        color="classe",
        title="Évolution de l'occupation du sol - Boeny (interactif)",
        markers=True,
    )
    fig.write_html("evolution_classes_interactif.html")
    print("Graphique interactif exporté : evolution_classes_interactif.html")
