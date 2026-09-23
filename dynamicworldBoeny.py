"""
Extraction de composites annuels Dynamic World pour la région Boeny (Madagascar)
=================================================================================

Objectif
--------
Pour chaque année de 2016 à 2025 :
  1. Filtrer la collection Dynamic World (Google/DYNAMICWORLD/V1) sur l'emprise
     de Boeny et sur l'année considérée ;
  2. Calculer un composite annuel "mode" (classe la plus fréquente par pixel) ;
  3. Fusionner certaines classes (herbe + arbustes/broussailles) ;
  4. Découper sur les limites de Boeny et exporter (GeoTIFF, Google Drive).

Prérequis
---------
    pip install earthengine-api

Avant la première utilisation :
    earthengine authenticate
et disposer d'un projet Google Cloud associé à Earth Engine (obligatoire
depuis 2022). Remplacez EE_PROJECT ci-dessous par l'ID de ce projet.

IMPORTANT — zone d'étude (AOI)
-------------------------------
Les régions administratives malgaches ("Faritra", créées en 2004, dont Boeny)
ne sont PAS présentes dans les jeux de limites administratives intégrés à
Earth Engine (ex. FAO GAUL), qui ne connaissent que les anciennes provinces
("Faritany"). Deux solutions :

  A) RECOMMANDÉ : téléchargez le contour de la région Boeny (ex. shapefile
     GADM niveau 1 pour Madagascar, ou couche officielle du BNGRC/INSTAT),
     importez-le comme asset Earth Engine (onglet "Assets" du Code Editor,
     ou `earthengine upload table`), puis renseignez son ID ci-dessous dans
     AOI_ASSET_ID.

  B) DÉPANNAGE RAPIDE : une bounding box approximative de Boeny est fournie
     par défaut (AOI_BBOX). Elle englobe largement la région mais n'en suit
     pas les limites administratives exactes — à ne pas utiliser pour les
     résultats finaux du mémoire, seulement pour tester le script.
"""

import ee

# ----------------------------------------------------------------------------
# 0. CONFIGURATION — à adapter
# ----------------------------------------------------------------------------

EE_PROJECT = "projet-sigd-01"  # ID de votre projet Google Cloud/Earth Engine
AOI_ASSET_ID = "projects/projet-sigd-01/assets/boeny-boundary"  # contour réel de la région Boeny (importé dans les Assets)
AOI_BBOX = [
    45.5,
    -18.3,
    49.5,
    -15.0,
]  # dépannage : [lon_min, lat_min, lon_max, lat_max]

YEAR_START = 2016
YEAR_END = 2025  # inclus

EXPORT_SCALE = 10  # mètres, résolution native Dynamic World
EXPORT_FOLDER = "boeny_dynamicworld"  # dossier créé dans votre Google Drive

# Fusion des classes Dynamic World :
#   0 water, 1 trees, 2 grass, 3 flooded_vegetation, 4 crops,
#   5 shrub_and_scrub, 6 built, 7 bare, 8 snow_and_ice
# -> nouvelle nomenclature (grass + shrub_and_scrub fusionnées) :
#   0 eau, 1 forêt, 2 végétation non forestière, 3 végétation inondée,
#   4 cultures, 5 bâti, 6 sol nu
REMAP_FROM = [0, 1, 2, 3, 4, 5, 6, 7, 8]
REMAP_TO = [
    0,
    1,
    2,
    3,
    4,
    2,
    5,
    6,
    6,
]  # 8 (neige/glace) rattachée à "sol nu" (cas marginal, non attendu à Boeny)

CLASS_NAMES = {
    0: "eau",
    1: "foret",
    2: "vegetation_non_forestiere",
    3: "vegetation_inondee",
    4: "cultures",
    5: "bati",
    6: "sol_nu",
}


# ----------------------------------------------------------------------------
# 1. INITIALISATION
# ----------------------------------------------------------------------------


def init_earth_engine():
    try:
        ee.Initialize(project=EE_PROJECT)
    except Exception:
        ee.Authenticate()
        ee.Initialize(project=EE_PROJECT)


def get_aoi():
    """Retourne la géométrie de la zone d'étude (Boeny)."""
    if AOI_ASSET_ID:
        return ee.FeatureCollection(AOI_ASSET_ID)
    print(
        "ATTENTION : AOI_ASSET_ID non défini — utilisation de la bounding box "
        "approximative (AOI_BBOX). À remplacer par le contour officiel de "
        "Boeny avant l'analyse finale."
    )
    geom = ee.Geometry.Rectangle(AOI_BBOX)
    return ee.FeatureCollection([ee.Feature(geom)])


# ----------------------------------------------------------------------------
# 2. COMPOSITE ANNUEL DYNAMIC WORLD
# ----------------------------------------------------------------------------


def annual_mode_composite(aoi_geometry, year):
    """Composite annuel 'mode' (classe la plus fréquente par pixel) sur l'AOI."""
    start = ee.Date.fromYMD(year, 1, 1)
    end = start.advance(1, "year")

    dw = (
        ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1")
        .filterBounds(aoi_geometry)
        .filterDate(start, end)
        .select("label")
    )

    mode_composite = dw.reduce(ee.Reducer.mode()).rename("label")

    # Fusion des classes selon REMAP_FROM -> REMAP_TO
    remapped = mode_composite.remap(REMAP_FROM, REMAP_TO).rename("classe")

    return remapped.clip(aoi_geometry).set("year", year)


# ----------------------------------------------------------------------------
# 3. EXPORT
# ----------------------------------------------------------------------------


def export_year(image, aoi_geometry, year):
    task = ee.batch.Export.image.toDrive(
        image=image,
        description=f"boeny_dw_{year}",
        folder=EXPORT_FOLDER,
        fileNamePrefix=f"boeny_occupation_sol_{year}",
        region=aoi_geometry.geometry(),
        scale=EXPORT_SCALE,
        crs="EPSG:32738",  # UTM 38S, adapté à Madagascar/Boeny
        maxPixels=1e13,
    )
    task.start()
    print(
        f"Export lancé pour {year} -> Google Drive/{EXPORT_FOLDER}/boeny_occupation_sol_{year}.tif"
    )
    return task


# ----------------------------------------------------------------------------
# 4. SCRIPT PRINCIPAL
# ----------------------------------------------------------------------------


def main():
    init_earth_engine()
    aoi = get_aoi()

    tasks = []
    for year in range(YEAR_START, YEAR_END + 1):
        composite = annual_mode_composite(aoi, year)
        tasks.append(export_year(composite, aoi, year))

    print(f"\n{len(tasks)} exports lancés (années {YEAR_START}-{YEAR_END}).")
    print(
        "Suivez leur avancement dans l'onglet 'Tasks' du Code Editor GEE, "
        "ou via task.status() en Python."
    )


if __name__ == "__main__":
    main()
