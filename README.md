# Itinéraire Convoi

Calcule un itinéraire qui reste sur le réseau convoi exceptionnel d'une catégorie
(Cat 1 / ATE, 48 t, 72 t, 94 t, 120 t+), puis l'ouvre dans Google Maps ou l'exporte
en GPX, KML et CSV (Google My Maps).

Données réseau : carte de Sven — https://sven85000-wq.github.io/reseau-convoi2/
Fond de carte : © contributeurs OpenStreetMap.

## Fichiers

| Fichier | Rôle |
|---|---|
| `index.html` | L'application (une seule page) |
| `graph.bin` | Réseau routier compacté en graphe (gzip), généré par `tools/build_graph.py` |
| `communes.txt` | Communes françaises avec code postal et coordonnées |
| `manifest.webmanifest`, `icon-*.png` | Installation sur l'écran d'accueil |
| `sw.js` | Fonctionnement sans réseau (garde l'appli et les données en cache) |

## Mettre à jour le réseau

Quand Sven publie un nouveau fichier :

```sh
python3 tools/build_graph.py "reseau_convois_france_2026(2).json" graph.bin
```

Puis incrémenter `VERSION` dans `sw.js` pour que les téléphones rechargent les données.

## Limites

Le calcul ignore les sens uniques, les hauteurs sous ouvrage, les ponts et les travaux.
Il ne remplace ni l'arrêté d'autorisation ni la reconnaissance d'itinéraire.
Les « raccords » affichés en rouge relient deux tronçons qui ne se touchent pas dans
les données (moins de 1,5 km) : à vérifier.
