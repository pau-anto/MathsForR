# MathsForR

> Diagonalisation, interpolation polynomiale, SVD et ACP implémentés en R — Projet trimestriel ESGI 4IABD 2026

## Contenu du projet

### 1. Diagonalisation et systèmes linéaires
- Résolution de suites récurrentes par diagonalisation
- Méthode générique pour k suites récurrentes d'ordre 1
- Cas non diagonalisable : décomposition de Jordan

### 2. Interpolation de données
- Génération aléatoire de nuages de points 2D
- Interpolation par matrice de Vandermonde
- Interpolation par matrice de Newton
- Tableau des différences divisées
- Application sur un jeu de données réel avec prévisions

### 3. Décomposition en Valeurs Singulières (SVD) & Pseudo-inverse
- Implémentation from scratch de la SVD
- SVD réduite et calcul de la pseudo-inverse
- Résolution de systèmes linéaires via pseudo-inverse

### 4. Analyse en Composantes Principales (ACP)
- Exploration du jeu de données Europe (25 pays, 9 variables)
- ACP complète : valeurs propres, cercle de corrélation, contributions
- Gestion des individus atypiques (Luxembourg)
- Interprétation des axes principaux

## Technologies

- **Langage :** R
- **IDE :** RStudio ou VSCode
- **Librairies :** `ggplot2`, `FactoMineR`, `factoextra`

## Structure du projet
```
MathsForR/
├── 1_diagonalisation/
│   ├── suites_recursives.R
│   └── jordanisation.R
├── 2_interpolation/
│   ├── vandermonde.R
│   ├── newton.R
│   └── differences_divisees.R
├── 3_svd/
│   ├── svd_implementation.R
│   └── pseudo_inverse.R
├── 4_acp/
│   ├── europe_data.R
│   └── acp_analyse.R
|   └── README.md
├── data/
│   └── europe.csv
├── rapport/
│   └── rapport.pdf
└── README.md
```

## 📝 Consignes

- Tous les résultats sont arrondis à 10⁻² près
- Chaque traitement est accompagné de son script R
- Les résultats sont présentés avec tracés et tableaux de calcul

---

*Projet académique — ESGI — Mathématiques avancées pour le big data avec R — Avril 2026*
