# 🌌 Exploration Numérique : Spirale d'Ulam, Topologie Fractale et Chaos Quantique

Ce projet est un laboratoire numérique explorant les ponts entre l'arithmétique discrète, la géométrie fractale et le chaos quantique.

## 📁 Structure du Projet

- **`modele.py`** : Script principal exécutant le pipeline unifié (spirale d'Ulam, dimension $D_0$ et opérateur spectral creux).
- **`test_scaling.py`** : Programme de vérification de la convergence de la dimension fractale $D_0$ selon la taille de la grille ($max\_k$).
- **`test_sensitivity.py`** : Test de robustesse et de sensibilité spectrale face aux variations des paramètres (zéros de Riemann, dimension $N$, bande passante).

## 🚀 Utilisation et Validation

Installez les dépendances :
```bash
pip install numpy matplotlib mpmath scipy
## 🔬 Validation Rigoureuse et Contrôle Nul

Pour s'assurer que les propriétés spectrales observées reflètent une structure spécifique (et non un artéfact de l'opérateur de Galerkin ou de l'unfolding polynomial), le projet intègre un **test de contrôle nul automatique** (`test_null_control_export.py`).

### Protocole expérimental
1. **Cas Réel :** Calcul du spectre et de la statistique de Kolmogorov-Smirnov (KS) à partir des vrais zéros de la fonction Zeta de Riemann.
2. **Contrôle Aléatoire (Monte-Carlo) :** Répétition de l'expérience avec des fréquences aléatoires uniformes de même ordre de grandeur sur plusieurs essais.
3. **Export :** Les métriques sont automatiquement consignées dans un rapport structuré `null_control_report.csv`.

### Exemple de Rapport de Synthèse

| Type d'Entrée | Essai (Trial) | Statistique KS | p-valeur | Espacement Moyen |
| :--- | :---: | :---: | :---: | :---: |
| **Riemann_Zeta** | 1 | 0.0412 | 3.12e-01 | 1.0002 |
| **Null_Control_Random** | 1 | 0.1245 | 1.05e-02 | 0.9987 |
| **Null_Control_Random** | 2 | 0.1189 | 1.42e-02 | 1.0011 |
| **Null_Control_Random** | 3 | 0.1310 | 8.90e-03 | 0.9994 |

*Interprétation :* La divergence observable entre les métriques des zéros de Riemann et celles des contrôles aléatoires permet de quantifier la sensibilité du pipeline face au bruit.
