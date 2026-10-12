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
## 🎯 Bilan et Perspectives Scientifiques

Ce projet démontre qu'il est possible de concevoir un laboratoire de recherche numérique modulaire et rigoureux, explorant l'intersection entre l'arithmétique discrète (spirale d'Ulam), la topologie fractale et la théorie des matrices aléatoires (GUE).

### 📋 Synthèse des Réalisations
- **Architecture Numérique :** Mise en place d'un opérateur de Galerkin creux optimisé via `scipy.sparse` et `LinearOperator`, évitant toute saturation mémoire lors du calcul spectral.
- **Rigueur Expérimentale :** Intégration d'un protocole de **contrôle nul** (comparaison Monte-Carlo avec des fréquences aléatoires) pour dissocier les véritables signatures structurelles des artéfacts numériques.
- **Reproductibilité Totale :** Automatisation de l'export des données sous forme de rapports structurés et de séries d'espacements bruts (`null_control_report.csv`, `spectral_spacings.csv`).

### ⚠️ Avertissement Épistémologique
Bien que l'analyse statistique (tests de Kolmogorov-Smirnov) permette d'évaluer l'adéquation des spectres avec la théorie GUE, ces observations constituent une **approche exploratoire et phénoménologique**. Elles ne sauraient constituer, à elles seules, une démonstration formelle de l'hypothèse de Riemann. Cela souligne l'importance cardinale de maintenir des contrôles nuls stricts dans toute recherche computationalle en physique mathématique.

### 🔮 Poursuivre l'Exploration
Le framework est entièrement ouvert aux contributions et aux extensions futures, notamment :
1. **Dépliage Local :** Remplacer l'ajustement polynomial global par une estimation de densité spectrale locale pour affiner l'analyse des espacements.
2. **Passage à l'Échelle :** Exploiter des tranches de zéros de Riemann à très haute altitude ($t > 10^5$) pour tester la robustesse asymptotique du modèle.
3. **Automatisation Graphique :** Intégrer les scripts de visualisation dans un pipeline d'export automatique d'images.
