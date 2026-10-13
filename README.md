# Ulam, von Mangoldt & GUE Spectral Models

Ce dépôt contient le cadre de recherche formel, les scripts de calcul et les analyses statistiques visant à séparer rigoureusement les illusions géométriques (spirale d'Ulam) des véritables structures arithmétiques et spectrales.

## 🏛️ Résumé des Travaux
1. **Démystification de la Spirale d'Ulam :** Démonstration que la forte densité de nombres premiers du polynôme d'Euler (C = 41) relève d'un crible algébrique pur dicté par son discriminant ($\Delta = -163$), la spirale n'agissant que comme un projecteur visuel de congruences.
2. **Opérateur Logarithmique de von Mangoldt :** Construction d'un opérateur hermitien de Galerkin fondé sur l'échelle logarithmique et la fonction $\Lambda(n)$.
3. **Validation par Contrôle Nul et Spécificité :** 
   - Comparaison par Monte-Carlo (essais de permutations aléatoires des poids).
   - Test de Kolmogorov-Smirnov (KS) : Écart de $z = -4.49$ à $-5.24$.
   - Analyse de l'Inverse Participation Ratio (IPR) : Écart de $z = -7.74$, prouvant une structure de phase non aléatoire.
   - Audit de spécificité avec la fonction de Möbius $\mu(n)$ ($z = +0.60$, confirmant l'absence de faux positifs).

## 📂 Structure du Dépôt
- `arithmetic_spectral_audit.py` : Script complet d'audit (von Mangoldt, Möbius, IPR, KS, Monte-Carlo).
- `generate_pdf.py` : Script de génération du mémoire de synthèse PDF officiel.
- `memoire_synthese_ulam_von_mangoldt_v2.pdf` : Le mémoire de synthèse scientifique complet et mis en page.
- `comparison_chart.png` : Graphique comparatif de l'audit de spécificité.

## 🚀 Utilisation
Exécutez l'audit complet avec Python :
```bash
python arithmetic_spectral_audit.py
