import os
import csv
import numpy as np
import matplotlib.pyplot as plt

def main():
    print("=" * 70)
    print(" DIAGNOSTIC ET VISUALISATION SPECTRALE")
    print("=" * 70)
    
    filename = "spectral_spacings.csv"
    riemann_spacings = []
    random_spacings = []
    
    # Vérification et chargement du fichier CSV
    if os.path.exists(filename):
        print(f"[OK] Fichier '{filename}' trouvé. Chargement des données...")
        with open(filename, mode="r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    val = float(row["Normalized_Spacing"])
                    if row["Type"] == "Riemann_Zeta":
                        riemann_spacings.append(val)
                    elif row["Type"] == "Null_Control_Random" and int(row["Trial"]) == 1:
                        random_spacings.append(val)
                except (ValueError, KeyError):
                    continue
    else:
        print(f"[ATTENTION] '{filename}' est introuvable.")
        print("-> Génération de données de secours (synthétiques) pour tester le graphique...")
        rng = np.random.default_rng(42)
        riemann_spacings = rng.rayleigh(scale=0.8, size=200)
        random_spacings = rng.normal(loc=1.0, scale=0.3, size=200)

    riemann_spacings = np.array(riemann_spacings)
    random_spacings = np.array(random_spacings)

    # Grille pour la courbe théorique GUE (Wigner-Dyson)
    s_vals = np.linspace(0, 3.5, 300)
    p_gue = (32.0 / (np.pi**2)) * (s_vals**2) * np.exp(-(4.0 / np.pi) * (s_vals**2))

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    # 1. Zéros de Riemann vs GUE
    axes[0].hist(riemann_spacings, bins=25, density=True, color='indigo', alpha=0.6, edgecolor='black', label="Zéros de Riemann")
    axes[0].plot(s_vals, p_gue, color='crimson', linewidth=2.5, label="Théorie GUE")
    axes[0].set_title("Espacements : Zéros de Riemann vs GUE")
    axes[0].set_xlabel("Espacement normalisé (s)")
    axes[0].set_ylabel("Densité de probabilité P(s)")
    axes[0].grid(True, linestyle=":", alpha=0.6)
    axes[0].legend()

    # 2. Contrôle Aléatoire vs GUE
    axes[1].hist(random_spacings, bins=25, density=True, color='forestgreen', alpha=0.6, edgecolor='black', label="Contrôle Aléatoire")
    axes[1].plot(s_vals, p_gue, color='crimson', linewidth=2.5, label="Théorie GUE")
    axes[1].set_title("Espacements : Contrôle Aléatoire vs GUE")
    axes[1].set_xlabel("Espacement normalisé (s)")
    axes[1].grid(True, linestyle=":", alpha=0.6)
    axes[1].legend()

    plt.tight_layout()
    output_fig = "comparaison_gue_spectre.png"
    plt.savefig(output_fig, dpi=300)
    print(f"[SUCCÈS] Graphique sauvegardé sous '{output_fig}'.")
    
    try:
        plt.show()
    except Exception:
        print("[INFO] Mode sans interface graphique détecté (l'image a bien été enregistrée dans le dossier).")

if __name__ == "__main__":
    main()
