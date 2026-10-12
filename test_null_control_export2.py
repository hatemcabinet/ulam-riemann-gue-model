
import csv
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from scipy.stats import kstest
from scipy.special import erf
from mpmath import zetazero

# ============================================================
# PARAMETRES
# ============================================================

NUM_ZEROS = 200
N = 300
BANDWIDTH = 20
NUM_TRIALS = 5
SEED = 42

REPORT_CSV = "null_control_report.csv"
SPACINGS_CSV = "spectral_spacings.csv"


# ============================================================
# CDF : APPROXIMATION DE WIGNER POUR LES ESPACEMENTS GUE
# ============================================================

def gue_cdf(s):
    s = np.asarray(s, dtype=float)
    x = np.maximum(s, 0.0)
    a = 4.0 / np.pi
    result = erf(np.sqrt(a) * x) - (4.0 * x / np.pi) * np.exp(-a * x**2)
    return np.where(s < 0, 0.0, result)


# ============================================================
# CONSTRUCTION DE L'OPERATEUR
# ============================================================

def construire_operateur(zeros_arr, N=300, bandwidth=20):
    zeros_arr = np.asarray(zeros_arr, dtype=float)

    if len(zeros_arr) == 0 or not np.all(np.isfinite(zeros_arr)):
        raise ValueError("Les fréquences doivent être finies et non vides.")

    indices = np.arange(1, N + 1, dtype=float)[:, None]
    weights = 1.0 / np.sqrt(np.abs(zeros_arr) + 1e-5)

    sin_mat = np.sin(indices * zeros_arr[None, :])
    weighted_sin = sin_mat * weights[None, :]

    # Matrice de Gram réelle symétrique :
    # Q_prime = sin_mat @ diag(weights) @ sin_mat.T
    Q_prime = sin_mat @ weighted_sin.T
    Q_prime = 0.5 * (Q_prime + Q_prime.T)

    # Terme complémentaire réel symétrique à bande limitée.
    rows, cols, values = [], [], []

    diag_vals = ((np.arange(1, N + 1) * np.pi) ** 2) / (2.0 * N)

    for i in range(N):
        rows.append(i)
        cols.append(i)
        values.append(diag_vals[i])

        for j in range(max(0, i - bandwidth),
                       min(N, i + bandwidth + 1)):
            if i != j:
                val = 0.05 * np.cos(i - j) / (1.0 + abs(i - j))
                if abs(val) > 1e-4:
                    rows.append(i)
                    cols.append(j)
                    values.append(val)

    Q_inf = sp.csr_matrix(
        (values, (rows, cols)), shape=(N, N), dtype=float
    )

    Q_inf = 0.5 * (Q_inf + Q_inf.T)

    # L'opérateur est explicitement réel et symétrique.
    # Cela évite la conjugaison élément par élément erronée
    # de la version initiale.
    def matvec(v):
        v = np.asarray(v, dtype=float)
        return Q_prime @ v + Q_inf @ v

    from scipy.sparse.linalg import LinearOperator

    H = LinearOperator(
        shape=(N, N),
        matvec=matvec,
        rmatvec=matvec,
        dtype=np.float64
    )

    return H


# ============================================================
# ESPACEMENTS SPECTRAUX
# ============================================================

def calculer_espacements(zeros_arr, N=300, bandwidth=20):
    H = construire_operateur(zeros_arr, N, bandwidth)

    k = min(50, N - 2)

    eigenvalues = eigsh(
        H,
        k=k,
        which="SA",
        tol=1e-6,
        return_eigenvectors=False
    )

    eigenvalues = np.sort(eigenvalues)

    if len(eigenvalues) < 5:
        raise ValueError("Nombre insuffisant de valeurs propres.")

    # Dépliement polynomial conservé pour reproduire l'expérience,
    # mais il devra être comparé à d'autres méthodes.
    ranks = np.arange(1, len(eigenvalues) + 1, dtype=float)
    coeffs = np.polyfit(eigenvalues, ranks, deg=2)
    unfolded = np.polyval(coeffs, eigenvalues)

    spacings = np.diff(unfolded)

    # Les espacements doivent être strictement positifs.
    spacings = spacings[spacings > 0]

    if len(spacings) < 5:
        raise ValueError("Trop peu d'espacements positifs.")

    mean_spacing = np.mean(spacings)
    if not np.isfinite(mean_spacing) or mean_spacing <= 0:
        raise ValueError("Normalisation des espacements impossible.")

    return spacings / mean_spacing


# ============================================================
# EXPERIENCE ET CONTROLES NULS
# ============================================================

def main():
    print("=" * 70)
    print("CONTROLE NUL — ANALYSE SPECTRALE")
    print("=" * 70)

    print("[1] Calcul des zéros de Riemann...")
    zeros_zeta = np.array(
        [float(zetazero(i).imag) for i in range(1, NUM_ZEROS + 1)]
    )

    rng = np.random.default_rng(SEED)

    cases = [("Riemann_Zeta", 1, zeros_zeta)]

    print(f"[2] Génération de {NUM_TRIALS} contrôles aléatoires...")

    for trial in range(1, NUM_TRIALS + 1):
        zeros_rand = rng.uniform(
            zeros_zeta.min(),
            zeros_zeta.max(),
            size=NUM_ZEROS
        )
        zeros_rand.sort()
        cases.append(("Null_Control_Random", trial, zeros_rand))

    reports = []
    spacing_rows = []

    for case_type, trial, frequencies in cases:
        print(f"[*] Analyse : {case_type}, essai {trial}")

        spacings = calculer_espacements(
            frequencies, N=N, bandwidth=BANDWIDTH
        )

        ks_stat, p_value = kstest(spacings, gue_cdf)

        reports.append({
            "Type": case_type,
            "Trial": trial,
            "N_Zeros": len(frequencies),
            "N_Spacings": len(spacings),
            "KS_Statistic": round(float(ks_stat), 6),
            "P_Value": f"{p_value:.6g}",
            "Mean_Spacing": round(float(np.mean(spacings)), 6)
        })

        for i, spacing in enumerate(spacings, start=1):
            spacing_rows.append({
                "Type": case_type,
                "Trial": trial,
                "Spacing_Index": i,
                "Normalized_Spacing": float(spacing)
            })

        print(
            f"    KS={ks_stat:.5f} | "
            f"p={p_value:.5g} | "
            f"espacements={len(spacings)}"
        )

    with open(REPORT_CSV, "w", newline="", encoding="utf-8") as f:
        fields = [
            "Type", "Trial", "N_Zeros", "N_Spacings",
            "KS_Statistic", "P_Value", "Mean_Spacing"
        ]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(reports)

    with open(SPACINGS_CSV, "w", newline="", encoding="utf-8") as f:
        fields = [
            "Type", "Trial", "Spacing_Index", "Normalized_Spacing"
        ]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(spacing_rows)

    print("=" * 70)
    print(f"Rapport : {REPORT_CSV}")
    print(f"Espacements : {SPACINGS_CSV}")
    print("=" * 70)
    print("Attention : les p-values KS sont exploratoires, pas des preuves.")


if __name__ == "__main__":
    main()
