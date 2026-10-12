import matplotlib.pyplot as plt
import numpy as np
from mpmath import zetazero
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator, eigsh
from scipy.special import erf
from scipy.stats import kstest

def est_premier(n):
    if n < 2: return False
    if n in (2, 3): return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0: return False
        i += 6
    return True

def generer_trous_premiers(max_k=15):
    nombres_couverts = set()
    for k in range(max_k + 1):
        nombres_couverts.add(4*k**2 - 2*k + 41)
        nombres_couverts.add(4*k**2 + 2*k + 41)
        nombres_couverts.add(4*k**2 + 1)
        nombres_couverts.add(4*k**2 + 8*k + 5)
        nombres_couverts.add(4*k**2 - 4*k + 2)
        nombres_couverts.add(4*k**2 + 6*k + 2)

    max_val = (2 * max_k + 1) ** 2
    coords_x, coords_y = [], []
    x, y = 0, 0
    dx, dy = 1, 0
    current_val = 1

    if current_val not in nombres_couverts and est_premier(current_val):
        coords_x.append(x); coords_y.append(y)

    step_len = 1
    while current_val < max_val:
        for _ in range(2):
            for _ in range(step_len):
                if current_val >= max_val: break
                x += dx; y += dy; current_val += 1
                if current_val not in nombres_couverts and est_premier(current_val):
                    coords_x.append(x); coords_y.append(y)
            dx, dy = -dy, dx
            if current_val >= max_val: break
        step_len += 1

    return np.column_stack((coords_x, coords_y))

def calculer_dimension_fractale(points):
    if len(points) < 10: return 0.0
    min_vals = np.min(points, axis=0)
    points_shifted = points - min_vals
    max_range = np.max(points_shifted)
    if max_range == 0: return 0.0
    sizes = np.unique(np.logspace(0, np.log10(max_range / 3), num=8, dtype=int))
    sizes = sizes[sizes > 0]

    epsilons, counts = [], []
    for s in sizes:
        boxed = np.floor(points_shifted / s)
        unique_boxes = len(np.unique(boxed, axis=0))
        epsilons.append(s); counts.append(unique_boxes)

    epsilons = np.array(epsilons, dtype=float)
    counts = np.array(counts, dtype=float)
    coeffs = np.polyfit(np.log(1.0 / epsilons), np.log(counts), 1)
    return coeffs[0]

def analyser_spectre_quantique_rapide(num_zeros=200, N=300, bandwidth=20):
    # Réduction de num_zeros et N pour éviter le gel du processeur
    zeros_t = [float(zetazero(n).imag) for n in range(1, num_zeros + 1)]
    zeros_arr = np.array(zeros_t)
    weights = 1.0 / np.sqrt(zeros_arr)
    indices = np.arange(1, N + 1)[:, np.newaxis]
    sin_mat = np.sin(indices * zeros_arr)
    weighted_sin = sin_mat * weights
    diag_vals = ((np.arange(1, N + 1) * np.pi)**2) / (2.0 * N)
    row_idx, col_idx, data_vals = [], [], []

    for i in range(N):
        row_idx.append(i); col_idx.append(i); data_vals.append(diag_vals[i])
        for j in range(max(0, i - bandwidth), min(N, i + bandwidth + 1)):
            if i != j:
                val = 0.05 * np.cos(i - j) / (1.0 + abs(i - j))
                if abs(val) > 1e-4:
                    row_idx.append(i); col_idx.append(j); data_vals.append(val)

    Q_inf_sparse = sp.csr_matrix((data_vals, (row_idx, col_idx)), shape=(N, N), dtype=complex)
    def matvec(v):
        term_q_prime = sin_mat @ (weighted_sin.T @ v)
        term_q_inf = Q_inf_sparse @ v
        total = term_q_prime + term_q_inf
        return 0.5 * (total + np.conj(total))

    H_op = LinearOperator((N, N), matvec=matvec, dtype=complex)
    k_evals = min(50, N - 2)  # Réduction du nombre de valeurs propres cherchées
    evals, _ = eigsh(H_op, k=k_evals, which='SA', tol=1e-4)
    evals = np.sort(evals)
    indices_eval = np.arange(1, len(evals) + 1)
    poly = np.polyfit(evals, indices_eval, deg=2)
    unfolded = np.polyval(poly, evals)
    spacings = np.diff(unfolded)
    spacings = spacings / np.mean(spacings)

    return evals, spacings

def gue_cdf(s):
    if s < 0: return 0.0
    a = 4.0 / np.pi
    return float(erf(np.sqrt(a) * s) - (4.0 * s / np.pi) * np.exp(-a * s**2))

vectorized_gue_cdf = np.vectorize(gue_cdf)

if __name__ == "__main__":
    print("=" * 70)
    print(" MODÈLE OPTIMISÉ : EXÉCUTION RAPIDE")
    print("=" * 70)

    points_trous = generer_trous_premiers(max_k=15)
    dim_fractale = calculer_dimension_fractale(points_trous)
    print(f" -> Dimension fractale (D_0) des trous premiers : {dim_fractale:.4f}")

    evals, spacings = analyser_spectre_quantique_rapide(num_zeros=200, N=300, bandwidth=20)
    print(f" -> Espacement moyen normalisé : {np.mean(spacings):.4f}")

    ks_stat, p_value = kstest(spacings, vectorized_gue_cdf)
    print(f" -> Test de Kolmogorov-Smirnov : Statistique KS = {ks_stat:.4f}, p-valeur = {p_value:.4e}")

    print("=" * 70)
    print(" [SUCCÈS] Exécution terminée en quelques secondes.")
    print("=" * 70)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].scatter(points_trous[:, 0], points_trous[:, 1], color='indigo', s=12, alpha=0.7)
    axes[0].set_title(f"Distribution Spatiale (D_0: {dim_fractale:.3f})")
    axes[0].grid(True, linestyle=":", alpha=0.6)

    axes[1].hist(spacings, bins=25, density=True, color='crimson', alpha=0.5, edgecolor='black', label="Simulé (N=300)")
    s_vals = np.linspace(0, max(max(spacings), 3.0), 200)
    p_gue = (32.0 / (np.pi**2)) * (s_vals**2) * np.exp(-(4.0 / np.pi) * (s_vals**2))
    axes[1].plot(s_vals, p_gue, color='blue', linewidth=2.5, label="Théorie Wigner-Dyson (GUE)")

    axes[1].set_title(f"Espacement Spectral vs GUE\n(N=300 | KS Stat: {ks_stat:.3f})")
    axes[1].set_xlabel("Espacement normalisé (s)")
    axes[1].set_ylabel("Densité de probabilité P(s)")
    axes[1].grid(True, linestyle=":", alpha=0.6)
    axes[1].legend()

    plt.tight_layout()
    plt.savefig("modele_rapide.png", dpi=300)
    plt.show()
