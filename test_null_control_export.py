import csv
import numpy as np
from mpmath import zetazero
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator, eigsh
from scipy.special import erf
from scipy.stats import kstest

def calculer_statistiques_spectrales(zeros_arr, N=300, bandwidth=20):
    weights = 1.0 / np.sqrt(np.abs(zeros_arr) + 1e-5)
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
    k_evals = min(50, N - 2)
    evals, _ = eigsh(H_op, k=k_evals, which='SA', tol=1e-3)
    evals = np.sort(evals)
    
    indices_eval = np.arange(1, len(evals) + 1)
    poly = np.polyfit(evals, indices_eval, deg=2)
    unfolded = np.polyval(poly, evals)
    spacings = np.diff(unfolded)
    spacings = spacings / np.mean(spacings)
    return spacings

def gue_cdf(s):
    if s < 0: return 0.0
    a = 4.0 / np.pi
    return float(erf(np.sqrt(a) * s) - (4.0 * s / np.pi) * np.exp(-a * s**2))

if __name__ == "__main__":
    print("=" * 70)
    print(" AUTOMATISATION DU CONTRÔLE NUL & EXPORT DES RÉSULTATS")
    print("=" * 70)
    
    num_zeros = 200
    N = 300
    bandwidth = 20
    num_trials = 5  # Nombre de répétitions pour le contrôle aléatoire
    
    vectorized_gue_cdf = np.vectorize(gue_cdf)
    csv_filename = "null_control_report.csv"
    
    results = []
    
    # 1. Évaluation avec les zéros de Riemann réels
    print("[*] Analyse du cas réel : Zéros de Riemann...")
    zeros_zeta = np.array([float(zetazero(n).imag) for n in range(1, num_zeros + 1)])
    spacings_zeta = calculer_statistiques_spectrales(zeros_zeta, N, bandwidth)
    ks_zeta, p_zeta = kstest(spacings_zeta, vectorized_gue_cdf)
    
    results.append({
        "Type": "Riemann_Zeta",
        "Trial": 1,
        "KS_Statistic": round(ks_zeta, 5),
        "P_Value": f"{p_zeta:.4e}",
        "Mean_Spacing": round(np.mean(spacings_zeta), 4)
    })
    
    # 2. Évaluation avec les contrôles nuls (fréquences aléatoires sur plusieurs essais)
    print(f"[*] Analyse des contrôles nuls ({num_trials} essais aléatoires)...")
    np.random.seed(42)
    for trial in range(1, num_trials + 1):
        zeros_rand = np.random.uniform(min(zeros_zeta), max(zeros_zeta), size=num_zeros)
        zeros_rand.sort()
        spacings_rand = calculer_statistiques_spectrales(zeros_rand, N, bandwidth)
        ks_rand, p_rand = kstest(spacings_rand, vectorized_gue_cdf)
        
        results.append({
            "Type": "Null_Control_Random",
            "Trial": trial,
            "KS_Statistic": round(ks_rand, 5),
            "P_Value": f"{p_rand:.4e}",
            "Mean_Spacing": round(np.mean(spacings_rand), 4)
        })
        print(f"    - Essai {trial}/{num_trials} terminé (KS = {ks_rand:.4f})")
        
    # 3. Exportation des résultats au format CSV
    print(f"[*] Exportation des données vers {csv_filename}...")
    with open(csv_filename, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["Type", "Trial", "KS_Statistic", "P_Value", "Mean_Spacing"])
        writer.writeheader()
        writer.writerows(results)
        
    print("=" * 70)
    print(f" [SUCCÈS] Rapport généré et enregistré dans '{csv_filename}'.")
    print("=" * 70)
