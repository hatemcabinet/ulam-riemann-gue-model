import numpy as np
from scipy.stats import kstest
from scipy.integrate import quad

def gue_cdf(s):
    if s < 0:
        return 0.0
    val, _ = quad(lambda x: (32.0 / (np.pi**2)) * (x**2) * np.exp(-(4.0 / np.pi) * (x**2)), 0, s)
    return val

vectorized_gue_cdf = np.vectorize(gue_cdf)

def compute_von_mangoldt(max_val):
    lam = np.zeros(max_val + 1)
    for n in range(2, max_val + 1):
        temp = n
        p = 2
        factors = []
        while p * p <= temp:
            while temp % p == 0:
                factors.append(p)
                temp //= p
            p += 1
        if temp > 1:
            factors.append(temp)
        if len(set(factors)) == 1:
            lam[n] = np.log(factors[0])
    return lam

def compute_mobius(max_val):
    mu_exact = np.zeros(max_val + 1, dtype=int)
    mu_exact[1] = 1
    for i in range(1, max_val + 1):
        if mu_exact[i] != 0:
            for j in range(2 * i, max_val + 1, i):
                mu_exact[j] -= mu_exact[i]
    return mu_exact

def build_metrics(weights, N=150, max_n=3000):
    n_vals = np.arange(1, len(weights) + 1)
    indices = np.arange(1, N + 1)[:, None]
    scaling = np.log(max_n + 1)

    log_n = np.log(n_vals + 1e-8)
    cos_mat = np.cos(indices * log_n / scaling)

    K = (cos_mat * weights) @ cos_mat.T
    np.fill_diagonal(K, K.diagonal() + (np.arange(1, N + 1) * np.pi)**2 / (2.0 * N))

    H = 0.5 * (K + K.T)
    evals, evecs = np.linalg.eigh(H)
    
    # Spacings for KS test
    poly = np.polyfit(evals, np.arange(1, len(evals) + 1), deg=2)
    unfolded = np.polyval(poly, evals)
    spacings = np.diff(unfolded)
    mean_s = np.mean(spacings)
    if mean_s > 0:
        spacings = spacings / mean_s
        
    # IPR calculation
    iprs = np.sum(evecs**4, axis=0)
    return spacings, np.mean(iprs)

def run_full_audit():
    print("======================================================================")
    print(" AUDIT COMPLET : VON MANGOLDT vs MÖBIUS (KS & IPR vs CONTRÔLE NUL)")
    print("======================================================================")

    N = 150
    max_n = 3000

    lam = compute_von_mangoldt(max_n)
    weights_lam = lam[1:] / np.sqrt(np.arange(1, max_n + 1) + 1e-8)
    sp_lam, ipr_lam = build_metrics(weights_lam, N=N, max_n=max_n)
    ks_lam = kstest(sp_lam, vectorized_gue_cdf).statistic

    mu = compute_mobius(max_n)
    weights_mu = mu[1:] / np.sqrt(np.arange(1, max_n + 1) + 1e-8)
    sp_mu, ipr_mu = build_metrics(weights_mu, N=N, max_n=max_n)
    ks_mu = kstest(sp_mu, vectorized_gue_cdf).statistic

    rng = np.random.default_rng(42)
    ks_ctrls, ipr_ctrls = [], []
    
    for _ in range(30):
        shuffled = rng.permutation(weights_lam)
        sp_c, ipr_c = build_metrics(shuffled, N=N, max_n=max_n)
        ks_ctrls.append(kstest(sp_c, vectorized_gue_cdf).statistic)
        ipr_ctrls.append(ipr_c)

    mean_ks_c, std_ks_c = np.mean(ks_ctrls), np.std(ks_ctrls)
    mean_ipr_c, std_ipr_c = np.mean(ipr_ctrls), np.std(ipr_ctrls)

    z_ks_lam = (ks_lam - mean_ks_c) / std_ks_c
    z_ipr_lam = (ipr_lam - mean_ipr_c) / std_ipr_c
    
    # Pour Möbius vs contrôle
    z_ks_mu = (ks_mu - mean_ks_c) / std_ks_c

    print(f" [Lambda(n)] KS = {ks_lam:.4f} (z = {z_ks_lam:.2f}) | IPR = {ipr_lam:.4f} (z = {z_ipr_lam:.2f})")
    print(f" [Möbius mu(n)] KS = {ks_mu:.4f} (z = {z_ks_mu:.2f}) | IPR = {ipr_mu:.4f}")
    print(f" [Contrôle Nul] KS Moyen = {mean_ks_c:.4f} (±{std_ks_c:.4f}) | IPR Moyen = {mean_ipr_c:.4f} (±{std_ipr_c:.4f})")
    print("======================================================================")

if __name__ == "__main__":
    run_full_audit()
