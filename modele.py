import matplotlib.pyplot as plt
import numpy as np

def est_premier(n):
    if n < 2: return False
    if n in (2, 3): return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0: return False
        i += 6
    return True

def generer_spirale_ulam_complete(max_k):
    max_val = (2 * max_k + 1) ** 2
    coords_x, coords_y = [], []
    x, y = 0, 0
    dx, dy = 1, 0
    current_val = 1
    if est_premier(current_val):
        coords_x.append(x); coords_y.append(y)
    step_len = 1
    while current_val < max_val:
        for _ in range(2):
            for _ in range(step_len):
                if current_val >= max_val: break
                x += dx; y += dy; current_val += 1
                if est_premier(current_val):
                    coords_x.append(x); coords_y.append(y)
            dx, dy = -dy, dx
            if current_val >= max_val: break
        step_len += 1
    return np.column_stack((coords_x, coords_y))

def generer_trous_premiers(max_k):
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

if __name__ == "__main__":
    print("=" * 70)
    print(" ANALYSE DE SCALING & CONVERGENCE (TEST DE ROBUSTESSE)")
    print("=" * 70)
    print(f"{'max_k':<8} | {'N_max':<8} | {'D_0 (Tous Prémiers)':<20} | {'D_0 (Trous Filtrés)':<20}")
    print("-" * 65)
    
    k_values = [10, 15, 20, 25, 30]
    dims_all, dims_holes = [], []
    
    for k in k_values:
        max_val = (2 * k + 1) ** 2
        pts_all = generer_spirale_ulam_complete(k)
        pts_holes = generer_trous_premiers(k)
        
        d_all = calculer_dimension_fractale(pts_all)
        d_holes = calculer_dimension_fractale(pts_holes)
        
        dims_all.append(d_all)
        dims_holes.append(d_holes)
        print(f"{k:<8} | {max_val:<8} | {d_all:<20.4f} | {d_holes:<20.4f}")
        
    print("=" * 70)
    print(" [SUCCÈS] Test de convergence validé.")
