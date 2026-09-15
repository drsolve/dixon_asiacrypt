"""Shared Dixon Step 1/4 complexity models used by Figures 6(d), 7, and 8.

The functions return base-2 logarithms of operation-count upper bounds.  This
module is intentionally output-free so notebooks and the vector exporter can
reuse one authoritative implementation.
"""

from math import comb, inf, isfinite, lgamma, log, log2

import matplotlib.pyplot as plt

DEFAULT_OMEGA = 2.37

def log2_add_exp(a, b):
    if not isfinite(a):
        return b
    if not isfinite(b):
        return a
    if a < b:
        a, b = b, a
    return a + log2(1.0 + 2.0 ** (b - a))

def log2_factorial(n):
    if n <= 1:
        return 0.0
    return lgamma(n + 1.0) / log(2.0)

def log2_soft_fft_multiply_from_degree_log2(log2_degree_bound):
    if not (log2_degree_bound > 0.0) or not isfinite(log2_degree_bound):
        return 0.0
    result = log2_degree_bound
    if log2_degree_bound > 1.0:
        result += log2(log2_degree_bound)
    if log2_degree_bound > 4.0:
        result += log2(log2(log2_degree_bound))
    return result

def log2_dense_monomial_count_upper(degree, num_vars):
    if degree <= 0 or num_vars <= 0:
        return 0.0
    return log2(comb(num_vars + degree, num_vars))

def log2_binomial_upper(n, k):
    if k < 0 or k > n:
        return float("-inf")
    return (lgamma(n + 1.0) - lgamma(k + 1.0) - lgamma((n - k) + 1.0)) / log(2.0)

def dixon_size(degrees):
    n0 = len(degrees)
    if n0 == 0:
        return 0

    slopes = sorted([max(d - 1, 0) for d in degrees], reverse=True)

    a = [0] * n0
    current_h = 0
    for i in range(n0 - 1):
        current_h += slopes[i]
        a[i + 1] = current_h

    D = [0] * (n0 + 1)
    D[0] = 1

    for k in range(1, n0 + 1):
        total = 0
        for i in range(1, k + 1):
            row_idx = i - 1
            col_idx = k - 1
            upper = a[row_idx] + 1
            lower = col_idx - row_idx + 1
            m_val = comb(upper, lower) if 0 <= lower <= upper else 0
            term = m_val * D[i - 1]
            if (k - i) % 2 == 1:
                term = -term
            total += term
        D[k] = total

    return D[n0]

def dixon_complexity_step4(degrees, ambient_n, omega=DEFAULT_OMEGA):
    size = dixon_size(degrees)
    if size <= 0:
        return 0.0

    m = len(degrees)
    total_deg = sum(degrees)

    if m == ambient_n + 1:
        d_factor = 1.0
    elif m == ambient_n:
        d_factor = float(total_deg)
    elif m < ambient_n:
        d_factor = float((total_deg + 1) ** (ambient_n - m + 1))
    else:
        return 0.0

    if d_factor <= 0.0:
        return 0.0
    return log2(d_factor) + omega * log2(float(size))

def build_sorted_degree_prefix(degrees):
    sorted_deg = sorted([max(d, 0) for d in degrees], reverse=True)
    prefix = [0]
    for d in sorted_deg:
        prefix.append(prefix[-1] + d)
    return sorted_deg, prefix

def prefix_degree_bound(prefix, count):
    if count <= 0:
        return 0
    return max(prefix[count] - count, 0)

def log2_kronecker_degree_from_bounds(bounds):
    log2_weight = 0.0
    log2_degree = float("-inf")
    for bound in bounds:
        if bound > 0:
            log2_degree = log2_add_exp(log2_degree, log2(float(bound)) + log2_weight)
        log2_weight += log2(float(bound) + 1.0)
    return 0.0 if (not isfinite(log2_degree) or log2_degree < 0.0) else log2_degree

def step1_row_highest_active_var_index(row_idx, num_elim_vars, num_parameter_vars):
    total_vars = 2 * num_elim_vars + num_parameter_vars
    if total_vars <= 0:
        return -1
    if num_parameter_vars > 0:
        return total_vars - 1
    if row_idx <= 0:
        return num_elim_vars - 1 if num_elim_vars > 0 else -1
    if row_idx <= num_elim_vars:
        return num_elim_vars + row_idx - 1
    return -1

def step1_entry_degree_log2(row_idx, col_degree, var_weight_log2, num_elim_vars, num_parameter_vars):
    highest_idx = step1_row_highest_active_var_index(row_idx, num_elim_vars, num_parameter_vars)
    effective_degree = col_degree
    if row_idx > 0:
        effective_degree -= 1
    if highest_idx < 0 or effective_degree <= 0:
        return float("-inf")
    return log2(float(effective_degree)) + var_weight_log2[highest_idx]

def recursive_degree_surrogate(degrees, num_elim):
    if num_elim <= 0 or len(degrees) < num_elim:
        return None
    selected = sorted([max(d, 1) for d in degrees])[:num_elim]
    if not selected:
        return None
    return list(reversed(selected))

def equal_degree_report(m, d, s=1, omega=DEFAULT_OMEGA, q=257, char=None):
    if char is None:
        char = q

    num_polys = m + 1
    num_all_vars = m + s
    degrees = [d] * num_polys
    total_degree_sum = sum(degrees)

    det_size = dixon_size(degrees)
    det_size_log2 = log2(float(det_size)) if det_size > 0 else 0.0

    step1_det_total_degree = total_degree_sum - m
    step1_var_count = 2 * m + s

    _, prefix = build_sorted_degree_prefix(degrees)
    var_degree_bounds = [0] * step1_var_count
    step1_partial_degree_bound = 0

    for i in range(m):
        count_orig = i + 1
        count_dual = m - i
        orig_bound = prefix_degree_bound(prefix, count_orig)
        dual_bound = prefix_degree_bound(prefix, count_dual)
        var_degree_bounds[i] = orig_bound
        var_degree_bounds[m + i] = dual_bound
        step1_partial_degree_bound = max(step1_partial_degree_bound, orig_bound, dual_bound)

    for i in range(s):
        param_bound = max(step1_det_total_degree, 0)
        var_degree_bounds[2 * m + i] = param_bound
        step1_partial_degree_bound = max(step1_partial_degree_bound, param_bound)

    var_weight_log2 = [0.0] * step1_var_count
    log2_grid_points = 0.0
    log2_tensor_sum = float("-inf")
    for i, bound in enumerate(var_degree_bounds):
        var_weight_log2[i] = log2_grid_points
        log2_grid_points += log2(float(bound) + 1.0)
        interp_part = (
            log2_soft_fft_multiply_from_degree_log2(log2(float(bound)))
            if bound > 0 else 0.0
        )
        log2_tensor_sum = log2_add_exp(log2_tensor_sum, interp_part)
    if not isfinite(log2_tensor_sum) or log2_tensor_sum < 0.0:
        log2_tensor_sum = 0.0

    log2_kronecker = log2_kronecker_degree_from_bounds(var_degree_bounds)
    log2_fft = log2_soft_fft_multiply_from_degree_log2(log2_kronecker)

    direct_naive_factorial = log2_factorial(num_polys)
    direct_naive = direct_naive_factorial + log2_fft
    direct_dense_la = omega * log2(float(num_polys)) if num_polys > 1 else 0.0
    direct_dense = direct_dense_la + log2_fft
    direct_best = min(direct_naive, direct_dense)
    direct_entry_monomials = log2_dense_monomial_count_upper(d, num_all_vars)
    direct_m2_log2 = 2.0 * det_size_log2 if det_size > 0 else 0.0
    direct_mpoly_mul_proxy = direct_m2_log2 + direct_entry_monomials
    direct_mpoly = direct_naive_factorial + direct_mpoly_mul_proxy
    direct_mpoly_split = (
        (log2(float(num_polys)) + float(num_polys)) if num_polys > 0 else 0.0
    ) + direct_mpoly_mul_proxy
    bareiss = (3.0 * log2(float(num_polys)) if num_polys > 1 else 0.0) + 4.0 * det_size_log2

    log2_col_sum = float("-inf")
    for col in range(num_polys):
        col_degree = degrees[col]
        col_log2_degree = float("-inf")
        for row in range(num_polys):
            entry_log2 = step1_entry_degree_log2(row, col_degree, var_weight_log2, m, s)
            if not isfinite(col_log2_degree) or entry_log2 > col_log2_degree:
                col_log2_degree = entry_log2
        log2_col_sum = log2_add_exp(log2_col_sum, col_log2_degree)

    log2_row_sum = float("-inf")
    for row in range(num_polys):
        row_log2_degree = float("-inf")
        for col in range(num_polys):
            entry_log2 = step1_entry_degree_log2(row, degrees[col], var_weight_log2, m, s)
            if not isfinite(row_log2_degree) or entry_log2 > row_log2_degree:
                row_log2_degree = entry_log2
        log2_row_sum = log2_add_exp(log2_row_sum, row_log2_degree)

    log2_col_avg = log2_col_sum - log2(float(num_polys)) if num_polys > 0 else 0.0
    log2_row_avg = log2_row_sum - log2(float(num_polys)) if num_polys > 0 else 0.0
    log2_s = min(log2_row_avg, log2_col_avg)
    if not isfinite(log2_s) or log2_s < 0.0:
        log2_s = 0.0
    hnf = direct_dense_la + log2_s

    B = total_degree_sum
    log2_small_det = direct_dense_la
    log2_shared_entry_eval = log2(float(num_polys)) + log2(float(B) + 1.0) if num_polys > 0 else 0.0
    log2_matrix_assembly = 2.0 * log2(float(num_polys)) if num_polys > 1 else 0.0
    log2_entry_eval = log2_add_exp(log2_shared_entry_eval, log2_matrix_assembly)
    log2_L = log2_add_exp(log2_small_det, log2_entry_eval)

    ordinary_probe = log2_grid_points + log2_L
    ordinary_tensor = log2_grid_points + log2_tensor_sum
    ordinary = log2_add_exp(ordinary_probe, ordinary_tensor)

    sparse_T = 2.0 * det_size_log2 + s * log2(float(B) + 1.0) if det_size > 0 else 0.0
    log2_q = log2(float(q)) if q > 1 else 0.0
    log2_logq = log2(log2_q) if log2_q > 0.0 else 0.0
    sparse = log2_add_exp(sparse_T + log2_L + log2_logq, sparse_T + 2.0 * log2_logq)

    recursive_seq = recursive_degree_surrogate(degrees, m)
    if recursive_seq is None:
        step12_recursive = inf
        step12_standard = inf
        recursive_m1 = None
    else:
        recursive_m1 = recursive_seq[0]
        step12_recursive = (
            2.0 * log2(float(recursive_m1))
            + 3.0 * log2_factorial(m)
            + 3.0 * sum(log2(float(di)) for di in recursive_seq[1:])
        )
        step12_standard = log2_factorial(m) + 4.0 * sum(log2(float(di)) for di in recursive_seq)
    step1_standard = step12_standard

    theorem_char_ok = char >= step1_partial_degree_bound
    theorem_expected = sparse + log2(4.0 / 3.0) if theorem_char_ok else inf

    if q <= 1 or sparse_T <= 0.0 or step1_partial_degree_bound <= 0:
        base_prob = 1.0
        base_expected = sparse
    else:
        collision = (
            float(step1_partial_degree_bound)
            * (2.0 ** sparse_T)
            * ((2.0 ** sparse_T) - 1.0)
        ) / (2.0 * (q - 1.0))
        base_prob = max(0.0, 1.0 - collision)
        base_expected = sparse - log2(base_prob) if base_prob > 0.0 else inf

    bezout = 1
    for deg in degrees:
        bezout *= max(deg, 0)

    step4_matrix_la = omega * det_size_log2 if det_size > 1 else 0.0

    if s <= 0:
        step4_hnf_degree_density = 0.0
        step4_hnf = step4_matrix_la
        step4_ordinary_grid_points = 0.0
        step4_ordinary_tensor_sum = 0.0
        step4_ordinary_probe = step4_matrix_la
        step4_ordinary_tensor = 0.0
        step4_ordinary = step4_matrix_la
        step4_sparse_T = 0.0
        step4_sparse = step4_matrix_la
    else:
        step4_hnf_degree_density = log2_kronecker_degree_from_bounds([B] * s)
        step4_hnf = step4_matrix_la + step4_hnf_degree_density

        step4_ordinary_grid_points = s * log2(bezout + 1)
        step4_ordinary_tensor_sum = (
            log2(float(s)) + log2_soft_fft_multiply_from_degree_log2(log2(bezout))
            if bezout > 0 else 0.0
        )
        step4_ordinary_probe = step4_ordinary_grid_points + step4_matrix_la
        step4_ordinary_tensor = step4_ordinary_grid_points + step4_ordinary_tensor_sum
        step4_ordinary = log2_add_exp(step4_ordinary_probe, step4_ordinary_tensor)

        step4_sparse_T = log2(comb(bezout + s, s))
        step4_sparse = log2_add_exp(
            step4_sparse_T + step4_matrix_la + log2_logq,
            step4_sparse_T + 2.0 * log2_logq,
        )

    step1_candidates = {
        "direct Kronecker": direct_best,
        "direct multivariate": direct_mpoly,
        "Bareiss": bareiss,
        "ordinary interpolation": ordinary,
        "Kronecker + HNF": hnf,
        "sparse interpolation": sparse,
        "standard Dixon construction": step1_standard,
        "recursive block Dixon construction": step12_recursive,
    }
    step1_best_method, step1_best = min(step1_candidates.items(), key=lambda item: item[1])

    step4_candidates = {
        "HNF" if s == 1 else "Kronecker + HNF": step4_hnf,
        "ordinary interpolation": step4_ordinary,
        "sparse interpolation": step4_sparse,
    }
    step4_best_method, step4 = min(step4_candidates.items(), key=lambda item: item[1])
    overall = max(step1_best, step4)

    gb_n = num_all_vars if num_all_vars > 0 else m
    if gb_n < 1:
        gb_n = 1
    grobner_dreg = gb_n * max(d - 1, 0) + 1
    grobner = omega * log2_binomial_upper(gb_n + grobner_dreg, gb_n)
    if not isfinite(grobner) or grobner < 0.0:
        grobner = 0.0

    if s == 1:
        fglm = (log2(float(num_all_vars)) if num_all_vars > 0 else 0.0) + omega * log2(float(bezout))
        if not isfinite(fglm) or fglm < 0.0:
            fglm = 0.0
    else:
        fglm = None

    total_direct = log2_add_exp(direct_best, step4)
    total_direct_mpoly = log2_add_exp(direct_mpoly, step4)
    total_bareiss = log2_add_exp(bareiss, step4)
    total_ordinary = log2_add_exp(ordinary, step4)
    total_hnf = log2_add_exp(hnf, step4)
    total_sparse = log2_add_exp(sparse, step4)
    total_standard = log2_add_exp(step1_standard, step4)
    total_recursive = log2_add_exp(step12_recursive, step4)

    return {
        "m": m,
        "s": s,
        "n0": num_polys,
        "n_all": num_all_vars,
        "d": d,
        "omega": omega,
        "M": det_size,
        "log2_M": det_size_log2,
        "B": B,
        "D": step1_partial_degree_bound,
        "log2_kronecker": log2_kronecker,
        "log2_grid_points": log2_grid_points,
        "log2_tensor_sum": log2_tensor_sum,
        "log2_L": log2_L,
        "step1_direct": direct_best,
        "step1_direct_naive": direct_naive,
        "step1_direct_dense": direct_dense,
        "step1_direct_mpoly": direct_mpoly,
        "step1_direct_mpoly_split": direct_mpoly_split,
        "step1_bareiss": bareiss,
        "step1_direct_mpoly_mul_proxy": direct_mpoly_mul_proxy,
        "step1_ordinary": ordinary,
        "step1_ordinary_probe": ordinary_probe,
        "step1_ordinary_tensor": ordinary_tensor,
        "step1_hnf": hnf,
        "step1_sparse": sparse,
        "step1_standard": step1_standard,
        "step12_recursive": step12_recursive,
        "step12_standard": step12_standard,
        "step12_recursive_m1": recursive_m1,
        "step12_recursive_seq": ",".join(str(x) for x in recursive_seq) if recursive_seq else None,
        "step1_best_method": step1_best_method,
        "step1_best": step1_best,
        "step1_sparse_theorem_expected": theorem_expected,
        "step1_sparse_base_expected": base_expected,
        "step4": step4,
        "step4_best_method": step4_best_method,
        "step4_hnf": step4_hnf,
        "step4_hnf_degree_density": step4_hnf_degree_density,
        "step4_ordinary": step4_ordinary,
        "step4_ordinary_grid_points": step4_ordinary_grid_points,
        "step4_ordinary_tensor_sum": step4_ordinary_tensor_sum,
        "step4_ordinary_probe": step4_ordinary_probe,
        "step4_ordinary_tensor": step4_ordinary_tensor,
        "step4_sparse": step4_sparse,
        "step4_sparse_T": step4_sparse_T,
        "overall": overall,
        "grobner_dreg": grobner_dreg,
        "grobner": grobner,
        "fglm": fglm,
        "total_direct": total_direct,
        "total_direct_mpoly": total_direct_mpoly,
        "total_bareiss": total_bareiss,
        "total_direct_mpoly_split": log2_add_exp(direct_mpoly_split, step4),
        "total_ordinary": total_ordinary,
        "total_hnf": total_hnf,
        "total_sparse": total_sparse,
        "total_standard": total_standard,
        "total_recursive": total_recursive,
        "theorem_char_ok": theorem_char_ok,
        "base_prob": base_prob,
    }

def degree_sweep_data(m, d_min, d_max, s=1, omega=DEFAULT_OMEGA, q=257, char=None):
    return [
        equal_degree_report(m=m, d=d, s=s, omega=omega, q=q, char=char)
        for d in range(d_min, d_max + 1)
    ]

def variable_sweep_data(d, m_min, m_max, s=1, omega=DEFAULT_OMEGA, q=257, char=None):
    return [
        equal_degree_report(m=m, d=d, s=s, omega=omega, q=q, char=char)
        for m in range(m_min, m_max + 1)
    ]

def _make_plot(rows, x_key, y_specs, title, xlabel, ylabel):
    fig, ax = plt.subplots(figsize=(10, 6))
    for key, label, style in y_specs:
        filtered = [row for row in rows if row.get(key) is not None]
        if not filtered:
            continue
        xs = [row[x_key] for row in filtered]
        ys = [row[key] for row in filtered]
        ax.plot(xs, ys, style, linewidth=2, label=label)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    return fig, ax

def plot_step1_step4_vs_degree(rows, show=True):
    step4_hnf_label = "Step 4 HNF" if rows and all(row.get("s") == 1 for row in rows) else "Step 4 Kronecker+HNF"
    omega = rows[0].get("omega") if rows else None
    omega_txt = f", omega={omega:.2f}" if omega is not None else ""
    fig, ax = _make_plot(
        rows,
        "d",
        [
            #("step1_direct_mpoly", "Step 1 direct expansion", "-o"),
            ("step1_direct_mpoly_split", "Step 1 cached expansion", "--o"),
            #("step1_bareiss", "Step 1 Bareiss", "-."),
            #("step1_hnf", "Step 1 HNF", "-^"),
            #("step1_ordinary", "Step 1 ordinary interpolation", "-s"),
            ("step1_sparse", "Step 1 sparse interpolation", "-d"),
            #("step1_standard", "Step 1 standard Dixon", "--h"),
            ("step12_recursive", "Step 1/2 FDixon", "-p"),
            ("step4_hnf", step4_hnf_label, "-x"),
            ("step4_ordinary", "Step 4 ordinary interpolation", "-+"),
            ("step4_sparse", "Step 4 sparse interpolation", "-*"),
            ("grobner", "Groebner basis", "-v"),
            ("fglm", "FGLM", "-*"),
        ],
        f"Equal-degree Dixon model: complexity vs degree (n={rows[0].get('m') + rows[0].get('s') if rows else '?'}{omega_txt})",
        "common degree d",
        "log2 complexity estimate",
    )
    if show:
        plt.show()
    return fig, ax

def plot_step1_step4_vs_variables(rows, show=True):
    step4_hnf_label = "Step 4 HNF" if rows and all(row.get("s") == 1 for row in rows) else "Step 4 Kronecker+HNF"
    omega = rows[0].get("omega") if rows else None
    omega_txt = f", omega={omega:.2f}" if omega is not None else ""
    fig, ax = _make_plot(
        rows,
        "n_all",
        [
            #("step1_direct_mpoly", "Step 1 direct expansion", "-o"),
            ("step1_direct_mpoly_split", "Step 1 cached expansion", "--o"),
            #("step1_bareiss", "Step 1 Bareiss", "-."),
            #("step1_hnf", "Step 1 HNF", "-^"),
            #("step1_ordinary", "Step 1 ordinary interpolation", "-s"),
            ("step1_sparse", "Step 1 sparse interpolation", "-d"),
            #("step1_standard", "Step 1 standard Dixon", "--h"),
            ("step12_recursive", "Step 1/2 FDixon", "-p"),
            ("step4_hnf", step4_hnf_label, "-x"),
            ("step4_ordinary", "Step 4 ordinary interpolation", "-+"),
            ("step4_sparse", "Step 4 sparse interpolation", "-*"),
            ("grobner", "Groebner basis", "-v"),
            ("fglm", "FGLM", "-*"),
        ],
        f"Equal-degree Dixon model: complexity vs variable count (d={rows[0].get('d') if rows else '?'}{omega_txt})",
        "variable count n",
        "log2 complexity estimate",
    )
    if show:
        plt.show()
    return fig, ax
