load("../../drsolve/drsolve_sage_interface.sage")


def xhash_c2s2_m1_3(dixon_path="../../drsolve/drsolve", live_output=True):
    set_dixon_path(dixon_path)

    field = 4611686018427388039
    R = PolynomialRing(ZZ, names=["x0", "x1", "x2", "x3", "x4", "x5"])
    x0, x1, x2, x3, x4, x5 = R.gens()
    ideal = [
        (x2^3, -3*x0^3 + 4*x0^2*x1 + 8*x0^2 - 3*x0*x1^2 - 8*x0*x1 - x0 + 2*x1^3 - x1^2 + 3*x1),
        (x3^3, 2*x0^3 - x0^2*x1 + x0^2 + 2*x0*x1^2 - 2*x0*x1 - 8*x0 - 3*x1^3 - x1^2 - 6*x1 - 2),
        (x4^3, 7*x0^3 - x0^2*x1 - 7*x0^2 + 6*x0*x1^2 - 4*x0*x1 + 8*x0 + 2*x1^3 - 6*x1^2 - 7*x1 - 1),
        (x5^3, 2*x0^3 + 4*x0^2*x1 + x0*x1^2 + 7*x0*x1 - 2*x0 + 7*x1^3 - 6*x1^2 + 2*x1 + 1),
    ]

    f0 = -3*x0^3 + 4*x0^2*x1 + 8*x0^2 - 3*x0*x1^2 - 8*x0*x1 - x0 + 2*x1^3 - x1^2 + 3*x1 - x2^3
    f1 = 2*x0^3 - x0^2*x1 + x0^2 + 2*x0*x1^2 - 2*x0*x1 - 8*x0 - 3*x1^3 - x1^2 - 6*x1 - x3^3 - 2
    f2 = 7*x0^3 - x0^2*x1 - 7*x0^2 + 6*x0*x1^2 - 4*x0*x1 + 8*x0 + 2*x1^3 - 6*x1^2 - 7*x1 - x4^3 - 1
    f3 = 2*x0^3 + 4*x0^2*x1 + x0*x1^2 + 7*x0*x1 - 2*x0 + 7*x1^3 - 6*x1^2 + 2*x1 - x5^3 + 1
    f4 = 4*x2^3 + 3*x2^2*x3 + 4*x2^2 + 6*x2*x3^2 + 2*x2*x3 - 6*x2 + 7*x3^3 - 6*x3^2 + 5*x3 + 8*x4 + 6*x5 - 4
    f5 = 4*x2^3 - 3*x2^2*x3 + x2^2 - 8*x2*x3^2 + 3*x2*x3 - 2*x2 - 3*x3^3 + 3*x3^2 - 3*x3 - 8*x4 + 8*x5 - 1

    common = dict(
        field_size=field,
        live_output=live_output,
        verbosity=1,
    )
    subres = dict(common, resultant_method="subres")

    print("=== xhash elimination using Sage interface ===")
    total_start = walltime()
    r1 = DixonIdeal([f5, f3], ideal, [x5], **common)
    r2 = DixonIdeal([r1, f2], ideal, [x4], **common)
    r3 = DixonIdeal([r2, f1], ideal, [x3], **common)
    r4 = DixonIdeal([r3, f0], ideal, [x2], **common)
    l1 = DixonIdeal([f4, f3], ideal, [x5], **common)
    l2 = DixonIdeal([l1, f2], ideal, [x4], **common)
    l3 = DixonIdeal([l2, f1], ideal, [x3], **common)
    l4 = DixonIdeal([l3, f0], ideal, [x2], **common)
    r = DixonRes([r4, l4], [x1], **subres)
    print("Total solve time: %.6f s" % walltime(total_start))
    return r


if __name__ == "__main__":
    xhash_c2s2_m1_3()
