load("../../drsolve/drsolve_sage_interface.sage")


def vision_2r_attack(dixon_path="../../drsolve/drsolve", live_output=True):
    set_dixon_path(dixon_path)

    R0.<u> = GF(2)[]
    K.<z8> = GF(2^8, modulus=u^8 + u^4 + u^3 + u^2 + 1)
    field = K
    R = PolynomialRing(K, names=["x0", "x1", "x2", "x3", "x4", "x5", "x6"])
    x0, x1, x2, x3, x4, x5, x6 = R.gens()

    p0 = (z8^7 + z8^6 + z8^5 + 1)*x1^4*x3^4 + (z8^4 + z8^3 + z8^2 + z8 + 1)*x2^4*x3^4 + (z8^6 + 1)*x1^4*x3^2 + (z8^7 + z8^6 + z8^5 + z8^3 + z8^2 + z8 + 1)*x2^4*x3^2 + (z8^6 + 1)*x1^2*x3^4 + (z8^7 + z8^6 + z8^5 + z8^3 + z8^2 + z8 + 1)*x2^2*x3^4 + (z8^7 + z8^4 + z8^3 + z8^2 + 1)*x1^4*x3 + (z8^6 + z8^4 + z8^3 + z8^2 + 1)*x2^4*x3 + (z8^7 + z8^4 + z8^3 + z8^2 + 1)*x1*x3^4 + (z8^6 + z8^4 + z8^3 + z8^2 + 1)*x2*x3^4 + (z8^7 + z8^4 + z8^2)*x1^2*x3^2 + (z8^7 + z8^6 + z8^4 + z8^3 + z8^2 + z8)*x2^2*x3^2 + (z8^4 + z8^3 + 1)*x3^4 + (z8^7 + z8^6 + z8^4 + z8^3 + z8)*x1^2*x3 + (z8^7 + z8^5 + z8^4 + z8^2 + z8 + 1)*x2^2*x3 + (z8^7 + z8^6 + z8^4 + z8^3 + z8)*x1*x3^2 + (z8^7 + z8^5 + z8^4 + z8^2 + z8 + 1)*x2*x3^2 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^2)*x1*x3 + (z8^7 + z8^3 + z8^2 + z8)*x2*x3 + (z8^6 + z8^5 + z8^3 + z8)*x3^2 + (z8^6 + z8^3 + z8^2 + z8 + 1)*x3 + 1
    p1 = (z8^5 + z8^4 + z8^3 + z8^2 + z8)*x1^4*x4^4 + (z8^7 + z8^6)*x2^4*x4^4 + (z8^7 + z8^6 + z8 + 1)*x1^4*x4^2 + (z8^6 + z8^5 + z8^3 + z8^2 + 1)*x2^4*x4^2 + (z8^7 + z8^6 + z8 + 1)*x1^2*x4^4 + (z8^6 + z8^5 + z8^3 + z8^2 + 1)*x2^2*x4^4 + (z8^7 + z8^5 + z8^4 + z8^3 + z8)*x1^4*x4 + (z8^6 + z8^5 + z8^4 + z8^3 + z8)*x2^4*x4 + (z8^7 + z8^5 + z8^4 + z8^3 + z8)*x1*x4^4 + (z8^6 + z8^5 + z8^4 + z8^3 + z8)*x2*x4^4 + (z8^7 + z8^5 + 1)*x1^2*x4^2 + (z8^7 + z8^6 + z8^5 + z8^3 + z8 + 1)*x2^2*x4^2 + (z8^6 + z8^5 + z8 + 1)*x4^4 + (z8^6 + z8^5 + z8^4 + z8 + 1)*x1^2*x4 + (z8^4 + z8^3 + z8^2 + z8)*x2^2*x4 + (z8^6 + z8^5 + z8^4 + z8 + 1)*x1*x4^2 + (z8^4 + z8^3 + z8^2 + z8)*x2*x4^2 + x1*x4 + (z8^6 + z8^5 + z8^4 + z8^3 + z8 + 1)*x2*x4 + (z8^6 + z8^5 + z8^4 + z8^2)*x4^2 + (z8^5 + z8^4 + z8^2)*x4 + 1
    p2 = z8*x0*x1 + z8*x1 + 1
    p3 = (z8^2 + z8)*x0*x2 + (z8^7 + z8^5 + z8^3 + z8 + 1)*x2 + 1
    p4 = z8*x3*x5 + (z8 + 1)*x4*x5 + (z8^5 + z8^4 + z8^2)*x5 + 1
    p5 = (z8^2 + z8)*x3*x6 + (z8^2 + z8 + 1)*x4*x6 + (z8^6 + z8^3 + z8^2 + z8)*x6 + 1
    p6 = (z8^7 + z8^6 + 1)*x5^4 + (z8^5 + z8^3 + z8^2 + z8 + 1)*x6^4 + (z8^5 + z8^3 + 1)*x5^2 + (z8^7 + z8^5 + z8^4 + z8 + 1)*x6^2 + (z8^6 + z8^5 + z8^4 + z8^3 + z8^2 + 1)*x5 + (z8^7 + z8^6 + z8^3 + z8^2 + 1)*x6 + (z8^6 + z8^4 + z8)

    common = dict(
        field_size=field,
        live_output=live_output,
        verbosity=1,
        time=True,
    )
    subres = dict(common, resultant_method="subres")
    dixon = dict(common, resultant_method="dixon")

    print("=== Vision Attack using Sage interface  ===")
    total_start = walltime()
    r1 = DixonRes([p2, p3], [x0], **subres)
    r2 = DixonRes([r1, p0, p1], [x1, x2], **dixon)
    s1 = DixonRes([p6, p5], [x6], **subres)
    s2 = DixonRes([s1, p4], [x5], **subres)
    d = DixonRes([r2, s2], [x3], **subres)
    print("Total solve time: %.6f s" % walltime(total_start))
    return d


if __name__ == "__main__":
    vision_2r_attack()
