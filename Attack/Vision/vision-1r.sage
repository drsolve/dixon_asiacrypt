load("../../drsolve/drsolve_sage_interface.sage")


def vision_1r_attack(dixon_path="../../drsolve/drsolve", live_output=True):
    set_dixon_path(dixon_path)

    R0.<u> = GF(2)[]
    K.<z8> = GF(2^8, modulus=u^8 + u^4 + u^3 + u^2 + 1)
    field = K
    R = PolynomialRing(K, names=["x0", "x1", "x2"])
    x0, x1, x2 = R.gens()

    p0 = z8*x0*x1 + (z8^7 + z8^5 + z8^3 + z8)*x1 + 1
    p1 = (z8^2 + z8)*x0*x2 + (z8^6 + z8^5 + z8^4 + z8^3 + z8 + 1)*x2 + 1
    p2 = (z8^6 + z8^5 + z8^3)*x1^4 + (z8^6 + z8^4 + z8^3 + z8^2)*x2^4 + (z8^6 + z8^4 + z8^3 + z8^2)*x1^2 + (z8^6 + z8^5 + z8^4 + z8)*x2^2 + (z8^7 + z8^6 + z8^4)*x1 + (z8^7 + z8^5 + z8^4 + z8^3)*x2 + (z8^5 + z8^4 + z8^3 + 1)

    common = dict(
        field_size=field,
        live_output=live_output,
        verbosity=1,
        time=True,
    )
    subres = dict(common, resultant_method="subres")

    print("=== Vision Attack using Sage interface ===")
    total_start = walltime()
    r1 = DixonRes([p0, p1], [x0], **subres)
    d = DixonRes([r1, p2], [x1], **subres)
    print("Total solve time: %.6f s" % walltime(total_start))
    return d


if __name__ == "__main__":
    vision_1r_attack()
