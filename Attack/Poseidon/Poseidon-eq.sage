from dataclasses import dataclass
from pathlib import Path
from time import perf_counter, process_time

from sage.all import *

SCRIPT_DIR = Path(globals().get("__file__", "attack/Poseidon.sage")).resolve().parent
load(str(SCRIPT_DIR / "poseidon.sage"))


@dataclass
class PoseidonCICOConfig:
    prime: int
    t: int
    capacity: int
    digest: int
    R_F: int
    R_P: int
    alpha: int = 3
    target_digest: object = None
    result_prefix: str = "./experiments/poseidon"

    @property
    def rate(self):
        return self.t - self.capacity

    @property
    def linearizable_rounds(self):
        return self.rate - self.digest

    def validate(self):
        assert self.capacity >= 1, "Capacity must be at least 1."
        assert self.rate >= self.digest, (
            f"Rate ({self.rate}) must be >= digest ({self.digest})."
        )
        assert gcd(self.alpha, self.prime - 1) == 1, (
            f"Exponent {self.alpha} does not define a permutation over GF({self.prime})."
        )
        if self.target_digest is not None:
            assert len(self.target_digest) == self.digest, (
                f"Target digest length ({len(self.target_digest)}) must match digest ({self.digest})."
            )

    def to_summary(self):
        return {
            "prime": self.prime,
            "state_size": self.t,
            "capacity": self.capacity,
            "rate": self.rate,
            "digest": self.digest,
            "linearizable_rounds": self.linearizable_rounds,
            "R_F": self.R_F,
            "R_P": self.R_P,
            "total_rounds": self.R_F + self.R_P,
            "alpha": self.alpha,
            "target_digest": self.target_digest,
        }


class PoseidonEquationGenerator:
    """
    Generate the direct input/output CICO system for Poseidon.

    This file follows the modeling choice described in `dixon_appendix.tex`:
    we only keep the unknown digest inputs as variables and treat the remaining
    rate coordinates and the capacity coordinates as fixed zeros.
    """

    def __init__(self, config):
        self.config = config
        self.config.validate()

    def build_system(self):
        if self.config.target_digest is None:
            self.config.target_digest = self._build_demo_target_digest()

        field = GF(self.config.prime)
        ring = PolynomialRing(field, "x", self.config.digest)
        input_vars = list(ring.gens())
        padded_input = input_vars + [field(0)] * self.config.linearizable_rounds

        poseidon = Poseidon(
            prime=self.config.prime,
            R_F=self.config.R_F,
            R_P=self.config.R_P,
            t=self.config.t,
            alpha=self.config.alpha,
        )
        symbolic_output = self._compute_symbolic_output(poseidon, padded_input)

        system = []
        for i in range(self.config.digest):
            system.append(symbolic_output[i] - self.config.target_digest[i])
        return system

    def _build_demo_target_digest(self):
        field = GF(self.config.prime)
        poseidon = Poseidon(
            prime=self.config.prime,
            R_F=self.config.R_F,
            R_P=self.config.R_P,
            t=self.config.t,
            alpha=self.config.alpha,
        )
        input_unknowns = [field(i) for i in range(self.config.digest)]
        input_padding = [field(0)] * self.config.linearizable_rounds
        input_capacity = [field(0)] * self.config.capacity
        full_input = input_unknowns + input_padding + input_capacity
        full_output = poseidon(full_input)
        return list(full_output[: self.config.digest])

    def build_demo_vector(self):
        field = GF(self.config.prime)
        input_unknowns = [field(i) for i in range(self.config.digest)]
        input_padding = [field(0)] * self.config.linearizable_rounds
        input_capacity = [field(0)] * self.config.capacity
        full_input = input_unknowns + input_padding + input_capacity
        poseidon = Poseidon(
            prime=self.config.prime,
            R_F=self.config.R_F,
            R_P=self.config.R_P,
            t=self.config.t,
            alpha=self.config.alpha,
        )
        full_output = poseidon(full_input)
        return {
            "input_unknowns": input_unknowns,
            "input_padding": input_padding,
            "input_capacity": input_capacity,
            "full_input": full_input,
            "full_output": full_output,
            "target_digest": list(full_output[: self.config.digest]),
        }

    def _compute_symbolic_output(self, poseidon, padded_input):
        t = self.config.t
        R_F = poseidon.R_F
        R_P = poseidon.R_P
        alpha = self.config.alpha

        state = list(padded_input) + [0] * self.config.capacity
        state = list(poseidon.MDS_matrix_field * vector(state))

        num_rounds = R_F + R_P - 1
        full_round_prefix = (R_F - 1) // 2
        for round_index in range(num_rounds):
            state = [
                state[i] + poseidon.round_constants_field[round_index][i]
                for i in range(t)
            ]
            if round_index < full_round_prefix or round_index >= full_round_prefix + R_P:
                state = [entry ** alpha for entry in state]
            else:
                state[0] = state[0] ** alpha
            state = list(poseidon.MDS_matrix_field * vector(state))
        return state

    def summarize_system(self, system, build_time):
        degrees = [poly.degree() for poly in system]
        histogram = {d: degrees.count(d) for d in range(min(degrees), max(degrees) + 1)}
        return {
            "equation_type": "direct_high_degree",
            "R_F": self.config.R_F,
            "R_P": self.config.R_P,
            "state_size": self.config.t,
            "capacity": self.config.capacity,
            "rate": self.config.rate,
            "digest": self.config.digest,
            "linearizable_rounds": self.config.linearizable_rounds,
            "target_digest": self.config.target_digest,
            "macaulay_bound": 1 + sum(poly.degree() - 1 for poly in system),
            "build_time": build_time,
            "num_equations": len(system),
            "num_variables": len(system[0].parent().gens()) if system else 0,
            "max_degree": max(degrees) if degrees else 0,
            "degree_histogram": histogram if degrees else {},
            "num_coefficients": sum(len(poly.coefficients()) for poly in system),
        }

    def write_outputs(self, system, summary):
        prefix = Path(self.config.result_prefix)
        prefix.parent.mkdir(parents=True, exist_ok=True)
        system_path = Path(str(prefix) + "_system.txt")
        summary_path = Path(str(prefix) + "_summary.txt")

        with system_path.open("w") as handle:
            for index, poly in enumerate(system):
                suffix = "\n" if index == len(system) - 1 else ",\n"
                handle.write(f"{poly}{suffix}")

        with summary_path.open("w") as handle:
            for key, value in summary.items():
                handle.write(f"{key}: {value}\n")

        return system_path, summary_path

    def run(self):
        total_start = perf_counter()
        start = process_time()
        system = self.build_system()
        build_time = process_time() - start
        summary = self.summarize_system(system, build_time)
        summary["total_time"] = perf_counter() - total_start
        system_path, summary_path = self.write_outputs(system, summary)
        return system, summary, system_path, summary_path


class ExperimentStarter:
    """
    Backward-compatible wrapper around `PoseidonEquationGenerator`.
    """

    def __init__(self, result_path, t, capacity, digest, target_digest=None, alpha=3):
        self.result_path = result_path
        self.t = t
        self.capacity = capacity
        self.digest = digest
        self.target_digest = target_digest
        self.alpha = alpha

    def __call__(self, primitive_name, prime, R_F, R_P):
        if primitive_name != "poseidon":
            raise ValueError(f"No primitive with name {primitive_name} defined.")

        config = PoseidonCICOConfig(
            prime=prime,
            t=self.t,
            capacity=self.capacity,
            digest=self.digest,
            R_F=R_F,
            R_P=R_P,
            alpha=self.alpha,
            target_digest=self.target_digest,
            result_prefix=self.result_path.rstrip("_"),
        )
        generator = PoseidonEquationGenerator(config)
        system, _, _, _ = generator.run()
        return system


def print_poseidon_demo(config, demo_vector):
    print("=" * 70)
    print("Poseidon Parameter Configuration")
    print("=" * 70)
    for key, value in config.to_summary().items():
        print(f"{key}: {value}")
    print("=" * 70)
    print("\nTest vector:")
    print(f"  Unknown input ({config.digest} elements): {demo_vector['input_unknowns']}")
    print(f"  Padding ({config.linearizable_rounds} elements): {demo_vector['input_padding']}")
    print(f"  Capacity ({config.capacity} elements): {demo_vector['input_capacity']}")
    print(f"  Full input: {demo_vector['full_input']}")
    print(f"  Full output: {demo_vector['full_output']}")
    print(f"  Target digest: {demo_vector['target_digest']}")


def run_demo():
    total_start = perf_counter()
    set_verbose(2)

    config = PoseidonCICOConfig(
        prime=4611686018427388073,
        t=6,
        capacity=1,
        digest=3,
        R_F=4,
        R_P=0,
        alpha=3,
        result_prefix="./experiments/poseidon_t6_c1_d3_RF4_RP0_direct_high_degree",
    )

    generator = PoseidonEquationGenerator(config)
    demo_vector = generator.build_demo_vector()
    config.target_digest = demo_vector["target_digest"]
    print_poseidon_demo(config, demo_vector)

    system, summary, system_path, summary_path = generator.run()

    print("\nPolynomial system constructed.")
    print(f"  Equations: {summary['num_equations']}")
    print(f"  Variables: {summary['num_variables']}")
    print(f"  Maximum degree: {summary['max_degree']}")
    print(f"  Build time: {summary['build_time']:.6f} s")
    print(f"  Total time: {summary['total_time']:.6f} s")
    print(f"  System file: {system_path}")
    print(f"  Summary file: {summary_path}")

    print("\nFirst equations:")
    for poly in system[: min(3, len(system))]:
        poly_str = str(poly)
        if len(poly_str) > 1000:
            print(f"  {poly_str[:120]}...")
        else:
            print(f"  {poly}")

    print(f"\nEnd-to-end demo time: {perf_counter() - total_start:.6f} s")


if __name__ == "__main__":
    run_demo()
