# Generate the system described in Liu et al. "Modelling Ciphers with Overdefined Systems of Quadratic Equations: Application to Friday, Vision, RAIN and Biscuit"

import itertools
from dataclasses import dataclass
from time import perf_counter

from sage.all import *

@dataclass
class VisionMITMConfig:
    prime: int = 2
    field_size: int = 8
    num_rounds: int = 1
    state_size: int = 2
    rate: int = 1
    security_level: int = 0

    @property
    def capacity(self):
        return self.state_size - self.rate

    def validate(self):
        assert self.rate >= 1, "Rate must be positive."
        assert self.capacity >= 0, "Capacity must be non-negative."
        assert self.capacity * self.field_size >= 2 * self.security_level, (
            f"Capacity size {self.capacity * self.field_size} cannot provide the desired security level {self.security_level}"
        )


class VisionMITMGenerator:
    """
    Generate the MITM polynomial system for Vision with `s=2`, matching the
    equation style described in `dixon_appendix.tex`.

    This refactor keeps the original algebraic model but makes the data flow
    explicit: affine layers, constants, permutation evaluation, trace recovery,
    and exported equation systems now live in one object.
    """

    def __init__(self, config):
        self.config = config
        self.config.validate()
        self.field = GF(self.config.prime ** self.config.field_size)

    def create_affine_layers(self):
        n = self.config.field_size
        field = self.field

        coeffs = [field(0)] * (n + 1)
        for i in range(3):
            coeffs[i] = field.random_element()
        matrix_affine = matrix(field, n, n, lambda i, j: coeffs[(i - j + n) % n] ** (2 ** j))
        while matrix_affine.rank() != n:
            for i in range(3):
                coeffs[i] = field.random_element()
            matrix_affine = matrix(field, n, n, lambda i, j: coeffs[(i - j + n) % n] ** (2 ** j))

        basis_vector = [field(0)] * n
        basis_vector[0] = field(1)
        affine_inverse = matrix_affine.solve_right(vector(basis_vector))
        return coeffs, affine_inverse

    def create_constants(self):
        field = self.field
        n = self.config.field_size
        s = self.config.state_size
        num_rounds = self.config.num_rounds

        g = field.gen()
        vandermonde = matrix([[g ** (i * j) for j in range(2 * s)] for i in range(s)]).echelon_form()
        mds = vandermonde[:, s:].transpose()
        round_constants = matrix(field, s, 2 * num_rounds + 1, lambda i, j: field.random_element())
        affine_coeffs, affine_inverse = self.create_affine_layers()
        return mds, round_constants.transpose(), affine_coeffs, affine_inverse

    def affine_layer(self, coeffs, state):
        return sum(coeffs[i] * state ** (2 ** i) for i in range(len(coeffs)))

    def evaluate_permutation(self, mds, constants, affine_coeffs, affine_inverse, input_state, verbose=False):
        output = matrix(self.field, self.config.state_size, 1, input_state)
        for j in range(self.config.state_size):
            output[j, 0] += constants[0][j]

        for round_index in range(self.config.num_rounds):
            for j in range(self.config.state_size):
                output[j, 0] = output[j, 0] ** -1
                output[j, 0] = self.affine_layer(affine_inverse, output[j, 0])
            if verbose:
                print(f"State after round {round_index} inverse/B^-1:\n{output}")

            output = mds * output
            for j in range(self.config.state_size):
                output[j, 0] += constants[2 * round_index + 1][j]

            for j in range(self.config.state_size):
                output[j, 0] = output[j, 0] ** -1
                output[j, 0] = self.affine_layer(affine_coeffs, output[j, 0])
            if verbose:
                print(f"State after round {round_index} inverse/B:\n{output}")

            output = mds * output
            for j in range(self.config.state_size):
                output[j, 0] += constants[2 * round_index + 2][j]

        return output

    def build_system(self, mds, constants, affine_coeffs, affine_inverse, output_state):
        s = self.config.state_size
        r = self.config.rate
        num_rounds = self.config.num_rounds
        ring_size = 2 * s * (num_rounds - 1) + s + r
        ring = PolynomialRing(self.field, ring_size, "x")
        x = ring.gens()

        var_counter = 0
        w_layers = [[[x[i]] for i in range(r)] + [[0] for _ in range(self.config.capacity)]]
        var_counter += r
        for i in range(self.config.capacity):
            w_layers[0][r + i][0] = self.affine_layer(affine_inverse, constants[0][r + i] ** -1)

        z_layers = [[list(x[var_counter + i : var_counter + i + 1]) for i in range(s)]]
        var_counter += s
        for _ in range(1, num_rounds - 1):
            w_layers.append([list(x[var_counter + i : var_counter + i + 1]) for i in range(s)])
            var_counter += s
            z_layers.append([list(x[var_counter + i : var_counter + i + 1]) for i in range(s)])
            var_counter += s
        w_layers.append([list(x[var_counter + i : var_counter + i + 1]) for i in range(s)])
        var_counter += s
        z_layers.append([list(x[var_counter + i : var_counter + i + 1]) for i in range(s)])

        output_values = list(itertools.chain.from_iterable(list(output_state.transpose())))
        equations = []

        for round_index in range(1, num_rounds):
            for j in range(s):
                backward = sum(affine_coeffs[k] * w_layers[round_index][j][0] ** (2 ** k) for k in range(3))
                forward = sum(
                    mds[j][k]
                    * sum(affine_coeffs[h] * z_layers[round_index - 1][k][0] ** (2 ** h) for h in range(3))
                    for k in range(s)
                )
                equations.append((forward + constants[2 * round_index][j]) * backward + 1)

        for j in range(s):
            sum_mds = sum(mds[j][k] * w_layers[0][k][0] for k in range(s)) + constants[1][j]
            equations.append(z_layers[0][j][0] * sum_mds + 1)

        for round_index in range(1, num_rounds):
            for j in range(s):
                sum_mds = sum(mds[j][k] * w_layers[round_index][k][0] for k in range(s)) + constants[2 * round_index + 1][j]
                equations.append(z_layers[round_index][j][0] * sum_mds + 1)

        for j in range(r):
            forward = sum(
                mds[j][k]
                * sum(affine_coeffs[h] * z_layers[num_rounds - 1][k][0] ** (2 ** h) for h in range(3))
                for k in range(s)
            )
            equations.append(output_values[j] + forward + constants[2 * num_rounds][j])

        return equations

    def build_test_trace(self, mds, constants, affine_coeffs, affine_inverse, input_state):
        s = self.config.state_size
        r = self.config.rate
        num_rounds = self.config.num_rounds

        trace = []
        output = matrix(self.field, s, 1, input_state)
        for j in range(s):
            output[j, 0] += constants[0][j]

        for round_index in range(num_rounds):
            for j in range(s):
                output[j, 0] = output[j, 0] ** -1
                output[j, 0] = self.affine_layer(affine_inverse, output[j, 0])

            for j in range(s):
                if round_index > 0:
                    trace.extend([output[j, 0], output[j, 0] ** 2, output[j, 0] ** 4])
                elif j < r:
                    trace.append(output[j, 0])

            output = mds * output
            for j in range(s):
                output[j, 0] += constants[2 * round_index + 1][j]

            for j in range(s):
                output[j, 0] = output[j, 0] ** -1
                trace.extend([output[j, 0], output[j, 0] ** 2, output[j, 0] ** 4])
                output[j, 0] = self.affine_layer(affine_coeffs, output[j, 0])

            output = mds * output
            for j in range(s):
                output[j, 0] += constants[2 * round_index + 2][j]

        return trace

    def verify_system(self, equations, mds, constants, affine_coeffs, affine_inverse, input_state):
        test_value = self.field.random_element()
        assert self.affine_layer(affine_inverse, self.affine_layer(affine_coeffs, test_value)) == test_value, (
            "Affine layer inverse is not correct."
        )
        trace = self.build_test_trace(mds, constants, affine_coeffs, affine_inverse, input_state)
        for poly in equations:
            assert poly(trace) == 0, f"Polynomial {poly} does not vanish on the trace."
        return True

    def summarize(self, equations):
        degrees = [poly.degree() for poly in equations]
        return {
            "field": self.field,
            "num_rounds": self.config.num_rounds,
            "state_size": self.config.state_size,
            "rate": self.config.rate,
            "capacity": self.config.capacity,
            "num_equations": len(equations),
            "num_variables": len(equations[0].parent().gens()) if equations else 0,
            "max_degree": max(degrees) if degrees else 0,
            "degree_histogram": {d: degrees.count(d) for d in range(min(degrees), max(degrees) + 1)} if degrees else {},
        }

    def export_system(self, equations, path):
        with open(path, "w") as handle:
            for index, poly in enumerate(equations):
                suffix = "\n" if index == len(equations) - 1 else ",\n"
                handle.write(f"{poly}{suffix}")


def run_demo():
    total_start = perf_counter()
    config = VisionMITMConfig()
    generator = VisionMITMGenerator(config)
    mds, constants, affine_coeffs, affine_inverse = generator.create_constants()

    print("Field polynomial:", generator.field.polynomial())
    input_state = [generator.field.random_element()] * config.rate + [0] * config.capacity
    output_state = generator.evaluate_permutation(
        mds, constants, affine_coeffs, affine_inverse, input_state, verbose=True
    )
    equations = generator.build_system(mds, constants, affine_coeffs, affine_inverse, output_state)
    summary = generator.summarize(equations)

    print("Vision MITM System Summary:")
    for key, value in summary.items():
        print(f"  {key}: {value}")

    generator.export_system(equations, f"./system_{config.num_rounds}.txt")
    print(f"  total_time: {perf_counter() - total_start:.6f} s")


if __name__ == "__main__":
    run_demo()
