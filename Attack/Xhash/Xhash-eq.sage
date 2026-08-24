from dataclasses import dataclass
from hashlib import shake_256
from time import perf_counter

from sage.all import *


@dataclass
class XHashConfig:
    p: int = 65537
    state_size: int = 3
    rate: int = 1
    capacity: int = 2
    alpha: int = 3
    mds_first_row: object = None
    num_constants: object = None


class XHashPolynomialSystem:
    """
    XHash system generator.

    The core polynomial formulas are kept unchanged. This refactor only
    standardizes configuration, summaries, and export helpers, and adds an
    export mode tailored to `drsolve_sage_interface`.
    """
    
    def __init__(self, p=65537, state_size=3, rate=1, capacity=2, alpha=3,
                 mds_first_row=None, num_constants=None):
        """
        Initialize XHash polynomial system
        
        Args:
            p: Prime field characteristic
            state_size: State size (t)
            rate: Rate (r) 
            capacity: Capacity (c)
            alpha: S-box exponent
            mds_first_row: First row of MDS matrix (if None, use default)
            num_constants: Number of constants to generate
        """
        self.p = p
        self.state_size = state_size
        self.rate = rate
        self.capacity = capacity
        self.alpha = alpha
        
        # Initialize field and polynomial ring
        self.Fp = GF(p)
        self.max_vars = 200
        self.R = PolynomialRing(self.Fp, self.max_vars, 'x', order='lex')
        self.gens = self.R.gens()
        self.counter = 0
        
        # Generate MDS matrix
        if mds_first_row is None:
            if state_size == 3:
                mds_first_row = [7, 23, 8]
            elif state_size == 4:
                mds_first_row = [7, 23, 8, 26]
            elif state_size == 12:
                mds_first_row = [7, 23, 8, 26, 13, 10, 9, 7, 6, 22, 21, 8]
            else:
                mds_first_row = list(range(1, state_size + 1))
        
        self.mds = self._generate_mds_matrix([self.Fp(x) for x in mds_first_row])
        
        # Generate constants
        if num_constants is None:
            num_constants = 10 * state_size  # Default: enough for several rounds
        self.constants = self._generate_round_constants(num_constants)
        self.constant_index = 0
    
    def _rotate(self, lst):
        """Rotate list one position to the right"""
        return lst[-1:] + lst[:-1]
    
    def _generate_mds_matrix(self, first_row):
        """Generate circulant MDS matrix"""
        matrix = [first_row]
        size = len(first_row)
        for _ in range(1, size):
            next_row = self._rotate(matrix[-1])
            matrix.append(next_row)
        return Matrix(self.Fp, matrix)
    
    def _generate_round_constants(self, num_constants):
        """Generate round constants using SHAKE256"""
        base_string = f"XHash({self.p},{self.state_size},{self.capacity},{self.alpha})"
        hash_result = shake_256(base_string.encode()).digest(9 * 2 * num_constants)
        
        constants = []
        for i in range(0, len(hash_result), 9):
            chunk = hash_result[i:i + 9]
            int_value = int.from_bytes(chunk, byteorder='little')
            reduced_value = int_value % self.p
            constants.append(self.Fp(reduced_value))
        return constants

    @classmethod
    def from_config(cls, config):
        return cls(
            p=config.p,
            state_size=config.state_size,
            rate=config.rate,
            capacity=config.capacity,
            alpha=config.alpha,
            mds_first_row=config.mds_first_row,
            num_constants=config.num_constants,
        )
    
    def new_var(self, n):
        """Generate n new variables"""
        if self.counter + n > self.max_vars:
            raise ValueError(f"Too many variables needed: {self.counter + n}")
        selected_vars = self.gens[self.counter:self.counter + n]
        self.counter += n
        return list(selected_vars)
    
    def reset_variables(self):
        """Reset variable counter"""
        self.counter = 0
        self.constant_index = 0
    
    def get_constants(self, n):
        """Get next n constants"""
        result = self.constants[self.constant_index:self.constant_index + n]
        self.constant_index += n
        return result
    
    def _add_ideal_relation(self, ideal, lhs, rhs):
        """Add an ideal generator as a Sage-object relation [lhs, rhs]."""
        if ideal is not None:
            ideal.append([lhs, rhs])
    
    # Basic operations (following original naming pattern)
    def ADC(self, x):
        """Add constants"""
        constants = self.get_constants(len(x))
        return [a + constants[i] for i, a in enumerate(x)]
    
    def M(self, x):
        """MDS matrix multiplication"""
        return list(vector(x) * self.mds)
    
    def SF(self, x):
        """S-box forward: x^alpha"""
        return [a ** self.alpha for a in x]
    
    def SB(self, x):
        """S-box backward: x^alpha (to be used in reverse equations)"""
        return [a ** self.alpha for a in x]
    
    def SP3(self, x, block_size=3, power=None):
        """Extension field S-box with configurable block size and power"""
        if power is None:
            power = self.alpha
            
        if block_size == 3 and power == 7:
            return self._sp3_block3_power7(x)
        elif block_size == 3 and power == 3:
            return self._sp3_block3_power3(x)
        elif block_size == 2 and power == 7:
            return self._sp2_block2_power7(x)
        elif block_size == 2 and power == 3:
            return self._sp2_block2_power3(x)
        else:
            # Fallback: element-wise power
            return [xi ** power for xi in x]
    
    def _sp3_block3_power7(self, x):
        """F_p^3 with 7th power (exact formula from paper)"""
        if len(x) < 3:
            x = x + [self.Fp(0)] * (3 - len(x))
        x0, x1, x2 = x[0], x[1], x[2]
        
        y0 = (x0**7 + 35*x0**4*x1**3 + 21*x0**2*x1**5 + 7*x0*x1**6 + x1**7 +
              42*x0**5*x1*x2 + 140*x0**3*x1**3*x2 + 105*x0**2*x1**4*x2 + 42*x0*x1**5*x2 + 14*x1**6*x2 +
              105*x0**4*x1*x2**2 + 210*x0**3*x1**2*x2**2 + 210*x0**2*x1**3*x2**2 + 210*x0*x1**4*x2**2 + 42*x1**5*x2**2 +
              35*x0**4*x2**3 + 140*x0**3*x1*x2**3 + 420*x0**2*x1**2*x2**3 + 280*x0*x1**3*x2**3 + 105*x1**4*x2**3 +
              70*x0**3*x2**4 + 210*x0**2*x1*x2**4 + 315*x0*x1**2*x2**4 + 140*x1**3*x2**4 +
              63*x0**2*x2**5 + 168*x0*x1*x2**5 + 105*x1**2*x2**5 +
              35*x0*x2**6 + 49*x1*x2**6 +
              9*x2**7)
        
        y1 = (7*x0**6*x1 + 35*x0**4*x1**3 + 35*x0**3*x1**4 + 21*x0**2*x1**5 + 14*x0*x1**6 + 2*x1**7 +
              42*x0**5*x1*x2 + 105*x0**4*x1**2*x2 + 140*x0**3*x1**3*x2 + 210*x0**2*x1**4*x2 + 84*x0*x1**5*x2 + 21*x1**6*x2 +
              21*x0**5*x2**2 + 105*x0**4*x1*x2**2 + 420*x0**3*x1**2*x2**2 + 420*x0**2*x1**3*x2**2 + 315*x0*x1**4*x2**2 + 84*x1**5*x2**2 +
              70*x0**4*x2**3 + 280*x0**3*x1*x2**3 + 630*x0**2*x1**2*x2**3 + 560*x0*x1**3*x2**3 + 175*x1**4*x2**3 +
              105*x0**3*x2**4 + 420*x0**2*x1*x2**4 + 525*x0*x1**2*x2**4 + 245*x1**3*x2**4 +
              105*x0**2*x2**5 + 294*x0*x1*x2**5 + 189*x1**2*x2**5 +
              63*x0*x2**6 + 84*x1*x2**6 +
              16*x2**7)
        
        y2 = (21*x0**5*x1**2 + 35*x0**3*x1**4 + 21*x0**2*x1**5 + 7*x0*x1**6 + 2*x1**7 +
              7*x0**6*x2 + 105*x0**4*x1**2*x2 + 140*x0**3*x1**3*x2 + 105*x0**2*x1**4*x2 + 84*x0*x1**5*x2 + 14*x1**6*x2 +
              21*x0**5*x2**2 + 105*x0**4*x1*x2**2 + 210*x0**3*x1**2*x2**2 + 420*x0**2*x1**3*x2**2 + 210*x0*x1**4*x2**2 + 63*x1**5*x2**2 +
              35*x0**4*x2**3 + 280*x0**3*x1*x2**3 + 420*x0**2*x1**2*x2**3 + 420*x0*x1**3*x2**3 + 140*x1**4*x2**3 +
              70*x0**3*x2**4 + 315*x0**2*x1*x2**4 + 420*x0*x1**2*x2**4 + 175*x1**3*x2**4 +
              84*x0**2*x2**5 + 210*x0*x1*x2**5 + 147*x1**2*x2**5 +
              49*x0*x2**6 + 63*x1*x2**6 +
              12*x2**7)
        
        return [y0, y1, y2]
    
    def _sp3_block3_power3(self, x):
        """F_p^3 with 3rd power (original SP3 logic)"""
        if len(x) < 3:
            x = x + [self.Fp(0)] * (3 - len(x))
        x0, x1, x2 = x[0], x[1], x[2]
        
        y0 = (3*x0**3 + 16*x0**2*x1 - 6*x0**2 - 15*x0*x1**2 + 3*x0*x1*x2 + 6*x0*x1 + 
              16*x0*x2**2 - 5*x0*x2 - 13*x0 - 16*x1**3 + 8*x1**2*x2 - 18*x1**2 - 
              9*x1*x2**2 + 7*x1*x2 + 8*x1 - 13*x2**3 + 9*x2**2 - 18*x2 - 4)
        
        y1 = (-4*x0**3 - 4*x0**2*x1 - 15*x0**2*x2 + 10*x0**2 + 15*x0*x1**2 - 12*x0*x1*x2 + 
              18*x0*x1 - x0*x2**2 + 4*x0*x2 + 2*x0 + 11*x1**3 + 6*x1**2*x2 - 10*x1**2 - 
              4*x1*x2**2 + 8*x1*x2 + 16*x1 + 12*x2**3 - 18*x2**2 + 7)
        
        y2 = (16*x0**3 - 11*x0**2*x1 + 15*x0**2*x2 + 3*x0**2 + 4*x0*x1**2 - 15*x0*x1*x2 - 
              10*x0*x1 + 14*x0*x2**2 - 5*x0*x2 - 10*x0 + 2*x1**3 - 5*x1**2*x2 - 16*x1**2 - 
              8*x1*x2**2 + x1 + 15*x2**3 + 3*x2**2 - 11*x2 - 2)
        
        return [y0, y1, y2]
    
    def _sp2_block2_power7(self, x):
        """F_p^2 with 7th power"""
        if len(x) < 2:
            x = x + [self.Fp(0)] * (2 - len(x))
        x0, x1 = x[0], x[1]
        
        y0 = 2*x0^7 - 3*x0^6*x1 - 4*x0^6 + 7*x0^5*x1^2 + x0^5*x1 + x0^5 + 7*x0^4*x1^3 - x0^4*x1^2 + 3*x0^4*x1 + 5*x0^4 - x0^3*x1^4 - 3*x0^3*x1^3 - x0^3*x1^2 + 6*x0^3*x1 - 5*x0^3 - 6*x0^2*x1^5 + 6*x0^2*x1^4 - 6*x0^2*x1^3 - 6*x0^2*x1^2 + 7*x0^2*x1 + 8*x0^2 + 7*x0*x1^6 + 7*x0*x1^5 + 8*x0*x1^4 - 7*x0*x1^3 - 4*x0*x1^2 - 5*x0*x1 - 5*x0 + 4*x1^7 + 2*x1^5 + x1^4 + x1^3 - 7*x1^2 + 5*x1
        y1 = 3*x0^7 - x0^6*x1 - 6*x0^6 - 3*x0^5*x1^2 + 8*x0^5*x1 + x0^5 + 6*x0^4*x1^2 + 3*x0^4*x1 + 8*x0^4 - 2*x0^3*x1^4 + 4*x0^3*x1^3 + x0^3*x1^2 - 5*x0^3*x1 + x0^3 + 4*x0^2*x1^5 - 2*x0^2*x1^4 + 8*x0^2*x1^3 + x0^2*x1^2 - 3*x0^2 + 3*x0*x1^6 + 3*x0*x1^5 - 3*x0*x1^4 + x0*x1^3 - 2*x0*x1^2 - 2*x0*x1 + 2*x0 + 3*x1^7 + 2*x1^6 - x1^5 + 7*x1^4 - 6*x1^3 - 7*x1^2 - x1 + 5
        
        return [y0, y1]
    
    def _sp2_block2_power3(self, x):
        """F_p^2 with 3rd power"""
        if len(x) < 2:
            x = x + [self.Fp(0)] * (2 - len(x))
        x0, x1 = x[0], x[1]
        
        y0 = 5*x0^3 + 4*x0^2*x1 + 4*x0^2 + x0*x1^2 - 8*x0*x1 - 8*x0 - 5*x1^3 - 8*x1^2 - 4*x1 + 6
        y1 = 6*x0^3 + x0^2*x1 - 5*x0^2 - 2*x0*x1^2 + 4*x0*x1 - 3*x0 - x1^3 - 5*x1^2 + 8*x1 + 1

        return [y0, y1]
    
    # Round functions with ideal string generation
    def F_step(self, x, equations, ideal=None):
        """F step: x -> ADC -> M -> SF -> y"""
        x1 = self.ADC(x)
        x2 = self.M(x1)
        x3 = self.SF(x2)
        y = self.new_var(len(x))
        
        # Add equations and ideal strings
        for i in range(len(x)):
            equations.append(y[i] - x3[i])
            self._add_ideal_relation(ideal, y[i], x3[i])
        
        return y
    
    def B_step(self, x, equations, ideal=None):
        """B step: x -> M -> ADC -> SB -> y (using reverse equation)"""
        x1 = self.M(x)
        x2 = self.ADC(x1)
        y = self.new_var(len(x))
        
        # Use reverse equation: SB(y) = x2, i.e., y^alpha = x2
        for i in range(len(x)):
            equations.append(y[i] ** self.alpha - x2[i])
            self._add_ideal_relation(ideal, y[i] ** self.alpha, x2[i])
        
        return y
    
    def P3_step(self, x, equations, ideal=None, block_size=3, power=None):
        """P3 step: x -> ADC -> SP3 -> y with configurable block size and power"""
        if power is None:
            power = self.alpha
            
        x1 = self.ADC(x)
        
        # Apply SP3 to appropriate number of elements based on block size
        if len(x) >= block_size:
            x2_partial = self.SP3(x1[:block_size], block_size, power)
            x2 = x2_partial + x1[block_size:]  # Keep remaining elements unchanged
        else:
            x2 = self.SP3(x1, block_size, power)  # SP3 will handle padding
            
        y = self.new_var(len(x))
        
        # Map back to original size and add equations/ideal strings
        for i in range(len(x)):
            if i < len(x2):
                equations.append(y[i] - x2[i])
                self._add_ideal_relation(ideal, y[i], x2[i])
            else:
                equations.append(y[i] - x[i])  # Identity for extra elements
                self._add_ideal_relation(ideal, y[i], x[i])
        
        return y
    
    def MC_step(self, x, equations, ideal=None):
        """MC step: x -> M -> ADC -> y"""
        x1 = self.M(x)
        x2 = self.ADC(x1)
        y = self.new_var(len(x))
        
        for i in range(len(x)):
            equations.append(y[i] - x2[i])
            self._add_ideal_relation(ideal, y[i], x2[i])
        
        return y
    
    def F0_step(self, x, equations, ideal=None):
        """F0 step: x -> ADC -> M -> y (no S-box)"""
        x1 = self.ADC(x)
        x2 = self.M(x1)
        y = self.new_var(len(x))
        
        for i in range(len(x)):
            equations.append(x2[i] - y[i])
            self._add_ideal_relation(ideal, y[i], x2[i])
        
        return y
    
    def FB_step(self, x, equations, ideal=None):
        """FB step: x -> SF -> M -> ADC -> SB^-1 -> y"""
        x1 = self.SF(x)
        x2 = self.M(x1)
        x3 = self.ADC(x2)
        y = self.new_var(len(x))
        
        # Reverse equation for SB: SF(y) = x3
        y1 = self.SF(y)
        for i in range(len(x)):
            equations.append(x3[i] - y1[i])
            self._add_ideal_relation(ideal, y[i] ** self.alpha, x3[i])
        
        return y
    
    def F0A_step(self, x, equations, ideal=None):
        """F0A step: complex multi-step operation"""
        x1 = self.ADC(x)
        if len(x1) >= 3:
            x2 = self.SP3(x1[:3]) + x1[3:]  # Apply SP3 to first 3 elements
        else:
            x2 = x1
        x3 = self.ADC(x2)
        x4 = self.M(x3)
        x5 = self.SF(x4)
        x6 = self.M(x5)
        x7 = self.ADC(x6)
        y = self.new_var(len(x))
        
        # Reverse equation
        y1 = self.SF(y)
        for i in range(len(x)):
            equations.append(x7[i] - y1[i])
            self._add_ideal_relation(ideal, y[i] ** self.alpha, x7[i])
        
        return y
    
    def F0B_step(self, x, equations, ideal=None):
        """F0B step"""
        x1 = self.ADC(x)
        x2 = self.M(x1)
        x3 = self.SF(x2)
        x4 = self.M(x3)
        x5 = self.ADC(x4)
        y = self.new_var(len(x))
        
        # Reverse equation
        y1 = self.SF(y)
        for i in range(len(x)):
            equations.append(x5[i] - y1[i])
            self._add_ideal_relation(ideal, y[i] ** self.alpha, x5[i])
        
        return y
    
    def P3F_step(self, x, equations, ideal=None, block_size=3, power=None):
        """P3F step with configurable block size and power"""
        if power is None:
            power = self.alpha
            
        x1 = self.ADC(x)
        if len(x) >= block_size:
            x2_partial = self.SP3(x1[:block_size], block_size, power)
            x2 = x2_partial + x1[block_size:]
        else:
            x2 = self.SP3(x1, block_size, power)
        x3 = self.ADC(x2)
        x4 = self.M(x3)
        y = self.new_var(len(x))
        
        for i in range(len(x)):
            equations.append(x4[i] - y[i])
            self._add_ideal_relation(ideal, y[i], x4[i])
        
        return y
    
    def P3MC_step(self, x, equations, ideal=None, block_size=3, power=None):
        """P3MC step with configurable block size and power"""
        if power is None:
            power = self.alpha
            
        x1 = self.ADC(x)
        if len(x) >= block_size:
            x2_partial = self.SP3(x1[:block_size], block_size, power)
            x2 = x2_partial + x1[block_size:]
        else:
            x2 = self.SP3(x1, block_size, power)
        x3 = self.M(x2)
        x4 = self.ADC(x3)
        y = self.new_var(len(x))
        
        for i in range(len(x)):
            equations.append(x4[i] - y[i])
            self._add_ideal_relation(ideal, y[i], x4[i])
        
        return y
    
    # Predefined system patterns with ideal generation
    def generate_simple_system(self, pattern="F-B", num_rounds=2, block_size=3, power=None, ideal=None):
        """Generate simple system with given pattern"""
        if power is None:
            power = self.alpha
            
        self.reset_variables()
        equations = []
        if ideal is None:
            ideal = []
        
        # Input state
        xx = list(self.new_var(self.rate)) + [self.Fp(0)] * self.capacity
        
        if pattern == "F-B":
            for _ in range(num_rounds):
                xx = self.F_step(xx, equations, ideal)
                xx = self.B_step(xx, equations, ideal)
        elif pattern == "F-B-P3":
            for _ in range(num_rounds):
                xx = self.F_step(xx, equations, ideal)
                xx = self.B_step(xx, equations, ideal)
                xx = self.P3_step(xx, equations, ideal, block_size, power)
            xx = self.MC_step(xx, equations, ideal)
        elif pattern == "F0B-F0A":
            xx = self.F0B_step(xx, equations, ideal)
            for _ in range(num_rounds):
                xx = self.F0A_step(xx, equations, ideal)
        elif pattern == "F0B-P3F":
            for _ in range(num_rounds):
                xx = self.F0B_step(xx, equations, ideal)
                xx = self.P3F_step(xx, equations, ideal, block_size, power)
        elif pattern == "P3F-FB-P3MC":
            for _ in range(num_rounds):
                xx = self.P3F_step(xx, equations, ideal, block_size, power)
                xx = self.FB_step(xx, equations, ideal)
            xx = self.P3MC_step(xx, equations, ideal, block_size, power)
        
        return equations, xx, ideal
    
    def print_system_info(self, equations, name="System", ideal=None):
        """Print system information"""
        print(f"{name} Info:")
        print(f"  Field: GF({self.p})")
        print(f"  State size: {self.state_size}, Rate: {self.rate}, Capacity: {self.capacity}")
        print(f"  S-box: x^{self.alpha}")
        print(f"  Equations: {len(equations)}")
        print(f"  Variables: {self.counter}")
        
        if ideal:
            print(f"  Ideal generators: {len(ideal)}")
        
        if equations:
            degrees = []
            for eq in equations:
                try:
                    degrees.append(eq.degree())
                except:
                    degrees.append(1)
            if degrees:
                avg_degree = float(sum(degrees)/len(degrees))  # Convert to Python float
                print(f"  Degrees: min={min(degrees)}, max={max(degrees)}, avg={avg_degree:.1f}")

    def summarize(self, equations, ideal=None, pattern=None, num_rounds=None, block_size=None, power=None):
        degrees = []
        for eq in equations:
            try:
                degrees.append(eq.degree())
            except Exception:
                degrees.append(1)
        return {
            "field": f"GF({self.p})",
            "state_size": self.state_size,
            "rate": self.rate,
            "capacity": self.capacity,
            "alpha": self.alpha,
            "pattern": pattern,
            "num_rounds": num_rounds,
            "block_size": block_size,
            "power": power,
            "num_equations": len(equations),
            "num_variables": self.counter,
            "num_ideal_generators": len(ideal) if ideal else 0,
            "min_degree": min(degrees) if degrees else 0,
            "max_degree": max(degrees) if degrees else 0,
        }

    def export_system(self, equations, format="sage", name="system", ideal=None):
        """Export system in plain Sage/drsolve-friendly formats."""
        if format == "sage":
            print(f"# {name}")
            print(f"F = GF({self.p})")
            print(f"R = PolynomialRing(F, {self.counter}, 'x')")
            print("equations = [")
            for eq in equations:
                print(f"    {eq},")
            print("]")
        elif format in ("plain", "string", "equations"):
            for i, eq in enumerate(equations):
                suffix = "" if i == len(equations) - 1 else ","
                print(f"{eq}{suffix}")
        elif format in ("ideal", "ideal_string"):
            if ideal:
                for i, relation in enumerate(ideal):
                    lhs, rhs = relation
                    suffix = "" if i == len(ideal) - 1 else ","
                    print(f"{lhs} = {rhs}{suffix}")
            else:
                print("No ideal relations generated")
        elif format == "drsolve_sage_interface":
            print(f"# {name} for drsolve_sage_interface")
            print(f"field_size = {self.p}")
            print("ideal = [")
            for lhs, rhs in ideal or []:
                print(f"    [{lhs}, {rhs}],")
            print("]")
            print("equations = [")
            for eq in equations:
                print(f"    {eq},")
            print("]")
        else:
            raise ValueError(f"Unsupported export format: {format}")

def test_xhash_systems_with_ideal():
    """Small test/demo with ideal generation."""
    total_start = perf_counter()
    print("XHash Polynomial System Generator with Ideal Relations")
    print("=" * 50)
    
    config = XHashConfig(p=17, state_size=3, rate=1, capacity=2, alpha=7)
    
    print(f"\nTesting configuration: {config}")
    system = XHashPolynomialSystem.from_config(config)
    
    # Test F0B-P3F pattern
    pattern = "F0B-P3F"
    print(f"\nPattern: {pattern}")
    
    ideal_relations = []
    equations, final_state, ideal_relations = system.generate_simple_system(
        pattern, num_rounds=2, block_size=3, power=7, ideal=ideal_relations
    )
    
    system.print_system_info(equations, f"Config-{pattern}", ideal_relations)
    
    print(f"\nIdeal Export:")
    system.export_system(equations, format="ideal", name=f"Config-{pattern}", ideal=ideal_relations)
    
    print(f"\nFirst few ideal relations:")
    for i, relation in enumerate(ideal_relations[:10]):
        lhs, rhs = relation
        print(f"  {i+1}: {lhs} = {rhs}")

    print(f"\nTotal time: {perf_counter() - total_start:.6f} s")

def run_demo():
    total_start = perf_counter()
    print("XHash Polynomial System Generator with Ideal Relations")
    print("=" * 50)
    
    config = XHashConfig(p=17, state_size=4, rate=2, capacity=2, alpha=3)
    
    print(f"\nTesting configuration: {config}")
    system = XHashPolynomialSystem.from_config(config)
    
    pattern = "F0B-P3F"
    print(f"\nPattern: {pattern}")
    
    ideal_relations = []
    equations, final_state, ideal_relations = system.generate_simple_system(
        pattern, num_rounds=2, block_size=2, power=3, ideal=ideal_relations
    )
    
    system.print_system_info(equations, f"Config-{pattern}", ideal_relations)
    summary = system.summarize(
        equations,
        ideal=ideal_relations,
        pattern=pattern,
        num_rounds=2,
        block_size=2,
        power=3,
    )
    print("\nSummary:")
    for key, value in summary.items():
        print(f"  {key}: {value}")
    
    print(f"\nIdeal Export:")
    system.export_system(equations, format="ideal", name=f"Config-{pattern}", ideal=ideal_relations)
    print(f"\nDrsolve Sage Interface Export:")
    system.export_system(equations, format="drsolve_sage_interface", name=f"Config-{pattern}", ideal=ideal_relations)
    print(f"\nPlain Equation Export:")
    system.export_system(equations, format="plain", name=f"Config-{pattern}")
    system.export_system(equations, format="sage", name=f"Config-{pattern}")
    print(f"\nTotal time: {perf_counter() - total_start:.6f} s")


if __name__ == "__main__":
    run_demo()
