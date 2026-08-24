# Code and Documentation for the ASIACRYPT 2026 Paper "Efficient Polynomial System Solving via Dixon Resultants: Applications to AO Primitives"

---
## 1. TeX Source and figures
This directory includes the complete LaTeX source and figures in the paper.

---

## 2. drsolve
The `drsolve` directory contains a Dixon resultant computation tool implemented in C.
- **Dependencies:** FLINT library (required), PML library (included)  
- **Compilation:** After installing the FLINT library, compile the tool by running:
```bash
./configure
make
sudo make install
```
This produces the executable used in our experiments.
A statically linked binary is also provided for convenience. In case of compatibility issues, please recompile from source.
Further details are provided in `drsolve/README.md`, or can be accessed by running:
```bash
./drsolve
```

---

## 3. Attack

The `Attack` directory contains experiment scripts and data files.

- **`.dr` files:** Input data files to be processed by the compiled Dixon tool.
  Example usage (run the following command inside `Attack/Poseidon`, assuming `drsolve` is available in your `PATH`):

```bash
drsolve poseidon3-0-3.dr
```

- **`.sage` interface:** The SageMath interface is provided via `drsolve_sage_interface.sage`, which serves as a wrapper for calling `drsolve`. Proper execution requires that the `drsolve` binary is available.

---

## 4. OtherTools
The `OtherTools` directory contains various experimental utilities:
- SageMath scripts for auxiliary testing and verification
- A Dixon resultant implementation based on Magma, used for cross-validation and comparison.

These tools are intended for research experimentation and validation.

---

## Additional note: Correspondence Between Paper Figures/Tables and Source Files

- **Section 3.3, Figure 1** and **Appendix E, Figure 6** → `Complexity_Comparision.ipynb`
- **Section 4.2, Figure 2** → `benchmark_DixonMagmaMsolve.ipynb`
- **Section 5, Figures 3/4/5** → `AO_Complexity.ipynb`
- **Appendix C, Table 6** → `Bound_Comparion.ipynb`

To verify that the Dixon matrix sizes reported in Table 7 match the theoretical predictions:

```bash
drsolve --test 1
```

For reproducing **Appendix G, Table 7**, run:

```bash
drsolve --test 2
```
---

## AI Tool Disclosure
We used generative AI tools, including Claude (Opus) and Codex, to assist in drafting portions of the codebase. All AI-assisted code was carefully reviewed, tested, and validated by the authors. The authors take full responsibility for the correctness of the implementation and all experimental results reported in the paper.