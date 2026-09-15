## Code and Documentation for the ASIACRYPT 2026 Paper "Efficient Polynomial System Solving via Dixon Resultants: Applications to AO Primitives"

---
### 1. TeX Source and figures
This directory includes the complete LaTeX source and figures in the paper.

---

### 2. drsolve
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

### 3. Attack

The `Attack` directory contains experiment scripts and data files.

- **`.dr` files:** Input data files to be processed by the compiled Dixon tool.
  Example usage (run the following command inside `Attack/Poseidon`, assuming `drsolve` is available in your `PATH`):

```bash
drsolve poseidon3-0-3.dr
```

- **`.sage` interface:** The SageMath interface is provided via `drsolve_sage_interface.sage`, which serves as a wrapper for calling `drsolve`. Proper execution requires that the `drsolve` binary is available.

---

### 4. OtherTools
`OtherTools` contains the code and recorded inputs used to reproduce the
paper's numerical tables and plots. Generated figures are written to
`figures/`; notebooks do not store cell output in the repository.

- `export_vector_figures.py` regenerates all 14 numerical plots as vector PDFs
  without rerunning the solver benchmarks.
- `Complexity_Comparision.ipynb` generates Figure 1 and every panel of Figure 6,
  including the max-over-stages model in Figure 6(d).
- `dixon_stage_models.py` is the shared Step 1/4 implementation used by Figures
  6(d), 7, and 8; `comp_step1&4.ipynb` is the corresponding plotting notebook.
- `AO_Complexity.ipynb` contains only the models and final plots for the
  Poseidon2, Vision, and XHash applications (Figures 3--5).
- `benchmark_DixonMagmaMsolve.ipynb`, `fermat_sage.py`, and the two recorded data
  files reproduce Figure 2 without overwriting the recorded measurements.
- `Bound_Comparion.ipynb` generates the values in Table 6.
- `dixon_magma/` contains the Magma implementation used for cross-validation.

With Python 3, NumPy, and Matplotlib installed, run from the repository root:

```bash
python3 OtherTools/export_vector_figures.py
```

The benchmark itself additionally requires SageMath, Magma, msolve, Fermat,
and a compiled `drsolve`; those programs are not invoked by the command above.

---

### Additional note: Correspondence Between Paper Figures/Tables and Source Files

- **Section 3.3, Figure 1** and **Appendix E, Figure 6** → `Complexity_Comparision.ipynb`
- **Section 4.2, Figure 2** → `benchmark_DixonMagmaMsolve.ipynb`, `fermat_sage.py`, and the recorded timing data
- **Section 5, Figures 3/4/5** → `AO_Complexity.ipynb`
- **Appendix C, Table 6** → `Bound_Comparion.ipynb`
- **Appendix H, Figures 7/8** → `comp_step1&4.ipynb` and `dixon_stage_models.py`

To verify that the Dixon matrix sizes reported in Table 6 match the theoretical predictions:

```bash
drsolve --test 1
```

For reproducing **Appendix G, Table 7**, run:

```bash
drsolve --test 2
```
---

### AI Tool Disclosure
We used generative AI tools, including Claude and ChatGPT, to assist in drafting portions of the codebase. All AI-assisted code was carefully reviewed, tested, and validated by the authors. The authors take full responsibility for the correctness of the implementation and all experimental results reported in the paper.
