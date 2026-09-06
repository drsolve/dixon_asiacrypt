# Bundled TeX dependencies

These unmodified dependencies were used to compile the revised manuscript.
The generated `.sty` files and the LNCS class/bibliography style are placed in
the parent directory so that the standard build command works locally.

| Package | Files used | Source | License |
|---|---|---|---|
| Springer LLNCS, 2025/02/25 v2.26 | `llncs.cls`, `splncs04.bst` | [CTAN distribution](https://ctan.math.illinois.edu/macros/latex/contrib/llncs.zip) | CC BY 4.0 or later; copyright Springer |
| autobreak | `autobreak.sty`, generated with `latex autobreak.ins` | [CTAN distribution](https://mirrors.ctan.org/macros/latex/contrib/autobreak.zip) | See included source distribution |
| algorithms | `algorithm.sty`, `algorithmic.sty`, generated with `latex algorithms.ins` | [CTAN distribution](https://ctan.math.illinois.edu/macros/latex/contrib/algorithms.zip) | LGPL; see `COPYING` in the included archive |

The original distribution archives, including source and license notices,
are included in this directory. No dependency source was modified.
Other packages used in `main.tex` come from the installed TeX distribution.
