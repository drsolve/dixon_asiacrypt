# Dixon Resultants: English Beamer Talk

34 slides: 30 main slides, 4 backup slides. Suggested duration: 30–40 minutes.

## Files
- `dixon_english.tex`: editable 16:9 Beamer source, including speaker notes.
- `speaker_notes.md`: notes and suggested pacing.
- `assets/`: two original benchmark plots from the uploaded paper.
- `verify_example.py`: exact verification of the teaching examples, using the Python standard library only.

## Compile
Use XeLaTeX with Beamer, TikZ, fontspec, Latin Modern, and standard AMS packages. In Overleaf, choose XeLaTeX as the compiler.

```sh
xelatex dixon_english.tex
xelatex dixon_english.tex
```

Notes are hidden by default. To show them, replace `\setbeameroption{hide notes}` with `\setbeameroption{show notes}`.

## Check the worked examples

```sh
python3 verify_example.py
```

The script checks the cancellation identity, divided differences, the coefficient matrix, its determinant, kernel vectors, the rational solutions, every solution over F_17, and the false-converse example.

## Content
The main worked example computes the actual Dixon resultant of x+y-z, x-y-1, x²+y²-5. Separate slides explain exact polynomial division and why a common root produces a matrix kernel vector. A two-equation counterexample explains why resultant roots need verification.

Experimental figures and timings come from the uploaded full version of the paper; they were not rerun. The talk preserves the distinction between support size and rank, heuristic and proven degree bounds, and symbolic costs and measured times. The teaching examples are explanatory additions, not new research contributions.
