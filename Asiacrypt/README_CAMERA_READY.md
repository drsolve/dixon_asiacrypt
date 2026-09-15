# ASIACRYPT 2026 Camera-Ready Package

The camera-ready paper is built from `main.tex` with the LNCS class supplied in this directory.

```sh
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

The source uses T1-encoded Latin Modern, a Type 1 implementation of the Computer Modern design, to avoid bitmap Type 3 fonts. Verify the result with `pdfinfo main.pdf` and `pdffonts main.pdf` as recommended by IACR.

The generated `main.pdf` has 37 pages. The paper body, acknowledgements, and disclosure occupy 30 pages; the bibliography begins on page 30 and is excluded from the conference's 30-page limit. The appendices are available in the full version on the Cryptology ePrint Archive (Paper 2026/1969).

The upload archive `asiacrypt2026-camera-ready.zip` contains only the source files, required figures and styles, bibliography, and generated PDF needed for the proceedings submission.
