Dear llncs user,

The files in this directory belong to the LaTeX2e package
for Springer's Lecture Notes in Computer Science (LNCS) and
other proceedings book series.

It consists of the following files:

  readme.txt         this file

  history.txt        the version history of the package

  llncs.cls          the LaTeX2e document class

  samplepaper.tex    a sample paper
  fig1.eps           a figure used in the sample paper

  llncsdoc.pdf       the documentation of the class (PDF version)

  splncs04.bst       current LNCS BibTeX style with alphabetic sorting


Upload your paper's final version
This page is intended for authors to upload a final version of a paper submitted to an IACR conference or transactions. In the form below, you must submit all LaTeX source files, as well as all graphics and BibTeX files and style files that are required to produce the final PDF file. You must also submit the PDF file that was created from these source LaTeX file(s)

The proceedings will be published in Springer's Lecture Notes in Computer Science (LNCS). Authors must prepare their camera-ready version in LaTeX2e. You should obtain the LNCS LaTeX2e class file (llncs.cls), which can be obtained off Springer's site at:

http://www.springer.de/comp/lncs/authors.html

Click on the link "Information for LNCS Authors", which will take you to a page with a zip file that includes the class file and documentation on how to use that class file. You should ignore the information about the Springer's copyright form and you must not sign Springer's LNCS copyright form.

The page limit for your camera-ready paper is 30 pages total, including all appendices but not including bibliography.

ePrint:
@misc{cryptoeprint:2026/1969,
      author = {Haohai Suo and Jiamin Cui},
      title = {Efficient Polynomial System Solving via Dixon Resultants: Applications to {AO} Primitives},
      howpublished = {Cryptology {ePrint} Archive, Paper 2026/1969},
      year = {2026},
      url = {https://eprint.iacr.org/2026/1969}
}



Recommended submission style
Electronic submissions to Asiacrypt 2026 must be in Portable Document Format (PDF) and follow the standard LNCS guidelines. The submission should preferably use Type 1 fonts rather than Type 3 fonts, which usually look fuzzy and ugly when viewed on screen.

The following procedure is recommended for generating submissions.

Preparing the LaTeX file
To follow the standard LNCS guidelines, you obtain the llncs package and use \documentclass{llncs} at the beginning of your LaTeX file. You should not use any other command to set the margin and/or change the font. This LaTeX style will be used for the preproceedings.

Generating a PDF file with pdflatex
After using the above declaration, assuming that your paper is stored in the file paper.tex, it suffices to type the command: $ pdflatex paper. This generates a file paper.pdf ready for submission.

There are other, more complex, procedures to generate such PDF files. These alternative procedures are not recommended. If, for some reason, an alternative procedure is used, the resulting PDF file should be verified using the following commands:

$ pdfinfo paper.pdf
$ pdffonts paper.pdf
The above two commands respectively print general information (including paper size) and font information.

Including graphics
To insert graphics into your PDF file, there are two different options:

generate the graphics using a text description within LaTeX; OR
include an externally generated graphics file
For the first option, authors should consider the PGF package. It can be used by including \usepackage{pgf} in the LaTeX file.

To use externally generated graphics, a convenient method relies on the following package: \usepackage{graphicx,color}. With this package, a PDF file (drawing.pdf) can be included using \includegraphics{drawing}.

Authors should make sure that their externally generated graphics PDF files have a correct bounding box specification. A set of various cryptography related graphics source codes can be found at https://www.iacr.org/authors/tikz/.