# English speaker notes

30 main slides and 4 backup slides. Suggested duration: 30–40 minutes.

Suggested pacing: slides 1–6, 5 minutes; slides 7–16, 12–15 minutes; slides 17–23, 7–8 minutes; slides 24–30, 7–9 minutes.

## Slide 1

Plan for 30--40 minutes for the main 30 slides; the four backup slides are for discussion. This talk is based on the uploaded paper. We focus on what Dixon elimination computes, why it is valid, and how the paper improves its analysis and implementation.

## Slide 2

The first half builds the mechanism carefully through one example. The second half explains the paper's contributions and the evidence about when Dixon is useful.

## Slide 3

Fixing some input and output coordinates produces a constrained-input constrained-output problem. Power maps and mixing layers give polynomial relations in the input and internal states. The diagram suppresses details of any particular permutation.

## Slide 4

Here n means the number of equations, and n minus one is the number of eliminated variables. With one retained parameter the original system has equally many variables and equations. More retained parameters are useful for intermediate elimination.

## Slide 5

Gröbner bases also support elimination. Dixon offers a different construction, not an exclusive capability. Its cost profile and treatment of retained variables as parameters can be attractive in certain regimes.

## Slide 6

We will unpack these steps on a system with a square, generically nonsingular matrix. After understanding that clean case, we return to rectangular and rank-deficient matrices.

## Slide 7

Elementary substitution also solves this system, giving an independent check. We use the full Dixon procedure so every matrix entry and every factor is visible. The equation degrees are one, one, and two.

## Slide 8

Each column is one input polynomial. The rows evaluate the same polynomials under progressively more replacements. This pattern forces predictable factors in the determinant.

## Slide 9

The signs follow from later row minus earlier row, divided by original variable minus auxiliary variable. The third row uses original C2. We never divide numerically by zero; the difference quotients are polynomials.

## Slide 10

The coefficient of x times u is two, giving entry (2,2). The coefficient of the constant row monomial times v is one minus z. The resulting matrix has entries in Q[z].

## Slide 11

The matrix is nonsingular over Q(z) although singular at the relevant parameter values. This closes the calculation. Back-substitution uses the first two equations and then checks the third.

## Slide 12

This is a polynomial identity, not an assumption about numerical auxiliary values. The factor theorem applies because replacing x by u makes the determinant identically zero. Since x minus u and y minus v are distinct relatively prime polynomial factors, their product divides the determinant.

## Slide 13

The cancellation is in the polynomial ring of auxiliary variables. Even if u later equals x0 at some numerical point, x0 minus u is not the zero polynomial. Once Delta is the zero polynomial, every auxiliary-monomial coefficient vanishes. This creates the nonzero left-kernel vector.

## Slide 14

This multiplication is the key intuition. A nonlinear solution supplies values for the row monomials, and those values annihilate the coefficient matrix. The determinant detects the resulting linear dependence.

## Slide 15

These equations are inconsistent for every t because their difference is one. The determinant nonetheless has the root zero, where the x-dependent leading terms disappear. This illustrates the false converse without relying on a complicated extra-factor example.

## Slide 16

The proof for the generically nonsingular example does not justify arbitrary minors of singular matrices. Rank extraction is linear algebra, while validity as an elimination polynomial is an additional algebraic condition. Backup slide 31 gives a small example of the distinction.

## Slide 17

Every monomial has an exponent vector, which is a lattice point. The construction restricts sums of exponents. Reversing their order gives prefix constraints. Our worked example had degrees one, one, two and a three-by-three matrix; this separate quadratic example illustrates the generic bound five.

## Slide 18

The new support bound is attained generically under a characteristic hypothesis, but this is not a statement about full rank. The bar lengths compare mathematical bounds. They do not show a measured sevenfold matrix reduction or speedup.

## Slide 19

HNF means Hermite normal form. The costs count base-field arithmetic. The paper compares several determinant models because no method is uniformly best. We retain all symbolic stages, especially since the idealized exponent omega equals two does not generally justify ignoring construction.

## Slide 20

The row and column degree sums provide cheap proxies. The experiment measures how often the determinant degree exceeds the Bézout bound. It does not measure all extraneous factors, and we retain the one observed failure.

## Slide 21

Efficient polynomial-matrix arithmetic is a major reason the implementation improves over earlier Dixon software. We have not rerun the benchmarks or independently certified identical end-to-end timing boundaries across solvers.

## Slide 22

Read the curve qualitatively. In this tested family of dense random square systems, DRSolve has the favorable curve at high degree. This is evidence for a regime, not a guarantee about every structured system of a similar degree.

## Slide 23

This complementary result is important. Practical Gröbner-basis optimizations matter, and the paper does not establish Dixon as the best method for MQ systems. The two plots motivate method selection, not a universal ranking.

## Slide 24

We shift from generic systems to application patterns. The distinguishing structures are a manageable square system, a ladder of local state relations, and a triangular system that permits reduction. No complete round-function description is needed for this level of explanation.

## Slide 25

Both favorable and unfavorable instances are shown. Partial rounds can influence Dixon supports and ranks, as well as the degree reached by Gröbner-basis computation. Generic upper bounds do not fully predict the resulting practical costs.

## Slide 26

The two-variable restriction is essential. Two middle-state relations permit one final one-variable resultant. This does not claim that two equations can eliminate an arbitrarily large middle block. Output-degree assumptions remain part of the estimate.

## Slide 27

A normal form modulo triangular relations may still contain internal variables. Successive resultant steps are required. With several retained inputs, the final relations must then be merged, and intermediate degree growth may dominate.

## Slide 28

Separate the symbolic upper-bound comparison from measured time. The CICO-2 hybrid example is much slower than direct elimination despite the small generic theoretical margin. Constants, structure and intermediate work matter.

## Slide 29

Distinguish the theorem, the heuristic and the experiment. Future opportunities include exploiting rank and support structure, proving stronger factor removal methods, and improving models for specific AO systems.

## Slide 30

Return to the core mechanism: nonlinear solutions imply linear dependence, so a determinant supplies a necessary relation in the retained variables. The contribution makes this route more tightly estimated and more practical. The next four slides are optional backup material.

## Slide 31

In the nonsingular generic case a nonzero kernel forces determinant zero. A generically singular matrix already has a kernel everywhere, so more is needed. This example shows why arbitrary minor extraction is not a complete elimination correctness argument.

## Slide 32

No real-number ordering or geometry is needed. The same identities reduce to the finite field. Characteristic seventeen avoids degeneration of the factor four and permits division by two. The verification script exhausts all 4913 triples.

## Slide 33

This connects the hand calculation to the mixed-degree result. The smallest degree does not enter this nested-support formula. The two largest degree bounds, minus one, determine the prefix budgets.

## Slide 34

These are reading pointers. All experimental values are from the uploaded paper. The hand-worked examples are explanatory additions, not claimed as new research contributions.