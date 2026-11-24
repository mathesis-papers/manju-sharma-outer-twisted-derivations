# Outer twisted derivations of cyclotomic integer rings

This package contains revision 2 of **Manju–Sharma's outer-derivation conjecture: a correction from Noferini–Williams Smith forms**.

Creators, in order: Haobo Ma, ChronoAI Pte Ltd, Singapore; Wenlin Zhang, National University of Singapore and ChronoAI Pte Ltd, Singapore. License: Creative Commons Attribution 4.0 International, `cc-by-4.0`.

## Mathematical content

For distinct automorphisms `sigma(zeta)=zeta^u` and `tau(zeta)=zeta^v` on `O=Z[zeta_n]`, put `h=n/gcd(n,|v-u|)`. The outer module is `O/(zeta^v-zeta^u)O`. Its additive group is `(Z/ell Z)^(phi(n)/phi(h))` if `h` is a power of a prime `ell`, and zero otherwise. The arithmetic classification is a cyclotomic specialization of Noferini–Williams's Smith-form theorems, with an elementary norm proof in the note.

For `n=p^k`, the denominator accompanying the unscaled classical adjugate is `p^(p^e)`, with `e=v_p(|v-u|)`; the least positive integer clearing the inverse matrix is `p`. The Gaussian counterexample disproves the printed criteria, while universal prime-power outer existence remains true. The coverage table includes all six conjectures in Sections 4–5 of arXiv:2412.03507v2, Conjectures 3.12–3.13 of arXiv:2412.03500v3, and the related quadratic sufficient-condition boundary.

## Build the note

All paths are relative to the extracted directory. The complete main source is `note.tex`; its bibliography is included in that file.

Use Tectonic 0.17.0:

```sh
tectonic --untrusted note.tex
```

Alternatively, with a LaTeX installation providing the packages declared in `note.tex`:

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

No BibTeX step is needed because the bibliography is inline. `metadata.json` contains the title, abstract, creator order, license, revision, problem identifiers and a pinned formal-source link.

## Exact finite checks

Use Python 3.8 or newer; no third-party packages are needed:

```sh
python3 verify.py --output verification-rebuilt.json
cmp verification.json verification-rebuilt.json
sha256sum -c MANIFEST.sha256
```

The checker uses integer arithmetic, exact rational fractions, fraction-free determinants and modular ranks. It checks both quadratic laws on all four ordered basis pairs, the Gaussian adjugate vectors, the conductor-nine valuation-zero example, and every unordered distinct unit pair for the conductors

`3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 24, 25, 27, 30, 32, 36, 45, 48, 60`.

There are exactly 1,338 such pairs; the largest degree is 24. For each prime-power order `h`, a determinant `ell^d` together with nullity `d` modulo `ell` forces exactly `d` nontrivial Smith factors, all equal to `ell`. A determinant of one gives the zero quotient. These finite checks complement the general written proofs; the finite list is not a classification by enumeration.

## Formal statement coverage

At commit `1ff680e6cc902e88cfbba1b684826388e2ac8d37`:

- [ManjuSharmaConjectureFourFour.lean](https://github.com/the-omega-institute/trureturing/blob/1ff680e6cc902e88cfbba1b684826388e2ac8d37/D5/S3/Arith/TwistedDerivations/ManjuSharmaConjectureFourFour.lean) contains `gaussianEquiv` for the pointed Gaussian/cyclotomic integer-ring identification; `derivation_twisted` for the full twisted law; `derivation_not_inner` for non-innerness; `printed_test` for the integral printed vector; and `result` for negation of the printed Conjecture 4.4.
- [ManjuSharmaConjectureFiveTwo.lean](https://github.com/the-omega-institute/trureturing/blob/1ff680e6cc902e88cfbba1b684826388e2ac8d37/D5/S3/Arith/TwistedDerivations/ManjuSharmaConjectureFiveTwo.lean) defines `outerClass` in the actual derivation quotient; `result` negates the printed conjunction using the Gaussian class. The separate universal existence assertion remains true.

The arbitrary-conductor quotient and group classification, the related quadratic criterion and the all-pairs conductor classification have written proofs in `note.tex`. `claims-and-sources.json` records the mathematical source/statement mapping. The source archive includes no dependencies beyond the files listed in `MANIFEST.sha256`.
