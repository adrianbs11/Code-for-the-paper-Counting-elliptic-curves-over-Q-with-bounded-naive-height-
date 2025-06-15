# Code for the paper *Counting elliptic curves over Q with bounded naive height*

* Adrian Barquero-Sanchez
* Daniel Mora-Mora
    
This repository contains the SageMath worksheets with the code that was implemented for the computations that were used in the paper *Counting elliptic curves over Q with bounded naive height*

🔧 Auxiliary Functions

* j_inv(A, B): Computes the j-invariant of the elliptic curve y^2 = x^3 + Ax + B.

* a_j(j): Computes the rational constant a(j) appearing in the parametrization of curves with fixed j-invariant.

* N_a(a): Computes the auxiliary constant NaNa​ from the factorization of aa.

* disc(A, B): Computes the discriminant of the curve y^2= x^3 + Ax + B.

* prime_condition(A, B): Checks if the pair (A,B) satisfies the "primality condition" used to select canonical representatives.

* generalized_naive_height(A, B, alpha, beta): Computes the generalized naive height Hα,β of a curve.

* calibrated_naive_height(A, B): Computes the calibrated naive height of a curve.

* uncalibrated_naive_height(A, B): Computes the uncalibrated naive height Hncal of a curve.

📐 Exact and Asymptotic Formulas

* generalized_exact_count_formula_for_all_elliptic_curves(X, alpha, beta): Computes the exact count of elliptic curves with height <= X using Theorem 4.5.

* generalized_asymptotic_estimate_for_all_isomorphism_classes(X, alpha, beta): Computes the asymptotic estimate for the number of isomorphism classes using Theorem 4.5.

* generalized_count_formula_elliptic_curves_with_j(X, j, alpha, beta): Computes the exact count of curves with fixed j-invariant using Theorem 4.9.

* generalized_asymptotic_formula_elliptic_curves_with_j(X, j, alpha, beta): Computes the asymptotic count of curves with fixed j-invariant using Theorem 4.9.

* generalized_exact_count_formula_cm_curves(X, alpha, beta): Computes the exact number of CM elliptic curves with height ≤ X, using Theorem 4.10.

* K_constant_cm(alpha, beta): Computes the constant K(α,β) used in the CM asymptotic formula from Theorem 4.10.

* generalized_asymptotic_estimate_isomorphism_classes_cm_curves(X, alpha, beta): Computes the asymptotic estimate for the number of CM curves with height ≤ X, using Theorem 4.10.

🔍 Exhaustive Search Functions

* generalized_exact_count_all_curves_by_exhaustive_search(X, alpha, beta): Counts all (A,B) ∈ Z^2 with nonzero discriminant and height ≤X≤X.

* generalized_exact_count_all_isomorphism_classes_by_exhaustive_search(X, alpha, beta): Same as above, but only counts ℚ-isomorphism class representatives (satisfying the primality condition).

* generalized_exact_count_curves_with_j_by_exhaustive_search(X, j, alpha, beta): Counts curves with fixed j-invariant and bounded height by brute-force search.

* generalized_exact_count_cm_curves_by_exhaustive_search(X, alpha, beta): Counts CM curves with height ≤X≤X via exhaustive search.

* generalized_exact_count_isomorphism_classes_cm_curves_by_exhaustive_search(X, alpha, beta): Same as above, but only counts CM curves that satisfy the isomorphism class condition.

📐 Parametrization Functions

* parametrized_point_with_j(j, m): Returns the coefficients A(j,m),B(j,m) of a curve with fixed j-invariant using the parametrization from Theorem 2.12 or 4.6.

* parametrized_curve_with_j(j, m): Returns the elliptic curve E_{A(j,m),B(j,m)}​ with fixed j-invariant and integer parameter m.

* curve_with_minimal_height(j): Returns the elliptic curve with fixed j-invariant and minimal generalized naive height.

* minimal_generalized_naive_height_for_j(j, alpha, beta): Computes the minimal generalized naive height among curves with given j-invariant and fixed parameters α,β.
