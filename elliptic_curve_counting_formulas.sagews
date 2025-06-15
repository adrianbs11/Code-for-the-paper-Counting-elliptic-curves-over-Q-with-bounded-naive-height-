# This script contains functions used to compute both the exact and asymptotic formulas presented in the paper.
#
# In addition, it includes functions that perform exhaustive computer searches.
# These are used to verify the correctness of all the exact formulas established in the paper.
#
# All code is intended to support and illustrate the theoretical results and numerical data reported in the article.


# Auxiliary functions: these compute constants and intermediate values required in the main counting formulas of the paper.


def j_inv(A, B):
    # Computes the j-invariant of the elliptic curve E_{A,B} : y^2 = x^3 + A*x + B.
#
#     Parameters:
#         A (rational): Coefficient A of the elliptic curve.
#         B (rational): Coefficient B of the elliptic curve.
#
#     Returns:
#         rational: The j-invariant of E_{A,B}.
    return 1728*4*A^3 / (4*A^3 + 27*B^2)


def a_j(j):
    # Computes the number a(j) = 4(1728 - j)/27j defined in equation (2.12) and in Lemma 4.2 of the paper
#
#     Parameter:
#        j (rational).
#
#     Returns:
#         rational: The value of a(j).
    return 4*(1728 - j)/(27*j)


def N_a(a):
# Computes the rational number N_a defined for a ≠ 0 in equation (2.11) and also Proposition 3.1 of the paper.
#
#     Parameters:
#         a (rational): A nonzero rational number.
#
#     Returns:
#         rational: The corresponding value N_a.
    N = 1
    factorization = factor(a)
    for f in factorization:
        if f[1] >= 0:
            N *= f[0] ** ceil(f[1] / 2)
        else:
            N *= f[0] ** ceil(f[1] / 3)
    return N


def disc(A, B):
    # Computes the discriminant of the elliptic curve E_{A,B} : y^2 = x^3 + A*x + B.
# 
#     Parameters:
#         A (rational): Coefficient A of the elliptic curve.
#         B (rational): Coefficient B of the elliptic curve.
# 
#     Returns:
#         rational: The discriminant Δ(E_{A,B}) = -16*(4*A^3 + 27*B^2).
    return -16 * (4 * A^3 + 27 * B^2)


def prime_condition(A, B):
    # Checks the "primality condition" for the pair of integers (A, B).
# 
#     This condition requires that there is no prime number p such that:
#         - p^4 divides A, and simultaneously
#         - p^6 divides B.
# 
#     If A = 0, then the condition reduces to checking that no prime p satisfies p^6 | B.
# 
#     Parameters:
#         A (int): The coefficient A in the elliptic curve y^2 = x^3 + A*x + B.
#         B (int): The coefficient B in the elliptic curve y^2 = x^3 + A*x + B.
# 
#     Returns:
#         bool: True if the pair (A, B) satisfies the primality condition. False otherwise.
    if A == 0:
        # Only need to check if any prime p satisfies p^6 | B
        B_factors = list(ZZ(B).factor())
        for p, e in B_factors:
            if e >= 6:
                return False
        return True
    else:
        # Check if any prime p satisfies p^4 | A and p^6 | B
        A_factors = list(ZZ(A).factor())
        for p, e in A_factors:
            if e >= 4 and (p^6).divides(B):
                return False
        return True



def generalized_naive_height(A, B, alpha, beta):
    # Computes the generalized naive height H_{alpha, beta}(E_{A,B}) of the elliptic curve E_{A,B}.
# 
#     The generalized naive height is defined as:
#         H_{alpha, beta}(E_{A,B}) := max(alpha*|A|^3, beta*|B|^2)
# 
#     Parameters:
#         A (rational): Coefficient A of the elliptic curve.
#         B (rational): Coefficient B of the elliptic curve.
# 
#     Returns:
#         rational: The generalized naive height of E_{A,B}.
    return max(alpha * abs(A^3), beta * B^2)


def calibrated_naive_height(A, B):
    # Computes the calibrated naive height H^{cal}(E_{A,B}) of the elliptic curve E_{A,B}.
# 
#     The naive height is defined as:
#         H^{cal}(E_{A,B}) := max(4*|A|^3, 27*|B|^2)
# 
#     Parameters:
#         A (rational): Coefficient A of the elliptic curve.
#         B (rational): Coefficient B of the elliptic curve.
# 
#     Returns:
#         rational: The calibrated naive height of E_{A,B}.
    return max(4 * abs(A^3), 27 * B^2)


def uncalibrated_naive_height(A, B):
    # Computes the uncalibrated naive height H^{ncal}(E_{A,B}) of the elliptic curve E_{A,B}.
# 
#     The naive height is defined as:
#         H^{ncal}(E_{A,B}) := max(|A|^3, |B|^2)
# 
#     Parameters:
#         A (rational): Coefficient A of the elliptic curve.
#         B (rational): Coefficient B of the elliptic curve.
# 
#     Returns:
#         rational: The calibrated naive height of E_{A,B}.
    return max(abs(A^3), B^2)


# Generalized counting formulas for the height H_{alpha, beta}
# 
# The following functions compute the different exact and asymptotic formulas from Theorems 4.5, 4.6, 4.9 and 4.10
# 
# We also include formulas that carry exhaustive searches for curves that serve to check the correctnes of the values
# from our formulas


def generalized_exact_count_formula_for_all_elliptic_curves(X, alpha, beta):
    # Computes the exact count for the family of all elliptic curves from directly from the formula
#     appearing in Theorem 4.5
# 
#     Computes the quantity:
#         (2 * floor((X**(1/3)) / (alpha**(1/3))) + 1) *
#         (2 * floor((X**(1/2)) / (beta**(1/2))) + 1)
#         - 2 * floor(c(alpha, beta) * X**(1/6)) - 1
# 
#     where:
#         c(alpha, beta) = min(1 / (3**0.5 * alpha**(1/6)), 1 / (2**(1/3) * beta**(1/6)))
# 
#     Parameters:
#         X (float): A positive real number.
#         alpha (float): A positive real number.
#         beta (float): A positive real number.
# 
#     Returns:
#         int: The computed integer value.
    if X <= 0 or alpha <= 0 or beta <= 0:
        raise ValueError("X, alpha, and beta must all be positive real numbers.")

    term1 = 2 * floor((X ** (1/3)) / (alpha ** (1/3))) + 1
    term2 = 2 * floor((X ** (1/2)) / (beta ** (1/2))) + 1

    c_alpha_beta = min(
        1 / (math.sqrt(3) * alpha ** (1/6)),
        1 / ((2 ** (1/3)) * beta ** (1/6))
    )

    term3 = 2 * floor(c_alpha_beta * X ** (1/6))

    result = term1 * term2 - term3 - 1
    return result



def generalized_asymptotic_estimate_for_all_isomorphism_classes(X, alpha, beta):
    # Computes the main term in the asymptotic formula for the family of
#     representatives for isomorphism classes of elliptic curves from Theorem 4.5
# 
#     Computes the quantity:
#         (4 / (alpha**(1/3) * beta**(1/2) * zeta(10))) * X**(5/6)
# 
#     where zeta(10) is the value of the Riemann zeta function at 10.
# 
#     Parameters:
#         X (float): A positive real number.
#         alpha (float): A positive real number.
#         beta (float): A positive real number.
# 
#     Returns:
#         float: The computed value.
    if X <= 0 or alpha <= 0 or beta <= 0:
        raise ValueError("X, alpha, and beta must all be positive real numbers.")

    zeta_10 = float(zeta(10))
    denominator = alpha ** (1/3) * beta ** (1/2) * zeta_10
    result = (4 / denominator) * X ** (5/6)
    return numerical_approx( result )


def generalized_count_formula_elliptic_curves_with_j(X, j, alpha, beta):
    # Computes the exact number of elliptic curves E_{A,B} with j-invariant equal to j
#     and generalized naive height H_{alpha, beta}(E_{A,B}) <= X by using the formulas
#     from Theorem 4.9.
# 
#     The function considers three cases:
#       - If j = 0: returns 2 * floor(X^{1/2} / beta^{1/2})
#       - If j = 1728: returns 2 * floor(X^{1/3} / alpha^{1/3})
#       - Otherwise: returns 2 * floor(c(j; H_{alpha,beta}) * X^{1/6})
# 
#     where:
#       c(j; H_{alpha,beta}) = min{
#         |a(j)|^{1/2} / (N_{a(j)} * alpha^{1/6}),
#         |a(j)|^{1/3} / (N_{a(j)} * beta^{1/6})
#       }
# 
#     Parameters:
#         X (float): A nonnegative real number.
#         j (rational): The j-invariant.
#         alpha (float): A strictly positive real number.
#         beta (float): A strictly positive real number.
# 
#     Returns:
#         int: The number of elliptic curves with j-invariant j and height ≤ X.
    if X < 0 or alpha <= 0 or beta <= 0:
        raise ValueError("X must be nonnegative, alpha and beta must be strictly positive.")

    if j == 0:
        return 2 * floor((X ** 0.5) / (beta ** 0.5))

    elif j == 1728:
        return 2 * floor((X ** (1/3)) / (alpha ** (1/3)))

    else:
        a = a_j(j)
        N = N_a(a)
        abs_a = abs(a)

        term1 = abs_a ** (1/2) / (N * alpha ** (1/6))
        term2 = abs_a ** (1/3) / (N * beta ** (1/6))
        c = min(term1, term2)

        return 2 * floor(c * X ** (1/6))


def generalized_asymptotic_formula_elliptic_curves_with_j(X, j, alpha, beta):
    # Computes the main term in the asymptotic formula for the number of isomorphism class
#     representatives of elliptic curves E_{A,B} with j-invariant equal to j and
#     generalized naive height H_{alpha, beta}(E_{A,B}) <= X by using the formulas
#     from Theorem 4.9.
# 
#     The function considers three cases:
#       - If j = 0: returns (2 / (beta^{1/2} * zeta(6))) * X^{1/2}
#       - If j = 1728: returns (2 / (alpha^{1/3} * zeta(4))) * X^{1/3}
#       - Otherwise: returns (2 / zeta(2)) * c(j; H_{alpha,beta}) * X^{1/6}
# 
#     where:
#       c(j; H_{alpha,beta}) = min{
#         |a(j)|^{1/2} / (N_{a(j)} * alpha^{1/6}),
#         |a(j)|^{1/3} / (N_{a(j)} * beta^{1/6})
#       }
# 
#     Parameters:
#         X (float): A nonnegative real number.
#         j (rational): The j-invariant.
#         alpha (float): A strictly positive real number.
#         beta (float): A strictly positive real number.
# 
#     Returns:
#         float: The value of the main term in the asymptotic estimate.
    if X < 0 or alpha <= 0 or beta <= 0:
        raise ValueError("X must be nonnegative, alpha and beta must be strictly positive.")

    if j == 0:
        zeta6 = float(zeta(6))
        return (2 / (beta ** 0.5 * zeta6)) * X ** 0.5

    elif j == 1728:
        zeta4 = float(zeta(4))
        return (2 / (alpha ** (1/3) * zeta4)) * X ** (1/3)

    else:
        a = a_j(j)
        N = N_a(a)
        abs_a = abs(a)

        term1 = abs_a ** (1/2) / (N * alpha ** (1/6))
        term2 = abs_a ** (1/3) / (N * beta ** (1/6))
        c = min(term1, term2)

        zeta2 = float(zeta(2))
        return numerical_approx( (2 / zeta2) * c * X ** (1/6) )


def generalized_exact_count_formula_cm_curves(X, alpha, beta):
    # Computes the exact number of elliptic curves E_{A,B} with complex multiplication (CM)
#     and generalized naive height H_{alpha, beta}(E_{A,B}) <= X.
# 
#     The count is computed using the formula from Theorem 4.10 of the paper:
# 
#         2 * floor(X^{1/2} / beta^{1/2}) +
#         2 * floor(X^{1/3} / alpha^{1/3}) +
#         sum_{j in J_cm, j ≠ 0,1728} 2 * floor(c(j; H_{alpha,beta}) * X^{1/6})
# 
#     where:
#         J_cm = {
#             0, 1728, -3375, 8000, -32768, 54000, 287496,
#             -12288000, 16581375, -884736, -884736000,
#             -147197952000, -262537412640768000
#         }
# 
#     and:
#         c(j; H_{alpha,beta}) = min(
#             |a(j)|^{1/2} / (N_{a(j)} * alpha^{1/6}),
#             |a(j)|^{1/3} / (N_{a(j)} * beta^{1/6})
#         )
# 
#     Parameters:
#         X (float): A nonnegative real number.
#         alpha (float): A strictly positive real number.
#         beta (float): A strictly positive real number.
# 
#     Returns:
#         int: The number of CM elliptic curves with height <= X.
    if X < 0 or alpha <= 0 or beta <= 0:
        raise ValueError("X must be nonnegative, alpha and beta must be strictly positive.")

    # List of CM j-invariants
    j_cm = [
        0, 1728, -3375, 8000, -32768, 54000, 287496,
        -12288000, 16581375, -884736, -884736000,
        -147197952000, -262537412640768000
    ]

    total = 0

    # j = 0 term
    total += 2 * math.floor(X ** 0.5 / beta ** 0.5)

    # j = 1728 term
    total += 2 * math.floor(X ** (1/3) / alpha ** (1/3))

    # j ≠ 0, 1728 terms
    for j in j_cm:
        if j == 0 or j == 1728:
            continue
        a = a_j(j)
        N = N_a(a)
        abs_a = abs(a)

        term1 = abs_a ** (1/2) / (N * alpha ** (1/6))
        term2 = abs_a ** (1/3) / (N * beta ** (1/6))
        c = min(term1, term2)

        total += 2 * math.floor(c * X ** (1/6))

    return total


def K_constant_cm(alpha, beta):
    # Computes the constant K(alpha, beta) defined in Theorem 4.10 as:
# 
#         K(alpha, beta) := sum_{j in J_cm, j ≠ 0, 1728} c(j; H_{alpha,beta})
# 
#     where:
#         c(j; H_{alpha,beta}) = min(
#             |a(j)|^{1/2} / (N_{a(j)} * alpha^{1/6}),
#             |a(j)|^{1/3} / (N_{a(j)} * beta^{1/6})
#         )
# 
#     Parameters:
#         alpha (float): A strictly positive real number.
#         beta (float): A strictly positive real number.
# 
#     Returns:
#         float: The value of K(alpha, beta).
    if alpha <= 0 or beta <= 0:
        raise ValueError("alpha and beta must be strictly positive.")

    j_cm = [
        -3375, 8000, -32768, 54000, 287496,
        -12288000, 16581375, -884736, -884736000,
        -147197952000, -262537412640768000
    ]

    total = 0.0
    for j in j_cm:
        a = a_j(j)
        N = N_a(a)
        abs_a = abs(a)

        term1 = abs_a ** (1/2) / (N * alpha ** (1/6))
        term2 = abs_a ** (1/3) / (N * beta ** (1/6))
        c = min(term1, term2)

        total += c

    return numerical_approx(total)



def generalized_asymptotic_estimate_isomorphism_classes_cm_curves(X, alpha, beta):
    # Computes an asymptotic estimate for the number of isomorphism class representatives
#     of elliptic curves E_{A,B} with complex multiplication (CM) and
#     generalized naive height H_{alpha, beta}(E_{A,B}) <= X.
# 
#     The estimate is based on the main terms from Theorem 4.10:
# 
#         (2 / (beta^{1/2} * zeta(6))) * X^{1/2} +
#         (2 / (alpha^{1/3} * zeta(4))) * X^{1/3} +
#         (2 / zeta(2)) * K(alpha, beta) * X^{1/6}
# 
#     where K(alpha, beta) is computed as in the function K_constant_cm.
# 
#     Parameters:
#         X (float): A nonnegative real number.
#         alpha (float): A strictly positive real number.
#         beta (float): A strictly positive real number.
# 
#     Returns:
#         float: The value of the asymptotic estimate.
    if X < 0 or alpha <= 0 or beta <= 0:
        raise ValueError("X must be nonnegative, alpha and beta must be strictly positive.")

    zeta2 = float(zeta(2))
    zeta4 = float(zeta(4))
    zeta6 = float(zeta(6))

    term1 = (2 / (beta ** 0.5 * zeta6)) * X ** 0.5
    term2 = (2 / (alpha ** (1/3) * zeta4)) * X ** (1/3)
    K = K_constant_cm(alpha, beta)
    term3 = (2 / zeta2) * K * X ** (1/6)

    return numerical_approx(term1 + term2 + term3)


# The following functions perform exhaustive searches over all elliptic curves in the different
# families considered. They can be used to verify the correctness of the exact formulas
# presented in the paper by confirming that they yield matching values.


def generalized_exact_count_all_curves_by_exhaustive_search(X, alpha, beta):
    # Computes the number of integer pairs (A, B) such that:
# 
#         - |A| <= X^{1/3} / alpha^{1/3}
#         - |B| <= X^{1/2} / beta^{1/2}
#         - The discriminant of the curve E_{A,B} is nonzero: disc(A, B) != 0
# 
#     This is done by exhaustive search and corresponds to counting all
#     elliptic curves E_{A,B} over Z with generalized naive height
#     H_{alpha, beta}(E_{A,B}) <= X.
# 
#     Parameters:
#         X (float): A nonnegative real number.
#         alpha (float): A strictly positive real number.
#         beta (float): A strictly positive real number.
# 
#     Returns:
#         int: The total number of valid (A, B) pairs.
    if X < 0 or alpha <= 0 or beta <= 0:
        raise ValueError("X must be nonnegative, and alpha and beta must be strictly positive.")

    # Compute bounds for A and B
    bound_A = floor(X ** (1/3) / alpha ** (1/3))
    bound_B = floor(X ** (1/2) / beta ** (1/2))

    count = 0
    for A in range(-bound_A, bound_A + 1):
        for B in range(-bound_B, bound_B + 1):
            if disc(A, B) != 0:
                count += 1

    return count


def generalized_exact_count_all_isomorphism_classes_by_exhaustive_search(X, alpha, beta):
    # Counts the number of integer pairs (A, B) in Z^2 such that:
#         - |A| ≤ X^{1/3} / alpha^{1/3},
#         - |B| ≤ X^{1/2} / beta^{1/2},
#         - The discriminant Δ(E_{A,B}) ≠ 0,
#         - The pair (A, B) satisfies the primality condition:
#               There is no prime p such that p^4 | A and p^6 | B.
# 
#     This corresponds to counting isomorphism class representatives of elliptic curves
#     under the generalized naive height H_{alpha, beta}(E_{A,B}) ≤ X.
# 
#     Parameters:
#         X (float): A positive real number.
#         alpha (float): A positive real number.
#         beta (float): A positive real number.
# 
#     Returns:
#         int: The number of valid pairs (A, B) satisfying all conditions.
    if X <= 0 or alpha <= 0 or beta <= 0:
        raise ValueError("X, alpha, and beta must be positive real numbers.")

    bound_A = floor(X ** (1/3) / alpha ** (1/3))
    bound_B = floor(X ** (1/2) / beta ** (1/2))

    count = 0
    for A in range(-bound_A, bound_A + 1):
        for B in range(-bound_B, bound_B + 1):
            if disc(A, B) != 0 and prime_condition(A, B):
                count += 1

    return count



def generalized_exact_count_curves_with_j_by_exhaustive_search(X, j, alpha, beta):
    # Computes the number of elliptic curves E_{A,B} over Z with:
#         - Generalized naive height H_{alpha, beta}(E_{A,B}) <= X
#         - Nonzero discriminant
#         - j-invariant equal to the given rational number j
# 
#     Parameters:
#         X (float): A nonnegative real number.
#         j (rational): The prescribed j-invariant.
#         alpha (float): A strictly positive real number.
#         beta (float): A strictly positive real number.
# 
#     Returns:
#         int: The number of curves E_{A,B} with j(E_{A,B}) == j and height <= X.
    if X < 0 or alpha <= 0 or beta <= 0:
        raise ValueError("X must be nonnegative, alpha and beta must be strictly positive.")

    bound_A = floor(X ** (1/3) / alpha ** (1/3))
    bound_B = floor(X ** (1/2) / beta ** (1/2))

    count = 0
    for A in range(-bound_A, bound_A + 1):
        for B in range(-bound_B, bound_B + 1):
            if disc(A, B) != 0 and j_inv(A, B) == j:
                count += 1

    return count


def generalized_exact_count_cm_curves_by_exhaustive_search(X, alpha, beta):
    # Computes the number of CM elliptic curves E_{A,B} over Z such that:
#         - H_{alpha, beta}(E_{A,B}) <= X
#         - disc(A, B) != 0
#         - j-invariant of E_{A,B} lies in the list of 13 CM j-invariants
# 
#     Parameters:
#         X (float): A nonnegative real number.
#         alpha (float): A strictly positive real number.
#         beta (float): A strictly positive real number.
# 
#     Returns:
#         int: The number of CM elliptic curves with height <= X.
    if X < 0 or alpha <= 0 or beta <= 0:
        raise ValueError("X must be nonnegative, alpha and beta must be strictly positive.")

    j_cm = {
        0, 1728, -3375, 8000, -32768, 54000, 287496,
        -12288000, 16581375, -884736, -884736000,
        -147197952000, -262537412640768000
    }

    bound_A = floor(X ** (1/3) / alpha ** (1/3))
    bound_B = floor(X ** (1/2) / beta ** (1/2))

    count = 0
    for A in range(-bound_A, bound_A + 1):
        for B in range(-bound_B, bound_B + 1):
            if disc(A, B) != 0:
                j_val = j_inv(A, B)
                if j_val in j_cm:
                    count += 1

    return count


def generalized_exact_count_isomorphism_classes_cm_curves_by_exhaustive_search(X, alpha, beta):
    # Computes the number of CM elliptic curves E_{A,B} over Z such that:
#         - H_{alpha, beta}(E_{A,B}) <= X
#         - disc(A, B) != 0
#         - j-invariant of E_{A,B} lies in the list of 13 CM j-invariants
# 
#     Parameters:
#         X (float): A nonnegative real number.
#         alpha (float): A strictly positive real number.
#         beta (float): A strictly positive real number.
# 
#     Returns:
#         int: The number of CM elliptic curves with height <= X.
    if X < 0 or alpha <= 0 or beta <= 0:
        raise ValueError("X must be nonnegative, alpha and beta must be strictly positive.")

    j_cm = {
        0, 1728, -3375, 8000, -32768, 54000, 287496,
        -12288000, 16581375, -884736, -884736000,
        -147197952000, -262537412640768000
    }

    bound_A = floor(X ** (1/3) / alpha ** (1/3))
    bound_B = floor(X ** (1/2) / beta ** (1/2))

    count = 0
    for A in range(-bound_A, bound_A + 1):
        for B in range(-bound_B, bound_B + 1):
            if disc(A, B) != 0:
                j_val = j_inv(A, B)
                if j_val in j_cm and prime_condition(A, B):
                    count += 1

    return count


# The following functions deal with the parametrization for curves with a
# fixed j-invariant from Theorems 2.12 and 4.6.

def parametrized_point_with_j(j, m):
    # Computes the coefficients (A, B) of the elliptic curve E_{A,B} = E_{A(j,m), B(j,m)} defined
#     by the parametrization in Theorems 2.12 and 4.6 of the paper, for a fixed j-invariant j
#     not equal to 0 or 1728 and parameter m in Z.
# 
#     This is used to generate elliptic curves with a prescribed j-invariant j.
# 
#     Parameters:
#         j (rational): The j-invariant of the elliptic curve.
#         m (integer): A nonzero integer parameter.
# 
#     Returns:
#         list: A list [A(j,m), B(j,m)] of the coefficients of the elliptic curve.
    a = a_j(j)
    Ajm = N_a(a)^2 * m^2 / a
    Bjm = N_a(a)^3 * m^3 / a
    return [Ajm, Bjm]


def parametrized_curve_with_j(j, m):
    # Returns the elliptic curve E_{A(j,m), B(j,m)} with prescribed j-invariant j and parameter m,
#     as described in Theorems 2.12 and 4.6.
# 
#     This curve is of the form:
#         y^2 = x^3 + A(j,m)x + B(j,m)
# 
#     Parameters:
#         j (rational): The j-invariant of the curve.
#         m (rational): A nonzero rational parameter.
# 
#     Returns:
#         EllipticCurve: The curve defined over QQ with Weierstrass coefficients [0, 0, 0, A(j,m), B(j,m)].
    A, B = parametrized_point_with_j(j, m)
    return EllipticCurve(QQ, [0, 0, 0, A, B])


def curve_with_minimal_height(j):
    # Returns the elliptic curve E_{A,B} over QQ with j-invariant j that has minimal
#     calibrated naive height.
# 
#     Special cases:
#         - If j = 0, the curve y^2 = x^3 + 1 has minimal height.
#         - If j = 1728, the curve y^2 = x^3 + x has minimal height.
# 
#     General case:
#         For j ≠ 0, 1728, the curve is constructed from the parametrized point (A(j,1), B(j,1))
#         as in Theorems 2.12 and 4.6 of the paper.
# 
#     Parameter:
#         j (rational): The j-invariant of the desired elliptic curve.
# 
#     Returns:
#         EllipticCurve: An elliptic curve over QQ with j-invariant j and minimal height.
    if j == 0:
        return EllipticCurve(QQ, [0, 1])
    elif j == 1728:
        return EllipticCurve(QQ, [1, 0])
    else:
        A, B = parametrized_point_with_j(j, 1)
        return EllipticCurve(QQ, [A, B])



def minimal_generalized_naive_height_for_j(j, alpha, beta):
    # Computes the generalized naive height H_{alpha, beta}(E_{A,B}) of the elliptic curve
#     with j-invariant j and minimal calibrated naive height, using the parametrization
#     E_{A(j,1), B(j,1)} from Theorems 2.12 and 4.6.
# 
#     Special cases:
#         - If j = 0, the minimal curve is y^2 = x^3 + 1 and the height is beta.
#         - If j = 1728, the minimal curve is y^2 = x^3 + x and the height is alpha.
# 
#     General case:
#         For j ≠ 0, 1728, the curve is constructed from the parametrized point (A(j,1), B(j,1)),
#         and the height is computed as:
#             H_{alpha, beta}(E_{A,B}) := max(alpha * |A|^3, beta * |B|^2)
# 
#     Parameters:
#         j (rational): The j-invariant of the elliptic curve.
#         alpha (float or real): Positive scaling constant for the A-term in the height.
#         beta (float or real): Positive scaling constant for the B-term in the height.
# 
#     Returns:
#         real: The value of the generalized naive height H_{alpha, beta}(E_{A,B}).
    if j == 0:
        return beta
    elif j == 1728:
        return alpha
    else:
        A, B = parametrized_point_with_j(j, 1)
        return max(alpha * abs(A)^3, beta * B^2)











