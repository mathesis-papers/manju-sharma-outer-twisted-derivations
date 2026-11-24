import argparse
from fractions import Fraction
from functools import lru_cache
import itertools
import json
import math
from pathlib import Path


def trim(coefficients):
    result = list(coefficients)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def divide_exact(numerator, denominator):
    remainder = trim(numerator)
    divisor = trim(denominator)
    quotient = [0] * max(1, len(remainder) - len(divisor) + 1)
    while len(remainder) >= len(divisor) and remainder != [0]:
        offset = len(remainder) - len(divisor)
        assert remainder[-1] % divisor[-1] == 0
        factor = remainder[-1] // divisor[-1]
        quotient[offset] = factor
        for index, coefficient in enumerate(divisor):
            remainder[offset + index] -= factor * coefficient
        remainder = trim(remainder)
    assert remainder == [0]
    return tuple(trim(quotient))


@lru_cache(maxsize=128)
def cyclotomic(conductor):
    polynomial = [-1] + [0] * (conductor - 1) + [1]
    for divisor in range(1, conductor):
        if conductor % divisor == 0:
            polynomial = divide_exact(polynomial, cyclotomic(divisor))
    return tuple(polynomial)


def reduce_polynomial(coefficients, modulus):
    remainder = list(coefficients)
    degree = len(modulus) - 1
    remainder += [0] * max(0, degree - len(remainder))
    for index in range(len(remainder) - 1, degree - 1, -1):
        factor = remainder[index]
        for position, coefficient in enumerate(modulus):
            remainder[index - degree + position] -= factor * coefficient
    return tuple(remainder[:degree])


def difference_matrix(conductor, first, second):
    modulus = cyclotomic(conductor)
    degree = len(modulus) - 1
    columns = []
    for basis_index in range(degree):
        polynomial = [0] * (max(first, second) + basis_index + 1)
        polynomial[second + basis_index] += 1
        polynomial[first + basis_index] -= 1
        columns.append(reduce_polynomial(polynomial, modulus))
    return [list(row) for row in zip(*columns)]


def determinant(matrix):
    work = [row[:] for row in matrix]
    degree = len(work)
    previous_pivot = 1
    sign = 1
    for pivot_index in range(degree - 1):
        if work[pivot_index][pivot_index] == 0:
            replacement = next((index for index in range(pivot_index + 1, degree)
                                if work[index][pivot_index] != 0), None)
            if replacement is None:
                return 0
            work[pivot_index], work[replacement] = work[replacement], work[pivot_index]
            sign *= -1
        pivot = work[pivot_index][pivot_index]
        for row_index in range(pivot_index + 1, degree):
            for column_index in range(pivot_index + 1, degree):
                numerator = (work[row_index][column_index] * pivot
                             - work[row_index][pivot_index] * work[pivot_index][column_index])
                assert numerator % previous_pivot == 0
                work[row_index][column_index] = numerator // previous_pivot
            work[row_index][pivot_index] = 0
        previous_pivot = pivot
    return sign * work[-1][-1]


def prime_power_base(number):
    for candidate in range(2, math.isqrt(number) + 2):
        if number % candidate:
            continue
        remaining = number
        while remaining % candidate == 0:
            remaining //= candidate
        return candidate if remaining == 1 else None
    return number if number > 1 else None


def rank_modulo(matrix, prime):
    work = [[entry % prime for entry in row] for row in matrix]
    rank = 0
    for column in range(len(work)):
        pivot = next((row for row in range(rank, len(work)) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, prime)
        work[rank] = [entry * inverse % prime for entry in work[rank]]
        for row in range(rank + 1, len(work)):
            coefficient = work[row][column]
            work[row] = [(entry - coefficient * pivot_entry) % prime
                         for entry, pivot_entry in zip(work[row], work[rank])]
        rank += 1
    return rank


def solve_exact(matrix, coordinates):
    work = [[Fraction(entry) for entry in row] + [Fraction(value)]
            for row, value in zip(matrix, coordinates)]
    degree = len(work)
    for column in range(degree):
        pivot = next(row for row in range(column, degree) if work[row][column])
        work[column], work[pivot] = work[pivot], work[column]
        coefficient = work[column][column]
        work[column] = [entry / coefficient for entry in work[column]]
        for row in range(degree):
            if row != column:
                coefficient = work[row][column]
                work[row] = [entry - coefficient * pivot_entry
                             for entry, pivot_entry in zip(work[row], work[column])]
    return [row[-1] for row in work]


def totient(number):
    return sum(math.gcd(index, number) == 1 for index in range(1, number + 1))


def multiply_gaussian(first, second):
    real_first, imaginary_first = first
    real_second, imaginary_second = second
    return (real_first * real_second - imaginary_first * imaginary_second,
            real_first * imaginary_second + imaginary_first * real_second)


def add_gaussian(first, second):
    return first[0] + second[0], first[1] + second[1]


def gaussian_derivation(value):
    return value[1], 0


def check_gaussian():
    basis = [(1, 0), (0, 1)]
    for first, second in itertools.product(basis, repeat=2):
        conjugate_second = second[0], -second[1]
        right = add_gaussian(multiply_gaussian(gaussian_derivation(first), conjugate_second),
                             multiply_gaussian(first, gaussian_derivation(second)))
        assert gaussian_derivation(multiply_gaussian(first, second)) == right
    matrix = difference_matrix(4, 1, 3)
    assert matrix == [[0, 2], [-2, 0]]
    adjugate = [[matrix[1][1], -matrix[0][1]], [-matrix[1][0], matrix[0][0]]]
    coordinates = [1, 0]
    adjugate_column = [sum(value * coordinate for value, coordinate in zip(row, coordinates))
                       for row in adjugate]
    printed = [Fraction(value, 2) for value in adjugate_column]
    corrected = [Fraction(value, determinant(matrix)) for value in adjugate_column]
    assert printed == [0, 1] and corrected == [0, Fraction(1, 2)]
    return {"matrix": matrix, "determinant": determinant(matrix), "C": coordinates,
            "adjugate": adjugate, "printed_vector": list(map(str, printed)),
            "corrected_vector": list(map(str, corrected)),
            "twisted_law_basis_pairs": 4,
            "basis_reduction": "The twisted law is Z-bilinear in the two inputs.",
            "non_innerness": "The unique multiplier solving 1=-2*i*gamma is i/2, outside Z[i]."}


def check_determinants():
    conductors = [3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 15, 16,
                  18, 20, 21, 24, 25, 27, 30, 32, 36, 45, 48, 60]
    results = []
    for conductor in conductors:
        units = [index for index in range(1, conductor) if math.gcd(index, conductor) == 1]
        pairs = 0
        orders = set()
        determinants = set()
        group_types = set()
        for first, second in itertools.combinations(units, 2):
            order = conductor // math.gcd(conductor, second - first)
            degree_ratio, remainder = divmod(totient(conductor), totient(order))
            assert remainder == 0
            base = sum(cyclotomic(order))
            assert base == (prime_power_base(order) or 1)
            expected = base ** degree_ratio
            matrix = difference_matrix(conductor, first, second)
            actual = determinant(matrix)
            assert actual == expected and actual > 0
            if base > 1:
                nullity = len(matrix) - rank_modulo(matrix, base)
                assert nullity == degree_ratio
                group_types.add((base, nullity))
            else:
                group_types.add((1, 0))
            orders.add(order)
            determinants.add(actual)
            pairs += 1
        all_pairs_outer = 1 not in determinants
        expected_all = bool(prime_power_base(conductor) or
                            (conductor % 2 == 0 and (conductor // 2) % 2 == 1
                             and prime_power_base(conductor // 2)))
        assert all_pairs_outer == expected_all
        results.append({"conductor": conductor, "degree": len(cyclotomic(conductor)) - 1,
                        "unordered_twist_pairs": pairs, "root_orders": sorted(orders),
                        "positive_determinants": sorted(determinants),
                        "outer_group_types": [{"prime": prime, "dimension": dimension}
                                              for prime, dimension in sorted(group_types)],
                        "all_pairs_have_nonzero_outer_classes": all_pairs_outer})
    return results


def check_exponent_zero():
    matrix = difference_matrix(9, 1, 2)
    coordinates = [1] + [0] * (len(matrix) - 1)
    inverse_vector = solve_exact(matrix, coordinates)
    assert determinant(matrix) == 3
    assert any(entry.denominator > 1 for entry in inverse_vector)
    return {"conductor": 9, "twist_exponents": [1, 2], "valuation": 0,
            "C": coordinates, "inverse_vector": list(map(str, inverse_vector)),
            "printed_vector": list(map(str, [3 * entry for entry in inverse_vector])),
            "corrected_vector_nonintegral": True,
            "coefficient_divisibility_by_p_to_valuation": True}


def check_half_integer_quadratic():
    def multiply(first, second):
        real_first, root_first = first
        real_second, root_second = second
        return (real_first * real_second - root_first * root_second,
                real_first * root_second + root_first * real_second + root_first * root_second)

    basis = [(1, 0), (0, 1)]
    for first, second in itertools.product(basis, repeat=2):
        conjugate_second = second[0] + second[1], -second[1]
        right = add_gaussian(multiply(gaussian_derivation(first), conjugate_second),
                             multiply(first, gaussian_derivation(second)))
        assert gaussian_derivation(multiply(first, second)) == right
    matrix = difference_matrix(6, 1, 5)
    assert matrix == [[1, 2], [-2, -1]]
    inverse_vector = solve_exact(matrix, [1, 0])
    assert inverse_vector == [Fraction(-1, 3), Fraction(2, 3)]
    return {"radicand": -3, "generator": "theta=(1+sqrt(-3))/2", "conductor": 6,
            "D_theta": [1, 0], "beta_even": True, "corrected_divisibility": False,
            "matrix": matrix, "determinant": determinant(matrix),
            "inverse_vector": list(map(str, inverse_vector)), "twisted_law_basis_pairs": 4}


def main():
    parser = argparse.ArgumentParser(description="Exact cyclotomic multiplication-matrix checks.")
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    determinants = check_determinants()
    result = {"schema_version": 1, "arithmetic": "integer arithmetic and exact rational fractions",
              "scope": "Gaussian twisted law and adjugate vectors; the conductor-9 valuation-zero and conductor-6 half-integer quadratic examples; positive determinants and additive quotient groups for every unordered distinct unit pair in the stated finite conductor list.",
              "gaussian": check_gaussian(), "determinant_checks": determinants,
              "valuation_zero_example": check_exponent_zero(),
              "half_integer_quadratic_example": check_half_integer_quadratic(),
              "group_method": "For order h a prime power, determinant ell^d and nullity d modulo ell force exactly d nontrivial Smith factors, all equal to ell. A determinant of one gives the zero quotient.",
              "total_unordered_twist_pairs": sum(item["unordered_twist_pairs"] for item in determinants),
              "status": "all assertions satisfied"}
    content = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.output:
        arguments.output.write_text(content)
    print(content, end="")


if __name__ == "__main__":
    main()
