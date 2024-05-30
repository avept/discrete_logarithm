import numpy as np
from sympy import factorint, Matrix
import math

def factor_base(n):
    limit = 3.38 * np.exp((1/2) * (np.log(n) * np.log(np.log(n))) ** (1/2)) # Adjusted c value to reduce limit
    print("limit: ", limit)
    
    primes = [True] * (n+1)
    primes[0] = primes[1] = False
    
    p = 2
    while (p * p <= n):
        if primes[p] == True:
            for i in range(p * p, n+1, p):
                primes[i] = False
        p += 1
    
    factor_base = []
    for p in range(2, int(limit)+1):
        if primes[p]:
            factor_base.append(p)
    
    return factor_base

def is_system_solvable(a, b):
    A = np.array(a)
    B = np.array(b).reshape(-1, 1)
    rankA = np.linalg.matrix_rank(A)
    rankAB = np.linalg.matrix_rank(np.hstack((A, B)))
    
    return rankA == rankAB == len(A)

def isSquare(m):
    return all(len(row) == len(m) for row in m)

def solve_system(A, B, n):
    A = np.hstack((A, B.reshape(-1, 1))) 
    n_mod = n - 1
    num_rows, num_cols = A.shape
    
    for i in range(num_rows):
        pivot = A[i, i]
        pivot_row = i
        
        if np.gcd(pivot, n_mod) != 1:
            for row in range(i + 1, num_rows):
                if np.gcd(A[row, i], n_mod) == 1:
                    A[[i, row]] = A[[row, i]]
                    pivot = A[i, i]
                    pivot_row = row
                    break
        
        if np.gcd(pivot, n_mod) != 1:
            raise ValueError("Matrix is singular and cannot be solved under mod {}".format(n_mod))
        
        pivot_inv = pow(int(pivot), -1, n_mod)
        A[i] = (A[i] * pivot_inv) % n_mod
        
        for j in range(i + 1, num_rows):
            factor = A[j, i]
            A[j] = (A[j] - factor * A[i]) % n_mod
    
    solution = [0] * (num_cols - 1)
    for i in range(num_rows - 1, -1, -1):
        solution[i] = (A[i, -1] - sum(A[i, j] * solution[j] for j in range(i + 1, num_cols - 1))) % n_mod
        
    if 0 in solution:
        return None
    
    return solution

def linear_system(factor_base, a, n):
    system = [] # Ax = B
    right_part = []
    
    iterations_count = 0
    while True:
        random_number = np.random.randint(1, n)
        sub_logarithm = pow(a, random_number, n)
        
        log_decomposition = factorint(sub_logarithm)
        is_smooth_log = all(prime_number in factor_base for prime_number in log_decomposition.keys())
        
        if random_number in right_part:
            continue
        
        if not is_smooth_log:
            continue
        
        equation = []
        for prime_number in factor_base:
            if prime_number in log_decomposition:
                equation.append(log_decomposition[prime_number] % (n - 1))
            else:
                equation.append(0)
                
        system.append(equation)
        right_part.append(random_number % (n - 1))
        iterations_count += 1
        
        # is_solvable = is_system_solvable(system, right_part)
        # if not is_solvable:
            # system.pop()
            # right_part.pop()
        # if len(system) == len(factor_base):
        A = np.array(system)
        B = np.array(right_part)
        try:
            # print("A: ", A)
            solution = solve_system(A, B, n)
            if solution is not None:
                print("correct A: ", A)
                print("correct B: ", B)
                return solution
        except ValueError as e:
            print(e)
            print("A: ", A)
            print("B: ", B)
            
            if iterations_count >= len(factor_base) + 10:
                system = []
                right_part = []
                iterations_count = 0
                # system.pop()
                # right_part.pop()
        
def evaluateLogarithm(factor_base, system_solution, a, b, n):
    while True:
        random_number = np.random.randint(1, n)
        number = (b * pow(a, random_number, n)) % n
        
        number_decomposition = factorint(number)
        is_smooth_number = all(prime_number in factor_base for prime_number in number_decomposition.keys())
        
        if not is_smooth_number:
            continue
        
        result = 0
        for i in range(len(factor_base)):
            if factor_base[i] in number_decomposition:
                result += (number_decomposition[factor_base[i]] * (system_solution[i] % (n-1)))
        
        return ((result - random_number) % (n-1))

def index_calculus(a, b, n):
    base = factor_base(n)
    solution = linear_system(base, a, n)
    print("solution: ", solution)
    result = evaluateLogarithm(base, solution, a, b, n)
    
    print("factor base: ", base)
    
    return result

if __name__ == "__main__":
    a = 10
    b = 17
    n = 47
                    
    print("result of index calculus: ", index_calculus(a, b, n))