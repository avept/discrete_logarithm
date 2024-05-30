import numpy as np
from sympy.ntheory import factorint, primerange
import argparse
import math
import time
import matplotlib.pyplot as plt

def factor_base(n):
    limit = int(3.38 * np.exp(0.5 * np.sqrt(np.log2(n) * np.log2(np.log2(n)))))
    print("limit: ", limit)
    
    factor_base = list(primerange(2, limit + 1))
    
    return factor_base

def is_system_solvable(a, b):
    A = np.array(a)
    B = np.array(b).reshape(-1, 1)
    rankA = np.linalg.matrix_rank(A)
    rankAB = np.linalg.matrix_rank(np.hstack((A, B)))
    
    return rankA == rankAB == len(A)

def isSquare(m):
    return all(len(row) == len(m) for row in m)

def solve_system2(A, B, n):
    sle = np.hstack((A, B.reshape(-1, 1))) 
    num_rows, num_cols = sle.shape

    processed = []
    for j in range(num_cols - 1):
        gcd_col_mod = [math.gcd(elem, (n - 1)) for elem in sle[:, j]]
        
        for i, gcd_val in enumerate(gcd_col_mod):
            if i in processed:
                continue
        
            if gcd_val == 1:
                inv_elem = pow(int(sle[i, j]), -1, (n - 1))
                processed.append(i)

                sle[i] = (sle[i] * inv_elem) % (n - 1)
                
                mask = np.arange(num_rows) != i
                sle[mask] = (sle[mask] - np.outer(sle[mask, j], sle[i])) % (n - 1)
                break
        
    solution = []
    for j in range(num_cols - 1):
        non_zero_indices = np.nonzero(sle[:, j])[0]
        if non_zero_indices.size > 0:
            solution.append(sle[non_zero_indices[0], -1])
        else:
            solution.append(0)
    
    return solution

def linear_system(factor_base, a, n):
    system = [] # Ax = B
    right_part = []
    
    i = 0
    while True:
        sub_logarithm = pow(a, i, n)
        
        log_decomposition = factorint(sub_logarithm)
        is_smooth_log = all(prime_number in factor_base for prime_number in log_decomposition.keys())
        
        if not is_smooth_log:
            i += 1
            continue
        
        equation = []
        for prime_number in factor_base:
            if prime_number in log_decomposition:
                equation.append(log_decomposition[prime_number])
            else:
                equation.append(0)
                
        system.append(equation)
        right_part.append(i)
        
        A = np.array(system)
        B = np.array(right_part)
        if len(system) >= len(factor_base) + 20:
            print("system: ", system)
            print("right part: ", right_part)
            solution = solve_system2(A, B, n)
            
            # if solution is not None:
            return solution
                
        i += 1
        
def evaluate_logarithm(factor_base, system_solution, a, b, n):
    i = 0
    while True:
        number = (b * pow(a, i, n)) % n
        
        number_decomposition = factorint(number)
        is_smooth_number = all(prime_number in factor_base for prime_number in number_decomposition.keys())
        
        if not is_smooth_number:
            i += 1
            continue
        
        result = 0
        for j in range(len(factor_base)):
            if factor_base[j] in number_decomposition:
                result = (result + (number_decomposition[factor_base[j]] * system_solution[j])) % (n-1)
        
        return ((result - i) % (n-1))

def index_calculus(a, b, n):
    base = factor_base(n)
    solution = linear_system(base, a, n)
    print("solution: ", solution)
    result = evaluate_logarithm(base, solution, a, b, n)
        
    return result

def time_vs_order_plot():
    digit_lengths = [3, 4, 5, 6, 7, 8]
    computation_times = [0.001, 0.002, 0.004, 0.014, 0.04, 0.12]

    plt.figure(figsize=(10, 6))
    plt.plot(digit_lengths, computation_times, marker='o', linestyle='-', color='b')

    plt.title('Computation Time of Discrete Logarithm vs Order of Number')
    plt.xlabel('Number of Digits')
    plt.ylabel('Computation Time (seconds)')
    plt.grid(True)

    plt.show()

if __name__ == "__main__":
    start_time = time.time()

    a = 179
    b = 97
    n = 191
                        
    print("result of index calculus: ", index_calculus(a, b, n))

    end_time = time.time()
    execution_time = end_time - start_time
    print("Execution time:", execution_time, "seconds")