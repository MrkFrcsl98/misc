# Miscellaneous Mathematical Concepts and Algorithms

This repository contains a collection of mathematical concepts and algorithms implemented in C++. The code covers various topics such as modular arithmetic, set theory, GCD calculation, Euclidean algorithm, matrix operations, and more.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Functions and Classes](#functions-and-classes)
  - [Mathematical Functions](#mathematical-functions)
  - [Set Operations](#set-operations)
  - [Matrix Operations](#matrix-operations)
  - [ASCII and Binary Conversions](#ascii-and-binary-conversions)
- [Usage](#usage)
- [Tests](#tests)

## Overview

This project provides implementations of various mathematical algorithms and concepts, including:
- Modular arithmetic
- Greatest Common Divisor (GCD)
- Extended Euclidean Algorithm
- Prime number computations
- Matrix operations
- Set operations
- ASCII and binary conversions

## Features

- Comprehensive mathematical functions
- Set and matrix operations
- Prime number generation
- Statistical functions (mean, median, mode, range)
- ASCII and binary conversion utilities
- Extensive inline documentation and comments

## Functions and Classes

### Mathematical Functions

1. **Modulo Class Calculation**
   - **Function:** `getModuloClass(int n, int m)`
   - **Description:** Returns the remainder of `n` divided by `m` if both are positive; otherwise, returns 0.

2. **Greatest Common Divisor (GCD)**
   - **Function:** `getGCD(int a, int b)`
   - **Description:** Computes the GCD of two integers using the Euclidean algorithm.

3. **Extended Euclidean Algorithm**
   - **Function:** `extendedEuclideanAlgorithm(int a, int m, int &x, int &y)`
   - **Description:** Computes the GCD of `a` and `m` and also finds integers `x` and `y` such that `ax + my = gcd(a, m)`.

4. **Coprime Check**
   - **Function:** `areCoprime(int a, int b)`
   - **Description:** Checks if two integers are coprime (GCD is 1).

5. **Prime Number Computation**
   - **Function:** `primeComputation(size_t threshold)`
   - **Description:** Generates an array of prime numbers up to a specified threshold.
   - **Function:** `sieveOfEratosthenes(size_t limit)`
   - **Description:** Computes prime numbers using the Sieve of Eratosthenes algorithm.

6. **Exponentiation**
   - **Function:** `raise(int b, int e)`
   - **Description:** Computes `b` raised to the power of `e`.

7. **Fibonacci Sequence**
   - **Function:** `fibonacciSequence()`
   - **Description:** Generates the Fibonacci sequence up to `N` elements.

8. **Even/Odd Check**
   - **Functions:** `isEven(T a)`, `isOdd(T a)`
   - **Description:** Checks if a number is even or odd.

9. **Statistics Calculations**
   - **Functions:** `getMode(const std::array<T, N> &arr)`, `getMedian(const std::array<T, N> _s)`, `getMean(const std::array<T, N> _s)`, `getRange(const std::array<T, N> _s)`
   - **Description:** Computes the mode, median, mean, and range of an array of numbers.

### Set Operations

1. **Set Class**
   - **Class:** `Set`
   - **Description:** Implements a set with operations like union, intersection, complement, and checks for equality and subsets.

### Matrix Operations

1. **Matrix Class**
   - **Class:** `Matrix`
   - **Description:** Implements a matrix with operations such as addition, transpose, and element distribution.

### ASCII and Binary Conversions

1. **Conversions**
   - **Functions:** `asciiToBinary(const std::string &input)`, `asciiToHex(const std::string &_d)`, `binaryToAscii(const std::string &binary)`, `hexToAscii(const std::string &hex)`, `binaryToHex(const std::string &binary)`, `hexToBinary(const std::string &hex)`
   - **Description:** Converts ASCII strings to binary/hexadecimal and vice versa.

## Usage

To use the functions and classes provided in this repository, include the relevant header files in your C++ project. The `main.cpp` file demonstrates how to utilize various functions and classes with example usage.

### Example Usage

```cpp
#include <iostream>
#include "math.hpp"

int main() {
    int n = 13, m = 5;
    std::cout << "Modulo Class of: 13 mod 5 = " << Math::getModuloClass(n, m) << std::endl;

    int a = 18, b = 5;
    std::cout << "GCD of 18 and 5 = " << Math::getGCD(a, b) << std::endl;

    return 0;
}
```

## Tests

The `main.cpp` file includes a `runTests` function within the `Matrices` namespace that performs various tests on the matrix operations. Additional tests for other functions can be added in a similar manner to ensure correctness.

### Running Tests

Compile and run the `main.cpp` file to execute the tests:

```sh
g++ -o main main.cpp
./main
```


This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

Feel free to reach out if you have any questions or need further assistance. Happy coding!
