```cpp
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>

namespace Math {
// Function to get the modulo class of n with respect to m
// Returns n % m if both n and m are positive; otherwise, returns 0.
static const int getModuloClass(const int n, const int m) {
  if (n > 0 && m > 0) {
    return n % m; // Return the remainder of n divided by m
  }
  return 0; // Return 0 if either n or m is non-positive
};

// Function to compute the Greatest Common Divisor (GCD) of two integers a and b
// Uses the Euclidean algorithm, which is efficient for this purpose.
// The 'noexcept' specifier indicates that this function does not throw exceptions.
static const int getGCD(const int a, const int b) noexcept {
  int dvd = a; // Dividend
  int dvs = b; // Divisor
  int rem = 0; // Remainder
  // Loop until the divisor becomes zero
  while (dvs != 0) {
    rem = dvd % dvs; // Calculate the remainder
    dvd = dvs;       // Update dividend to the current divisor
    dvs = rem;       // Update divisor to the remainder
  }
  return dvd; // Return the GCD
};

// Function to find the Extended Euclidean Algorithm for integers a and m
// This function computes the GCD of a and m, and also finds integers x and y
// such that ax + my = gcd(a, m).
// The 'noexcept' specifier indicates that this function does not throw exceptions.
static const int extendedEuclideanAlgorithm(int a, int m, int &x, int &y) noexcept {
  int x0 = 1; // Coefficient for a
  int y0 = 0; // Coefficient for m
  int x1 = 0; // Next coefficient for a
  int y1 = 1; // Next coefficient for m
  int x2 = 0; // Temporary variable for x
  int y2 = 0; // Temporary variable for y
  int r = 0;  // Remainder
  int q = 0;  // Quotient
  // Loop until m becomes zero
  while (m != 0) {
    r = a % m; // Calculate the remainder
    q = a / m; // Calculate the quotient
    a = m;     // Update a to m
    m = r;     // Update m to the remainder
    // Update coefficients using the previous values
    x2 = x0 - q * x1;
    y2 = y0 - q * y1;
    x0 = x1; // Move to the next coefficient for a
    y0 = y1; // Move to the next coefficient for m
    x1 = x2; // Update x1 to the new coefficient for a
    y1 = y2; // Update y1 to the new coefficient for m
  }
  x = x0;   // Set the output parameter x to the coefficient for a
  y = y0;   // Set the output parameter y to the coefficient for m
  return a; // Return the GCD
};

// Function to check if two integers a and b are coprime
// Returns true (1) if they are coprime (GCD is 1), otherwise returns false (0).
static const int areCoprime(const int a, const int b) noexcept {
  return getGCD(a, b) == 1; // Check if the GCD of a and b is 1
};

template <std::size_t N> static std::array<int, N> primeComputation(const std::size_t threshold) noexcept {
  std::array<int, N> _r; // Array to store prime numbers
  _r[0] = 2;             // The first prime number is 2
  std::size_t k = 1;     // Index for storing primes, starting from the second position

  // Iterate through odd numbers starting from 3 up to the threshold
  for (int i = 3; i <= threshold; i += 2) {
    if (k >= N)
      break;       // Stop if the array is full
    bool p = true; // Assume i is prime

    // Check for factors from 3 to the square root of i
    for (int j = 3; j <= std::sqrt(i); j += 2) {
      if (i % j == 0) { // If i is divisible by j, it's not prime
        p = false;
        break; // Exit the loop early
      }
    }

    // If i is prime, store it in the array
    if (p) {
      _r[k++] = i; // Store the prime number and increment the index
    }
  }

  return _r; // Return the array of prime numbers
};

// compute primes using Sieve of Eratosthenes...
constexpr std::vector<int> sieveOfEratosthenes(const std::size_t limit) {
  std::vector<bool> isPrime(limit + 1, true); // Initialize a boolean vector
  std::vector<int> primes;                    // Vector to store prime numbers

  isPrime[0] = isPrime[1] = false; // 0 and 1 are not prime numbers

  for (int p = 2; p * p <= limit; ++p) {
    if (isPrime[p]) {
      // Mark all multiples of p as non-prime
      for (int i = p * p; i <= limit; i += p) {
        isPrime[i] = false;
      }
    }
  }

  // Collect all prime numbers
  for (int p = 2; p <= limit; ++p) {
    if (isPrime[p]) {
      primes.push_back(p);
    }
  }

  return primes;
};

// Function to raise a base 'b' to the exponent 'e'
// This function computes b^e using a simple loop.
static constexpr int raise(const int b, const int e) noexcept {
  int i = b; // Initialize the result with the base
  int c = 0; // Counter for the exponent
  // Loop to multiply the base 'b' by itself 'e' times
  while (++c < e) {
    i *= b; // Multiply the current result by the base
  }
  return i; // Return the final result
}

// Template function to generate a Fibonacci sequence of size N
// The function fills an array with Fibonacci numbers up to the T-th-1 term.
template <std::size_t N> static constexpr std::array<int, N> fibonacciSequence(const std::size_t T) noexcept {
  std::array<int, N> v; // Array to hold the Fibonacci sequence
  int a = 0;            // First Fibonacci number
  int b = 1;            // Second Fibonacci number
  int c = a + b;        // Next Fibonacci number (initially 1)
  int r = 0;            // Index for filling the array
  v[0] = c;             // Store the first Fibonacci number in the array

  // Loop to generate Fibonacci numbers until the T-th term
  while (r < T - 1) {
    a = b;      // Update 'a' to the previous 'b'
    b = c;      // Update 'b' to the current Fibonacci number
    c = a + b;  // Calculate the next Fibonacci number
    v[++r] = c; // Store the new Fibonacci number in the array
  }
  return v; // Return the array containing the Fibonacci sequence
}

// Function to check if a number is even
// Returns true if 'a' is even, false otherwise.
template <typename T> static constexpr bool isEven(const T a) noexcept {
  return a % 2 == 0; // Check if 'a' is divisible by 2
}

// Function to check if a number is odd
// Returns true if 'a' is odd, false otherwise.
template <typename T> static constexpr bool isOdd(const T a) noexcept {
  return a % 2 != 0; // Check if 'a' is not divisible by 2
}

// Function to compute the mode of an array
// The mode is the value that appears most frequently in the array.
// T: Type of elements in the array
// N: Size of the array
// H: Maximum possible value in the array (used for the size of the frequency array)
template <typename T, std::size_t N, std::size_t H> static constexpr int getMode(const std::array<T, N> _s) noexcept {
  std::size_t M = 0;                        // Variable to store the mode (most frequent value)
  T R[H + 1];                               // Frequency array to count occurrences of each value
  std::fill(std::begin(R), std::end(R), 0); // Initialize frequency array to 0

  // Count occurrences of each value in the input array
  for (int i{0}; i < N; ++i) {
    R[_s[i]] += 1; // Increment the count for the value _s[i]
  }

  // Find the mode by checking the frequency array
  for (int C = 0; C < std::size(R); ++C) {
    if (M >= R[C]) {
      continue; // Skip if the current count is not greater than the current mode
    }
    M = C; // Update mode to the current index if it has a higher count
  }
  return M; // Return the mode
}

// Function to calculate the median of an array
template <typename T, std::size_t N> static constexpr int getMedian(const std::array<T, N> _s) noexcept {
  // Create a vector from the input array for sorting
  std::vector<T> O(_s.begin(), _s.end());

  // Sort the vector in ascending order
  std::sort(O.begin(), O.end());

  // Check if the number of elements is even or odd
  // If even, return the average of the two middle elements
  // If odd, return the middle element
  return N % 2 == 0 ? (O[(N / 2) - 1] + O[(N / 2)]) / 2 : O[N / 2];
}

// Function to calculate the mean (average) of an array
template <typename T, std::size_t N> static constexpr int getMean(const std::array<T, N> _s) noexcept {
  std::size_t r{0}; // Initialize a variable to hold the sum
  // Iterate through each element in the array and accumulate the sum
  for (const T i : _s) {
    r += i;
  }
  // Return the mean by dividing the total sum by the number of elements
  return r / _s.size();
}

// Function to calculate the range of an array
template <typename T, std::size_t N> static constexpr int getRange(const std::array<T, N> _s) noexcept {
  // Create a copy of the input array to sort
  std::array<T, N> r(_s);

  // Sort the copied array in ascending order
  std::sort(r.begin(), r.end());

  // Return the difference between the maximum and minimum values
  return r[r.size() - 1] - r[0];
};



}; // namespace Math

int main(int argc, char **argv) {

  {
    constexpr int M = 5;

    // Definition of congruence classes for modulo "5"
    std::array<std::array<int, 10>, M> classSet;
    classSet[0] = {-25, -20, -15, -10, -5, 0, 5, 10, 15, 20};
    classSet[1] = {-24, -19, -14, -9, -4, 1, 6, 11, 16, 21};
    classSet[2] = {-23, -18, -13, -8, -3, 2, 7, 12, 17, 22};
    classSet[3] = {-22, -17, -12, -7, -2, 3, 8, 13, 18, 23};
    classSet[4] = {-21, -16, -11, -6, -1, 4, 9, 14, 19, 24}; // n - 1

    // display congruence classes...
    std::cout << "Congruence Classes for modulo " << M << ": \n";
    for (int i = 0; i < classSet.size(); ++i) {
      std::cout << "Class " << i << " : {";
      for (int j = 0; j < classSet[i].size(); ++j) {
        std::cout << classSet[i][j] << (j < classSet[i].size() - 1 ? ", " : "");
      }
      std::cout << "}" << "\n";
    }

    // each class can have infinite n of numbers, represented as: n{...,-5,0,5,10,...}
    // here we defined 10 elements only.
    // everytime we operate with this modulo(5), the result of the operation: n mod 5
    // will always find its result within one of these congruence classes.
    // There are no duplicates within the sets, and when we do n mod m the result
    // will always be in one of the sets.
    // The classes start from 0 and go up to m-1, where m is the modulo value.
    // Let's define some operations: n = b mod M, where n is "13", b is ?, and M is 5,
    // now we find which class n belongs to by dividing n by m and saving the remainder,
    // when we work with modular arithmetic, we don't care about the quotient, we only care
    // about the remainder, which will give us the class for n.
    // For example, 13(n) mod 5(m) will give us 3(b), 3 in this case is the congruence class, and
    // if we look into the class 3(b), we will find 13(n).

    std::cout << "Modulo Class of: 13 mod " << M << " = " << Math::getModuloClass(13, M) << "\n";

    // now let's verify if n mod m = b, b = 4, n = 22, m = 5:

    std::cout << "Is 22 mod " << M << " = 4 ? " << std::boolalpha << (Math::getModuloClass(22, M) == 4) << "\n";

    // let's try with other values..
    std::cout << "Is 22 mod " << M << " = 3 ? " << std::boolalpha << (Math::getModuloClass(22, M) == 3) << "\n";
    std::cout << "Is 22 mod " << M << " = 2 ? " << std::boolalpha << (Math::getModuloClass(22, M) == 2) << "\n";
    std::cout << "Is 22 mod " << M << " = 1 ? " << std::boolalpha << (Math::getModuloClass(22, M) == 1) << "\n";
    std::cout << "Is 22 mod " << M << " = 0 ? " << std::boolalpha << (Math::getModuloClass(22, M) == 0) << "\n";
  }
  // Let's define another set of congruence classes for another modulo, this time using modulo 7.
  // And instead of populating the sets manually, we will do it programatically.
  {
    constexpr int M = 9;                        // Modulus
    constexpr int E = 15;                       // Maximum number of entries per class
    std::array<std::array<int, E>, M> classSet; // 2D array to hold congruence classes

    // Fill the classSet with congruence classes
    for (int i = 0; i < M; ++i) {   // Iterate over each congruence class
      for (int j = 0; j < E; ++j) { // Iterate over the entries in each class
        if (i == 0) {               // Special case for the first congruence class (0 mod M)
          if (j == 0) {
            classSet[i][0] = M * (-E / 2); // Initialize the first entry for class 0
            continue;                      // Skip to the next iteration
          }
        }
        if (j == 0) {                              // First entry for other classes
          classSet[i][j] = classSet[i - 1][0] + 1; // Set the first entry based on the previous class
          continue;                                // Skip to the next iteration
        }

        // Fill the rest of the entries for the current class
        classSet[i][j] = classSet[i][j - 1] + M; // Increment by M to get the next entry
      }
    }
    // display congruence classes...
    std::cout << "Congruence Classes for modulo " << M << ": \n";
    for (int i = 0; i < classSet.size(); ++i) {
      std::cout << "Class " << i << " : {";
      for (int j = 0; j < classSet[i].size(); ++j) {
        std::cout << classSet[i][j] << (j < classSet[i].size() - 1 ? ", " : "");
      }
      std::cout << "}" << "\n";
    }
  }

  // Modular arithmetic is mainly used for cryptography operations.
  // The most common type of set used in cryptography is the set of integers
  // modulo M, starting from 0 up to M-1. The 2 primary operations performed
  // on this type of sets are addition and multiplication.
  // An important aspect of these sets is the closure, which guarantees that
  // when multiplying or adding any two numbers from that set, the result will
  // always be within the set, this ensures that the operation never leaves
  // the set. Another important aspect is associativity, meaning that for any
  // * or + operations the order does not matter.
  // The set also includes a neutral element, which is '0': a+0 ≡ a (mod M).
  // Additionally, each element in the set has an additive inverse, denoted
  // as -a: a + (-a) ≡ 0 (mod M).
  // Addition, subtraction and multiplication always have an additive inverse,
  // but division not always.
  // For division, the additive inverse exists only if the result of GCD(a, m) = 1,
  // if not, then it has not an additive inverse.
  // These properties define a structure known as a RING in abstract algebra.
  // In a RING there are 4 primary operations, addition, subtraction, multiplication
  // and division. As said, the first 3 are straighforward, division is not always.
  {

    // Declare some variables
    int x = 0;   // Variable for storing inverses or intermediate results
    int y = 0;   // Variable for storing inverses or intermediate results
    int n = 0;   // Variable for the number to operate on
    int m = 0;   // Variable for the modulus
    int b = 0;   // Variable for results of operations
    int g = 0;   // Variable for GCD
    int r = 0;   // Variable for results of division
    int inv = 0; // Variable for storing the multiplicative inverse

    // Addition operation
    n = 13;              // This represents n
    m = 5;               // This is the modulo
    x = m - (n % m) % m; // Find the additive inverse of n mod m
    b = (x + n) % m;     // Result of (n + x) mod m

    std::cout << "(Addition): The additive inverse of " << n << " mod " << m << " = " << x << "\n";
    std::cout << "(Addition): Result: " << n << " + " << x << " = " << b << " mod " << m << "\n";

    // Subtraction operation
    n = 15;                // New value for n
    m = 6;                 // New modulus
    x = (m - (n % m)) % m; // Find the additive inverse of n mod m
    b = (n + x) % m;       // Result of (n - x) mod m

    std::cout << "(Subtraction): The additive inverse of " << n << " mod " << m << " = " << x << "\n";
    std::cout << "(Subtraction): Result: " << n << " - " << x << " = " << b << " mod " << m << "\n";

    // Multiplication operation
    n = 18; // New value for n
    m = 5;  // New modulus

    g = Math::getGCD(n, m); // Calculate GCD of n and m

    if (g == 1) {                                       // Check if n and m are coprime
      g = Math::extendedEuclideanAlgorithm(n, m, x, y); // Find x and y such that n*x + m*y = gcd(n, m)
      inv = (x % m + m) % m;                            // Ensure the inverse is positive
      std::cout << "(Multiplication): The multiplicative inverse of " << n << " mod " << m << " = " << inv << "\n";
      std::cout << "(Multiplication): Result: (" << n << "*" << inv << ") + (" << m << "*" << y << ") = " << g << " mod " << m << "\n";
    } else {
      std::cout << n << " and " << m << " are not coprime, so there is no multiplicative inverse.\n";
    }

    // Division operation
    n = 23; // Dividend
    b = 4;  // Divisor
    m = 9;  // Modulus

    g = Math::getGCD(b, m); // Calculate GCD of b and m

    if (g == 1) {                                       // Check if b and m are coprime
      g = Math::extendedEuclideanAlgorithm(b, m, x, y); // Find the multiplicative inverse of b
      inv = (x % m + m) % m;                            // Ensure the inverse is positive
      r = (n * inv) % m;                                // Result of (n / b) mod m

      std::cout << "(Division): The multiplicative inverse of " << b << " mod " << m << " = " << inv << "\n";
      std::cout << "(Division): Result: " << n << " * " << inv << " = " << r << " mod " << m << "\n";
    } else {
      std::cout << b << " and " << m << " are not coprime, so there is no multiplicative inverse.\n";
    }
  }

  // Prime Numbers computation, prime numbers are numbers that divide only by 1 or themselves,
  // for example, 3 is prime because it only divides by 3 and 1.
  constexpr std::size_t threshold = 60;                                  // prime value threshold, computation stops at this value
  constexpr std::size_t PMAX = 20;                                        // max entries for prime computation result
  std::array<int, PMAX> primes = Math::primeComputation<PMAX>(threshold); // compute primes and store max 100 values
  std::vector<int> sievePrimes = Math::sieveOfEratosthenes(threshold); // using sieve...
  std::cout << "Computed Prime Numbers(Standard): ";
  for (const auto i : primes) {
    std::cout << i << ", ";
  }
  std::cout << "\n";
  std::cout << "Computed Prime Numbers(Sieve):    ";
  for (const auto i : sievePrimes) {
    std::cout << i << ", ";
  }
  std::cout << "\n";

  constexpr int base = 2;
  constexpr int exp = 8;
  constexpr int raise2 = Math::raise(base, exp); // this is = 256, 2^8 = 256
  std::cout << "Raising " << base << " to the power of " << exp << ", " << base << "^" << exp << " = " << raise2 << "\n";

  constexpr std::size_t fst = 15;
  constexpr std::array<int, fst> fibSeq = Math::fibonacciSequence<fst>(fst);

  std::cout << "Fibonacci Sequence: {";
  for (const int i : fibSeq) {
    std::cout << i << ", ";
  }
  std::cout << "}\n";

  constexpr std::uint8_t N1{10}, N2{13};

  std::cout << "Is Even " << (int)N1 << "? " << std::boolalpha << Math::isEven<std::uint8_t>(N1) << "\n";
  std::cout << "Is Odd " << (int)N2 << "? " << std::boolalpha << Math::isOdd<std::uint8_t>(N2) << "\n";

  constexpr std::size_t MODE_SET_MAX = 10;
  constexpr std::array<int, MODE_SET_MAX> mode_set = {5, 10, 2, 5, 4, 5, 7, 5, 4, 9};
  constexpr int getMode = Math::getMode<int, MODE_SET_MAX, 10>(mode_set);

  std::cout << "Mode of Set{";
  for (const int i : mode_set) {
    std::cout << i << ", ";
  }
  std::cout << "} = " << getMode << "\n";

  constexpr std::array<int, 10> median_set = {5, 1, 9, 4, 3, 8, 2, 4, 7, 8}; // 1,2,3,4, 4,5 ,7,8,8,9
  constexpr int getMedian = Math::getMedian<int, 10>(median_set);

  std::cout << "Median from Set {";
  for (const int i : median_set) {
    std::cout << i << ", ";
  }
  std::cout << "} = " << getMedian << "\n";

  constexpr std::array<int, 10> mean_set = {15, 7, 8, 13, 1, 5, 3, 28, 9, 4};
  constexpr int getMean = Math::getMean(mean_set);
  std::cout << "Mean from Set {";
  for (const int i : mean_set)
    std::cout << i << ", ";
  std::cout << "} = " << getMean << "\n";

  constexpr std::array<int, 10> range_set = {5, 7, 8, 6, 4, 3, 2, 4, 11, 5};
  constexpr int getRange = Math::getRange(range_set);
  std::cout << "Range from Set {";
  for (const int i : range_set)
    std::cout << i << ", ";
  std::cout << "} = " << getRange << "\n";

  return 0;
}
```
