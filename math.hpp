#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <set>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <unordered_set>
#include <deque>
#include <list>
#include <stdexcept>
#include <cassert>


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
  x = 1;
  y = 0;              // Coefficients for a and m respectively
  int x1 = 0, y1 = 1; // Next coefficients
  while (m != 0) {
    int q = a / m; // Calculate the quotient
    int temp = a % m;
    a = m;
    m = temp; // Update a and m to the remainder
    temp = x;
    x = x1;
    x1 = temp - q * x1; // Update coefficients for x
    temp = y;
    y = y1;
    y1 = temp - q * y1; // Update coefficients for y
  }
  return a; // Return the GCD
}

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

template <std::size_t N> static constexpr std::array<int, N> fibonacciSequence() noexcept {
  std::array<int, N> v = {0, 1}; // Initialize the array with the first two Fibonacci numbers
  for (std::size_t i = 2; i < N; ++i) {
    v[i] = v[i - 1] + v[i - 2]; // Calculate the next Fibonacci number
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
template <typename T, std::size_t N, T H> static constexpr T getMode(const std::array<T, N> &arr) noexcept {
  std::array<std::size_t, static_cast<std::size_t>(H) + 1> frequency = {}; // Frequency array to count occurrences of each value
  // Count occurrences of each value in the input array
  for (const T &value : arr) {
    ++frequency[static_cast<std::size_t>(value)];
  }
  // Find the mode by checking the frequency array
  return std::distance(frequency.begin(), std::max_element(frequency.begin(), frequency.end()));
};

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

template <typename T, std::size_t N> class Set {
  // Internal array to hold the elements of the set
  std::array<T, N> _set;

public:
  // Default constructor
  explicit inline Set() noexcept {};

  // Constructor that initializes the set with a given array
  inline Set(const std::array<T, N> &_s) noexcept : _set(_s) {};

  // Returns the number of elements in the set (cardinality)
  inline constexpr std::size_t cardinality() const noexcept { return this->_set.size(); };

  // Checks if the current set is equal to another array
  template <std::size_t eN> inline constexpr bool equal(const std::array<T, eN> &_o) const noexcept {
    // If sizes are different, sets are not equal
    if (N != eN)
      return false;

    // Create copies of the sets to sort
    auto a1 = _set;
    auto a2 = _o;

    // Sort both sets for comparison
    std::sort(a1.begin(), a1.end());
    std::sort(a2.begin(), a2.end());

    // Compare sorted arrays
    return a1 == a2;
  };

  // Checks if the current set is a subset of another array
  template <typename eT, std::size_t eN> inline constexpr bool subsetOf(const std::array<eT, eN> &_o) const noexcept {
    // Create an unordered_set from the other array for efficient lookup
    std::unordered_set<eT> otherSet(_o.begin(), _o.end());

    // Check if each element in the current set exists in the other set
    return std::all_of(this->_set.begin(), this->_set.end(), [&otherSet](const T &e) { return otherSet.find(e) != otherSet.end(); });
  };

  // Checks if the current set is a proper subset of another array
  template <std::size_t eN> inline constexpr bool properSubsetOf(const std::array<T, eN> &_o) const noexcept {
    std::set<T> a1(this->_set.begin(), this->_set.end());
    std::set<T> a2(_o.begin(), _o.end());

    // Check if a1 is a subset of a2 and strictly smaller
    return std::includes(a2.begin(), a2.end(), a1.begin(), a1.end()) && a1.size() < a2.size();
  };

  // Allocates a new array to the set
  inline void alloc(const std::array<T, N> _s) noexcept {
    this->_set = _s; // Assign the new array to the internal set
  };

  // Union of the current set with another array
  // A∪B={x∣x∈A or x∈B}
  template <std::size_t aS> inline const std::set<T> unionOf(const std::array<T, aS> &_s) noexcept {
    std::set<T> s1(this->_set.begin(), this->_set.end()), s2(_s.begin(), _s.end());
    s1.insert(s2.begin(), s2.end()); // Insert elements from the second set
    return s1;                       // Return the union set
  };

  // Intersection of the current set with another array
  // A∩B={x∣x∈A and x∈B}
  template <std::size_t aS> inline const std::set<T> intersectionOf(const std::array<T, aS> &_s) noexcept {
    std::set<T> s1(this->_set.begin(), this->_set.end()), s2(_s.begin(), _s.end()), s3;
    // Find common elements between the two sets
    for (const T &e : s1) {
      if (s2.find(e) != s2.end()) {
        s3.insert(e); // Insert common element into the result set
      }
    }
    return s3; // Return the intersection set containing common elements
  };

  // Complement of the current set with respect to another array
  // A′={x∈U∣x∈/A}
  template <std::size_t aS> inline const std::set<T> complementOf(const std::array<T, aS> &_s) noexcept {
    std::set<T> s1(this->_set.begin(), this->_set.end()), s2(_s.begin(), _s.end()), s3;
    // Find elements in the second set that are not in the current set
    for (const T &e : s2) {
      if (s1.find(e) == s1.end()) {
        s3.insert(e); // Insert element into the result set if not found in the current set
      }
    }
    return s3; // Return the complement set
  };

  // Destructor to clean up the set
  ~Set() noexcept {
    this->_set = {}; // Clear the internal array (optional, as it will be automatically destroyed)
  };
};

// Converts a string of ASCII characters to a binary string representation.
const std::string asciiToBinary(const std::string &input) noexcept {
  // Lookup table for converting each ASCII character to its 8-bit binary representation.
  static const char *const lookup[256] = {
      "00000000", "00000001", "00000010", "00000011", "00000100", "00000101", "00000110", "00000111", "00001000", "00001001", "00001010", "00001011", "00001100", "00001101", "00001110", "00001111",
      "00010000", "00010001", "00010010", "00010011", "00010100", "00010101", "00010110", "00010111", "00011000", "00011001", "00011010", "00011011", "00011100", "00011101", "00011110", "00011111",
      "00100000", "00100001", "00100010", "00100011", "00100100", "00100101", "00100110", "00100111", "00101000", "00101001", "00101010", "00101011", "00101100", "00101101", "00101110", "00101111",
      "00110000", "00110001", "00110010", "00110011", "00110100", "00110101", "00110110", "00110111", "00111000", "00111001", "00111010", "00111011", "00111100", "00111101", "00111110", "00111111",
      "01000000", "01000001", "01000010", "01000011", "01000100", "01000101", "01000110", "01000111", "01001000", "01001001", "01001010", "01001011", "01001100", "01001101", "01001110", "01001111",
      "01010000", "01010001", "01010010", "01010011", "01010100", "01010101", "01010110", "01010111", "01011000", "01011001", "01011010", "01011011", "01011100", "01011101", "01011110", "01011111",
      "01100000", "01100001", "01100010", "01100011", "01100100", "01100101", "01100110", "01100111", "01101000", "01101001", "01101010", "01101011", "01101100", "01101101", "01101110", "01101111",
      "01110000", "01110001", "01110010", "01110011", "01110100", "01110101", "01110110", "01110111", "01111000", "01111001", "01111010", "01111011", "01111100", "01111101", "01111110", "01111111",
      "10000000", "10000001", "10000010", "10000011", "10000100", "10000101", "10000110", "10000111", "10001000", "10001001", "10001010", "10001011", "10001100", "10001101", "10001110", "10001111",
      "10010000", "10010001", "10010010", "10010011", "10010100", "10010101", "10010110", "10010111", "10011000", "10011001", "10011010", "10011011", "10011100", "10011101", "10011110", "10011111",
      "10100000", "10100001", "10100010", "10100011", "10100100", "10100101", "10100110", "10100111", "10101000", "10101001", "10101010", "10101011", "10101100", "10101101", "10101110", "10101111",
      "10110000", "10110001", "10110010", "10110011", "10110100", "10110101", "10110110", "10110111", "10111000", "10111001", "10111010", "10111011", "10111100", "10111101", "10111110", "10111111",
      "11000000", "11000001", "11000010", "11000011", "11000100", "11000101", "11000110", "11000111", "11001000", "11001001", "11001010", "11001011", "11001100", "11001101", "11001110", "11001111",
      "11010000", "11010001", "11010010", "11010011", "11010100", "11010101", "11010110", "11010111", "11011000", "11011001", "11011010", "11011011", "11011100", "11011101", "11011110", "11011111",
      "11100000", "11100001", "11100010", "11100011", "11100100", "11100101", "11100110", "11100111", "11101000", "11101001", "11101010", "11101011", "11101100", "11101101", "11101110", "11101111",
      "11110000", "11110001", "11110010", "11110011", "11110100", "11110101", "11110110", "11110111", "11111000", "11111001", "11111010", "11111011", "11111100", "11111101", "11111110", "11111111"};

  // Initialize an empty string to hold the binary result.
  std::string result;
  // Reserve space for the result string to optimize memory allocation.
  result.reserve(input.length() * 8); // Each character will be represented by 8 binary digits.

  // Iterate over each character in the input string.
  for (unsigned char character : input) {
    // Append the corresponding binary representation from the lookup table.
    result.append(lookup[character]);
  }

  // Return the final binary string representation.
  return result;
};

// Converts a string of ASCII characters to a hexadecimal string representation.
const std::string asciiToHex(const std::string &_d) noexcept {
  // String containing hexadecimal digits for conversion.
  static const char hexDigits[17] = "0123456789ABCDEF";
  // Initialize an empty string to hold the hexadecimal result.
  std::string _r;
  // Reserve space for the result string to optimize memory allocation.
  _r.reserve(_d.length() * 2); // Each character will be represented by 2 hex digits.

  // Iterate over each character in the input string.
  for (unsigned char E : _d) {
    // Get the high nibble (4 bits) of the character and convert it to a hex digit.
    _r += hexDigits[E >> 4]; // Shift right by 4 bits to get the high nibble.
    // Get the low nibble (4 bits) of the character and convert it to a hex digit.
    _r += hexDigits[E & 0x0F]; // Mask with 0x0F to get the low nibble.
  }

  // Return the final hexadecimal string representation.
  return _r;
};

// Converts a binary string representation to a string of ASCII characters.
const std::string binaryToAscii(const std::string &binary) noexcept {
  std::string result;
  // Reserve space for the result string to optimize memory allocation.
  result.reserve(binary.length() / 8); // Each ASCII character is represented by 8 bits.

  // Process each 8 bits (1 byte) in the binary string.
  for (size_t i = 0; i < binary.length(); i += 8) {
    // Convert the 8-bit binary substring to an ASCII character.
    std::string byteString = binary.substr(i, 8);
    char asciiChar = static_cast<char>(std::stoi(byteString, nullptr, 2)); // Convert binary to decimal.
    result += asciiChar;                                                   // Append the ASCII character to the result.
  }

  // Return the final ASCII string representation.
  return result;
}

// Converts a hexadecimal string representation to a string of ASCII characters.
const std::string hexToAscii(const std::string &hex) noexcept {
  std::string result;
  // Reserve space for the result string to optimize memory allocation.
  result.reserve(hex.length() / 2); // Each ASCII character is represented by 2 hex digits.

  // Process each pair of hex digits in the string.
  for (size_t i = 0; i < hex.length(); i += 2) {
    // Convert the 2 hex digits to an ASCII character.
    std::string byteString = hex.substr(i, 2);
    char asciiChar = static_cast<char>(std::stoi(byteString, nullptr, 16)); // Convert hex to decimal.
    result += asciiChar;                                                    // Append the ASCII character to the result.
  }

  // Return the final ASCII string representation.
  return result;
}

// Converts a binary string representation to a hexadecimal string representation.
const std::string binaryToHex(const std::string &binary) noexcept {
  static const char hexDigits[] = "0123456789ABCDEF";
  std::string result;
  // Reserve space for the result string to optimize memory allocation.
  result.reserve((binary.length() + 3) / 4); // Each 4 bits (1 nibble) is represented by 1 hex digit.

  // Process each 4 bits (1 nibble) in the binary string.
  for (size_t i = 0; i < binary.length(); i += 4) {
    // Convert the 4-bit binary substring to a hex digit.
    std::string nibbleString = binary.substr(i, 4);
    int hexValue = std::stoi(nibbleString, nullptr, 2); // Convert binary to decimal.
    result += hexDigits[hexValue];                      // Append the hex digit to the result.
  }

  // Return the final hexadecimal string representation.
  return result;
}

// Converts a hexadecimal string representation to a binary string representation.
const std::string hexToBinary(const std::string &hex) noexcept {
  std::string result;
  // Reserve space for the result string to optimize memory allocation.
  result.reserve(hex.length() * 4); // Each hex digit is represented by 4 bits.

  // Process each hex digit in the string.
  for (char hexChar : hex) {
    // Convert the hex digit to its binary representation.
    int hexValue = (hexChar >= '0' && hexChar <= '9') ? hexChar - '0' : (hexChar - 'A' + 10);
    // Append the 4-bit binary representation to the result.
    for (int i = 3; i >= 0; --i) {
      result += (hexValue & (1 << i)) ? '1' : '0'; // Check each bit and append '1' or '0'.
    }
  }

  // Return the final binary string representation.
  return result;
};

// Function to compute the factorial of a number 'n' down to 'r'
// Example: factorial(5, 2) computes 5 * 4 * 3
static constexpr std::size_t factorial(std::size_t n, std::size_t r) {
  // If n is 0 or r is greater than or equal to n, return 0
  if (n == 0 || r >= n)
    return 0;

  // Initialize result to n
  std::size_t result = n;

  // Multiply result by each number from n-1 down to r+1
  while (--n > r) {
    result *= n;
  }

  // Return the computed factorial
  return result;
}

// Function to compute the number of permutations (with or without repetition)
// If rep is true, compute n^r (permutations with repetition)
// If rep is false, compute n! / (n-r)! (permutations without repetition)
static constexpr std::size_t permutationFunction(std::size_t n, std::size_t r, bool rep) { return rep ? raise(n, r) : factorial(n, r); }

// Function to compute the number of combinations (with or without repetition)
// If rep is true, compute (n+r-1)! / (r! * (n-1)!)
// If rep is false, compute n! / (r! * (n-r)!)
static constexpr std::size_t combinationFunction(std::size_t n, std::size_t r, bool rep) {
  // If n or r is 0, return 0
  if (n == 0 || r == 0)
    return 0;

  // Compute the numerator
  std::size_t num = rep ? factorial(n + r - 1, 1) : factorial(n, n - r);

  // Compute the denominator
  std::size_t denom = factorial(r, 1) * (rep ? factorial(n - 1, 1) : 1);

  // Return the result of num / denom
  return num / denom;
}

// Function to compute the number of combinations (with or without repetition)
// This version includes an output parameter 'o' for additional information
static constexpr std::size_t combinationFunction(std::size_t n, std::size_t r, bool rep, unsigned int &o) {
  // If n is 0, n is less than r, or r is 0, return 0
  if (n == 0 || n < r || r == 0)
    return 0;

  // Compute the numerator
  std::size_t num = rep ? factorial(n + r - 1, 1) : factorial(n, n - r);

  // Compute the denominator
  std::size_t denom = factorial(r, 1) * (rep ? factorial(n - 1, 1) : 1);

  // Return the result of num / denom
  return num / denom;
};

namespace Matrices {

// Template structure to define a matrix with rows (rows) and columns (cols)
template <int Rows, int Cols>
struct MatrixDimensions {
  std::uint32_t rows{Rows}; // Number of rows
  std::uint32_t cols{Cols}; // Number of columns

  // Default constructor
  MatrixDimensions() = default;

  // Parameterized constructor to initialize matrix dimensions
  MatrixDimensions(std::uint32_t rows, std::uint32_t cols) : rows(Rows), cols(Cols) {}
};

// Template specialization to check if a type is a container (default is false)
template <typename T>
struct IsContainer : std::false_type {};

// Template specialization to check if a std::vector is a container (true)
template <typename T, typename Alloc>
struct IsContainer<std::vector<T, Alloc>> : std::true_type {};

// Template specialization to check if a std::deque is a container (true)
template <typename T, typename Alloc>
struct IsContainer<std::deque<T, Alloc>> : std::true_type {};

// Template specialization to check if a std::list is a container (true)
template <typename T, typename Alloc>
struct IsContainer<std::list<T, Alloc>> : std::true_type {};

// Template structure to define a matrix with elements of type T and fixed dimensions (rows, cols)
template <typename T, std::uint32_t Rows, std::uint32_t Cols>
struct Matrix {
  std::vector<std::deque<T>> elements{}; // 2D vector to store matrix elements
  MatrixDimensions<Rows, Cols> dimensions{Rows, Cols}; // Matrix dimensions

  // Default constructor
  Matrix() = default;

  // Parameterized constructor to initialize matrix elements and dimensions
  Matrix(std::vector<std::deque<T>> elements, MatrixDimensions<Rows, Cols> dimensions)
      : elements(std::move(elements)), dimensions(dimensions) {}

  // Assignment operator overload for Matrix
  template <typename U, std::uint32_t R, std::uint32_t C>
  Matrix<U, R, C>& operator=(const Matrix<U, R, C>& other) {
    std::cout << "Copy constructor overload\n";
    elements = other.elements; // Copy elements
    dimensions = other.dimensions; // Copy dimensions
    return *this;
  }

  // Equality operator overload to compare matrices
  bool operator==(const Matrix<T, Rows, Cols>& other) const {
    std::cout << "Equality comparison\n";
    // Compare dimensions first
    if (dimensions.rows != other.dimensions.rows || dimensions.cols != other.dimensions.cols) {
      return false;
    }
    // Compare elements
    for (std::uint32_t i = 0; i < dimensions.rows; ++i) {
      for (std::uint32_t j = 0; j < dimensions.cols; ++j) {
        if (elements[i][j] != other.elements[i][j]) {
          return false;
        }
      }
    }
    return true;
  }

  // Distribute elements from a vector into the matrix row-wise
  void distributeRowWise(const std::vector<T>& input) {
    reset(); // Reset the matrix
    std::size_t index = 0;
    for (std::uint32_t i = 0; i < dimensions.rows && index < input.size(); ++i) {
      for (std::uint32_t j = 0; j < dimensions.cols && index < input.size(); ++j) {
        elements[i][j] = input[index++];
      }
    }
  }

  // Distribute elements from a 2D vector into the matrix
  void distributeFrom2DVector(const std::vector<std::deque<T>>& input) {
    reset(); // Reset the matrix
    for (std::uint32_t i = 0; i < std::min(dimensions.rows, static_cast<std::uint32_t>(input.size())); ++i) {
      for (std::uint32_t j = 0; j < std::min(dimensions.cols, static_cast<std::uint32_t>(input[i].size())); ++j) {
        elements[i][j] = input[i][j];
      }
    }
  }

  // Set new dimensions for the matrix and resize the internal vectors
  template <typename SizeType>
  void setDimensions(SizeType newRows, SizeType newCols) {
    dimensions.rows = newRows;
    dimensions.cols = newCols;
    elements.resize(newRows);
    for (auto& row : elements) {
      row.resize(newCols);
    }
  }

  // Reset the matrix to its original dimensions
  void reset() {
    elements.resize(dimensions.rows);
    for (auto& row : elements) {
      row.resize(dimensions.cols);
    }
  }

  // Transpose the matrix (swap rows and columns)
  void transpose() {
    std::vector<std::deque<T>> transposedElements(dimensions.cols, std::deque<T>(dimensions.rows));
    for (std::uint32_t i = 0; i < dimensions.rows; ++i) {
      for (std::uint32_t j = 0; j < dimensions.cols; ++j) {
        transposedElements[j][i] = elements[i][j];
      }
    }
    elements.swap(transposedElements);
    std::swap(dimensions.rows, dimensions.cols);
  }

  // Addition operator overload for matrices
  template <typename U, std::uint32_t R, std::uint32_t C>
  Matrix<U, R, C> operator+(const Matrix<U, R, C>& other) const {
    if (other.dimensions.rows != dimensions.rows || other.dimensions.cols != dimensions.cols) {
      throw std::runtime_error("Matrices do not match!");
    }
    Matrix<U, R, C> result;
    result.setDimensions(dimensions.rows, dimensions.cols);
    for (std::uint32_t i = 0; i < dimensions.rows; ++i) {
      for (std::uint32_t j = 0; j < dimensions.cols; ++j) {
        result.elements[i][j] = elements[i][j] + other.elements[i][j];
      }
    }
    return result;
  }

  // Get the size of the matrix (total number of elements)
  std::size_t size() const {
    return dimensions.rows * dimensions.cols;
  }

  void print() const {
    for (std::uint32_t i = 0; i < dimensions.rows; ++i) {
        for (std::uint32_t j = 0; j < dimensions.cols; ++j) {
            std::cout << elements[i][j] << " "; // Print each element followed by a space
        }
        std::cout << std::endl; // New line after each row
    }
};

};

static void runTests() {
    // Test 1: Default Constructor
    Matrix<int, 2, 3> matrix1;
    assert(matrix1.size() == 6); // 2 rows * 3 cols
    std::cout << "Test 1 (Default Constructor): Passed\n";

    // Test 2: Parameterized Constructor
    Matrix<int, 2, 3> matrix2({{1, 2, 3}, {4, 5, 6}}, MatrixDimensions<2, 3>());
    assert(matrix2.elements[0][0] == 1);
    assert(matrix2.elements[1][2] == 6);
    std::cout << "Test 2 (Parameterized Constructor): Passed\n";

    // Test 3: Assignment Operator
    Matrix<int, 2, 3> matrix3;
    matrix3 = matrix2;
    assert(matrix3 == matrix2);
    std::cout << "Test 3 (Assignment Operator): Passed\n";

    // Test 4: Equality Operator
    assert(matrix2 == matrix3);
    std::cout << "Test 4 (Equality Operator): Passed\n";

    // Test 5: Distribute Row Wise
    std::vector<int> input1 = {7, 8, 9, 10, 11, 12};
    matrix1.distributeRowWise(input1);
    assert(matrix1.elements[0][0] == 7);
    assert(matrix1.elements[1][2] == 12);
    std::cout << "Test 5 (Distribute Row Wise): Passed\n";

    // Test 6: Distribute From 2D Vector
    std::vector<std::deque<int>> input2 = {{13, 14, 15}, {16, 17, 18}};
    matrix1.distributeFrom2DVector(input2);
    assert(matrix1.elements[0][0] == 13);
    assert(matrix1.elements[1][1] == 17);
    std::cout << "Test 6 (Distribute From 2D Vector): Passed\n";

    // Test 7: Set Dimensions
    matrix1.setDimensions(3, 2);
    assert(matrix1.dimensions.rows == 3);
    assert(matrix1.dimensions.cols == 2);
    std::cout << "Test 7 (Set Dimensions): Passed\n";

    // Test 8: Reset
    matrix1.reset();
    assert(matrix1.elements.size() == 3);
    assert(matrix1.elements[0].size() == 2);
    std::cout << "Test 8 (Reset): Passed\n";

    // Test 9: Transpose
    matrix1.distributeRowWise(input1);
    matrix1.transpose();
    assert(matrix1.elements[0][0] == 7); // After transpose, this should be the first element
    assert(matrix1.elements[1][0] == 8); // This should be the second element
    std::cout << "Test 9 (Transpose): Passed\n";

    // Test 10: Addition Operator
    Matrix<int, 2, 3> matrix4({{1, 1, 1}, {1, 1, 1}}, MatrixDimensions<2, 3>());
    Matrix<int, 2, 3> result = matrix2 + matrix4;
    assert(result.elements[0][0] == 2);
    assert(result.elements[1][2] == 7);
    std::cout << "Test 10 (Addition Operator): Passed\n";

    // Test 11: Size
    assert(matrix2.size() == 6);
    std::cout << "Test 11 (Size): Passed\n";

    // Test 12: Exception on Addition of Different Dimensions
    Matrix<int, 2, 2> matrix5({{1, 2}, {3, 4}}, MatrixDimensions<2, 2>());
    try {
        matrix2 + matrix5; // This should throw an exception
        std::cout << "Test 12 (Exception on Addition of Different Dimensions): Failed\n";
    } catch (const std::runtime_error&) {
        std::cout << "Test 12 (Exception on Addition of Different Dimensions): Passed\n";
    }
}

} // namespace Matrices

}; // namespace Math
