#include "math.hpp"

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
  constexpr std::size_t threshold = 60;                                   // prime value threshold, computation stops at this value
  constexpr std::size_t PMAX = 20;                                        // max entries for prime computation result
  std::array<int, PMAX> primes = Math::primeComputation<PMAX>(threshold); // compute primes and store max 100 values
  std::vector<int> sievePrimes = Math::sieveOfEratosthenes(threshold);    // using sieve...
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
  constexpr std::array<int, fst> fibSeq = Math::fibonacciSequence<fst>();

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

  std::array<int, 5> setA({1, 2, 3, 4, 5});
  std::array<int, setA.size()> setB({1, 2, 3, 4, 5});

  Math::Set<int, setA.size()> newSet(setA);
  std::cout << "SetA cardinality: " << newSet.cardinality() << "\n";
  std::cout << "SetA == SetB? " << std::boolalpha << newSet.equal(setB) << "\n";
  std::cout << "SetA subset of SetB? " << std::boolalpha << newSet.subsetOf(setB) << "\n";
  std::cout << "SetA proper Subset of SetB? " << std::boolalpha << newSet.properSubsetOf(setB) << "\n";

  setA[1] = 1;
  setA[2] = 1;
  newSet.alloc(setA);
  std::cout << "SetA == SetB? " << std::boolalpha << newSet.equal(setB) << "\n";
  std::cout << "SetA subset of SetB? " << std::boolalpha << newSet.subsetOf(setB) << "\n";
  std::cout << "SetA proper Subset of SetB? " << std::boolalpha << newSet.properSubsetOf(setB) << "\n";

  std::set<int> unionSet = newSet.unionOf(setB);
  std::set<int> intersectionSet = newSet.intersectionOf(setB);
  std::set<int> complementSet = newSet.complementOf(setB);

  std::cout << "Union of: {";
  for (const auto e : setA)
    std::cout << e << ", ";
  std::cout << "} and {";
  for (const auto e : setB)
    std::cout << e << ", ";
  std::cout << "} => {";
  for (const auto e : unionSet)
    std::cout << e << ", ";
  std::cout << "}\n";

  std::cout << "Intersection of: {";
  for (const auto e : setA)
    std::cout << e << ", ";
  std::cout << "} and {";
  for (const auto e : setB)
    std::cout << e << ", ";
  std::cout << "} => {";
  for (const auto e : intersectionSet)
    std::cout << e << ", ";
  std::cout << "}\n";

  std::cout << "Complement of: {";
  for (const auto e : setA)
    std::cout << e << ", ";
  std::cout << "} from U {";
  for (const auto e : setB)
    std::cout << e << ", ";
  std::cout << "} => {";
  for (const auto e : complementSet)
    std::cout << e << ", ";
  std::cout << "}\n";

  std::string asciiText("AB");
  std::string ascii_to_binary(Math::asciiToBinary(asciiText));
  std::cout << asciiText << " in binary = " << ascii_to_binary << "\n";

  std::string ascii_to_hex(Math::asciiToHex(asciiText));
  std::cout << asciiText << " in Hex = " << ascii_to_hex << "\n";

  std::string binary_to_ascii(Math::binaryToAscii(ascii_to_binary));
  std::string hex_to_ascii(Math::hexToAscii(ascii_to_hex));

  std::cout << ascii_to_binary << " to ASCII = " << binary_to_ascii << "\n";
  std::cout << ascii_to_hex << " to ASCII = " << hex_to_ascii << "\n";

  std::string binary_to_hex(Math::binaryToHex(ascii_to_binary));
  std::string hex_to_binary(Math::hexToBinary(ascii_to_hex));

  std::cout << ascii_to_binary << " to Hex = " << binary_to_hex << "\n";
  std::cout << ascii_to_hex << " to binary = " << hex_to_binary << "\n";

  constexpr std::size_t f1{5}, f2{2};
  constexpr std::size_t factorization = Math::factorial(f1, f2);
  std::cout << "Factorial of " << f1 << ", " << f2 << " = " << factorization << "\n";

  constexpr bool REP = true;
  constexpr std::size_t permutation = Math::permutationFunction(f1, f2, REP);
  std::cout << "Permutation of " << f1 << ", " << f2 << " = " << permutation << "\n";

  constexpr std::size_t combination = Math::combinationFunction(f1, f2, REP);
  std::cout << "Combination of " << f1 << ", " << f2 << " = " << combination << "\n";

  Math::Matrices::runTests();

  return 0;
}
