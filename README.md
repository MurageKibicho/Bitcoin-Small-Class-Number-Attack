# Bitcoin-Small-Class-Number-Attack
Satoshi chose an elliptic curve with a bizarre prime number. It splits *unusually* well over √-3. 
The √-3 is special because it has primes and unique factorizations.
Math guys would say √-3 "has a class number of 1".

There are nine special square roots of imaginary numbers, called Heegner numbers, and math researchers have long speculated the existence of an attack built around these numbers.
Trump defunded my math PhD, now I identify as an internet math autist while I figure things out. I might as well explote this attack while I'm funemployed.

### Different Levels of Understanding Explanation

#### For Absolute Beginners (no math background, non-technical)
We use imaginary numbers to attack imaginary internet money.

#### Super Technical (moderate math background, technical)
We *attempt* an attack on secp256k1's field characteristic over both of its unique factorization domains, the regular integers and the ring of Eisenstein integers.
This is experimental mathematics and our tests involve:
1. Index calculus to find an Eisenstein integer, somewhat similar to the fundamental solution of the Pell equation as we saw here.
2. A Xedni calculus approach that leverages the fact that Hyperbolas encode an elliptic curve's addition property.
3. Division polynomials encode elliptic curve information and permit multiplication. We can exploit the factorization of p-1.

Here's our thesis: if we f around then we shall find out. All we need are the field characteristic's logarithms to get started.

## Getting Started
We provide a C library that runs best on CPU, not GPU. It's tes# ted on Linux. One can run locally, on Google Colab or a cloud provider like Runpod.

### Running Locally
This library requires LibGMP, FLINT and MPFR. 
We provide a Linux Bash script for a fresh installation.

### Runpod Tutorial

First install m4 
```
sudo apt-get update
sudo apt-get install -y m4
```

Then follow the 

### Running on Googe Colab






