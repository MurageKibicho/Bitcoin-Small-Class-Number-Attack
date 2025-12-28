# Bitcoin-Small-Class-Number-Attack
Satoshi chose an elliptic curve with a bizarre prime number. It splits *unusually* well over √-3. 

The √-3 is special because it has primes and unique factorizations and math guys would say "it has a class number of 1".

There are nine special square roots of imaginary numbers, called Heegner numbers, and math researchers have speculated on the existence of an attack built around these numbers.

I might as well explore this attack while I'm funemployed since Trump defunded my math PhD. 

### Different Levels of Understanding

#### For Absolute Beginners (no math background, non-technical)
We use imaginary numbers to attack imaginary internet money.

#### Super Technical (decent math background, technical)
We *attempt* an attack on secp256k1's field characteristic over both of its unique factorization domains, the regular integers and the ring of Eisenstein integers.
This is experimental mathematics and our tests involve:
1. Index calculus to find an Eisenstein integer, somewhat similar to the fundamental solution of the Pell equation as we saw here.
2. A Xedni calculus approach that leverages the fact that Hyperbolas encode an elliptic curve's addition property.
3. Division polynomials encode elliptic curve information and permit multiplication. We can exploit the factorization of p-1.

Here's our thesis: if we f around then we shall find out. All we need are the field characteristic's logarithms to get started.

## Getting Started
We provide a C library that runs best on CPU, not GPU. It's tested on Linux. One can run locally, on Google Colab or a cloud provider like Runpod.

### Running on Google Colab
Here's the Google Colab notebook.

### Running Locally or on Runpod
This library requires LibGMP, FLINT and MPFR. 
We provide a Bash script for a fresh installation.

1. The `Download.sh` is located in the folder `Step 0: Getting Started`. You can skip installation if you already installed FLINT 3.4.0.
2. Clone the repo.
```
git clone https://github.com/MurageKibicho/Bitcoin-Small-Class-Number-Attack.git
```
3.  cd into the folder.
```
cd 'Step 0: Getting Started'
```

4. Ensure m4 is available
```
sudo apt-get update
sudo apt-get install -y m4
```

5. Make Download.sh executable, and run Download.sh
```
chmod +x Download.sh

./Download.sh
```

6. Run sample program
   
I installed Flint in home. Modify this command to match your build

```
clear && gcc main.c -o m.o   -I/home/Bitcoin-Small-Class-Number-Attack/build/local/include   -L/home/Bitcoin-Small-Class-Number-Attack/build/local/lib   -lflint -lmpfr -lgmp -lm   -Wl,-rpath,/home/Bitcoin-Small-Class-Number-Attack/build/local/lib   && ./m.o
```






