#define STB_DS_IMPLEMENTATION
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <dirent.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpzi.h>
#include <flint/fmpq.h>
#include <flint/fmpz_factor.h>
#include "stb_ds.h"
#define fmp_printf(label, var)printf("%s: ", label);fmpz_print(var);printf("\n");
#define MAX_FILENAME_LENGTH 2056
#define INFINITE_LOOP_CHECK_SUM 1000
#define MPF_DEFAULT_PRECISION 512
//Run: clear && gcc Relations.c -o m.o -lm -lgmp -lmpfr -lflint && ./m.o  

typedef struct relation_struct *Relation;
typedef struct discrete_log_system_struct *Bitcoin;
typedef struct eisenstein_integer_struct *fmpzE;
struct eisenstein_integer_struct
{
	fmpz_t a;
	fmpz_t b;
	fmpz_t norm;
	int exponent;
};
struct discrete_log_system_struct
{
	fmpz_t prime;
	fmpz_t rootOfUnity0;
	fmpz_t rootOfUnity1;
	fmpz_t T;
	fmpz_t V;
	fmpz_t heegnerNumber;
	fmpz_t heegnerNumberRoot;
	fmpz_t ell;
};
struct relation_struct
{
	int aSign;
	int reconstructionSuccess;
	bool validRelation;
	int posVCountA;
	int negVCountA;
	int posVCountB;
	int negVCountB;
	int bestInt;
	int bestComplex;
	fmpz_t k;
	fmpz_t rhs1;
	fmpz_t a;
	fmpz_t b;
	
	fmpz_factor_t aNum;
	fmpz_factor_t aDen;
	fmpz_factor_t bNum;
	fmpz_factor_t bDen;
	fmpzE *aEis;
	fmpzE *bEis;
};


fmpzE fmpzE_init()
{
	fmpzE e = malloc(sizeof(struct eisenstein_integer_struct));
	e->exponent = 0;
	fmpz_init(e->a);
	fmpz_init(e->b);
	fmpz_init(e->norm);
	return e;
}

void fmpzE_clear(fmpzE e)
{
	fmpz_clear(e->a);
	fmpz_clear(e->b);
	fmpz_clear(e->norm);
	free(e);
}

void fmpzE_add(fmpzE r, fmpzE a, fmpzE b)
{
	fmpz_add(r->a, a->a, b->a);
	fmpz_add(r->b, a->b, b->b);
}

void fmpzE_sub(fmpzE r, fmpzE a, fmpzE b)
{
	fmpz_sub(r->a, a->a, b->a);
	fmpz_sub(r->b, a->b, b->b);
}

void fmpzE_mul(fmpzE r, fmpzE a, fmpzE b)
{
	fmpz_t temp0;fmpz_init(temp0);
	fmpz_mul(r->a, a->a, b->a);
	fmpz_mul(temp0, a->b, b->b);
	fmpz_sub(r->a,r->a,temp0);
	
	
	fmpz_set(r->b, temp0);
	fmpz_neg(r->b, r->b);	
	fmpz_mul(temp0, a->a, b->b);
	fmpz_add(r->b,r->b,temp0);
	fmpz_mul(temp0, a->b, b->a);
	fmpz_add(r->b,r->b,temp0);
	fmpz_clear(temp0);
}

void fmpzE_conjugate(fmpzE result, fmpzE e)
{
	fmpz_sub(result->a, e->a, e->b);
	fmpz_neg(result->b, e->b);
}

void fmpzE_norm(fmpz_t norm, fmpz_t a, fmpz_t b)
{
	fmpz_t temp0;
	fmpz_init(temp0);
	fmpz_mul(temp0, a,a);//a^2
	fmpz_mul(norm, b,b);//b^2
	fmpz_add(norm, norm, temp0);
	
	fmpz_mul(temp0, a, b);//ab
	fmpz_sub(norm, norm, temp0);
	fmpz_clear(temp0);
}

void fmpz_set_mpf_round(fmpz_t result, mpf_t x)
{
	mpf_t f_tmp,frac;
	mpf_init(f_tmp);mpf_init(frac);
	
	mpf_floor(f_tmp, x);
	mpf_sub(frac, x, f_tmp);
	if(mpf_cmp_d(frac, 0.5) < 0) 
	{
		fmpz_set_mpf(result, f_tmp); 
	}
	else
	{
		mpf_ceil(f_tmp, x); 
		fmpz_set_mpf(result, f_tmp);
	}
	mpf_clear(f_tmp);mpf_clear(frac);
	
}

void fmpzE_set_fmpzE(fmpzE result, fmpzE e)
{
	fmpz_set(result->a, e->a);
	fmpz_set(result->b, e->b);
	fmpz_set(result->norm, e->norm);
}
void fmpzE_divide(fmpzE r, fmpzE a, fmpzE b)
{
	fmpzE bConjugate = fmpzE_init();
	fmpzE num = fmpzE_init();
	fmpzE_conjugate(bConjugate, b);
	
	fmpzE_mul(num, a, bConjugate);
	fmpzE_norm(b->norm, b->a, b->b);
	
	mpf_t f_NumA, f_NumB, f_norm;
	mpf_init(f_NumA);mpf_init(f_NumB);mpf_init(f_norm);

	fmpz_get_mpf(f_NumA, num->a);
	fmpz_get_mpf(f_NumB, num->b);
	fmpz_get_mpf(f_norm, b->norm);

	mpf_div(f_NumA, f_NumA, f_norm);
	mpf_div(f_NumB, f_NumB, f_norm);

	
	fmpz_set_mpf_round(r->a, f_NumA);
	fmpz_set_mpf_round(r->b, f_NumB);
	
	mpf_clear(f_NumA);mpf_clear(f_NumB);mpf_clear(f_norm);
	fmpzE_clear(num);
	fmpzE_clear(bConjugate);
}

void fmpzE_mod(fmpzE r, fmpzE a, fmpzE b)
{
	fmpzE q = fmpzE_init();
	fmpzE temp = fmpzE_init();
	
	fmpzE_divide(q, a, b);  
	fmpzE_mul(temp, q, b);   
	fmpzE_sub(r, a, temp);  

	fmpzE_clear(q);
	fmpzE_clear(temp);
}

void fmpzE_gcd(fmpzE result, fmpzE alpha, fmpzE beta)
{
	fmpzE a = fmpzE_init();
	fmpzE b = fmpzE_init();
	fmpzE r = fmpzE_init();
	fmpzE zero = fmpzE_init();
	fmpz_set(a->a, alpha->a);fmpz_set(a->b, alpha->b);
	fmpz_set(b->a, beta->a);fmpz_set(b->b, beta->b);
	fmpz_zero(zero->a);fmpz_zero(zero->b);
	int infiniteLoopCheck = 0;
	while(!(fmpz_equal(b->a, zero->a) && fmpz_equal(b->b, zero->b)))
	{
		//r = a % b
		fmpzE_mod(r, a, b);

		//a = b
		fmpz_set(a->a, b->a);fmpz_set(a->b, b->b);

		//b = r
		fmpz_set(b->a, r->a);fmpz_set(b->b, r->b);
    		
    		if(infiniteLoopCheck >= INFINITE_LOOP_CHECK_SUM)
		{
			printf("Infinite loop in fmpzE_gcd\n");
			assert(infiniteLoopCheck < INFINITE_LOOP_CHECK_SUM);
		}
		infiniteLoopCheck += 1;
    	
    	}

	//result = a
	fmpz_set(result->a, a->a);
	fmpz_set(result->b, a->b);

	fmpzE_clear(a);
	fmpzE_clear(b);
	fmpzE_clear(r);
	fmpzE_clear(zero);
}

void fmpzE_powm_fmpzE(fmpzE result, fmpzE alpha, fmpz_t exponent, fmpzE modulo)
{
	//Case 0^0 = 0
	if(fmpz_is_zero(alpha->a) && fmpz_is_zero(alpha->b) && fmpz_is_zero(exponent))
	{
		fmpz_set_ui(result->a, 0);
		fmpz_set_ui(result->b, 0);
	}
	else
	{
		fmpzE temp0 = fmpzE_init();
		//Initialize result to 1
		fmpz_set_ui(result->a, 1);fmpz_set_ui(result->b, 0);
		int size = fmpz_sizeinbase(exponent, 2);
		for(int i = size - 1; i >= 0; i--)
		{
			fmpzE_mul(temp0, result, result);	
			fmpzE_set_fmpzE(result, temp0);
			
			fmpzE_mod(temp0, result, modulo);	
			fmpzE_set_fmpzE(result, temp0);
			if(fmpz_tstbit(exponent, i))
			{
				fmpzE_mul(temp0, result, alpha);	
				fmpzE_set_fmpzE(result, temp0);
			
				fmpzE_mod(temp0, result, modulo);	
				fmpzE_set_fmpzE(result, temp0);
			}
		}
		fmpzE_clear(temp0);
	}
}


void fmpzE_GCDSplit(fmpz_t prime, fmpz_t heegnerNumber, fmpzE complexFactor)
{
	fmpz_t temp2;
	fmpz_init(temp2);
	int foundSquareRoot = fmpz_sqrtmod(temp2, heegnerNumber, prime);
	assert(foundSquareRoot == 1);
	fmpzE pE = fmpzE_init();
	fmpzE gcdTemp = fmpzE_init();
	fmpzE divTemp = fmpzE_init();
	fmpz_set(pE->a, prime);fmpz_set_ui(pE->b, 0);
	fmpz_add_ui(gcdTemp->a, temp2, 1);fmpz_set_ui(gcdTemp->b, 2);

	//Set complex factor to gcd result
	fmpz_set(pE->a, prime);fmpz_set_ui(pE->b, 0);
	fmpzE_gcd(complexFactor, gcdTemp, pE);
	
	fmpzE_clear(pE);fmpzE_clear(gcdTemp);fmpzE_clear(divTemp);fmpz_clear(temp2);
}
 
void fmpzE_print(fmpzE e)
{
	printf("(");fmpz_print(e->a);printf(" + ");fmpz_print(e->b);printf(")\n");
}

void fmpzE_printPrimeFactors(fmpzE *primeFactors)
{
	for(size_t i = 0; i < arrlen(primeFactors); i++)
	{
		printf("(");fmpz_print(primeFactors[i]->a);printf(" + ");fmpz_print(primeFactors[i]->b);printf(") ^ %d, ",primeFactors[i]->exponent);	
	}
}

void fmpzE_freePrimeFactors(fmpzE *primeFactors)
{
	for(size_t i = 0; i < arrlen(primeFactors); i++)
	{
		fmpzE_clear(primeFactors[i]);
	}
	arrfree(primeFactors);
}

fmpzE *PrimeFactorization(fmpzE alpha)
{
	fmpzE *primeFactors = NULL;
	fmpzE pE = fmpzE_init();
	fmpzE gcdTemp = fmpzE_init();
	fmpzE divTemp = fmpzE_init();
	fmpzE alphaCopy = fmpzE_init();
	fmpzE modResult = fmpzE_init();
	fmpz_t mod, neg3,neg3Root,halfP;
	fmpz_factor_t normFactors;fmpz_factor_init(normFactors);
	fmpz_init(halfP);fmpz_init(mod);fmpz_init(neg3);fmpz_init(neg3Root);
	fmpz_set_si(neg3, -3);
	fmpzE_norm(alpha->norm, alpha->a, alpha->b);
	
	fmpz_factor(normFactors, alpha->norm);
	//Set alphaCopy
	fmpz_set(alphaCopy->a, alpha->a);fmpz_set(alphaCopy->b, alpha->b);fmpz_set(alphaCopy->norm, alpha->norm);
	//printf("Norm: ");PrintFactors(normFactors);
	for(slong i = 0; i < normFactors->num; i++)
	{
		fmpzE complexFactor = fmpzE_init();fmpzE complexConjugate = fmpzE_init();
		complexFactor->exponent = 0;complexConjugate->exponent = 0;
		fmpz_set(pE->a, &normFactors->p[i]);fmpz_set_ui(pE->b, 0);
		if(fmpz_cmp_ui(&normFactors->p[i], 3) == 0)
		{
			//if p[i] = 3, 1 - w is a factor
			fmpz_set_ui(complexFactor->a, 1);
			fmpz_set_si(complexFactor->b, -1);
			complexFactor->exponent = normFactors->exp[i];	
			for(int j = 0; j < normFactors->exp[i]; j++)
			{
				fmpzE_divide(divTemp, alphaCopy, complexFactor);
				fmpz_set(alphaCopy->a,divTemp->a);fmpz_set(alphaCopy->b,divTemp->b);	
			}
		}
		else
		{
			int mod3 = fmpz_mod_ui(mod, &normFactors->p[i], 3);
			assert(mod3 != 0);
			if(mod3 == 2)
			{
				//Exponent must be even, if odd sth is wrong
				assert(normFactors->exp[i] % 2 == 0);
				//if p[i] % 3 == 2, remove two factors of p[i]
				fmpz_set(complexFactor->a, &normFactors->p[i]);
				fmpz_set_ui(complexFactor->b, 0);
				complexFactor->exponent = normFactors->exp[i] / 2;
				for(int j = 0; j < normFactors->exp[i]; j += 2)
				{
					fmpzE_divide(divTemp, alphaCopy, pE);
					fmpz_set(alphaCopy->a,divTemp->a);fmpz_set(alphaCopy->b,divTemp->b);			
				}
			}
			else
			{
				int foundSquareRoot = fmpz_sqrtmod(neg3Root, neg3, &normFactors->p[i]);
				//fmpz_sub(neg3Root,&normFactors->p[i],neg3Root);		
				//fmpz_neg(neg3Root, neg3Root);
				//printf("neg3Root: ");fmpz_print(neg3Root);printf(": ");fmpz_print(&normFactors->p[i]);printf("\n");

				assert(foundSquareRoot == 1);
				fmpz_add_si(gcdTemp->a, neg3Root, 1);fmpz_set_si(gcdTemp->b, 2);
				
				//Set complex factor to gcd result
				fmpzE_gcd(complexFactor, gcdTemp, pE);
				//Find Complex conjugate
				fmpzE_conjugate(complexConjugate, complexFactor);
				fmpzE_norm(complexFactor->norm,complexFactor->a,complexFactor->b);
				for(int j = 0; j < normFactors->exp[i]; j++)
				{
					//Test if complex factor divides alphaCopy
					fmpzE_mod(modResult, alphaCopy, complexFactor);
					if(fmpz_cmp_ui(modResult->a, 0) == 0 && fmpz_cmp_ui(modResult->b, 0) == 0)
					{
						fmpzE_divide(divTemp, alphaCopy, complexFactor);
						fmpz_set(alphaCopy->a,divTemp->a);fmpz_set(alphaCopy->b,divTemp->b);			
						complexFactor->exponent += 1;
					}
					else
					{
						fmpzE_divide(divTemp, alphaCopy, complexConjugate);
						fmpz_set(alphaCopy->a,divTemp->a);fmpz_set(alphaCopy->b,divTemp->b);
						complexConjugate->exponent += 1;
					}
				}
				
			}
		}
		if(complexConjugate->exponent > 0){arrput(primeFactors, complexConjugate);}else{fmpzE_clear(complexConjugate);}
		if(complexFactor->exponent > 0){arrput(primeFactors, complexFactor);}else{fmpzE_clear(complexFactor);}
	}
	//Push alphaCopy as final factor if not unit
	fmpzE_norm(alphaCopy->norm, alphaCopy->a, alphaCopy->b);
	alphaCopy->exponent = 1;
	arrput(primeFactors, alphaCopy);

	fmpz_clear(halfP);fmpz_clear(mod);fmpz_clear(neg3);fmpz_clear(neg3Root);
	fmpz_factor_clear(normFactors);
	fmpzE_clear(modResult);fmpzE_clear(pE);fmpzE_clear(gcdTemp);fmpzE_clear(divTemp);
	return primeFactors;	
}


Bitcoin CreateBitcoin()
{
	char *primeNumberHexadecimal = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F";	
	char *factor4Base10 = "205115282021455665897114700593932402728804164701536103180137503955397371";	
	Bitcoin bitcoin = malloc(sizeof(struct discrete_log_system_struct));
	fmpz_init(bitcoin->prime);
	fmpz_init(bitcoin->rootOfUnity0);	
	fmpz_init(bitcoin->rootOfUnity1);	
	fmpz_init(bitcoin->T);	
	fmpz_init(bitcoin->V);	
	fmpz_init(bitcoin->heegnerNumber);	
	fmpz_init(bitcoin->heegnerNumberRoot);	
	fmpz_init(bitcoin->ell);
	
	fmpz_set_str(bitcoin->prime, primeNumberHexadecimal, 16);
	fmpz_set_si(bitcoin->heegnerNumber, -3);
	fmpz_set_str(bitcoin->rootOfUnity0, "55594575648329892869085402983802832744385952214688224221778511981742606582254",10);
	fmpz_set_str(bitcoin->rootOfUnity1, "60197513588986302554485582024885075108884032450952339817679072026166228089408",10);
	fmpz_set_str(bitcoin->T, "367917413016453100223835821029139468249",10);
	fmpz_set_str(bitcoin->V, "303414439467246543595250775667605759171",10);
	fmpz_set_str(bitcoin->ell,factor4Base10, 10);
	return bitcoin;	
}

void DestroyBitcoin(Bitcoin bitcoin)
{
	fmpz_clear(bitcoin->prime);
	fmpz_clear(bitcoin->rootOfUnity0);	
	fmpz_clear(bitcoin->rootOfUnity1);	
	fmpz_clear(bitcoin->T);	
	fmpz_clear(bitcoin->V);	
	fmpz_clear(bitcoin->heegnerNumber);	
	fmpz_clear(bitcoin->heegnerNumberRoot);	
	fmpz_clear(bitcoin->ell);
	free(bitcoin);
}

void PrintFactors(fmpz_factor_t factors)
{
	for(slong i = 0; i < factors->num; i++)
	{
		int size = fmpz_sizeinbase(&factors->p[i], 2);
		printf("(|%d| ", size);
		fmpz_print(&factors->p[i]); printf(" ^ ");
		fmpz_print(&factors->exp[i]); printf("), ");
	}
	printf("\n");
}

Relation CreateRelation()
{
	Relation relation = malloc(sizeof(struct relation_struct));
	fmpz_init(relation->k);
	fmpz_init(relation->rhs1);
	fmpz_init(relation->a);
	fmpz_init(relation->b);
	fmpz_factor_init(relation->aNum);
	fmpz_factor_init(relation->aDen);
	fmpz_factor_init(relation->bNum);
	fmpz_factor_init(relation->bDen);
	relation->aEis  = NULL;
	relation->bEis  = NULL;
	relation-> aSign = 0;
	relation->reconstructionSuccess = 0;
	relation->validRelation = false;
	relation->posVCountA = 0;
	relation->negVCountA = 0;
	relation->posVCountB = 0;
	relation->negVCountB = 0;	
	return relation;
}

void PrintRelation(Relation relation)
{
	printf("Rational (%d)\n",relation->aSign);
	fmp_printf("a: ", relation->a);
	fmp_printf("b: ", relation->b);
	printf("posVCountA:(%d)\n", relation->posVCountA);
	printf("negVCountA:(%d)\n", relation->negVCountA);
	printf("posVCountB:(%d)\n", relation->posVCountB);
	printf("negVCountB:(%d)\n", relation->negVCountB);

	
	fmp_printf("k: ", relation->k);
	fmp_printf("rhs: ", relation->rhs1);
	printf("aNumerator Int :");PrintFactors(relation->aNum);
	printf("aDenominator Int :");PrintFactors(relation->aDen);
	printf("a Eisenstein :");fmpzE_printPrimeFactors(relation->aEis);printf("\n");
		
	printf("bNumerator Int :");PrintFactors(relation->bNum);
	printf("bDenominator Int :");PrintFactors(relation->bDen);
	printf("b Eisenstein :");fmpzE_printPrimeFactors(relation->bEis);printf("\n");
	
}


void ReconstructNumberWithoutV(fmpz_t result, fmpz_factor_t numerator, fmpz_factor_t denominator,fmpzE *eisenstein, fmpz_t rootOfUnity, fmpz_t prime)
{
	fmpz_t temp0, temp1,denom;
	fmpz_init(temp0);fmpz_init(temp1);fmpz_init(denom);
	fmpz_set_ui(temp0, 1);
	for(size_t i = 0; i < arrlen(eisenstein); i++)
	{
		// factor = a_i + b_i * S mod p
		fmpz_mul(temp1, eisenstein[i]->b, rootOfUnity);
		fmpz_add(temp1, temp1, eisenstein[i]->a);
		fmpz_mod(temp1, temp1, prime);

		// power = factor^{e_i} mod p
		fmpz_powm_ui(temp1, temp1, eisenstein[i]->exponent, prime);

		// rhs = rhs * power mod p
		fmpz_mul(temp0, temp0, temp1);
		fmpz_mod(temp0, temp0, prime);
	}

	//fmpz_mul(temp0, temp0, V);
	//fmpz_mod(temp0, temp0, prime);
	
	for(slong i = numerator->num-1; i >= 0; i--)
	{
		fmpz_powm_ui(temp1, &numerator->p[i], numerator->exp[i], prime);
		fmpz_mul(temp0, temp0, temp1);
		fmpz_mod(temp0, temp0, prime);
	}
	fmpz_mod(temp0, temp0, prime);
	
	fmpz_set_ui(denom, 1);
	for(slong i = denominator->num-1; i >= 0; i--)
	{
		fmpz_powm_ui(temp1, &denominator->p[i], denominator->exp[i], prime);
		fmpz_mul(denom, denom, temp1);
		fmpz_mod(denom, denom, prime);
	}
	fmpz_invmod(denom, denom, prime);
	fmpz_mul(temp0, temp0, denom);
	fmpz_mod(temp0, temp0, prime);
	fmpz_set(result, temp0);
	fmpz_clear(temp0);fmpz_clear(temp1);fmpz_clear(denom);
}

bool TestRelation(Bitcoin bitcoin, Relation relation)
{
	bool result = true;
	//Reconstruct a/b mod p
	fmpz_t temp0,temp1,temp2,recon,VA,VB;
	fmpz_init(temp0);fmpz_init(temp1);fmpz_init(temp2);fmpz_init(recon);fmpz_init(VA);fmpz_init(VB);
	//Set V values
	fmpz_neg(VA, bitcoin->V);
	fmpz_powm_ui(VA, VA, relation->negVCountA, bitcoin->prime);
	fmpz_powm_ui(temp0, bitcoin->V, relation->posVCountA, bitcoin->prime);
	fmpz_mul(VA,VA,temp0);fmpz_mod(VA, VA, bitcoin->prime);
	
	fmpz_neg(VB, bitcoin->V);
	fmpz_powm_ui(VB, VB, relation->negVCountB, bitcoin->prime);
	fmpz_powm_ui(temp0, bitcoin->V, relation->posVCountB, bitcoin->prime);
	fmpz_mul(VB,VB,temp0);fmpz_mod(VB, VB, bitcoin->prime);
	
	//Ensure a/b matches target
	fmpz_invmod(temp2, relation->b, bitcoin->prime);
	fmpz_mul_si(temp1, relation->a, relation->aSign);
	
	fmpz_mul(temp2, temp2, temp1);
	fmpz_mod(temp2, temp2, bitcoin->prime);
	
	//Reconstruct a and b without V
	ReconstructNumberWithoutV(temp0, relation->aNum, relation->aDen, relation->aEis, bitcoin->rootOfUnity1, bitcoin->prime);		
	ReconstructNumberWithoutV(temp1, relation->bNum, relation->bDen, relation->bEis, bitcoin->rootOfUnity1, bitcoin->prime);								

	fmpz_mul(temp0,VA,temp0);fmpz_mod(temp0, temp0, bitcoin->prime);
	fmpz_mul(temp1,VB,temp1);fmpz_mod(temp1, temp1, bitcoin->prime);
	
	fmp_printf("a with V: ", temp0);
	fmp_printf("b with V: ", temp1);
	fmp_printf("a", relation->a);
	fmp_printf("b", relation->b);
		
	if(fmpz_cmp(temp2, relation->rhs1) != 0)
	{
		printf("a/b reconstruction fail\n");
		result = false;
	}
	if(fmpz_cmp(temp0, relation->a) != 0)
	{
		printf("a reconstruction fail\n");
		fmp_printf("a with V: ", temp0);
		fmp_printf("b with V: ", temp1);
		fmp_printf("a", relation->a);
		fmp_printf("b", relation->b);
		result = false;
	}
	if(fmpz_cmp(temp1, relation->b) != 0)
	{
		printf("b reconstruction fail\n");
		fmp_printf("a with V: ", temp0);
		fmp_printf("b with V: ", temp1);
		fmp_printf("a", relation->a);
		fmp_printf("b", relation->b);
		result = false;
	}
	
	
	fmpz_clear(temp0);fmpz_clear(temp1);fmpz_clear(temp2);fmpz_clear(recon);fmpz_clear(VA);fmpz_clear(VB);
	return result;	
}

void DestroyRelation(Relation relation)
{
	for(size_t i = 0; i < arrlen(relation->bEis); i++)
	{
		fmpzE_clear(relation->bEis[i]);
	}
	for(size_t i = 0; i < arrlen(relation->aEis); i++)
	{
		fmpzE_clear(relation->aEis[i]);
	}
	arrfree(relation->aEis);
	arrfree(relation->bEis);
	fmpz_factor_clear(relation->aNum);
	fmpz_factor_clear(relation->aDen);
	fmpz_factor_clear(relation->bNum);
	fmpz_factor_clear(relation->bDen);
	fmpz_clear(relation->a);
	fmpz_clear(relation->b);
	fmpz_clear(relation->k);
	fmpz_clear(relation->rhs1);
	free(relation);
}



size_t GetFileSize(char *fileName)
{
	FILE *fp = fopen(fileName, "rb");assert(fp != NULL);
	fseek(fp, 0L, SEEK_END);
	size_t currentFileSize = ftell(fp);rewind(fp);
	fclose(fp);
	return currentFileSize;
}
void SortFileNames(char **fileNames, int fileCount)
{
	char temporary[MAX_FILENAME_LENGTH] = {0};
	for(int k = 0; k < fileCount; k++)
	{
		for(int j = k + 1; j < fileCount; j++)
		{
			if(strcasecmp(fileNames[k], fileNames[j]) > 0)
			{
				strcpy(temporary, fileNames[k]);
				strcpy(fileNames[k], fileNames[j]);
				strcpy(fileNames[j], temporary);
				memset(temporary,0,MAX_FILENAME_LENGTH);
			}
		}
	}
}
char *CreateStorageDirectory(char *path)
{
	char outputFolder[MAX_FILENAME_LENGTH] = {0};
	snprintf(outputFolder, sizeof(outputFolder), "%s",path);
	DIR *directory = opendir(outputFolder);
	if(directory == NULL)
	{
		int outputFolderCreationSuccess = mkdir(outputFolder,0777);assert(outputFolderCreationSuccess == 0);
	}
	else
	{
		closedir(directory);
	}
	char *result = calloc(MAX_FILENAME_LENGTH, sizeof(char));
	for(int i = 0; i < MAX_FILENAME_LENGTH; i++)
	{
		result[i] = outputFolder[i];
	}
	return result;
}
int CountNumberOfFilesInDirectory(char *path)
{
	int result = 0;
	DIR *directory = opendir(path);assert(directory != NULL);
	struct dirent *directoryEntry;
	while((directoryEntry = readdir(directory)) != NULL)
	{
		if(directoryEntry->d_name[0] != '.')
		{
			int fileNameLength = (int)sizeof(directoryEntry->d_name);assert(fileNameLength <= MAX_FILENAME_LENGTH);result++;
		}
	}
	closedir(directory);
	return result;
}
char **GetFileNamesInDirectory(char *path, int *fileCount)
{	
	/*Count files and ensure length <= 256*/
	*fileCount = CountNumberOfFilesInDirectory(path);
	/*Array to store file names*/
	char **result = malloc(*fileCount * sizeof(char*));
	for(int i = 0; i < *fileCount; i++)
	{
		result[i] = calloc(MAX_FILENAME_LENGTH, sizeof(char));
	}
	/*Get file names*/
	DIR *directory = opendir(path);assert(directory != NULL);
	struct dirent *directoryEntry;
	int currentIndex = 0;
	while((directoryEntry = readdir(directory)) != NULL)
	{
		if(directoryEntry->d_name[0] != '.')
		{
			assert(strlen(directoryEntry->d_name) < MAX_FILENAME_LENGTH);
			for(int i = 0; i < strlen(directoryEntry->d_name); i++)
			{
				result[currentIndex][i] = directoryEntry->d_name[i];
			}
			currentIndex++;
		}
	}
	SortFileNames(result, *fileCount);
	closedir(directory);
	return result;
}
void PrintFileNames(int fileCount, char **fileNames)
{
	for(int i = 0; i < fileCount; i++)
	{
		printf("%3d %s\n", i, fileNames[i]);
	}
}

void DestroyFileNames(int fileCount, char **fileNames)
{
	for(int i = 0; i < fileCount; i++)
	{
		free(fileNames[i]);
	}
	free(fileNames);
}


void ReadBatchFile(char *filename, fmpz_t k, fmpz_t rhs1)
{
	FILE *in = fopen(filename, "r");
	assert(in != NULL);
	char line[2048];
	while(fgets(line, sizeof(line), in))
	{
		//Remove trailing newline if present
		line[strcspn(line, "\n")] = 0;
		//Split at comma
		char *comma = strchr(line, ',');
		if(!comma)
		{
			fprintf(stderr, "Malformed line: %s\n", line);
			continue;
		}
		*comma = '\0';
		char *kStr = line;
		char *rhs1Str = comma + 1;

		// Convert strings to fmpz
		fmpz_set_str(k, kStr, 10);
		fmpz_set_str(rhs1, rhs1Str, 10);

	}
	fclose(in);
}


int FindLargestFactorSize(fmpz_factor_t factors, int *largePrimeCount, int *largestIndex, int factorBound)
{
	int largest = 0;
	*largePrimeCount = 0;
	*largestIndex = 0;
	for(slong i = 0; i < factors->num; i++)
	{
		int size = fmpz_sizeinbase(&factors->p[i], 2);
		if(size > largest)
		{
			*largestIndex = i;
			largest = size;
		}
		if(size > factorBound)
		{
			*largePrimeCount += 1;
		}
		//fmpz_print(&factors->p[i]); printf(" ^ ");
		//fmpz_print(&factors->exp[i]); printf(")");
	}
	return largest;
}


void GaussReduce(fmpz_t x1, fmpz_t y1, fmpz_t x2, fmpz_t y2)
{
	//Shortest solution is stored in x1,y1
	fmpz_t mu, dot12, dot11, temp0,tmp2,tmp3;
	fmpz_init(mu);fmpz_init(tmp2);fmpz_init(tmp3); fmpz_init(dot12); fmpz_init(dot11); fmpz_init(temp0);
	int infiniteLoopCheck = 0;
	while(1)
	{
		//mu = round((v1·v2)/(v1·v1))
		fmpz_mul(dot12, x1, x2);fmpz_mul(temp0, y1, y2);fmpz_add(dot12, dot12, temp0);
		fmpz_mul(dot11, x1, x1);fmpz_mul(temp0, y1, y1);fmpz_add(dot11, dot11, temp0);
		fmpz_tdiv_q(mu, dot12, dot11);

		// v2 = v2 - mu * v1
		fmpz_mul(temp0, mu, x1); fmpz_sub(x2, x2, temp0);
		fmpz_mul(temp0, mu, y1); fmpz_sub(y2, y2, temp0);

		// ||v2|| < ||v1||, swap
		fmpz_mul(temp0, x1, x1);
		fmpz_mul(tmp2, y1, y1); fmpz_add(temp0, temp0, tmp2);
		fmpz_mul(tmp2, x2, x2);
		fmpz_mul(tmp3, y2, y2); fmpz_add(tmp2, tmp2, tmp3);

		if(fmpz_cmp(tmp2, temp0) < 0){fmpz_swap(x1, x2);fmpz_swap(y1, y2);}
		else{break;}
		infiniteLoopCheck += 1;
		if(infiniteLoopCheck > INFINITE_LOOP_CHECK_SUM)
		{
			printf("Gauss reduce may be in infinite loop\n");
			assert(infiniteLoopCheck < INFINITE_LOOP_CHECK_SUM);
		}	 
	}

	fmpz_clear(mu);fmpz_clear(tmp2);fmpz_clear(tmp3); fmpz_clear(dot12); fmpz_clear(dot11); fmpz_clear(temp0);
}

bool VerifyLatticeSolution(fmpz_t x, fmpz_t y, fmpz_t T, fmpz_t V, fmpz_t p)
{
	fmpz_t temp0,t;
	fmpz_init(temp0);fmpz_init(t);

	//temp0 = V*x - T*y
	fmpz_mul(temp0, V, x);
	fmpz_mul(t, T, y);
	fmpz_sub(temp0, temp0, t);
	fmpz_mod(temp0, temp0, p);

	bool ok = (fmpz_cmp_ui(temp0, 0) == 0);
	fmpz_clear(temp0);
	fmpz_clear(t);
	return ok;
}

int CV_DT_SmoothnessTest(fmpz_t c, fmpz_t d, fmpz_t T, fmpz_t V, fmpz_t r, fmpz_t cV_min_dT, fmpz_factor_t factors, int targetSmallerBits)
{
	fmpz_t temp0,t;
	fmpz_init(temp0);fmpz_init(t);

	//temp0 = (V*c - T*d)/r
	fmpz_mul(temp0, V, c);fmpz_mul(t, T, d);fmpz_sub(temp0, temp0, t);
	fmpz_set(cV_min_dT, temp0);
	fmpz_divexact(temp0, temp0, r);
	
	int proved = -1;
	int largePrimeCount = 0;
	int largestIndex = 0;
	int factorResult = fmpz_factor_smooth(factors, temp0,targetSmallerBits, proved);
	int largestBits = FindLargestFactorSize(factors, &largePrimeCount, &largestIndex, targetSmallerBits);
	//printf("(cV-dT) / r: (%d %d : %d) ", largestBits, targetSmallerBits,factorResult);fmpz_print(temp0);printf("\n");
	fmpz_clear(temp0);fmpz_clear(t);
	return largestBits;
}

int DotProductSearch(Bitcoin bitcoin, fmpz_factor_t factors, fmpz_factor_t bestFactors, int testRange, int targetSmallerBits, fmpz_t targetPrime, fmpz_t temp0, fmpz_t c,fmpz_t  d, fmpz_t x0, fmpz_t y0, fmpz_t  x1,fmpz_t  y1)
{
	fmpz_t dot, norm1_sq, k_num, k_floor, k_ceil,vx,vy,gcdtmp;
	fmpz_init(dot); fmpz_init(norm1_sq);fmpz_init(k_num);fmpz_init(k_floor);fmpz_init(k_ceil);fmpz_init(vx);fmpz_init(vy);fmpz_init(gcdtmp);
	
	fmpz_mul(dot, x0, x1);
	fmpz_addmul(dot, y0, y1);
	fmpz_mul(norm1_sq, x1, x1);
	fmpz_addmul(norm1_sq, y1, y1);
	fmpz_neg(k_num, dot);
	fmpz_fdiv_q(k_floor, k_num, norm1_sq);
    	fmpz_cdiv_q(k_ceil,  k_num, norm1_sq);
    	
	int smallest = 8000000;
    	for(slong delta = -testRange; delta <= testRange; delta+=1)
	{
		//v = x0 + k*x1
		fmpz_set(k_num, k_floor);
		fmpz_add_si(k_num, k_num, delta);
		fmpz_mul(vx, x1, k_num);
		fmpz_mul(vy, y1, k_num);
		fmpz_add(vx, vx, x0);
		fmpz_add(vy, vy, y0);
		//Ensure coprime vectors
		fmpz_gcd(gcdtmp, vx, vy); 
		if(VerifyLatticeSolution(vx, vy, bitcoin->T, bitcoin->V, targetPrime) && fmpz_cmp_ui(gcdtmp, 1) == 0)
		{
			int vbits = CV_DT_SmoothnessTest(vx, vy, bitcoin->T, bitcoin->V, targetPrime, temp0, factors, targetSmallerBits);
            		if(vbits < smallest && vbits > 2)
            		{
           			//printf("\n%d: Current Best:%d ,Target:%d, ",k, smallest, targetSmallerBits);fmpz_print(c);printf(", ");fmpz_print(d);printf("\n");
            			_fmpz_factor_fit_length(bestFactors, factors->num);
            			bestFactors->num = factors->num;
            			for(slong i = 0; i < factors->num; i++)
            			{
            				fmpz_set(&bestFactors->p[i], &factors->p[i]);
            				bestFactors->exp[i] = factors->exp[i];
            			}
            			smallest = vbits;
            			fmpz_set(c, vx);
            			fmpz_set(d, vy);
            			if(smallest <= targetSmallerBits){break;}
            			//printf("Current Best:%d ",smallest);fmpz_print(c);printf(", ");fmpz_print(d);printf("\n");
            		}
		}
		
		fmpz_set(k_num, k_ceil);
		fmpz_add_si(k_num, k_num, delta);
		fmpz_mul(vx, x1, k_num);
		fmpz_mul(vy, y1, k_num);
		fmpz_add(vx, vx, x0);
		fmpz_add(vy, vy, y0);
		//Ensure coprime vectors
		fmpz_gcd(gcdtmp, vx, vy); 
		if(VerifyLatticeSolution(vx, vy, bitcoin->T, bitcoin->V, targetPrime) && fmpz_cmp_ui(gcdtmp, 1) == 0)
		{
			int vbits = CV_DT_SmoothnessTest(vx, vy, bitcoin->T, bitcoin->V, targetPrime, temp0, factors, targetSmallerBits);
            		if(vbits < smallest && vbits > 2)
            		{
           			//printf("\n%d: Current Best:%d ,Target:%d, ",k, smallest, targetSmallerBits);fmpz_print(c);printf(", ");fmpz_print(d);printf("\n");
            			_fmpz_factor_fit_length(bestFactors, factors->num);
            			bestFactors->num = factors->num;
            			for(slong i = 0; i < factors->num; i++)
            			{
            				fmpz_set(&bestFactors->p[i], &factors->p[i]);
            				bestFactors->exp[i] = factors->exp[i];
            			}
            			smallest = vbits;
            			fmpz_set(c, vx);
            			fmpz_set(d, vy);
            			if(smallest <= targetSmallerBits){break;}
            			//printf("Current Best:%d ",smallest);fmpz_print(c);printf(", ");fmpz_print(d);printf("\n");
            		}
		}
	}
	fmpz_clear(dot);fmpz_clear(norm1_sq);fmpz_clear(k_num);fmpz_clear(k_floor);fmpz_clear(k_ceil);fmpz_clear(vx);fmpz_clear(vy);fmpz_clear(gcdtmp);
	return smallest;
}

int SmallSearch(Bitcoin bitcoin, fmpz_factor_t factors, fmpz_factor_t bestFactors, int testRange, int targetSmallerBits, fmpz_t targetPrime, fmpz_t temp0, fmpz_t c,fmpz_t  d, fmpz_t x0, fmpz_t y0, fmpz_t  x1,fmpz_t  y1)
{
	fmpz_t k, A, B,vx,vy,gcdtmp, tmp;fmpz_t klist[4];;
	fmpz_init(k);fmpz_init(A);fmpz_init(B);fmpz_init(tmp);fmpz_init(vx);fmpz_init(vy);fmpz_init(gcdtmp);
	for(int i = 0; i < 4; i++){fmpz_init(klist[i]);}
	int smallest = 800000;
	//Euclidean CVP 
	fmpz_mul(tmp, x0, x1);
	fmpz_addmul(tmp, y0, y1);          
	fmpz_neg(tmp, tmp);
	fmpz_mul(A, x1, x1);
	fmpz_addmul(A, y1, y1);            
	fmpz_fdiv_q(klist[0], tmp, A);
	// k = -x0/x1
	if(!fmpz_is_zero(x1)){fmpz_neg(tmp, x0);fmpz_fdiv_q(klist[1], tmp, x1);}else{fmpz_zero(klist[1]);}
	// k = -y0/y1
	if(!fmpz_is_zero(y1)){fmpz_neg(tmp, y0);fmpz_fdiv_q(klist[2], tmp, y1);}else fmpz_zero(klist[2]);
	//k = -(x0*V - y0*T)/(x1*V - y1*T)
	fmpz_mul(A, x0, bitcoin->V);
	fmpz_submul(A, y0, bitcoin->T);    
	fmpz_mul(B, x1, bitcoin->V);
	fmpz_submul(B, y1, bitcoin->T);    
	if(!fmpz_is_zero(B)){fmpz_neg(A, A);fmpz_fdiv_q(klist[3], A, B);}else fmpz_zero(klist[3]);
	
	//Test each k and small neighborhood
	for(int i = 0; i < 4; i++)
	{
		for(int delta = -testRange; delta <= testRange; delta++)
		{
			fmpz_set(k, klist[i]);
			fmpz_add_si(k, k, delta);
			fmpz_mul(vx, x1, k);
			fmpz_add(vx, vx, x0);
			fmpz_mul(vy, y1, k);
			fmpz_add(vy, vy, y0);

			fmpz_gcd(gcdtmp, vx, vy);
			if(fmpz_cmp_ui(gcdtmp, 1) != 0)continue;
			if(!VerifyLatticeSolution(vx, vy, bitcoin->T, bitcoin->V, targetPrime))continue;
			int vbits = CV_DT_SmoothnessTest(vx, vy,bitcoin->T, bitcoin->V,targetPrime, temp0,factors, targetSmallerBits);

			if (vbits < smallest && vbits > 2)
			{
				_fmpz_factor_fit_length(bestFactors, factors->num);
				bestFactors->num = factors->num;
				for(slong j = 0; j < factors->num; j++)
				{
					fmpz_set(&bestFactors->p[j], &factors->p[j]);
					bestFactors->exp[j] = factors->exp[j];
				}

				smallest = vbits;
				fmpz_set(c, vx);
				fmpz_set(d, vy);

				if(smallest <= targetSmallerBits)goto done_small;
			}
		}
	}
	done_small:
		for (int i = 0; i < 4; i++) fmpz_clear(klist[i]);
	fmpz_clear(k);fmpz_clear(A);fmpz_clear(B);fmpz_clear(tmp);fmpz_clear(vx);fmpz_clear(vy);fmpz_clear(gcdtmp);
	return smallest;
}

int MakePrimeSmaller(Bitcoin bitcoin, fmpz_t targetPrime, fmpz_t c, fmpz_t d,fmpz_factor_t bestFactors, int *largestInteger, int *largestComplex)
{
	int targetSmallerBits = 27;	
	//Find particular solution Vu − Tv = 1 and scale to targetPrime modulus
	fmpz_factor_t factors;fmpz_factor_init(factors);
	fmpz_t x0, y0,x1,y1,vx, vy, gcdtmp,cV_min_dT,temp0, unityRoot,k_num;
	fmpz_init(temp0);fmpz_init(cV_min_dT);fmpz_init(x0);fmpz_init(y0);fmpz_init(k_num);fmpz_init(x1);fmpz_init(y1); fmpz_init(gcdtmp);fmpz_init(vx);fmpz_init(vy);fmpz_init(unityRoot);
	fmpz_xgcd(gcdtmp, x0, y0, bitcoin->V, bitcoin->T);
	fmpz_neg(y0, y0);
	fmpz_mul(x0, x0, targetPrime);
	fmpz_mul(y0, y0, targetPrime);
	assert(VerifyLatticeSolution(x0, y0, bitcoin->T, bitcoin->V, targetPrime) == true);
	//Set homogeneous solution
	fmpz_set(x1, bitcoin->T);fmpz_set(y1, bitcoin->V);
	assert(VerifyLatticeSolution(x1, y1, bitcoin->T, bitcoin->V, targetPrime) == true);
	
	//Find shortest solution	
	GaussReduce(x0,y0,x1,y1);
	assert(VerifyLatticeSolution(x0, y0, bitcoin->T, bitcoin->V, targetPrime) == true);
	int testRange = 500;
	int smallest = 8000000;
	smallest = DotProductSearch(bitcoin, factors, bestFactors, testRange,targetSmallerBits, targetPrime, temp0, c,d, x0, y0, x1, y1);
	

	
	fmpz_factor_clear(factors);
	fmpz_clear(temp0);fmpz_clear(cV_min_dT);fmpz_clear(x0);fmpz_clear(k_num);fmpz_clear(y0);fmpz_clear(x1);fmpz_clear(y1);fmpz_clear(gcdtmp);fmpz_clear(vx);fmpz_clear(vy);fmpz_clear(unityRoot);
	*largestInteger = smallest;
	return smallest;
}

bool CompareFactorization_Norm(fmpz_t target, fmpz_factor_t denomFactors, fmpz_t c, fmpz_t d, fmpz_t V, fmpz_t S, fmpz_t p, fmpzE *primeFactors, int *largestFactor, int *findVSign)
{
	fmpz_t lhs, rhs, factor, power,mod;
	fmpz_init(lhs);fmpz_init(rhs);fmpz_init(factor);fmpz_init(power);fmpz_init(mod);
	

	// Compute LHS = cV - dT mod p
	fmpz_mul(lhs, d, S);
	fmpz_add(lhs, lhs, c);            // lhs = c + dS
	fmpz_mul(lhs, lhs, V);            // lhs = V(c + dS)
	fmpz_mod(lhs, lhs, p);

	// Compute RHS = product of (a_i + b_i S)^e_i mod p
	fmpz_set_ui(rhs, 1);
	*largestFactor = 0;
	for(size_t i = 0; i < arrlen(primeFactors); i++)
	{
		// factor = a_i + b_i * S mod p
		int sizeA = fmpz_sizeinbase(primeFactors[i]->a,2);
		int sizeB = fmpz_sizeinbase(primeFactors[i]->b,2);
		if(sizeA > *largestFactor){*largestFactor = sizeA;}
		if(sizeB > *largestFactor){*largestFactor = sizeB;}
		fmpz_mul(factor, primeFactors[i]->b, S);
		fmpz_add(factor, factor, primeFactors[i]->a);
		fmpz_mod(factor, factor, p);

		// power = factor^{e_i} mod p
		fmpz_powm_ui(power, factor, primeFactors[i]->exponent, p);

		// rhs = rhs * power mod p
		fmpz_mul(rhs, rhs, power);
		fmpz_mod(rhs, rhs, p);
	}

	fmpz_mul(rhs, rhs, V);
	fmpz_mod(rhs, rhs, p);
	
	bool result = fmpz_equal(lhs, rhs);
	if(result == true)
	{
		//Find sign
		fmpz_mod(mod, lhs, target);
		//If not 0 then we need to sub
		if(fmpz_cmp_ui(mod, 0) != 0)
		{
			fmpz_sub(mod,p,lhs);
			fmpz_mod(mod, mod, target);
			//If zero now then make V negative
			if(fmpz_cmp_ui(mod, 0) == 0)
			{
				*findVSign = -1;
			}
		}
		else
		{
			*findVSign = 1;
		}
		//printf("\nV(c+dS) mod p = "); fmpz_print(lhs); printf("\n");
		//printf("Product of factors mod p = "); fmpz_print(rhs); printf("\n");
		//printf("Match? %s\n", fmpz_equal(lhs, rhs) ? "YES" : "NO");
		//printf("Target mod p = "); fmpz_print(target); printf("\n");
		//printf("lhs mod target = "); fmpz_print(mod); printf("\n");
		
	}
	//fmpz_set_str(power, "14791611737913357511537", 10);
	//if(fmpz_cmp(target, power) == 0)
	{
		//fmpzE_printPrimeFactors(primeFactors);printf("\n");
		//exit(1);
	}

	fmpz_clear(lhs);fmpz_clear(rhs);fmpz_clear(factor);fmpz_clear(power);fmpz_clear(mod);
	return result;
}

void HandleFactors(bool numOrDen, Bitcoin bitcoin, Relation relation, fmpz_factor_t factors, int factorbaseSize)
{
	int largePrimeCount = 0;
	int largestIndex = 0;
	int largestFactorSize = FindLargestFactorSize(factors, &largePrimeCount, &largestIndex, factorbaseSize);
	int bestInt = 0;
	int bestComplex = 0;
	int findVSign = 0;
	fmpz_factor_t bestFactors;
	fmpz_t c,d;
	fmpzE cd = fmpzE_init();
	fmpz_init(c);fmpz_init(d);
	//printf("\nlargePrimeCount : %d\n",largePrimeCount);
	//PrintFactors(factors);
	for(slong j = factors->num-1; j >= 0; j--)
	{
		int size = fmpz_sizeinbase(&factors->p[j], 2);
		if(size > factorbaseSize)
		{
			fmpz_factor_init(bestFactors);
			int smallFound = MakePrimeSmaller(bitcoin, &factors->p[j], cd->a, cd->b, bestFactors, &bestInt, &bestComplex);

			fmpzE *cv_dt_Factors = PrimeFactorization(cd);
			findVSign = 0;
			bool factorizationTest = CompareFactorization_Norm(&factors->p[j], bestFactors, cd->a, cd->b, bitcoin->V, bitcoin->rootOfUnity1, bitcoin->prime,cv_dt_Factors, &bestComplex, &findVSign);		
			if(numOrDen == true)
			{
				for(slong k = bestFactors->num-1; k >= 0; k--)
				{
					_fmpz_factor_append(relation->aDen, &bestFactors->p[k], bestFactors->exp[k]);				
				}
				for(size_t k = 0; k < arrlen(cv_dt_Factors); k++)
				{
					fmpzE eisFactor = fmpzE_init();
					fmpz_set(eisFactor->a, cv_dt_Factors[k]->a);
					fmpz_set(eisFactor->b, cv_dt_Factors[k]->b);
					eisFactor->exponent = cv_dt_Factors[k]->exponent;
					arrput(relation->aEis, eisFactor);
				}
				if(findVSign == 1)
				{
					relation->posVCountA += 1;
				}
				else if(findVSign == -1)
				{
					relation->negVCountA += 1;
				}
			}
			else
			{
				for(slong k = bestFactors->num-1; k >= 0; k--)
				{
					_fmpz_factor_append(relation->bDen, &bestFactors->p[k], bestFactors->exp[k]);				
				}
				for(size_t k = 0; k < arrlen(cv_dt_Factors); k++)
				{
					fmpzE eisFactor = fmpzE_init();
					fmpz_set(eisFactor->a, cv_dt_Factors[k]->a);
					fmpz_set(eisFactor->b, cv_dt_Factors[k]->b);
					eisFactor->exponent = cv_dt_Factors[k]->exponent;
					arrput(relation->bEis, eisFactor);
				}
				if(findVSign == 1)
				{
					relation->posVCountB += 1;
				}
				else if(findVSign == -1)
				{
					relation->negVCountB += 1;
				}
			}
			//fmpzE_printPrimeFactors(cv_dt_Factors);printf("\n");
			//PrintFactors(bestFactors);
			//printf("\n|%d| with %d %d| VA:(%d %d), VB:(%d %d) \n",size, bestInt, bestComplex, relation->posVCountA, relation->negVCountA,relation->posVCountB, relation->negVCountB);
			relation->bestInt = bestInt;
			relation->bestComplex = bestComplex;
			fmpzE_freePrimeFactors(cv_dt_Factors);
			fmpz_factor_clear(bestFactors);
		}
		else
		{
			if(numOrDen == true)
			{
				_fmpz_factor_append(relation->aNum, &factors->p[j], factors->exp[j]);
			}
			else
			{
				_fmpz_factor_append(relation->bNum, &factors->p[j], factors->exp[j]);			
			}
		}
	}
	fmpz_clear(c);fmpz_clear(d);fmpzE_clear(cd);
}

void CustomMurage(Bitcoin bitcoin)
{
	fmpz_t temp0;fmpz_init(temp0);
	fmpz_factor_t splitFactors;fmpz_factor_init(splitFactors);
	//Add split factors
	fmpz_set_ui(temp0,3);_fmpz_factor_append(splitFactors, temp0, 1);
	fmpz_set_ui(temp0,53);_fmpz_factor_append(splitFactors, temp0, 1);
	fmpz_set_ui(temp0,127);_fmpz_factor_append(splitFactors, temp0, 1);
	fmpz_set_ui(temp0,622793);_fmpz_factor_append(splitFactors, temp0, 1);
	fmpz_set_ui(temp0,96590783);_fmpz_factor_append(splitFactors, temp0, 1);
	fmpz_set_ui(temp0,10405113293);_fmpz_factor_append(splitFactors, temp0, 1);
	fmpz_set_str(temp0,"126147778933237103171",10);_fmpz_factor_append(splitFactors, temp0, 1);
	fmpz_set_str(temp0,"39831306940520651376621602713",10);_fmpz_factor_append(splitFactors, temp0, 1);
	fmpz_set_str(temp0,"11600063",10);_fmpz_factor_append(splitFactors, temp0, -1);
	
	PrintFactors(splitFactors);
	fmpz_factor_clear(splitFactors);
	fmpz_clear(temp0);	
}

void GenerateStartKey(int bitCount, fmpz_t secretKey)
{
	srand(time(NULL) ^ getpid());
	fmpz_set_ui(secretKey, 0);
	//Always set start bit
	fmpz_setbit(secretKey, bitCount - 1);
	for(int i = bitCount - 2; i >= 0; i--)
	{
		int decide = rand() % 2;
		if(decide)
		{
			fmpz_setbit(secretKey, i);	
		}
	}
}

void ProcessRelation(Bitcoin bitcoin, Relation relation, fmpz_t temp0,fmpz_factor_t factorNum,fmpz_factor_t factorDen,fmpq_t rationalReconstruction)
{
	int proved = -1;
	int targetSmallerBits = 60;
	int factorbaseSize = 27;
	int validCount = 0; 
	
	relation->reconstructionSuccess = fmpq_reconstruct_fmpz(rationalReconstruction, relation->rhs1, bitcoin->prime);
	relation->validRelation = false;
	if(relation->reconstructionSuccess)
	{
		fmpz_set(relation->a, fmpq_numref(rationalReconstruction));
		fmpz_set(relation->b, fmpq_denref(rationalReconstruction));

		relation->aSign = 1;
		if(fmpz_cmp_ui(relation->a, 0) < 0)
		{
			//Find absolute value of numerator and sign
			relation->aSign = -1;
			fmpz_mul_si(relation->a,relation->a,-1);
		}
		
		//Factorize num and den
		int factorNumResult = fmpz_factor_smooth(factorNum, relation->a, targetSmallerBits, proved);
		int factorDenResult = fmpz_factor_smooth(factorDen, relation->b, targetSmallerBits, proved);
		HandleFactors(true, bitcoin, relation, factorNum, factorbaseSize);
		HandleFactors(false, bitcoin, relation, factorDen, factorbaseSize);
		//printf("\n(%ld / %ld)\n", i,  arrlen(relations));
		//PrintRelation(relation);
		relation->validRelation = TestRelation(bitcoin, relation);
		
	}
}

void FullPrintRelation(FILE *fp, Relation relation)
{
	fmpz_fprint(fp,relation->k);fprintf(fp,"\n");
	fmpz_fprint(fp,relation->rhs1);fprintf(fp,"\n");
	fprintf(fp,"%d\n",relation->aSign);
	fprintf(fp,"%d\n",relation->posVCountA);
	fprintf(fp,"%d\n",relation->negVCountA);
	fprintf(fp,"%d\n",relation->posVCountB);
	fprintf(fp,"%d\n",relation->negVCountB);
	fprintf(fp,"%d\n",relation->bestInt);
	fprintf(fp,"%d\n",relation->bestComplex);
	fprintf(fp,"%ld\n",relation->aNum->num);
	for(slong i = 0; i < relation->aNum->num; i++)
	{
		fmpz_fprint(fp,&relation->aNum->p[i]);fprintf(fp,"\n");
		fprintf(fp,"%ld\n",relation->aNum->exp[i]);
	}
	fprintf(fp,"%ld\n",relation->aDen->num);
	for(slong i = 0; i < relation->aDen->num; i++)
	{
		fmpz_fprint(fp,&relation->aDen->p[i]);fprintf(fp,"\n");
		fprintf(fp,"%ld\n",relation->aDen->exp[i]);
	}
	fprintf(fp,"%ld\n",relation->bNum->num);
	for(slong i = 0; i < relation->bNum->num; i++)
	{
		fmpz_fprint(fp,&relation->bNum->p[i]);fprintf(fp,"\n");
		fprintf(fp,"%ld\n",relation->bNum->exp[i]);
	}
	fprintf(fp,"%ld\n",relation->bDen->num);
	for(slong i = 0; i < relation->bDen->num; i++)
	{
		fmpz_fprint(fp,&relation->bDen->p[i]);fprintf(fp,"\n");
		fprintf(fp,"%ld\n",relation->bDen->exp[i]);
	}
	fprintf(fp,"%ld\n",arrlen(relation->aEis));
	for(size_t i = 0; i < arrlen(relation->aEis); i++)
	{
		fmpz_fprint(fp,relation->aEis[i]->a);fprintf(fp,",");fmpz_fprint(fp,relation->aEis[i]->b);fprintf(fp,"\n");
		fprintf(fp,"%d\n",relation->aEis[i]->exponent);
	}
	fprintf(fp,"%ld\n",arrlen(relation->bEis));
	for(size_t i = 0; i < arrlen(relation->bEis); i++)
	{
		fmpz_fprint(fp,relation->bEis[i]->a);fprintf(fp,",");fmpz_fprint(fp,relation->bEis[i]->b);fprintf(fp,"\n");
		fprintf(fp,"%d\n",relation->bEis[i]->exponent);
	}
	fprintf(fp,"\n");		
}

void CollectRelations()
{
	int bitCount = rand() % 256;
	int factorBaseSize = 27;
	if(bitCount == 0){bitCount = 1;}
	
	fmpz_t base, exponent,result;fmpz_init(base);fmpz_init(exponent);fmpz_init(result);
	fmpz_set_str(base, "70926458516305583538736597864263618897131839054605814437836363826034175889389", 10);
	fmpz_t temp0;fmpz_init(temp0);
	fmpz_factor_t factorNum;fmpz_factor_init(factorNum);
	fmpz_factor_t factorDen;fmpz_factor_init(factorDen);
	fmpq_t rationalReconstruction;fmpq_init(rationalReconstruction);
	GenerateStartKey(bitCount, exponent);
	
	Bitcoin bitcoin = CreateBitcoin();
	fmpz_powm(result, base, exponent, bitcoin->prime);
	
	fmp_printf("exponent", exponent);
	fmp_printf("res", result);
	
	int validRelationCount = 0;
	int totalRelationCount = 0;
	char *outputDirectory = "Relations";
	time_t now = time(NULL);struct tm *tm_now = localtime(&now);
	char timebuf[64];strftime(timebuf, sizeof(timebuf), "%Y%m%d_%H%M%S", tm_now);
	char kOffsetStr[256];fmpz_get_str(kOffsetStr, 10, exponent);
	char filename[512];
	snprintf(filename, sizeof(filename),"%s/%s_%s.txt",outputDirectory,timebuf,kOffsetStr);
	FILE *out = fopen(filename, "a");
	assert(out != NULL);
	
	for(int i = 0; i < 1; i++)
	{
		Relation relation = CreateRelation();
		fmpz_mul(result, result, base);
		fmpz_mod(result, result, bitcoin->prime);
		
		fmpz_add_ui(relation->k, exponent, i);
		fmpz_set(relation->rhs1, result);
		ProcessRelation(bitcoin, relation, temp0, factorNum, factorDen, rationalReconstruction);
		if(relation->reconstructionSuccess && relation->validRelation == true)
		{
			printf("\nFound %d %d| VA:(%d %d), VB:(%d %d) \n", relation->bestInt, relation->bestComplex, relation->posVCountA, relation->negVCountA,relation->posVCountB, relation->negVCountB);	
			if(relation->bestInt <= factorBaseSize && relation->bestComplex <= factorBaseSize)
			{
				validRelationCount += 1;
				printf("Status: Accepted(%d / %d)\n", validRelationCount, totalRelationCount);
				FullPrintRelation(out, relation);
			}
			else
			{
				printf("Status: Rejected\n");
			}
		}
		totalRelationCount += 1;
		DestroyRelation(relation);
	}
	
	fmpz_clear(base);fmpz_clear(exponent);fmpz_clear(result);
	fmpz_factor_clear(factorNum);
	fmpz_factor_clear(factorDen);
	fmpz_clear(temp0);
	fmpq_clear(rationalReconstruction);
	DestroyBitcoin(bitcoin);
	fclose(out);
}

int main()
{
	Bitcoin bitcoin = CreateBitcoin();
	CollectRelations();
	DestroyBitcoin(bitcoin);
	flint_cleanup();
	return 0;
}
