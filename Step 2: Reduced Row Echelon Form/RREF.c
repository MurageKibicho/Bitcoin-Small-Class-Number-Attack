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
/*This is WIP*/
//Run: clear && gcc RREF.c -o m.o -lm -lgmp -lmpfr -lflint && ./m.o
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
bool TestRelationRREF(Bitcoin bitcoin, Relation relation)
{
	bool result = true;
	//Reconstruct a/b mod p
	fmpz_t temp0,temp1,temp2,recon,VA,VB,base;
	fmpz_init(temp0);fmpz_init(temp1);fmpz_init(temp2);fmpz_init(recon);fmpz_init(VA);fmpz_init(VB);fmpz_init(base);
	//Set V values
	fmpz_neg(VA, bitcoin->V);
	fmpz_powm_ui(VA, VA, relation->negVCountA, bitcoin->prime);
	fmpz_powm_ui(temp0, bitcoin->V, relation->posVCountA, bitcoin->prime);
	fmpz_mul(VA,VA,temp0);fmpz_mod(VA, VA, bitcoin->prime);
	
	fmpz_neg(VB, bitcoin->V);
	fmpz_powm_ui(VB, VB, relation->negVCountB, bitcoin->prime);
	fmpz_powm_ui(temp0, bitcoin->V, relation->posVCountB, bitcoin->prime);
	fmpz_mul(VB,VB,temp0);fmpz_mod(VB, VB, bitcoin->prime);

	fmpz_set_str(base, "70926458516305583538736597864263618897131839054605814437836363826034175889389", 10);
	fmpz_powm(temp2, base, relation->k, bitcoin->prime);
	if(fmpz_cmp(temp2, relation->rhs1) != 0)
	{
		printf("base ^ k != relation->rhs1\n");
		fmp_printf("relation->rhs1", relation->rhs1);
		fmp_printf("temp2", temp2);
		result = false;
	}
	//Find rational reconstruction
	fmpq_t rationalReconstruction;fmpq_init(rationalReconstruction);
	relation->reconstructionSuccess = fmpq_reconstruct_fmpz(rationalReconstruction, relation->rhs1, bitcoin->prime);
	if(relation->reconstructionSuccess)
	{
		fmpz_set(relation->a, fmpq_numref(rationalReconstruction));
		fmpz_set(relation->b, fmpq_denref(rationalReconstruction));

		int aSign = 1;
		if(fmpz_cmp_ui(relation->a, 0) < 0)
		{
			//Find absolute value of numerator and sign
			aSign = -1;
			fmpz_mul_si(relation->a,relation->a,-1);
		}
		if(aSign != relation->aSign)
		{
			printf("Reconstruct sign mismatch\n");
			result = false;
		}
	}
	else
	{
		printf("fmpq_reconstruct_fmpz failed\n");
		result = false;
	}
	fmpq_clear(rationalReconstruction);
	
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
	
	//fmp_printf("a with V: ", temp0);
	//fmp_printf("b with V: ", temp1);
	//fmp_printf("a", relation->a);
	//fmp_printf("b", relation->b);
		
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
	
	
	fmpz_clear(temp0);fmpz_clear(temp1);fmpz_clear(temp2);fmpz_clear(recon);fmpz_clear(VA);fmpz_clear(VB);fmpz_clear(base);
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

int TestReadRelations(Bitcoin bitcoin)
{
	Relation *relations = NULL;
	char *directoryName = "Relations";
	int directoryFileCount = 0;
	printf("Opening folder `%s`\n", directoryName);
	fmpz_t temp0; fmpz_init(temp0);
	char **fileNames = GetFileNamesInDirectory(directoryName, &directoryFileCount);assert(fileNames != NULL);
	for(int j = 0; j < directoryFileCount; j++)
	{
		char fileName[MAX_FILENAME_LENGTH] = {0};
		snprintf(fileName, sizeof(fileName), "%s/%s",directoryName,fileNames[j]);
		printf("%s\n",fileName);
		{
			size_t fileSize = GetFileSize(fileName);
			if(fileSize > 0)
			{
				FILE* fp = fopen(fileName, "r");
				assert(fp != NULL);
				
				char line[2048];
				int currentStep = 0;
				slong factorCount = 0;
				slong exponent = 0;
				while(true)
				{
					//Read Relations in a loop
					Relation relation = CreateRelation(); 
					//Read k
					if(!fgets(line, sizeof(line), fp))
					{
						DestroyRelation(relation);  
						break;                  
					}
					line[strcspn(line, "\r\n")] = 0; // remove newline
					if(fmpz_set_str(relation->k, line, 10) != 0)  return 0;;
					//Read rhs1
					if(!fgets(line, sizeof(line), fp)) return 0;
					line[strcspn(line, "\r\n")] = 0;
					if(fmpz_set_str(relation->rhs1, line, 10) != 0)  return 0;;
					if(!fgets(line, sizeof(line), fp)) return 0; relation->aSign = atoi(line);
					if(!fgets(line, sizeof(line), fp)) return 0; relation->posVCountA = atoi(line);
					if(!fgets(line, sizeof(line), fp)) return 0; relation->negVCountA = atoi(line);
					if(!fgets(line, sizeof(line), fp)) return 0; relation->posVCountB = atoi(line);
					if(!fgets(line, sizeof(line), fp)) return 0; relation->negVCountB = atoi(line);
					if(!fgets(line, sizeof(line), fp)) return 0; relation->bestInt = atoi(line);
					if(!fgets(line, sizeof(line), fp)) return 0; relation->bestComplex = atoi(line);
					//Read aNum
					if(!fgets(line, sizeof(line), fp)) return 0; factorCount = atol(line);
					for(slong i = 0; i < factorCount; i++)
					{
						if(!fgets(line, sizeof(line), fp)) return 0;line[strcspn(line, "\r\n")] = 0;
						if(fmpz_set_str(temp0, line, 10) != 0) return 0;

						if(!fgets(line, sizeof(line), fp)) return 0;
						exponent = atol(line);
						_fmpz_factor_append(relation->aNum, temp0, exponent);				
					} 

					//Read aDen
					if(!fgets(line, sizeof(line), fp)) return 0; factorCount = atol(line);
					for(slong i = 0; i < factorCount; i++)
					{
						if(!fgets(line, sizeof(line), fp)) return 0;line[strcspn(line, "\r\n")] = 0;
						if(fmpz_set_str(temp0, line, 10) != 0) return 0;

						if(!fgets(line, sizeof(line), fp)) return 0;
						exponent = atol(line);
						_fmpz_factor_append(relation->aDen, temp0, exponent);				
					}    

					//Read bNum
					if(!fgets(line, sizeof(line), fp)) return 0; factorCount = atol(line);
					for(slong i = 0; i < factorCount; i++)
					{
						if(!fgets(line, sizeof(line), fp)) return 0;line[strcspn(line, "\r\n")] = 0;
						if(fmpz_set_str(temp0, line, 10) != 0) return 0;

						if(!fgets(line, sizeof(line), fp)) return 0;
						exponent = atol(line);
						_fmpz_factor_append(relation->bNum, temp0, exponent);				
					} 

					//Read bDen
					if(!fgets(line, sizeof(line), fp)) return 0; factorCount = atol(line);
					for(slong i = 0; i < factorCount; i++)
					{
						if(!fgets(line, sizeof(line), fp)) return 0;line[strcspn(line, "\r\n")] = 0;
						if(fmpz_set_str(temp0, line, 10) != 0) return 0;

						if(!fgets(line, sizeof(line), fp)) return 0;
						exponent = atol(line);
						_fmpz_factor_append(relation->bDen, temp0, exponent);				
					}    
					//Read aEis
					if(!fgets(line, sizeof(line), fp)) return 0; size_t lenA = atol(line);
					relation->aEis = NULL;
					for(size_t i = 0; i < lenA; i++)
					{
						fmpzE eis = fmpzE_init();
						if(!fgets(line, sizeof(line), fp)) return 0;
						line[strcspn(line, "\r\n")] = 0;
						//printf("SDS%sm\n", line);
						//a,b split line
						char *comma = strchr(line, ','); 
						if(!comma)
						{
							fprintf(stderr, "aEis Malformed line: %s\n", line);
							continue;
						}
						*comma = '\0';
						char *aStr = line;
						char *bStr = comma + 1;
						//printf("F%s,%s ", aStr, bStr);
						if(fmpz_set_str(eis->a, aStr, 10) != 0) return 0;
						if(fmpz_set_str(eis->b, bStr, 10) != 0) return 0;
						if(!fgets(line, sizeof(line), fp)) return 0;
						eis->exponent = atoi(line);
						arrput(relation->aEis, eis);
					}
					//Read bEis
					if(!fgets(line, sizeof(line), fp)) return 0; size_t lenB = atol(line);
					relation->bEis = NULL;
					for(size_t i = 0; i < lenB; i++)
					{
						fmpzE eis = fmpzE_init();
						if(!fgets(line, sizeof(line), fp)) return 0;
						line[strcspn(line, "\r\n")] = 0;
						//a,b split line
						char *comma = strchr(line, ','); 
						if(!comma)
						{
							fprintf(stderr, "bEis Malformed line: %s\n", line);
							continue;
						}
						*comma = '\0';
						char *aStr = line;
						char *bStr = comma + 1;
						//printf("F%s,%s ", aStr, bStr);
						if(fmpz_set_str(eis->a, aStr, 10) != 0) return 0;
						if(fmpz_set_str(eis->b, bStr, 10) != 0) return 0;
						if(!fgets(line, sizeof(line), fp)) return 0;
						eis->exponent = atoi(line);
						arrput(relation->bEis, eis);
					}
					//Skip new line at end
					if(!fgets(line, sizeof(line), fp)) return 0;			
					//PrintRelation(relation); 
					relation->validRelation = TestRelationRREF(bitcoin, relation);
					if(relation->validRelation)
					{
						arrput(relations, relation);
					}
					
				}
				
				fclose(fp);	
			}
			
		}
	}
	printf("Found %ld valid relations\n", arrlen(relations));
	fmpz_clear(temp0);
	DestroyFileNames(directoryFileCount, fileNames);	
	for(size_t i = 0; i < arrlen(relations); i++)
	{	
		DestroyRelation(relations[i]);
	}
	arrfree(relations);
	return 0;
}

int main()
{
	Bitcoin bitcoin = CreateBitcoin();
	TestReadRelations(bitcoin);
	DestroyBitcoin(bitcoin);
	flint_cleanup();
	return 0;
}
