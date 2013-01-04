/*
// This code is a modification of code written by Jan Wigginto to implement exact SNP test of
// Hardy-Weinberg Equilibrium. The modification consist of adjusting the function parameter list to
// allow the function to be called from R, and the addition of a loop to calculate multiple tests in one call.
// Modified by Ian Painter. The following is the original header comment:
//
//
// This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
// Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of 
// Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000  
//
// Written by Jan Wigginton
//
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void SNPHWE(int *length, int *obs_hets_in, int *obs_hom1_in, int *obs_hom2_in, double *hw_out)
{
  
   int index;
   
   double p_hwe = 0.0;
   int lngth = *length;
   
   

   
   for (index = 0; index < lngth; index++)
   {	   
	   
	      
	   
	   int obs_hets = obs_hets_in[index];
	   int obs_hom1 = obs_hom1_in[index];
	   int obs_hom2 = obs_hom2_in[index];
	   

	   int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
	   int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;
	
	  
	   
	   int rare_copies = 2 * obs_homr + obs_hets;
	   int genotypes   = obs_hets + obs_homc + obs_homr;
	
	   double * het_probs = (double *) malloc((size_t) (rare_copies + 1) * sizeof(double));
	   if (het_probs == NULL) 
	      {
	/*
	      printf("FATAL ERROR - SNP-HWE: Unable to allocate array for heterozygote probabilities" );
	      exit(EXIT_FAILURE); */

	      }
	   else
	   {
		   int i;
		   for (i = 0; i <= rare_copies; i++)
		      het_probs[i] = 0.0;
	
		   /* start at midpoint */
		   int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);
	
		   /* check to ensure that midpoint and rare alleles have same parity */
		   if ((rare_copies & 1) ^ (mid & 1))
		      mid++;
	
		   int curr_hets = mid;
		   int curr_homr = (rare_copies - mid) / 2;
		   int curr_homc = genotypes - curr_hets - curr_homr;
	
		   het_probs[mid] = 1.0;
		   double sum = het_probs[mid];
		   for (curr_hets = mid; curr_hets > 1; curr_hets -= 2)
		      {
		      het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0)
					       / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
		      sum += het_probs[curr_hets - 2];
	
		      /* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
		      curr_homr++;
		      curr_homc++;
		      }
	
		   curr_hets = mid;
		   curr_homr = (rare_copies - mid) / 2;
		   curr_homc = genotypes - curr_hets - curr_homr;
		   for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2)
		      {
		      het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc
					    /((curr_hets + 2.0) * (curr_hets + 1.0));
		      sum += het_probs[curr_hets + 2];
	
		      /* add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote */
		      curr_homr--;
		      curr_homc--;
		      }
	
		   for (i = 0; i <= rare_copies; i++)
		      het_probs[i] /= sum;
	
		   /* alternate p-value calculation for p_hi/p_lo
		   double p_hi = het_probs[obs_hets];
		   for (i = obs_hets + 1; i <= rare_copies; i++)
		     p_hi += het_probs[i];
		   
		   double p_lo = het_probs[obs_hets];
		   for (i = obs_hets - 1; i >= 0; i--)
		      p_lo += het_probs[i];
	
		   
		   double p_hi_lo = p_hi < p_lo ? 2.0 * p_hi : 2.0 * p_lo;
		   */
	
		   p_hwe = 0.0;
		   /*  p-value calculation for p_hwe  */
		   for (i = 0; i <= rare_copies; i++)
		   {
		      if (het_probs[i] > het_probs[obs_hets])
			 continue;
		      p_hwe += het_probs[i];
		   }
		   
		   p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;
		   
		   
		   hw_out[index] = p_hwe;
	   
		   free(het_probs);
		}
   }
   
   

   }
