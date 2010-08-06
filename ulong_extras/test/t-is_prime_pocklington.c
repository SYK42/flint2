/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2009 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
   int i, result;
   ulong count = 0;
   mp_limb_t d;
   mpz_t d_m;
   
   printf("is_prime_pocklington....");
   fflush(stdout);
   
   for (i = 0; i < 10000; i++) /* Test that primes pass the test */
   {
      mpz_init(d_m);

      do
      {
         d = n_randtest() | 1;
         mpz_set_ui(d_m, d);
         mpz_nextprime(d_m, d_m);
         d = mpz_get_ui(d_m);
      } while (mpz_size(d_m) > 1);

      result = n_is_prime_pocklington(d, 100);
      if (result == -1) count++;

      if (result == 0)
      {
         printf("FAIL:\n");
         printf("d = %lu is declared composite\n", d); 
         abort();
      }

      mpz_clear(d_m);
   }
         
   for (i = 0; i < 10000; i++) /* Test that composites do not pass */
   {
      mpz_init(d_m);

      do
      {
         d = n_randtest() | 1;
         if (d == 1UL) d++;
         mpz_set_ui(d_m, d);
      } while (mpz_probab_prime_p(d_m, 12));

      result = n_is_prime_pocklington(d, 100);
      if (result == -1) count++;

      if (result == 1)
      {
         printf("FAIL:\n");
         printf("d = %lu is declared prime\n", d); 
         abort();
      }

      mpz_clear(d_m);
   }

   if (count > 2000)
   {
      printf("FAIL: Pocklington-Lehmer failed too many times (%lu times)\n", count);
      abort();
   }

   printf("PASS\n");
   return 0;
}
