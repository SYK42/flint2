/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);

    flint_printf("reconstruct_fmpz....");
    fflush(stdout);


    int fail_counter_mqrr = 0;
    for (i = 0; i < 1000*flint_test_multiplier(); i++)
    {
        int result;
        fmpq_t x, y;
        fmpz_t mod;
        fmpz_t res;
        fmpz_t T;

        fmpq_init(x);
        fmpq_init(y);
        fmpz_init(mod);
        fmpz_init(res);
        fmpz_init(T);

        fmpq_randtest(x, state, 1000);

        /* Modulus m > 2*max(|n|,d)^2 */
        if (fmpz_cmpabs(&x->num, &x->den) >= 0)
            fmpz_mul(mod, &x->num, &x->num);
        else
            fmpz_mul(mod, &x->den, &x->den);
        fmpz_mul_2exp(mod, mod, 1);
        do fmpz_add_ui(mod, mod, 1);
        while (!fmpq_mod_fmpz(res, x, mod));

        result = fmpq_reconstruct_fmpz(y, res, mod);

        if (!result || !fmpq_equal(x, y))
        {
            flint_printf("FAIL: reconstruction failed\n");
            flint_printf("input = ");
            fmpq_print(x);
            flint_printf("\nmodulus = ");
            fmpz_print(mod);
            flint_printf("\nresidue = ");
            fmpz_print(res);
            flint_printf("\nreconstructed = ");
            fmpq_print(y);
            flint_printf("\nfmpq_reconstruct_fmpz return value = %d", result);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_set_ui(T, fmpz_sizeinbase(mod, 2));
        fmpz_mul_ui(T, T, 1024);
        result = fmpq_reconstruct_fmpz_mqrr(y, res, mod, T);
        if (result && !fmpq_equal(x, y))
        {
            fail_counter_mqrr ++;
        }

        fmpq_clear(x);
        fmpq_clear(y);
        fmpz_clear(mod);
        fmpz_clear(res);
        fmpz_clear(T);
    }

    // mqrr should only fail in less than 1% of cases
    if (fail_counter_mqrr * 100 >  1000*flint_test_multiplier())
    {
        flint_printf("FAIL: maximal quotient reconstruction failed %i out of %i times.", fail_counter_mqrr, 1000*flint_test_multiplier());
        fflush(stdout);
        flint_abort();
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}

