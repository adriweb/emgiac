/* Copyright (C) 2004  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. */

/******************************************************************/
/*                                                                */
/*              PARI header file (common to all versions)         */
/*                                                                */
/******************************************************************/
#ifdef STMT_START /* perl headers */
#  undef STMT_START
#endif
#ifdef STMT_END
#  undef STMT_END
#endif
/* STMT_START { statements; } STMT_END;
 * can be used as a single statement, as in
 * if (x) STMT_START { ... } STMT_END; else ...
 * [ avoid "dangling else" problem in macros ] */
#define STMT_START        do
#define STMT_END        while (0)
/*=====================================================================*/
/* pari_CATCH(numer) {
 *   recovery
 * } pari_TRY {
 *   code
 * } pari_ENDCATCH
 * will execute 'code', then 'recovery' if exception 'numer' is thrown
 * [ any exception if numer == CATCH_ALL ].
 * pari_RETRY = as pari_TRY, but execute 'recovery', then 'code' again [still catching] */

extern THREAD jmp_buf *iferr_env;
extern const long CATCH_ALL;

#define pari_CATCH2(var,err) {         \
  jmp_buf *var=iferr_env;    \
  jmp_buf __env;             \
  iferr_env = &__env;        \
  if (setjmp(*iferr_env))    \
  {                          \
    GEN __iferr_data = pari_err_last(); \
    iferr_env = var; \
    if (err!=CATCH_ALL && err_get_num(__iferr_data) != err) \
      pari_err(0, __iferr_data);

#define pari_CATCH2_reset(var) (iferr_env = var)
#define pari_ENDCATCH2(var) iferr_env = var; } }

#define pari_CATCH(err) pari_CATCH2(__iferr_old,err)
#define pari_RETRY } iferr_env = &__env; {
#define pari_TRY } else {
#define pari_CATCH_reset() pari_CATCH2_reset(__iferr_old)
#define pari_ENDCATCH pari_ENDCATCH2(__iferr_old)

extern const double LOG2, LOG10_2, LOG2_10;
#ifndef  PI
#  define PI (3.141592653589)
#endif

/* Common global variables: */
extern int new_galois_format, factor_add_primes, factor_proven;
extern ulong DEBUGFILES, DEBUGLEVEL, DEBUGMEM, precdl;
extern long DEBUGVAR;
extern ulong pari_mt_nbthreads;
extern THREAD GEN  bernzone;
extern GEN primetab;
extern GEN gen_m1,gen_1,gen_2,gen_m2,ghalf,gen_0,gnil,err_e_STACK;
extern THREAD VOLATILE int PARI_SIGINT_block, PARI_SIGINT_pending;

extern const long lontyp[];
extern void (*cb_pari_ask_confirm)(const char *);
extern int  (*cb_pari_whatnow)(PariOUT *out, const char *, int);
extern void (*cb_pari_sigint)(void);
extern int (*cb_pari_handle_exception)(long);
extern void (*cb_pari_pre_recover)(long);
extern void (*cb_pari_err_recover)(long);
extern const char *pari_library_path;

enum manage_var_t {
  manage_var_create,
  manage_var_delete,
  manage_var_init,
  manage_var_next,
  manage_var_max_avail,
  manage_var_pop
};

/* pari_init_opts */
enum {
  INIT_JMPm = 1,
  INIT_SIGm = 2,
  INIT_DFTm = 4,
  INIT_noPRIMEm = 8,
  INIT_noIMTm = 16
};

#ifndef HAS_EXP2
#  undef exp2
#  define exp2(x) (exp((double)(x)*LOG2))
#endif
#ifndef HAS_LOG2
#  undef log2
#  define log2(x) (log((double)(x))/LOG2)
#endif

#define ONLY_REM ((GEN*)0x1L)
#define ONLY_DIVIDES ((GEN*)0x2L)

#define NEXT_PRIME_VIADIFF(p,d) STMT_START { (p) += *(d)++; } STMT_END
#define PREC_PRIME_VIADIFF(p,d) STMT_START { (p) -= *--(d); } STMT_END
#define NEXT_PRIME_VIADIFF_CHECK(p,d)  STMT_START \
  { if (!*(d)) pari_err_MAXPRIME(0); NEXT_PRIME_VIADIFF(p,d); } STMT_END
