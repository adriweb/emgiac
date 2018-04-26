/* Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. */

/* output of get_nf and get_bnf */
enum {
  typ_NULL = 0,
  typ_POL,
  typ_Q,
  typ_NF,
  typ_BNF,
  typ_BNR,
  typ_ELL, /* elliptic curve */
  typ_QUA, /* quadclassunit  */
  typ_GAL, /* galoisinit     */
  typ_BID,
  typ_PRID,
  typ_MODPR,
  typ_RNF
};

/* idealtyp */
enum {
  id_PRINCIPAL = 0,
  id_PRIME,
  id_MAT
};

typedef struct {
  GEN x; /* defining polynomial (monic, integral) */
  GEN x0; /* original defining polynomial (integral) */
  GEN bas;  /* Z-basis of O_K (t_VEC of t_POL) */
  long r1; /* number of real places of K */
/* possibly NULL = irrelevant or not computed */
  GEN dK; /* disc(K) */
  GEN dKP; /* "primes" dividing disc(K) [if we have a composite in the list
              then the structure may not be correct] */
  GEN index; /* [O_K : Z[X]/(x)] */
  GEN unscale; /* x = C*x0(X / unscale), rational */
  GEN dx;   /* disc(x) */
  GEN basden; /* [nums(bas), dens(bas)] */
} nfbasic_t;

typedef struct {
  GEN T, dT; /* monic defining polynomial, disc(T) */
  GEN T0; /* ORIGINAL polynomial T0 */
  GEN unscale; /* T = C*T0(x / unscale), rational */
  GEN dK; /* field discriminant */
  GEN index; /* index of power basis in maximal order */
  GEN dTP, dTE; /* (possibly partial) factorization of dT, primes / exponents */
  GEN dKP, dKE; /* (possibly partial) factorization of dK, primes / exponents */
  GEN basis; /* Z-basis for maximal order */
} nfmaxord_t;

typedef struct {
  GEN x;
  GEN ro;   /* roots of x */
  long r1;
  GEN basden;
  long prec;
/* possibly -1 = irrelevant or not computed */
  long extraprec;
/* possibly NULL = irrelevant or not computed */
  GEN M;
  GEN G;
} nffp_t;

/* qfr3 / qfr5 */
struct qfr_data { GEN D, sqrtD, isqrtD; };

/* various flags for nf/bnf routines */
enum {
  nf_ORIG = 1,
  nf_GEN = 1,
  nf_ABSOLUTE = 2,
  nf_FORCE = 2,
  nf_ALL = 4,
  nf_GENMAT = 4,
  nf_INIT = 4,
  nf_RAW = 8,
  nf_RED = 8,
  nf_PARTIALFACT = 16,
  nf_ROUND2 = 64, /* obsolete */
  nf_ADDZK =  256,
  nf_GEN_IF_PRINCIPAL = 512
};

enum {
  rnf_REL = 1,
  rnf_COND = 2
};

/* LLL */
enum {
  LLL_KER  = 1, /* only kernel */
  LLL_IM   = 2, /* only image */
  LLL_ALL  = 4, /* kernel & image */
  LLL_GRAM       = 0x100,
  LLL_KEEP_FIRST = 0x200,
  LLL_INPLACE    = 0x400
};

/* HNF */
enum { hnf_MODID = 1, hnf_PART = 2, hnf_CENTER = 4 };

/* for minim */
enum {
  min_ALL       = 0,
  min_FIRST     = 1,
  min_PERF      = 2,
  min_VECSMALL  = 3,
  min_VECSMALL2 = 4
};
/* for fincke_pohst() */
typedef struct FP_chk_fun {
  GEN (*f)(void *,GEN);
  /* f_init allowed to permute the columns of u and r */
  GEN (*f_init)(struct FP_chk_fun*,GEN,GEN);
  GEN (*f_post)(struct FP_chk_fun*,GEN,GEN);
  void *data;
  long skipfirst;
} FP_chk_fun;

/* for ideallog / zlog */
typedef struct {
  GEN lists; /* lists[i] = */
  GEN ind;  /* ind[i] = start of vector */
  GEN P, e; /* finit part of conductor = prod P^e */
  GEN archp; /* archimedean part of conductor, in permutation form */
  long n;  /* total number of generators for all (O_K/P^e)^* and (O_K/f_oo) */
  GEN U; /* base change matrix from generators to bid.gen */
} zlog_S;

GEN fincke_pohst(GEN a,GEN BOUND,long stockmax,long PREC, FP_chk_fun *CHECK);
void remake_GM(GEN nf, nffp_t *F, long prec);
GEN nfbasic_to_nf(nfbasic_t *T, GEN ro, long prec);

void init_zlog_bid(zlog_S *S, GEN bid);
GEN  log_gen_arch(zlog_S *S, long index);
GEN  log_gen_pr(zlog_S *S, long index, GEN nf, long e);
GEN  zlog(GEN nf, GEN a, GEN sgn, zlog_S *S);

/* conversions basis / alg */

/* nf a genuine NF, x an nfelt (t_COL) or t_MAT whose columns represent nfelts.
 * Return the corresponding elements as t_POLs (implicitly mod nf.pol) */
#define coltoliftalg(nf,x) (gmul(gel((nf),7), (x)))
GEN    poltobasis(GEN nf,GEN x);
GEN    coltoalg(GEN nf,GEN x);

/* Other number fields routines */
GEN    archstar_full_rk(GEN x, GEN bas, GEN v, GEN gen);
GEN    check_and_build_cycgen(GEN bnf);
long   check_LIMC(long LIMC, long LIMCMAX);
GEN    checkbid_i(GEN bid);
GEN    checkbnf_i(GEN bnf);
GEN    checknf_i(GEN nf);
GEN    pow_ei_mod_p(GEN nf, long I, GEN n, GEN p);
GEN    galoisbig(GEN x, long prec);
GEN    get_arch_real(GEN nf,GEN x,GEN *emb,long prec);
GEN    get_bas_den(GEN bas);
void   nf_set_multable(GEN nf, GEN bas, GEN basden);
GEN    get_nfindex(GEN bas);
GEN    get_proj_modT(GEN basis, GEN T, GEN p);
GEN    get_roots(GEN x,long r1,long prec);
GEN    get_theta_abstorel(GEN T, GEN pol, GEN k);
GEN    idealsqrtn(GEN nf, GEN x, GEN gn, int strict);
GEN    init_unif_mod_fZ(GEN L);
GEN    init_units(GEN BNF);
GEN    make_integral(GEN nf, GEN L0, GEN f, GEN listpr);
GEN    maxord_i(GEN p, GEN f, long mf, GEN w, long flag);
GEN    nf_deg1_prime(GEN nf);
GEN    nfpol_to_Flx(GEN nf, GEN pol, ulong *ptp);
GEN    nfroots_split(GEN nf, GEN pol);
GEN    pidealprimeinv(GEN nf, GEN x);
GEN    primedec_apply_kummer(GEN nf,GEN pol,long e,GEN p);
GEN    prodid(GEN nf, GEN I);
GEN    rnfallbase(GEN nf, GEN *ppol, GEN *pD, GEN *pd, GEN *pfi);
GEN    rnf_basM(GEN rnf);
GEN    special_anti_uniformizer(GEN nf, GEN pr);
GEN    subgroupcondlist(GEN cyc, GEN bound, GEN listKer);
void   testprimes(GEN bnf, GEN bound);
GEN    to_Fp_simple(GEN nf, GEN x, GEN ffproj);
GEN    unif_mod_fZ(GEN pr, GEN F);
GEN    unnf_minus_x(GEN x);
GEN    ideallog_sgn(GEN nf, GEN x, GEN sgn, GEN bid);
GEN    zlog_units(GEN nf, GEN U, GEN sgnU, GEN bid);
GEN    zlog_units_noarch(GEN nf, GEN U, GEN bid);

/* Dedekind zeta */
GEN  zeta_get_limx(long r1, long r2, long bit);
long zeta_get_i0(long r1, long r2, long bit, GEN limx);
long zeta_get_N0(GEN C,  GEN limx);
