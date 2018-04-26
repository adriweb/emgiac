#define PARI_TUNE

#ifdef PARI_TUNE
extern long SQRI_KARATSUBA_LIMIT;
extern long MULII_KARATSUBA_LIMIT;
extern long MULRR_MULII_LIMIT;
extern long SQRI_FFT_LIMIT;
extern long MULII_FFT_LIMIT;
extern long Fp_POW_REDC_LIMIT;
extern long Fp_POW_BARRETT_LIMIT;
extern long INVMOD_GMP_LIMIT;
extern long DIVRR_GMP_LIMIT;
extern long Flx_MUL_KARATSUBA_LIMIT;
extern long Flx_SQR_KARATSUBA_LIMIT;
extern long Flx_MUL_HALFMULII_LIMIT;
extern long Flx_SQR_HALFSQRI_LIMIT;
extern long Flx_MUL_MULII_LIMIT;
extern long Flx_SQR_SQRI_LIMIT;
extern long Flx_MUL_MULII2_LIMIT;
extern long Flx_SQR_SQRI2_LIMIT;
extern long Flx_INVBARRETT_LIMIT;
extern long Flx_DIVREM_BARRETT_LIMIT;
extern long Flx_REM_BARRETT_LIMIT;
extern long Flx_BARRETT_LIMIT;
extern long Flx_HALFGCD_LIMIT;
extern long Flx_GCD_LIMIT;
extern long Flx_EXTGCD_LIMIT;
extern long FpX_INVBARRETT_LIMIT;
extern long FpX_DIVREM_BARRETT_LIMIT;
extern long FpX_REM_BARRETT_LIMIT;
extern long FpX_BARRETT_LIMIT;
extern long FpX_HALFGCD_LIMIT;
extern long FpX_GCD_LIMIT;
extern long FpX_EXTGCD_LIMIT;
extern long EXPNEWTON_LIMIT;
extern long INVNEWTON_LIMIT;
extern long LOGAGM_LIMIT;
extern long LOGAGMCX_LIMIT;
extern long AGM_ATAN_LIMIT;
extern long RgX_SQR_LIMIT;
extern long RgX_MUL_LIMIT;
#else
#  define SQRI_KARATSUBA_LIMIT     __SQRI_KARATSUBA_LIMIT
#  define MULII_KARATSUBA_LIMIT    __MULII_KARATSUBA_LIMIT
#  define MULRR_MULII_LIMIT        __MULRR_MULII_LIMIT
#  define SQRI_FFT_LIMIT           __SQRI_FFT_LIMIT
#  define MULII_FFT_LIMIT          __MULII_FFT_LIMIT
#  define Fp_POW_REDC_LIMIT        __Fp_POW_REDC_LIMIT
#  define Fp_POW_BARRETT_LIMIT     __Fp_POW_BARRETT_LIMIT
#  define INVMOD_GMP_LIMIT         __INVMOD_GMP_LIMIT
#  define DIVRR_GMP_LIMIT          __DIVRR_GMP_LIMIT
#  define EXPNEWTON_LIMIT          __EXPNEWTON_LIMIT
#  define INVNEWTON_LIMIT          __INVNEWTON_LIMIT
#  define LOGAGM_LIMIT             __LOGAGM_LIMIT
#  define LOGAGMCX_LIMIT           __LOGAGMCX_LIMIT
#  define AGM_ATAN_LIMIT           __AGM_ATAN_LIMIT
#  define Flx_MUL_KARATSUBA_LIMIT  __Flx_MUL_KARATSUBA_LIMIT
#  define Flx_SQR_KARATSUBA_LIMIT  __Flx_SQR_KARATSUBA_LIMIT
#  define Flx_MUL_HALFMULII_LIMIT  __Flx_MUL_HALFMULII_LIMIT
#  define Flx_SQR_HALFSQRI_LIMIT   __Flx_SQR_HALFSQRI_LIMIT
#  define Flx_MUL_MULII_LIMIT      __Flx_MUL_MULII_LIMIT
#  define Flx_SQR_SQRI_LIMIT       __Flx_SQR_SQRI_LIMIT
#  define Flx_MUL_MULII2_LIMIT     __Flx_MUL_MULII2_LIMIT
#  define Flx_SQR_SQRI2_LIMIT      __Flx_SQR_SQRI2_LIMIT
#  define Flx_INVBARRETT_LIMIT     __Flx_INVBARRETT_LIMIT
#  define Flx_DIVREM_BARRETT_LIMIT __Flx_DIVREM_BARRETT_LIMIT
#  define Flx_REM_BARRETT_LIMIT    __Flx_REM_BARRETT_LIMIT
#  define Flx_BARRETT_LIMIT        __Flx_BARRETT_LIMIT
#  define Flx_HALFGCD_LIMIT        __Flx_HALFGCD_LIMIT
#  define Flx_GCD_LIMIT            __Flx_GCD_LIMIT
#  define Flx_EXTGCD_LIMIT         __Flx_EXTGCD_LIMIT
#  define FpX_INVBARRETT_LIMIT     __FpX_INVBARRETT_LIMIT
#  define FpX_DIVREM_BARRETT_LIMIT __FpX_DIVREM_BARRETT_LIMIT
#  define FpX_REM_BARRETT_LIMIT    __FpX_REM_BARRETT_LIMIT
#  define FpX_BARRETT_LIMIT        __FpX_BARRETT_LIMIT
#  define FpX_HALFGCD_LIMIT        __FpX_HALFGCD_LIMIT
#  define FpX_GCD_LIMIT            __FpX_GCD_LIMIT
#  define FpX_EXTGCD_LIMIT         __FpX_EXTGCD_LIMIT
#  define RgX_SQR_LIMIT            __RgX_SQR_LIMIT
#  define RgX_MUL_LIMIT            __RgX_MUL_LIMIT
#endif
