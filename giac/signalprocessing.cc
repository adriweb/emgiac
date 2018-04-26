/*
 * signalprocessing.cc
 *
 * (c) 2018 by Luka Marohnić
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#include "signalprocessing.h"
#include "giacPCH.h"
#include "giac.h"
#include <sstream>
#include <string>

using namespace std;

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

gen _cross_correlation(const gen &g,GIAC_CONTEXT) {
  if (g.type==_STRNG && g.subtype==-1) return g;
  if (g.type!=_VECT || g.subtype!=_SEQ__VECT)
    return gentypeerr(contextptr);
  vecteur &args=*g._VECTptr;
  if (args.size()!=2 || args.front().type!=_VECT || args.back().type!=_VECT)
    return gensizeerr(contextptr);
  vecteur A=*args.front()._VECTptr,B=*args.back()._VECTptr;
  int m=A.size(),n=B.size(),N=n>m?n:m;
  int N2=1;
  while (N2<N) N2*=2;
  N=N2;
  while (int(A.size())<N || int(B.size())<N) {
    if (int(A.size())<N)
      A.push_back(0);
    if (int(B.size())<N)
      B.push_back(0);
  }
  A.resize(2*N-1);
  B.resize(2*N-1);
  vecteur a=*_fft(A,contextptr)._VECTptr,b=*_fft(B,contextptr)._VECTptr;
  vecteur cc_ffted=*_pointprod(makesequence(a,conj(b,contextptr)),contextptr)._VECTptr;
  vecteur cc=*_epsilon2zero(_ifft(cc_ffted,contextptr),contextptr)._VECTptr;
  reverse(cc.begin(),cc.begin()+N);
  reverse(cc.begin()+N,cc.end());
  return vecteur(cc.begin()+N-m,cc.end()-N+n);
}
static const char _cross_correlation_s []="cross_correlation";
static define_unary_function_eval (__cross_correlation,&_cross_correlation,_cross_correlation_s);
define_unary_function_ptr5(at_cross_correlation,alias_at_cross_correlation,&__cross_correlation,0,true)

gen _auto_correlation(const gen &g,GIAC_CONTEXT) {
  if (g.type==_STRNG) {
    if (g.subtype==-1) return g;
    return _auto_correlation(_readwav(g,contextptr),contextptr);
  }
  return _cross_correlation(makesequence(g,g),contextptr);
}
static const char _auto_correlation_s []="auto_correlation";
static define_unary_function_eval (__auto_correlation,&_auto_correlation,_auto_correlation_s);
define_unary_function_ptr5(at_auto_correlation,alias_at_auto_correlation,&__auto_correlation,0,true)

gen _lowpass(const gen &g,GIAC_CONTEXT) {
  if (g.type==_STRNG && g.subtype==-1) return g;
  if (g.type!=_VECT || g.subtype!=_SEQ__VECT || g._VECTptr->front().type!=_VECT)
    return gentypeerr(contextptr);
  vecteur &args=*g._VECTptr,&data=*args.front()._VECTptr;
  if (args.size()<2)
    return gensizeerr(contextptr);
  if (_evalf(args.at(1),contextptr).type!=_DOUBLE_)
    return gentypeerr(contextptr);
  double cutoff=_evalf(args.at(1),contextptr).DOUBLE_val();
  int samplerate=44100;
  if (args.size()>2) {
    if (!is_integer(args.at(2)) || (samplerate=args.at(2).val)<=0)
      return gentypeerr(contextptr);
  }
  double rc=0.15915494309/cutoff;
  double dt=1.0/samplerate;
	double alpha=dt/(rc+dt);
  int n=data.size();
	vecteur output(n,0);
	for (int i=1;i<n;++i) {
		gen &s=output[i-1];
		output[i]=s+alpha*(data[i]-s);
	}
  return output;
}
static const char _lowpass_s []="lowpass";
static define_unary_function_eval (__lowpass,&_lowpass,_lowpass_s);
define_unary_function_ptr5(at_lowpass,alias_at_lowpass,&__lowpass,0,true)

gen _highpass(const gen &g,GIAC_CONTEXT) {
  if (g.type==_STRNG && g.subtype==-1) return g;
  if (g.type!=_VECT || g.subtype!=_SEQ__VECT || g._VECTptr->front().type!=_VECT)
    return gentypeerr(contextptr);
  vecteur &args=*g._VECTptr,&data=*args.front()._VECTptr;
  if (args.size()<2)
    return gensizeerr(contextptr);
  if (_evalf(args.at(1),contextptr).type!=_DOUBLE_)
    return gentypeerr(contextptr);
  double cutoff=_evalf(args.at(1),contextptr).DOUBLE_val();
  int samplerate=44100;
  if (args.size()>2) {
    if (!is_integer(args.at(2)) || (samplerate=args.at(2).val)<=0)
      return gentypeerr(contextptr);
  }
  double rc=0.15915494309/cutoff;
  double dt=1.0/samplerate;
	double alpha=rc/(rc+dt);
  int n=data.size();
	vecteur output(n,0);
	for (int i=1;i<n;++i) {
		output[i]=alpha*(output[i-1]+data[i]-data[i-1]);
	}
  return output;
}
static const char _highpass_s []="highpass";
static define_unary_function_eval (__highpass,&_highpass,_highpass_s);
define_unary_function_ptr5(at_highpass,alias_at_highpass,&__highpass,0,true)

gen _convolution(const gen &g,GIAC_CONTEXT) {
  if (g.type==_STRNG && g.subtype==-1) return g;
  if (g.type!=_VECT || g.subtype!=_SEQ__VECT)
    return gentypeerr(contextptr);
  vecteur &args=*g._VECTptr;
  if (args.size()!=2 || args.front().type!=_VECT || args.back().type!=_VECT)
    return gensizeerr(contextptr);
  vecteur A=*args.front()._VECTptr,B=*args.back()._VECTptr;
  int lenA=A.size(),lenB=B.size(),l=lenA>lenB?lenA:lenB,l2=1;
  while (l2<l) l2*=2;
  int len=2*l2-1;
  A.resize(len);
  B.resize(len);
  vecteur a=*_fft(A,contextptr)._VECTptr,b=*_fft(B,contextptr)._VECTptr;
  vecteur cv=*_ifft(_pointprod(makesequence(a,b),contextptr),contextptr)._VECTptr;
  cv.resize(lenA+lenB-1);
  return cv;
}
static const char _convolution_s []="convolution";
static define_unary_function_eval (__convolution,&_convolution,_convolution_s);
define_unary_function_ptr5(at_convolution,alias_at_convolution,&__convolution,0,true)

vecteur apply_window_function(const gen &expr,const identificateur &k,const vecteur &data,int start,int len,GIAC_CONTEXT) {
  vecteur output(len);
  for (int j=0;j<len;++j) {
    output[j]=_evalf(subst(expr,k,gen((double)j),false,contextptr),contextptr)*data[start+j];
  }
  return output;
}

bool nivelate(vecteur &data,int k,const gen &b,const gen &val,const unary_function_ptr *comp,GIAC_CONTEXT) {
  gen r;
  if (has_i(data[k]) && !is_zero((r=_abs(data[k],contextptr)))) {
    if (_eval(symbolic(comp,makevecteur(r,b)),contextptr).val!=0) {
      data[k]=val*data[k]/r;
      return true;
    }
  }
  else {
    if (_eval(symbolic(comp,makevecteur(data[k],b)),contextptr).val!=0) {
      data[k]=val;
      return true;
    }
  }
  return false;
}

gen _threshold(const gen &g,GIAC_CONTEXT) {
  if (g.type==_STRNG && g.subtype==-1) return g;
  if (g.type!=_VECT || g.subtype!=_SEQ__VECT)
    return gentypeerr(contextptr);
  vecteur &args=*g._VECTptr;
  if (int(args.size())<2)
    return gensizeerr(contextptr);
  if (args.front().type!=_VECT)
    return gentypeerr(contextptr);
  vecteur &data=*args.front()._VECTptr;
  gen &bnd=args.at(1);
  int n=data.size();
  vecteur output=data;
  if (bnd.type==_VECT) {
    if (int(bnd._VECTptr->size())!=2)
      return gensizeerr(contextptr);
    gen lb=bnd._VECTptr->front(),ub=bnd._VECTptr->back(),lval,uval;
    if (lb.is_symb_of_sommet(at_equal)) {
      lval=_rhs(lb,contextptr);
      lb=_lhs(lb,contextptr);
    } else lval=lb;
    if (ub.is_symb_of_sommet(at_equal)) {
      uval=_rhs(ub,contextptr);
      ub=_lhs(ub,contextptr);
    } else uval=ub;
    for (int k=0;k<n;++k) {
      if (!nivelate(output,k,lb,lval,at_inferieur_strict,contextptr))
        nivelate(output,k,ub,uval,at_superieur_strict,contextptr);
    }
  } else {
    gen val;
    if (bnd.is_symb_of_sommet(at_equal)) {
      val=_rhs(bnd,contextptr);
      bnd=_lhs(bnd,contextptr);
    } else val=bnd;
    if (_evalf(bnd,contextptr).type!=_DOUBLE_)
      return gentypeerr(contextptr);
    gen comp=at_inferieur_strict,isabs;
    bool absolute=false;
    for (const_iterateur it=args.begin()+2;it!=args.end();++it) {
      if (*it==at_superieur_strict || *it==at_superieur_egal || *it==at_inferieur_egal)
        comp=*it;
      isabs=gen(1);
      if (*it==at_abs || (it->is_symb_of_sommet(at_equal) &&
                          it->_SYMBptr->feuille._VECTptr->front()==at_abs &&
                          (isabs=it->_SYMBptr->feuille._VECTptr->back()).type==_INT_ &&
                          isabs.subtype==_INT_BOOLEAN)) {
        if (has_i(data) || !is_strictly_positive(bnd,contextptr))
          return gentypeerr(contextptr);
        absolute=(bool)isabs.val;
      }
    }
    for (int k=0;k<n;++k) {
      if (absolute) {
        if (_eval(symbolic(comp._FUNCptr,makevecteur(_abs(data[k],contextptr),bnd)),contextptr).val!=0)
          output[k]=is_positive(data[k],contextptr)?val:-val;
      } else nivelate(output,k,bnd,val,comp._FUNCptr,contextptr);
    }
  }
  return _eval(output,contextptr);
}
static const char _threshold_s []="threshold";
static define_unary_function_eval (__threshold,&_threshold,_threshold_s);
define_unary_function_ptr5(at_threshold,alias_at_threshold,&__threshold,0,true)

bool parse_window_parameters(const gen &g,vecteur &data,int &start,int &len,double *alpha,GIAC_CONTEXT) {
  start=0;
  if (g.type==_VECT && g.subtype!=_SEQ__VECT) {
    data=*g._VECTptr;
    len=data.size();
    return true;
  }
  if (g.type!=_VECT || g.subtype!=_SEQ__VECT || g._VECTptr->size()>3 || g._VECTptr->front().type!=_VECT)
    return false;
  vecteur &args=*g._VECTptr;
  data=*args.front()._VECTptr;
  len=data.size();
  bool has_alpha;
  if (_evalf(args.at(1),contextptr).type==_DOUBLE_) {
    has_alpha=true;
    if (!alpha)
      return false;
    *alpha=_evalf(args.at(1),contextptr).DOUBLE_val();
  } else if (args.size()>2) return false;
  if (args.back().is_symb_of_sommet(at_interval)) {
    gen lh=_lhs(args.back(),contextptr),rh=_rhs(args.back(),contextptr);
    if (!lh.is_integer() || !rh.is_integer() || lh.val<0 || rh.val>=len)
      return false;
    start=lh.val;
    len=rh.val-start+1;
  } else if (args.size()>2 || !has_alpha) return false;
  return true;
}

gen _bartlett_hann_window(const gen &g,GIAC_CONTEXT) {
  if (g.type==_STRNG && g.subtype==-1) return g;
  vecteur data;
  int start,N;
  identificateur k(" k");
  if (!parse_window_parameters(g,data,start,N,0,contextptr))
    return gentypeerr(contextptr);
  double a=0.62,b=0.48,c=0.38;
  gen expr=a-b*_abs(k/(N-1)-fraction(1,2),contextptr)-c*cos(2*k*_IDNT_pi()/(N-1),contextptr);
  return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _bartlett_hann_window_s []="bartlett_hann_window";
static define_unary_function_eval (__bartlett_hann_window,&_bartlett_hann_window,_bartlett_hann_window_s);
define_unary_function_ptr5(at_bartlett_hann_window,alias_at_bartlett_hann_window,&__bartlett_hann_window,0,true)

gen _blackman_harris_window(const gen &g,GIAC_CONTEXT) {
  if (g.type==_STRNG && g.subtype==-1) return g;
  vecteur data;
  int start,N;
  identificateur k(" k");
  if (!parse_window_parameters(g,data,start,N,0,contextptr))
    return gentypeerr(contextptr);
  gen a(0.35875),b(0.48829),c(0.14128),d(0.01168);
  gen K=k*_IDNT_pi()/(N-1),expr=a-b*cos(2*K,contextptr)+c*cos(4*K,contextptr)-d*cos(6*K,contextptr);
  return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _blackman_harris_window_s []="blackman_harris_window";
static define_unary_function_eval (__blackman_harris_window,&_blackman_harris_window,_blackman_harris_window_s);
define_unary_function_ptr5(at_blackman_harris_window,alias_at_blackman_harris_window,&__blackman_harris_window,0,true)

gen _blackman_window(const gen &g,GIAC_CONTEXT) {
  if (g.type==_STRNG && g.subtype==-1) return g;
  double alpha=0.16;
  vecteur data;
  int start,N;
  identificateur k(" k");
  if (!parse_window_parameters(g,data,start,N,&alpha,contextptr) || alpha<=0)
    return gentypeerr(contextptr);
  gen K=k*_IDNT_pi()/(N-1),expr=(1-alpha)/2-cos(2*K,contextptr)/2+alpha*cos(4*K,contextptr)/2;
  return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _blackman_window_s []="blackman_window";
static define_unary_function_eval (__blackman_window,&_blackman_window,_blackman_window_s);
define_unary_function_ptr5(at_blackman_window,alias_at_blackman_window,&__blackman_window,0,true)

gen _bohman_window(const gen &g,GIAC_CONTEXT) {
  if (g.type==_STRNG && g.subtype==-1) return g;
  vecteur data;
  int start,N;
  identificateur k(" k");
  if (!parse_window_parameters(g,data,start,N,0,contextptr))
    return gentypeerr(contextptr);
  gen K=_abs(2*k/(N-1)-1,contextptr),expr=(1-K)*cos(_IDNT_pi()*K,contextptr)+sin(_IDNT_pi()*K,contextptr)/_IDNT_pi();
  return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _bohman_window_s []="bohman_window";
static define_unary_function_eval (__bohman_window,&_bohman_window,_bohman_window_s);
define_unary_function_ptr5(at_bohman_window,alias_at_bohman_window,&__bohman_window,0,true)

gen _cosine_window(const gen &g,GIAC_CONTEXT) {
  if (g.type==_STRNG && g.subtype==-1) return g;
  double alpha=1.0;
  vecteur data;
  int start,N;
  identificateur k(" k");
  if (!parse_window_parameters(g,data,start,N,&alpha,contextptr) || alpha<=0)
    return gentypeerr(contextptr);
  gen expr=exp(alpha*ln(sin(k*_IDNT_pi()/(N-1),contextptr),contextptr),contextptr);
  return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _cosine_window_s []="cosine_window";
static define_unary_function_eval (__cosine_window,&_cosine_window,_cosine_window_s);
define_unary_function_ptr5(at_cosine_window,alias_at_cosine_window,&__cosine_window,0,true)

gen _gaussian_window(const gen &g,GIAC_CONTEXT) {
  if (g.type==_STRNG && g.subtype==-1) return g;
  double alpha=0.1;
  vecteur data;
  int start,N;
  identificateur k(" k");
  if (!parse_window_parameters(g,data,start,N,&alpha,contextptr) || alpha<=0 || alpha>0.5)
    return gentypeerr(contextptr);
  gen c=(N-1)/2.0,expr=exp(-pow((k-c)/(alpha*c),2)/2,contextptr);
  return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _gaussian_window_s []="gaussian_window";
static define_unary_function_eval (__gaussian_window,&_gaussian_window,_gaussian_window_s);
define_unary_function_ptr5(at_gaussian_window,alias_at_gaussian_window,&__gaussian_window,0,true)

gen _hamming_window(const gen &g,GIAC_CONTEXT) {
  if (g.type==_STRNG && g.subtype==-1) return g;
  vecteur data;
  int start,N;
  identificateur k(" k");
  if (!parse_window_parameters(g,data,start,N,0,contextptr))
    return gentypeerr(contextptr);
  gen a(0.54),b(0.46),expr=a-b*cos(2*_IDNT_pi()*k/(N-1),contextptr);
  return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _hamming_window_s []="hamming_window";
static define_unary_function_eval (__hamming_window,&_hamming_window,_hamming_window_s);
define_unary_function_ptr5(at_hamming_window,alias_at_hamming_window,&__hamming_window,0,true)

gen _hann_poisson_window(const gen &g,GIAC_CONTEXT) {
  if (g.type==_STRNG && g.subtype==-1) return g;
  double alpha=1;
  vecteur data;
  int start,N;
  identificateur k(" k");
  if (!parse_window_parameters(g,data,start,N,&alpha,contextptr))
    return gentypeerr(contextptr);
  gen K=2*_IDNT_pi()*k/(N-1);
  gen expr=(1-cos(K,contextptr))*exp(-alpha*_abs(N-1-2*k,contextptr)/(N-1),contextptr)/2;
  return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _hann_poisson_window_s []="hann_poisson_window";
static define_unary_function_eval (__hann_poisson_window,&_hann_poisson_window,_hann_poisson_window_s);
define_unary_function_ptr5(at_hann_poisson_window,alias_at_hann_poisson_window,&__hann_poisson_window,0,true)

gen _hann_window(const gen &g,GIAC_CONTEXT) {
  if (g.type==_STRNG && g.subtype==-1) return g;
  vecteur data;
  int start,N;
  identificateur k(" k");
  if (!parse_window_parameters(g,data,start,N,0,contextptr))
    return gentypeerr(contextptr);
  gen expr=pow(sin(_IDNT_pi()*k/(N-1),contextptr),2);
  return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _hann_window_s []="hann_window";
static define_unary_function_eval (__hann_window,&_hann_window,_hann_window_s);
define_unary_function_ptr5(at_hann_window,alias_at_hann_window,&__hann_window,0,true)

gen _parzen_window(const gen &g,GIAC_CONTEXT) {
  if (g.type==_STRNG && g.subtype==-1) return g;
  vecteur data;
  int start,N;
  identificateur k(" k");
  if (!parse_window_parameters(g,data,start,N,0,contextptr))
    return gentypeerr(contextptr);
  gen K=1-2*k/(N-1),cond=symb_inferieur_egal(symbolic(at_abs,(N-1)/2.0-k),(N-1)/4.0);
  gen f1=1-6*pow(K,2)*(1-_abs(K,contextptr)),f2=2*pow(1-_abs(K,contextptr),3);
  gen expr=symbolic(at_when,makevecteur(cond,f1,f2));
  return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _parzen_window_s []="parzen_window";
static define_unary_function_eval (__parzen_window,&_parzen_window,_parzen_window_s);
define_unary_function_ptr5(at_parzen_window,alias_at_parzen_window,&__parzen_window,0,true)

gen _poisson_window(const gen &g,GIAC_CONTEXT) {
  if (g.type==_STRNG && g.subtype==-1) return g;
  double alpha=1.0;
  vecteur data;
  int start,N;
  identificateur k(" k");
  if (!parse_window_parameters(g,data,start,N,&alpha,contextptr))
    return gentypeerr(contextptr);
  gen expr=exp(-alpha*_abs(2*k/(N-1)-1,contextptr),contextptr);
  return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _poisson_window_s []="poisson_window";
static define_unary_function_eval (__poisson_window,&_poisson_window,_poisson_window_s);
define_unary_function_ptr5(at_poisson_window,alias_at_poisson_window,&__poisson_window,0,true)

gen _riemann_window(const gen &g,GIAC_CONTEXT) {
  if (g.type==_STRNG && g.subtype==-1) return g;
  vecteur data;
  int start,N;
  identificateur k(" k");
  if (!parse_window_parameters(g,data,start,N,0,contextptr))
    return gentypeerr(contextptr);
  gen K=(2*k/(N-1)-1)*_IDNT_pi(),cond=symbolic(at_same,makevecteur(k,(N-1)/2.0));
  gen expr=symbolic(at_when,makevecteur(cond,1,sin(K,contextptr)/K));
  return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _riemann_window_s []="riemann_window";
static define_unary_function_eval (__riemann_window,&_riemann_window,_riemann_window_s);
define_unary_function_ptr5(at_riemann_window,alias_at_riemann_window,&__riemann_window,0,true)

gen _triangle_window(const gen &g,GIAC_CONTEXT) {
  if (g.type==_STRNG && g.subtype==-1) return g;
  double L=0;
  vecteur data;
  int start,N;
  identificateur k(" k");
  if (!parse_window_parameters(g,data,start,N,&L,contextptr) || (L!=1 && L!=-1 && L!=0))
    return gentypeerr(contextptr);
  gen expr=1-_abs((2*k-N+1)/(N+L),contextptr);
  return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _triangle_window_s []="triangle_window";
static define_unary_function_eval (__triangle_window,&_triangle_window,_triangle_window_s);
define_unary_function_ptr5(at_triangle_window,alias_at_triangle_window,&__triangle_window,0,true)

gen _tukey_window(const gen &g,GIAC_CONTEXT) {
  if (g.type==_STRNG && g.subtype==-1) return g;
  double alpha=0.5;
  vecteur data;
  int start,N;
  identificateur k(" k");
  if (!parse_window_parameters(g,data,start,N,&alpha,contextptr) || alpha<0 || alpha>1)
    return gentypeerr(contextptr);
  double p=alpha*(N-1)/2.0,q=1-alpha/2;
  gen cond1=symb_inferieur_strict(k,p),cond2=symb_inferieur_egal(k,q*(N-1));
  gen f1=(1+cos(_IDNT_pi()*(k/p-1),contextptr))/2,f2=(1+cos(_IDNT_pi()*(k/p+1-2/alpha),contextptr))/2;
  gen expr=symbolic(at_piecewise,makevecteur(cond1,f1,cond2,1,f2));
  return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _tukey_window_s []="tukey_window";
static define_unary_function_eval (__tukey_window,&_tukey_window,_tukey_window_s);
define_unary_function_ptr5(at_tukey_window,alias_at_tukey_window,&__tukey_window,0,true)

gen _welch_window(const gen &g,GIAC_CONTEXT) {
  if (g.type==_STRNG && g.subtype==-1) return g;
  vecteur data;
  int start,N;
  identificateur k(" k");
  if (!parse_window_parameters(g,data,start,N,0,contextptr))
    return gentypeerr(contextptr);
  double p=(N-1)/2.0;
  gen expr=1-pow(1-k/p,2);
  return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _welch_window_s []="welch_window";
static define_unary_function_eval (__welch_window,&_welch_window,_welch_window_s);
define_unary_function_ptr5(at_welch_window,alias_at_welch_window,&__welch_window,0,true)

#ifndef NO_NAMESPACE_GIAC
}
#endif // ndef NO_NAMESPACE_GIAC
