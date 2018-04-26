// -*- mode:C++ ; compile-command: "emcc opengl.cc -I. -I.. -DHAVE_CONFIG_H -DIN_GIAC -DGIAC_GENERIC_CONSTANTS -DNO_STDEXCEPT -Os -s ALLOW_MEMORY_GROWTH=1 -s LEGACY_GL_EMULATION=1" -*-
/*
 *  Copyright (C) 2006,2014 B. Parisse, Institut Fourier, 38402 St Martin d'Heres
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "SDL/SDL.h"

#include <fstream>
#include "vector.h"
#include <algorithm>
#include <fcntl.h>
#include <cmath>
#include <time.h> // for nanosleep
#include <stdio.h>
#include <dirent.h>
#include <sys/stat.h> // auto-recovery function
#include "path.h"
#ifndef IN_GIAC
#include <giac/misc.h>
#else
#include "misc.h"
#endif

using namespace std;
using namespace giac;


  void fltk_point(int deltax,int deltay,int i0,int j0,int epaisseur_point,int type_point){
    fl_line_style(FL_SOLID,epaisseur_point-1,0);
    switch (type_point){
    case 1: // losange
      fl_line(deltax+i0-epaisseur_point,deltay+j0,deltax+i0,deltay+j0-epaisseur_point);
      fl_line(deltax+i0,deltay+j0-epaisseur_point,deltax+i0+epaisseur_point,deltay+j0);
      fl_line(deltax+i0-epaisseur_point,deltay+j0,deltax+i0,deltay+j0+epaisseur_point);
      fl_line(deltax+i0,deltay+j0+epaisseur_point,deltax+i0+epaisseur_point,deltay+j0);
      break;
    case 2: // croix verticale
      fl_line(deltax+i0,deltay+j0-epaisseur_point,deltax+i0,deltay+j0+epaisseur_point);
      fl_line(deltax+i0-epaisseur_point,deltay+j0,deltax+i0+epaisseur_point,deltay+j0);
      break;
    case 3: // carre
      fl_line(deltax+i0-epaisseur_point,deltay+j0-epaisseur_point,deltax+i0-epaisseur_point,deltay+j0+epaisseur_point);
      fl_line(deltax+i0+epaisseur_point,deltay+j0-epaisseur_point,deltax+i0+epaisseur_point,deltay+j0+epaisseur_point);
      fl_line(deltax+i0-epaisseur_point,deltay+j0-epaisseur_point,deltax+i0+epaisseur_point,deltay+j0-epaisseur_point);
      fl_line(deltax+i0-epaisseur_point,deltay+j0+epaisseur_point,deltax+i0+epaisseur_point,deltay+j0+epaisseur_point);
      break;
    case 5: // triangle
      fl_line(deltax+i0-epaisseur_point,deltay+j0,deltax+i0,deltay+j0-epaisseur_point);
      fl_line(deltax+i0,deltay+j0-epaisseur_point,deltax+i0+epaisseur_point,deltay+j0);
      fl_line(deltax+i0-epaisseur_point,deltay+j0,deltax+i0+epaisseur_point,deltay+j0);
      break;
    case 7: // point
      if (epaisseur_point>2)
	fl_arc(deltax+i0-(epaisseur_point-1),deltay+j0-(epaisseur_point-1),2*(epaisseur_point-1),2*(epaisseur_point-1),0,360);
      else
	fl_line(deltax+i0,deltay+j0,deltax+i0+1,deltay+j0);
      break;
    case 6: // etoile
      fl_line(deltax+i0-epaisseur_point,deltay+j0,deltax+i0+epaisseur_point,deltay+j0);
      // no break to add the following lines
    case 0: // 0 croix diagonale
      fl_line(deltax+i0-epaisseur_point,deltay+j0-epaisseur_point,deltax+i0+epaisseur_point,deltay+j0+epaisseur_point);
      fl_line(deltax+i0-epaisseur_point,deltay+j0+epaisseur_point,deltax+i0+epaisseur_point,deltay+j0-epaisseur_point);
      break;
    default: // 4 nothing drawn
      break;
    }
    fl_line_style(0);
  }

  inline void swapint(int & i0,int & i1){
    int tmp=i0;
    i0=i1;
    i1=tmp;
  }

  void check_fl_draw(const char * ch,int i0,int j0,int imin,int jmin,int di,int dj,int delta_i,int delta_j){
    /* int n=fl_size();
       if (j0>=jmin-n && j0<=jmin+dj+n) */
    // cerr << i0 << " " << j0 << endl;
    if (strlen(ch)>2000)
      fl_draw("String too long for display",i0+delta_i,j0+delta_j);
    else
      fl_draw(ch,i0+delta_i,j0+delta_j);
  }

  void check_fl_point(int i0,int j0,int imin,int jmin,int di,int dj,int delta_i,int delta_j){
    /* if (i0>=imin && i0<=imin+di && j0>=jmin && j0<=jmin+dj) */
      fl_point(i0+delta_i,j0+delta_j);
  }

  void check_fl_rect(int i0,int j0,int i1,int j1,int imin,int jmin,int di,int dj,int delta_i,int delta_j){
    /*    bool clipped=false;
    if (imin>i0 || jmin>j0){ 
      check_fl_line(i0,j0,i0+i1,j0,imin,jmin,di,dj,delta_i,delta_j);
      check_fl_line(i0,j0,i0,j0+j1,imin,jmin,di,dj,delta_i,delta_j);
      check_fl_line(i0,j0+j1,i0+i1,j0+j1,imin,jmin,di,dj,delta_i,delta_j);
      check_fl_line(i0+i1,j0,i0+i1,j0+j1,imin,jmin,di,dj,delta_i,delta_j);
    }
    else 
      fl_rect(i0+delta_i,j0+delta_j,min(i1,di),min(j1,dj));
    */
    fl_rect(i0+delta_i,j0+delta_j,i1,j1);
  }

  void check_fl_rectf(int x,int y,int w,int h,int imin,int jmin,int di,int dj,int delta_i,int delta_j){
    /*    if (x<imin){ 
      w -= (imin-x); // di -= (imin-x); 
      x=imin; 
    }
    if (y<jmin){ 
      h -= (jmin-y); // dj -= (jmin-y); 
      y=jmin; 
    }
    if (w<=0 || di<=0 || h<=0 || dj<=0) 
      return ;
    fl_rectf(x+delta_i,y+delta_j,min(w,di),min(h,dj));
    */
    fl_rectf(x+delta_i,y+delta_j,w,h);
  }

  // Calls fl_line, checking with bounds
  void check_fl_line(int i0,int j0,int i1,int j1,int imin,int jmin,int di,int dj,int delta_i,int delta_j){
    /* int imax=imin+di,jmax=jmin+dj;
    if (i0>i1){
      swapint(i0,i1);
      swapint(j0,j1);
    }
    if (i0>=imax || i1<=imin)
      return;
    if (i0!=i1){
      // Compute line slope
      double m=(j1-j0)/double(i1-i0);
      if (i0<imin){ // replace i0 by imin, compute corresp. j0
	j0 += int(floor((imin-i0)*m+.5));
	i0 = imin;
      }
      if (i1>imax){
	j1 += int(floor((imax-i1)*m+.5));
	i1 = imax;
      }
    }
    if (j0>j1){
      if (j1>=jmax || j0<=jmin)
	return;
      // Compute line slope
      double m=(i1-i0)/double(j1-j0);
      if (j0>jmax){
	i0 += int(floor((jmax-j0)*m+.5));
	j0 = jmax;
      }
      if (j1<jmin){
	i1 += int(floor((jmin-j1)*m+.5));
	j1 = jmin;
      }
    }
    else {
      if (j0>=jmax || j1<=jmin)
	return;
      if (j0!=j1){
	// Compute line slope
	double m=(i1-i0)/double(j1-j0);
	if (j0<jmin){
	  i0 += int(floor((jmin-j0)*m+.5));
	  j0 = jmin;
	}
	if (j1>jmax){
	  i1 += int(floor((jmax-j1)*m+.5));
	  j1 = jmax;
	}
      }
      } */
    fl_line(i0+delta_i,j0+delta_j,i1+delta_i,j1+delta_j);
  }

  int logplot_points=20;

  void checklog_fl_line(double i0,double j0,double i1,double j1,double deltax,double deltay,bool logx,bool logy,double window_xmin,double x_scale,double window_ymax,double y_scale){
    if (!logx && !logy){
      fl_line(round(i0+deltax),round(j0+deltay),round(i1+deltax),round(j1+deltay));
      return;
    }
    // interpolate (i0,j0)->(i1,j1)
    // if log enabled i0=(log10(x0)-window_xmin)*x_scale -> x0=pow10(i0/x_scale+window_xmin)
    // j0=(window_ymax-log10(y1))*y_scale -> y1=pow10(window_ymax-j0/y_scale)
    if (logx && logy){
      double prevx=i0,prevy=j0,curx,cury;
      double I0=pow10(i0/x_scale+window_xmin),I1=pow10(i1/x_scale+window_xmin),
	J0=pow10(window_ymax-j0/y_scale),J1=pow10(window_ymax-j1/y_scale);
      for (int i=1;i<logplot_points;++i){
	double t=double(i)/logplot_points;
	curx=(log10(I0+t*(I1-I0))-window_xmin)*x_scale;
	cury=(window_ymax-log10(J0+t*(J1-J0)))*y_scale;
	fl_line(round(prevx+deltax),round(prevy+deltay),round(curx+deltax),round(cury+deltay));
	prevx=curx; prevy=cury;
      }
      return;
    }
    if (logy){
      double prevx=i0,prevy=j0,curx,cury;
      double J0=pow10(window_ymax-j0/y_scale),J1=pow10(window_ymax-j1/y_scale);
      for (int i=1;i<logplot_points;++i){
	double t=double(i)/logplot_points;
	curx=i0+t*(i1-i0);
	cury=(window_ymax-log10(J0+t*(J1-J0)))*y_scale;
	fl_line(round(prevx+deltax),round(prevy+deltay),round(curx+deltax),round(cury+deltay));
	prevx=curx; prevy=cury;
      }
      return;
    }
    // logx
    double prevx=i0,prevy=j0,curx,cury;
    double I0=pow10(i0/x_scale+window_xmin),I1=pow10(i1/x_scale+window_xmin);
    for (int i=1;i<logplot_points;++i){
      double t=double(i)/logplot_points;
      curx=(log10(I0+t*(I1-I0))-window_xmin)*x_scale;
      cury=j0+t*(j1-j0);
      fl_line(round(prevx+deltax),round(prevy+deltay),round(curx+deltax),round(cury+deltay));
      prevx=curx; prevy=cury;
    }
  }

  void draw_legende(const vecteur & f,int i0,int j0,int labelpos,const Graph2d * iptr,int clip_x,int clip_y,int clip_w,int clip_h,int deltax,int deltay){
    if (f.empty() ||!iptr->show_names )
      return;
    string legendes;
    const context * contextptr = get_context(iptr);
    if (f[0].is_symb_of_sommet(at_curve)){
      gen & f0=f[0]._SYMBptr->feuille;
      if (f0.type==_VECT && !f0._VECTptr->empty()){
	gen & f1 = f0._VECTptr->front();
	if (f1.type==_VECT && f1._VECTptr->size()>4 && (!is_zero((*f1._VECTptr)[4]) || (iptr->show_names & 2)) ){
	  gen legende=f1._VECTptr->front();
	  gen var=(*f1._VECTptr)[1];
	  gen r=re(legende,contextptr),i=im(legende,contextptr),a,b;
	  if (var.type==_IDNT && is_linear_wrt(r,*var._IDNTptr,a,b,contextptr)){
	    i=subst(i,var,(var-b)/a,false,contextptr);
	    legendes=i.print(contextptr);
	  }
	  else
	    legendes=r.print(contextptr)+","+i.print(contextptr);
	  if (legendes.size()>18){
	    if (legendes.size()>30)
	      legendes="";
	    else
	      legendes=legendes.substr(0,16)+"...";
	  }
	}
      }
    }
    if (f.size()>2)
      legendes=gen2string(f[2])+(legendes.empty()?"":":")+legendes;
    if (legendes.empty())
      return;
    if (abs_calc_mode(contextptr)==38 && legendes.size()>1 && legendes[0]=='G')
      legendes=legendes.substr(1,legendes.size()-1);
    fl_font(cst_greek_translate(legendes),iptr->labelsize());
    int dx=3,dy=1;
    find_dxdy(legendes,labelpos,iptr->labelsize(),dx,dy);
    check_fl_draw(legendes.c_str(),iptr->x()+i0+dx,iptr->y()+j0+dy,clip_x,clip_y,clip_w,clip_h,deltax,deltay);
  }

  void petite_fleche(double i1,double j1,double dx,double dy,int deltax,int deltay,int width){
    double dxy=std::sqrt(dx*dx+dy*dy);
    if (dxy){
      dxy/=min(5,int(dxy/10))+width;
      dx/=dxy;
      dy/=dxy;
      double dxp=-dy,dyp=dx;
      dx*=std::sqrt(3.0);
      dy*=sqrt(3.0);
      fl_polygon(round(i1)+deltax,round(j1)+deltay,round(i1+dx+dxp)+deltax,round(j1+dy+dyp)+deltay,round(i1+dx-dxp)+deltax,round(j1+dy-dyp)+deltay);
    }
  }

  // helper for Graph2d::draw method
  // plot_i is the position we are drawing in plot_instructions
  // f is the vector of arguments and s is the function: we draw s(f)
  // legende_x is the x position in pixels of the parameter subwindow
  // parameter_y is the current y position for parameter drawing
  // x_scale and y_scale are the current scales
  // horizontal_pixels and vertical_pixels the size of the window
  void fltk_draw(Graph2d & Mon_image,int plot_i,const gen & g,double x_scale,double y_scale,int clip_x,int clip_y,int clip_w,int clip_h){
    int deltax=Mon_image.x(),deltay=Mon_image.y();
    History_Pack * hp =get_history_pack(&Mon_image);
    context * contextptr=hp?hp->contextptr:get_context(hp);
    if (g.type==_VECT){
      vecteur & v =*g._VECTptr;
      const_iterateur it=v.begin(),itend=v.end();
      for (;it!=itend;++it){
	fltk_draw(Mon_image,plot_i,*it,x_scale,y_scale,clip_x,clip_y,clip_w,clip_h);
	fl_line_style(0); // back to default line style
      } // end for it
    }
    if (g.type!=_SYMB)
      return;
    unary_function_ptr s=g._SYMBptr->sommet;
    if (s==at_animation){
      fltk_draw(Mon_image,plot_i,get_animation_pnt(g,Mon_image.animation_instructions_pos),x_scale,y_scale,clip_x,clip_y,clip_w,clip_h);
      return;
    }
    if (g._SYMBptr->feuille.type!=_VECT)
      return;
    vecteur f=*g._SYMBptr->feuille._VECTptr;
    int mxw=Mon_image.w(),myw=Mon_image.h()-(Mon_image.show_axes?(Mon_image.title.empty()?1:2):0)*Mon_image.labelsize();
    double i0,j0,i0save,j0save,i1,j1;
    int fs=f.size();
    if ((fs==4) && (s==at_parameter)){
      return ;
    }
    string the_legend;
    vecteur style(get_style(f,the_legend));
    int styles=style.size();
    // color
    int ensemble_attributs = style.front().val;
    bool hidden_name = false;
    if (style.front().type==_ZINT){
      ensemble_attributs = mpz_get_si(*style.front()._ZINTptr);
      hidden_name=true;
    }
    else
      hidden_name=ensemble_attributs<0;
    int width           =(ensemble_attributs & 0x00070000) >> 16; // 3 bits
    int epaisseur_point =(ensemble_attributs & 0x00380000) >> 19; // 3 bits
    int type_line       =(ensemble_attributs & 0x01c00000) >> 22; // 3 bits
    if (type_line>4)
      type_line=(type_line-4)<<8;
    int type_point      =(ensemble_attributs & 0x0e000000) >> 25; // 3 bits
    int labelpos        =(ensemble_attributs & 0x30000000) >> 28; // 2 bits
    bool fill_polygon   =(ensemble_attributs & 0x40000000) >> 30;
    int couleur         =(ensemble_attributs & 0x0000ffff);
    epaisseur_point += 2;
    std::pair<Fl_Image *,Fl_Image *> * texture = 0;
    for (int i=2;i<styles;++i){
      gen & attr=style[i];
      if (attr.type==_VECT && attr._VECTptr->size()<=3 ){
	gen attrv0=attr._VECTptr->front();
	if (attrv0.type==_INT_ && attrv0.val==_GL_TEXTURE){
	  gen attrv1=(*attr._VECTptr)[1];
	  if (attrv1.type==_VECT && attrv1._VECTptr->size()==2 && attrv1._VECTptr->front().type==_STRNG && is_undef(attrv1._VECTptr->back())){
	    // reload cached image
	    attrv1=attrv1._VECTptr->front();
	    std::map<std::string,std::pair<Fl_Image *,Fl_Image*> *>::iterator it,itend=texture2d_cache.end();
	    it=texture2d_cache.find(attrv1._STRNGptr->c_str());
	    if (it!=itend){
	      std::pair<Fl_Image *,Fl_Image*> * old= it->second;
	      delete old;
	      texture2d_cache.erase(it);
	    }
	  }
	  if (attrv1.type==_STRNG){
	    get_texture2d(*attrv1._STRNGptr,texture);
	  }
	  // set texture
	  continue;
	} // end attrv0 = gl_texture
      }
    }
    if (s==at_pnt){ 
      // f[0]=complex pnt or vector of complex pnts or symbolic
      // f[1] -> style 
      // f[2] optional=label
      gen point=f[0];
      if (point.type==_VECT && point.subtype==_POINT__VECT)
	return;
      if ( (f[0].type==_SYMB) && (f[0]._SYMBptr->sommet==at_curve) && (f[0]._SYMBptr->feuille.type==_VECT) && (f[0]._SYMBptr->feuille._VECTptr->size()) ){
	// Mon_image.show_mouse_on_object=false;
	point=f[0]._SYMBptr->feuille._VECTptr->back();
      }
      if (is_undef(point))
	return;
      if ( equalposcomp(Mon_image.selected,plot_i))
	fl_color(FL_BLUE);
      else
	xcas_color(couleur);
      fl_line_style(type_line,width+1,0); 
      if (point.type==_SYMB) {
	if (point._SYMBptr->sommet==at_hyperplan || point._SYMBptr->sommet==at_hypersphere)
	  return;
	/* fl_curve is a filled curved, not a Bezier curve!
	if (point._SYMBptr->sommet==at_Bezier && point._SYMBptr->feuille.type==_VECT && point._SYMBptr->feuille._VECTptr->size()==4){
	  vecteur v=*point._SYMBptr->feuille._VECTptr;
	  double i2,j2,i3,j3; gen e,f0,f1;
	  evalfdouble2reim(v[0],e,f0,f1,contextptr);
	  i0=f0._DOUBLE_val; j0=f1._DOUBLE_val;
	  evalfdouble2reim(v[1],e,f0,f1,contextptr);
	  i1=f0._DOUBLE_val; j1=f1._DOUBLE_val;
	  evalfdouble2reim(v[2],e,f0,f1,contextptr);
	  i2=f0._DOUBLE_val; j2=f1._DOUBLE_val;
	  evalfdouble2reim(v[3],e,f0,f1,contextptr);
	  i3=f0._DOUBLE_val; j3=f1._DOUBLE_val;
	  fl_push_matrix();
	  fl_translate(deltax,deltay);
	  fl_mult_matrix(x_scale,0,0,-y_scale,0,0);
	  fl_translate(-Mon_image.window_xmin,-Mon_image.window_ymax);
	  fl_curve(i0,j0,i1,j1,i2,j2,i3,j3);
	  fl_end_complex_polygon();
	  fl_pop_matrix(); // Restore initial matrix
	  return;
	}
	*/
	if (point._SYMBptr->sommet==at_cercle && !(Mon_image.display_mode & 0xc00)){
	  vecteur v=*point._SYMBptr->feuille._VECTptr;
	  gen diametre=remove_at_pnt(v[0]);
	  gen e1=diametre._VECTptr->front().evalf_double(1,contextptr),e2=diametre._VECTptr->back().evalf_double(1,contextptr);
	  gen centre=rdiv(e1+e2,2.0,contextptr);
	  gen e12=e2-e1;
	  double ex=evalf_double(re(e12,contextptr),1,contextptr)._DOUBLE_val,ey=evalf_double(im(e12,contextptr),1,contextptr)._DOUBLE_val;
	  if (!Mon_image.findij(centre,x_scale,y_scale,i0,j0,contextptr))
	    return;
	  gen diam=std::sqrt(ex*ex+ey*ey);
	  gen angle=std::atan2(ey,ex);
	  gen a1=v[1].evalf_double(1,contextptr),a2=v[2].evalf_double(1,contextptr);
	  if ( (diam.type==_DOUBLE_) && (a1.type==_DOUBLE_) && (a2.type==_DOUBLE_) ){
	    i1=diam._DOUBLE_val*x_scale/2.0;
	    j1=diam._DOUBLE_val*y_scale/2.0;
	    double a1d=a1._DOUBLE_val,a2d=a2._DOUBLE_val,angled=angle._DOUBLE_val;
	    bool changer_sens=a1d>a2d;
	    if (changer_sens){
	      double tmp=a1d;
	      a1d=a2d;
	      a2d=tmp;
	    }
	    if (fill_polygon){
	      if (v[1]==0 && v[2]==cst_two_pi)
		fl_pie(deltax+round(i0-i1),deltay+round(j0-j1),round(2*i1),round(2*j1),0,360);
	      else
		fl_pie(deltax+round(i0-i1),deltay+round(j0-j1),round(2*i1),round(2*j1),(angled+a1d)*180/M_PI,(angled+a2d)*180/M_PI);
	    }
	    else {
	      double anglei=(angled+a1d),anglef=(angled+a2d),anglem=(anglei+anglef)/2;
	      fl_arc(deltax+round(i0-i1),deltay+round(j0-j1),round(2*i1),round(2*j1),anglei*180/M_PI,anglef*180/M_PI);
	      if (v.size()>=4){ // if cercle has the optionnal 5th arg
		if (v[3]==2)
		  petite_fleche(i0+i1*std::cos(anglem),j0-j1*std::sin(anglem),-i1*std::sin(anglem),-j1*std::cos(anglem),deltax,deltay,width);
		else {
		  if (changer_sens)
		    petite_fleche(i0+i1*std::cos(anglei),j0-j1*std::sin(anglei),-i1*std::sin(anglei),-j1*std::cos(anglei),deltax,deltay,width);
		  else
		    petite_fleche(i0+i1*std::cos(anglef),j0-j1*std::sin(anglef),i1*std::sin(anglef),j1*std::cos(anglef),deltax,deltay,width);
		}
	      }
	    }
	    // Label a few degrees from the start angle, 
	    // FIXME should use labelpos
	    double anglel=angled+a1d+0.3;
	    if (v.size()>=4 && v[3]==2)
	      anglel=angled+(0.45*a1d+0.55*a2d);
	    i0=i0+i1*std::cos(anglel); 
	    j0=j0-j1*std::sin(anglel);
	    if (!hidden_name)
	      draw_legende(f,round(i0),round(j0),labelpos,&Mon_image,clip_x,clip_y,clip_w,clip_h,0,0);
	    return;
	  }
	} // end circle
	if (point._SYMBptr->sommet==at_pixon){
	  // pixon (i,j,color)
	  if (point._SYMBptr->feuille.type!=_VECT)
	    return;
	  vecteur &v=*point._SYMBptr->feuille._VECTptr;
	  if (v.size()<3 || v[0].type!=_INT_ || v[1].type!=_INT_ || v[2].type!=_INT_)
	    return;
	  int delta_i=v[0].val,delta_j=v[1].val;
	  xcas_color(v[2].val);
#ifdef IPAQ
	  if (delta_i>0 && delta_i<mxw && delta_j>0 && delta_j<myw)
	    check_fl_point(deltax+delta_i,deltay+delta_j,clip_x,clip_y,clip_w,clip_h,0,0);
#else
	  delta_i *= 2;
	  delta_j *= 2;
	  if (delta_i>0 && delta_i<mxw && delta_j>0 && delta_j<myw){
	    check_fl_point(deltax+delta_i,deltay+delta_j,clip_x,clip_y,clip_w,clip_h,0,0);
	    check_fl_point(deltax+delta_i,deltay+delta_j+1,clip_x,clip_y,clip_w,clip_h,0,0);
	    check_fl_point(deltax+delta_i+1,deltay+delta_j,clip_x,clip_y,clip_w,clip_h,0,0);
	    check_fl_point(deltax+delta_i+1,deltay+delta_j+1,clip_x,clip_y,clip_w,clip_h,0,0);
	  }
#endif
	  return;
	}
	if (point._SYMBptr->sommet==at_bitmap){
	  // bitmap(vector of int (1 per line)), 1st line, 1st col, [type]
	  if (point._SYMBptr->feuille.type!=_VECT)
	    return;
	  vecteur &v=*point._SYMBptr->feuille._VECTptr;
	  if (v.size()<3 || v[0].type!=_VECT || v[1].type!=_INT_ || v[2].type!=_INT_ )
	    return;
	  int delta_i=v[1].val,delta_j=v[2].val;
	  double xmin=Mon_image.window_xmin,ymin=Mon_image.window_ymin,xmax=Mon_image.window_xmax,ymax=Mon_image.window_ymax;
	  //gen psize=_Pictsize(0);
	  int bitmap_w=Mon_image.w()-int(Mon_image.ylegende*(Mon_image.show_axes?Mon_image.labelsize():0)),
	    bitmap_h=Mon_image.h()-(Mon_image.show_axes?((Mon_image.title.empty()?1:2)*Mon_image.labelsize()):0);
	  if (v.size()>8){
	    xmin=v[3]._DOUBLE_val;
	    xmax=v[4]._DOUBLE_val;
	    ymin=v[5]._DOUBLE_val;
	    ymax=v[6]._DOUBLE_val;
	    bitmap_w=v[7].val;
	    bitmap_h=v[8].val;
	  }
	  double bitmap_scalex=(xmax-xmin)/bitmap_w,scalex=(Mon_image.window_xmax-Mon_image.window_xmin)/(Mon_image.w()-int(Mon_image.ylegende*(Mon_image.show_axes?Mon_image.labelsize():0)));
	  double bitmap_scaley=(ymax-ymin)/bitmap_h,scaley=(Mon_image.window_ymax-Mon_image.window_ymin)/(Mon_image.h()-(Mon_image.show_axes?((Mon_image.title.empty()?1:2)*Mon_image.labelsize()):0));
	  double X,Y;
	  int ii,jj;
	  const_iterateur it=v[0]._VECTptr->begin(),itend=v[0]._VECTptr->end();
	  for (;it!=itend;++it,++delta_i){
	    if (it->type!=_INT_ && it->type!=_ZINT)
	      continue;
	    gen z=*it;
	    mpz_t zz,zr;
	    if (it->type==_INT_)
	      mpz_init_set_ui(zz,it->val);
	    else
	      mpz_init_set(zz,*it->_ZINTptr);
	    mpz_init(zr);
	    for (int j=delta_j;mpz_sgn(zz);++j){
	      mpz_tdiv_r_2exp (zr, zz, 1);
	      mpz_tdiv_q_2exp (zz, zz, 1);
	      if (mpz_sgn(zr)){
		X=xmin+j*bitmap_scalex;
		ii=int(0.5+(X-Mon_image.window_xmin)/scalex);
		Y=ymax-delta_i*bitmap_scaley;
		jj=int(0.5+(Mon_image.window_ymax-Y)/scaley);
		if (ii>0 && ii<mxw && jj>0 && jj<myw)
		  check_fl_point(deltax+ii,deltay+jj,clip_x,clip_y,clip_w,clip_h,0,0);
	      }
	    }
	    mpz_clear(zr);
	    mpz_clear(zz);
	  }
	  return;
	} // end bitmap
	if (point._SYMBptr->sommet==at_legende){
	  gen & f=point._SYMBptr->feuille;
	  if (f.type==_VECT && f._VECTptr->size()==3){
	    vecteur & fv=*f._VECTptr;
	    if (fv[0].type==_VECT && fv[0]._VECTptr->size()>=2 && fv[1].type==_STRNG && fv[2].type==_INT_){
	      vecteur & fvv=*fv[0]._VECTptr;
	      if (fvv[0].type==_INT_ && fvv[1].type==_INT_){
		fl_font(FL_HELVETICA,Mon_image.labelsize());
		xcas_color(fv[2].val);
		int dx=0,dy=0;
		string legendes(*fv[1]._STRNGptr);
		find_dxdy(legendes,labelpos,Mon_image.labelsize(),dx,dy);
		fl_draw(legendes.c_str(),deltax+fvv[0].val+dx,deltay+fvv[1].val+dy);
	      }
	    }
	  }
	}
      } // end point.type==_SYMB
      if (point.type!=_VECT || (point.type==_VECT && (point.subtype==_GROUP__VECT || point.subtype==_VECTOR__VECT) && point._VECTptr->size()==2 && is_zero(point._VECTptr->back()-point._VECTptr->front())) ){ // single point
	if (!Mon_image.findij((point.type==_VECT?point._VECTptr->front():point),x_scale,y_scale,i0,j0,contextptr))
	  return;
	if (i0>0 && i0<mxw && j0>0 && j0<myw)
	  fltk_point(deltax,deltay,round(i0),round(j0),epaisseur_point,type_point);
	if (!hidden_name)
	  draw_legende(f,round(i0),round(j0),labelpos,&Mon_image,clip_x,clip_y,clip_w,clip_h,0,0);
	return;
      }
      // path
      const_iterateur jt=point._VECTptr->begin(),jtend=point._VECTptr->end();
      if (jt==jtend)
	return;
      bool logx=Mon_image.display_mode & 0x400,logy=Mon_image.display_mode & 0x800;
      if (jtend-jt==2 && logx && point.subtype==_LINE__VECT){
	// find points with + coordinates ax+cx*t=x, ay+cy*t=y=
	gen a=*jt,ax=re(a,contextptr),ay=im(a,contextptr),b=*(jt+1),c=b-a,cx=re(c,contextptr),cy=im(c,contextptr);
	if (!is_zero(cx)){
	  // y=ay+cy/cx*(x-ax)
	  gen x=pow(10,Mon_image.window_xmin,contextptr);
	  gen y=ay+cy/cx*(x-ax);
	  Mon_image.findij(x+cst_i*y,x_scale,y_scale,i0,j0,contextptr);
	  for (int i=1;i<=logplot_points;++i){
	    x=pow(10,Mon_image.window_xmin+i*(Mon_image.window_xmax-Mon_image.window_xmin)/logplot_points,contextptr);
	    y=ay+cy/cx*(x-ax);
	    Mon_image.findij(x+cst_i*y,x_scale,y_scale,i1,j1,contextptr);
	    fl_line(round(i0+deltax),round(j0+deltay),round(i1+deltax),round(j1+deltay));
	    i0=i1;
	    j0=j1;
	  }
	  return;
	}
      }
      if (texture && jtend-jt>2){
	// use *jt and *(jt+2) for the rectangle texture
	Mon_image.findij(*jt,x_scale,y_scale,i0,j0,contextptr);
	if (!Mon_image.findij(*(jt+2),x_scale,y_scale,i1,j1,contextptr))
	  return;
	if (i0>i1)
	  std::swap(i0,i1);
	if (j0>j1)
	  std::swap(j0,j1);
	int tx=int(i0+.5)+deltax;
	int tw=int(i1-i0+.5);
	int ty=int(j0+.5)+deltay;
	int th=int(j1-j0+.5);
	if (texture->second && texture->second->w()==tw && texture->second->h()==th)
	  texture->second->draw(tx,ty,tw,th);
	else {
	  if (texture->second)
	    delete texture->second;
	  if (texture->first){
	    texture->second=texture->first->copy(tw,th);
	    texture->second->draw(tx,ty,tw,th);
	  }
	}
	return;
      }
      if (jt->type==_VECT)
	return;
      if ( (type_point || epaisseur_point>2) && type_line==0 && width==0){
	for (;jt!=jtend;++jt){
	  if (!Mon_image.findij(*jt,x_scale,y_scale,i0,j0,contextptr))
	    return;
	  if (i0>0 && i0<mxw && j0>0 && j0<myw)
	    fltk_point(deltax,deltay,round(i0),round(j0),epaisseur_point,type_point);
	}
	if (!hidden_name)
	  draw_legende(f,round(i0),round(j0),labelpos,&Mon_image,clip_x,clip_y,clip_w,clip_h,0,0);
	return;
      }
      // initial point
      if (!Mon_image.findij(*jt,x_scale,y_scale,i0,j0,contextptr))
	return;
      if (fill_polygon && *jt==*(jtend-1)){
	const_iterateur jtsave=jt;
	gen e,f0,f1;
	// Compute matrix for complex drawing
	fl_push_matrix();
	fl_translate(deltax,deltay);
	fl_mult_matrix(x_scale,0,0,-y_scale,0,0);
	fl_translate(-Mon_image.window_xmin,-Mon_image.window_ymax);
	fl_begin_complex_polygon();
	for (;jt!=jtend;++jt){
	  evalfdouble2reim(*jt,e,f0,f1,contextptr);
	  if ((f0.type==_DOUBLE_) && (f1.type==_DOUBLE_))
	    fl_vertex(f0._DOUBLE_val,f1._DOUBLE_val);
	}
	fl_end_complex_polygon();
	fl_pop_matrix(); // Restore initial matrix
	if (!width){
	  if (!hidden_name)
	    draw_legende(f,round(i0),round(j0),labelpos,&Mon_image,clip_x,clip_y,clip_w,clip_h,0,0);
	  return;
	}
	jt=jtsave;
	fl_line_style(type_line,width,0); 
	fl_color(epaisseur_point-2+(type_point<<3));
      }
      i0save=i0;
      j0save=j0;
      ++jt;
      if (jt==jtend){
	if (i0>0 && i0<mxw && j0>0 && j0<myw)
	  check_fl_point(deltax+round(i0),deltay+round(j0),clip_x,clip_y,clip_w,clip_h,0,0);
	if (!hidden_name)
	  draw_legende(f,round(i0),round(j0),labelpos,&Mon_image,clip_x,clip_y,clip_w,clip_h,0,0);
	return;
      }
      bool seghalfline=( point.subtype==_LINE__VECT || point.subtype==_HALFLINE__VECT ) && (point._VECTptr->size()==2);
      // rest of the path
      for (;;){
	if (!Mon_image.findij(*jt,x_scale,y_scale,i1,j1,contextptr))
	  return;
	if (!seghalfline){
	  checklog_fl_line(i0,j0,i1,j1,deltax,deltay,logx,logy,Mon_image.window_xmin,x_scale,Mon_image.window_ymax,y_scale);
	  if (point.subtype==_VECTOR__VECT){
	    double dx=i0-i1,dy=j0-j1;
	    petite_fleche(i1,j1,dx,dy,deltax,deltay,width);
	  }
	}
	++jt;
	if (jt==jtend){ // label of line at midpoint
	  if (point.subtype==_LINE__VECT){
	    i0=(6*i1-i0)/5-8;
	    j0=(6*j1-j0)/5-8;
	  }
	  else {
	    i0=(i0+i1)/2-8;
	    j0=(j0+j1)/2;
	  }
	  break;
	}
	i0=i1;
	j0=j1;
      }
      // check for a segment/halfline/line
      if ( seghalfline){
	double deltai=i1-i0save,adeltai=std::abs(deltai);
	double deltaj=j1-j0save,adeltaj=std::abs(deltaj);
	if (point.subtype==_LINE__VECT){
	  if (deltai==0)
	    checklog_fl_line(i1,0,i1,clip_h,deltax,deltay,logx,logy,Mon_image.window_xmin,x_scale,Mon_image.window_ymax,y_scale);
	  else {
	    if (deltaj==0)
	      checklog_fl_line(0,j1,clip_w,j1,deltax,deltay,Mon_image.display_mode & 0x400,Mon_image.display_mode & 0x800,Mon_image.window_xmin,x_scale,Mon_image.window_ymax,y_scale);
	    else {
	      // Find the intersections with the 4 rectangle segments
	      // Horizontal x=0 or w =i1+t*deltai: y=j1+t*deltaj
	      vector< complex<double> > pts;
	      double y0=j1-i1/deltai*deltaj;
	      if (y0>=0 && y0<=clip_h)
		pts.push_back(complex<double>(0.0,y0));
	      double yw=j1+(clip_w-i1)/deltai*deltaj;
	      if (yw>=0 && yw<=clip_h)
		pts.push_back(complex<double>(clip_w,yw));
	      // Vertical y=0 or h=j1+t*deltaj, x=i1+t*deltai
	      double x0=i1-j1/deltaj*deltai;
	      if (x0>0 && x0<=clip_w)
		pts.push_back(complex<double>(x0,0.0));
	      double xh=i1+(clip_h-j1)/deltaj*deltai;
	      if (xh>=0 && xh<=clip_w)
		pts.push_back(complex<double>(xh,clip_h));
	      if (pts.size()>=2)
		checklog_fl_line(pts[0].real(),pts[0].imag(),pts[1].real(),pts[1].imag(),deltax,deltay,Mon_image.display_mode & 0x400,Mon_image.display_mode & 0x800,Mon_image.window_xmin,x_scale,Mon_image.window_ymax,y_scale);
	    } // end else adeltai==0 , adeltaj==0
	  } // end else adeltai==0
	} // end LINE_VECT
	else {
	  double N=1;
	  if (adeltai){
	    N=clip_w/adeltai+1;
	    if (adeltaj)
	      N=max(N,clip_h/adeltaj+1);
	  }
	  else {
	    if (adeltaj)
	      N=clip_h/adeltaj+1;
	  }
	  N *= 2; // increase N since rounding might introduce too small clipping
	  while (fabs(N*deltai)>10000)
	    N /= 2;
	  while (fabs(N*deltaj)>10000)
	    N /= 2;
	  checklog_fl_line(i0save,j0save,i1+N*deltai,j1+N*deltaj,deltax,deltay,Mon_image.display_mode & 0x400,Mon_image.display_mode & 0x800,Mon_image.window_xmin,x_scale,Mon_image.window_ymax,y_scale);
	}
      } // end seghalfline
      if ( (point.subtype==_GROUP__VECT) && (point._VECTptr->size()==2))
	; // no legend for segment
      else {
	if (!hidden_name)
	  draw_legende(f,round(i0),round(j0),labelpos,&Mon_image,clip_x,clip_y,clip_w,clip_h,0,0);
      }
    } // end pnt subcase
    
  }

