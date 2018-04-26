# -*- mode:Makefile -*-
# giac in javascript with emcc
# -s EXPORTED_FUNCTIONS="['_caseval']": export function
# then from javascript console
# caseval = Module.cwrap('caseval', 'string', ['string']) type number (inc. pointer) or string
# caseval = Module.cwrap('_ZN4giac7casevalEPKc', 'string', ['string'])
# Module['noExitRuntime']=true
# --closure 1 run closure compiler, this is long
# --pre-js: add javascript code
# --compression lzma/lzma-native,lzma/lzma-decoder.js,LZMA.decompress
# lzma/lzma-native < file_to_compress > compressed_file
# -s options are in emscripten src/settings.js
#PREC = -s PRECISE_I32_MUL=1 #-DGIAC_GGB 
PREC = -s PRECISE_I64_MATH=1 #-DGIAC_GGB 
CXXFLAGS =  $(PREC) -I. -I.. -DHAVE_CONFIG_H -DIN_GIAC -DGIAC_GENERIC_CONSTANTS -DNO_STDEXCEPT -Oz -s ALLOW_MEMORY_GROWTH=1 -s ASSERTIONS=1 -s GL_UNSAFE_OPTS=0  # --bind -DEMCC_BIND # -g2 -s SAFE_HEAP=1 # -DEMCC_GLUT -s USE_SDL=2 # --llvm-opts 2 -v
CFLAGS =  $(PREC) -I. -I.. -v  # -pg
CPPFLAGS = $(CXXFLAGS)
LDJS2FLAGS = $(PREC)  -v -s EXPORTED_FUNCTIONS="['__ZN4giac7casevalEPKc','__ZN4giac13giac_rendererEPKc']" -s ALLOW_MEMORY_GROWTH=1 -s LEGACY_GL_EMULATION=1 -s ASSERTIONS=1 -s GL_UNSAFE_OPTS=0 --memory-init-file 0 #-s USE_SDL=2 #--closure 1
LDJSFLAGS = $(PREC) -O1 -v -s EXPORTED_FUNCTIONS="['__ZN4giac7casevalEPKc','__ZN4giac13giac_rendererEPKc','caseval']" -s ALLOW_MEMORY_GROWTH=1 -s LEGACY_GL_EMULATION=1 -s ASSERTIONS=1 -s GL_UNSAFE_OPTS=0 --memory-init-file 0 # --bind -DEMCC_BIND   # -g2 -s SAFE_HEAP=1  # -s USE_SDL=2 # -s USE_TYPED_ARRAYS=0 # -s TOTAL_MEMORY=24000000 
CXX=emcc
CC=emcc
GIACOBJS = sym2poly.o gausspol.o threaded.o maple.o ti89.o mathml.o moyal.o misc.o permu.o quater.o desolve.o input_parser.o symbolic.o index.o modpoly.o modfactor.o ezgcd.o derive.o solve.o intg.o intgab.o risch.o lin.o series.o subst.o vecteur.o sparse.o csturm.o tex.o global.o ifactor.o alg_ext.o gauss.o isom.o help.o plot.o plot3d.o rpn.o prog.o pari.o cocoa.o TmpLESystemSolver.o TmpFGLM.o unary.o usual.o identificateur.o gen.o input_lexer.o tinymt32.o opengl.o
LIBS = libpari.a libmpfi.a libmpfr.a libgmp.a --js-library time.js
giac.js:	$(GIACOBJS) 
	$(CXX) $(LDJSFLAGS) $(GIACOBJS) $(LIBS) -o giac.js # --preload-file doc/fr/keywords #  --closure 1 >& log
#	/bin/cp giacggb.js giac.js
#	/bin/rm giac.js.lz
#	lzma/lzma-native < giacggb.js > giacggb.js.compress
#	cp giacggb.js /shared/ggbjs/giac.js
#	scp giac.js giac.js.mem malherbe.ujf-grenoble.fr:public_html/giac
#	lzip giac.js
#	scp giac.js.lz malherbe:public_html/giac
#	lzip -d giac.js.lz
giac:	$(GIACOBJS) 
	$(CXX) $(LDJS2FLAGS) $(GIACOBJS) $(LIBS) -o giac.js #--closure 1 
#	/bin/cp giacggb.js giac.js
#	/bin/rm giac.js.lz
#	lzma/lzma-native < giacggb.js > giacggb.js.compress
#	cp giacggb.js /shared/ggbjs/giac.js
#	scp giac.js malherbe:public_html/giac
#	lzip giac.js
#	scp giac.js.lz malherbe:public_html/giac
#	lzip -d giac.js.lz
icas.html:	$(GIACOBJS) icas.o
	$(CXX) $(LDJSFLAGS) $(GIACOBJS) $(LIBS) icas.o -o icas.html
clean:
	rm -f $(GIACOBJS) test.o icas.o
.cc.o:
	$(CXX) $(CXXFLAGS) -c $<
.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

input_lexer.o: input_lexer.cc
	$(CXX) $(CXXFLAGSNOASM) -c input_lexer.cc
