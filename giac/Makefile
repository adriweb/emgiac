# -*- mode:Makefile -*-
# emscripten 1.37: add -s ASSERTIONS=1 in CXXFLAGS/LDJSFLAGS, -s WASM=1 -s EXTRA_EXPORTED_RUNTIME_METHODS='["ccall","cwrap"]' -s "BINARYEN_TRAP_MODE='clamp'" in LDJSFLAGS, comment html5.h in opengl.cc
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
CXXFLAGS = $(PREC)  -I. -I.. -DHAVE_CONFIG_H -DIN_GIAC -DGIAC_GENERIC_CONSTANTS -Oz -s ALLOW_MEMORY_GROWTH=1 -s GL_UNSAFE_OPTS=0 #-s ASSERTIONS=1 -fno-exceptions -DNO_STDEXCEPT # --bind -DEMCC_BIND --llvm-opts 2 -v
CXXFLAGSNOASM = $(PREC) -fno-exceptions -DNO_STDEXCEPT -I. -I.. -DHAVE_CONFIG_H -DIN_GIAC -DGIAC_GENERIC_CONSTANTS -O1 -s ALLOW_MEMORY_GROWTH=1 -s GL_UNSAFE_OPTS=0 
CFLAGS =  $(PREC) -I. -I.. -v  # -pg
CPPFLAGS = $(CXXFLAGS)
LDJS1FLAGS = $(PREC) -O1 -v -s EXPORTED_FUNCTIONS="['_caseval','__ZN4giac13giac_rendererEPKc']" -s ALLOW_MEMORY_GROWTH=1 -s LEGACY_GL_EMULATION=1 -s GL_UNSAFE_OPTS=0 --memory-init-file 0 #-s DISABLE_EXCEPTION_CATCHING=0 #--closure 1
LDJS2FLAGS = $(PREC) -O2 -v -s EXPORTED_FUNCTIONS="['_caseval','__ZN4giac13giac_rendererEPKc']" -s ALLOW_MEMORY_GROWTH=1 -s LEGACY_GL_EMULATION=1 -s GL_UNSAFE_OPTS=0 --memory-init-file 0 #--closure 1
# enabling DISABLE_EXCEPTION_CATCHING=0 compiles forever
LDJSFLAGS = $(PREC) -O3 -v -s EXPORTED_FUNCTIONS="['_caseval','__ZN4giac13giac_rendererEPKc']" -s TOTAL_MEMORY=67108864 -s LEGACY_GL_EMULATION=1 -s GL_UNSAFE_OPTS=0 --memory-init-file 0 #-s EXTRA_EXPORTED_RUNTIME_METHODS='["ccall","cwrap"]' -s DISABLE_EXCEPTION_CATCHING=0 # --bind -DEMCC_BIND -s USE_TYPED_ARRAYS=0 #  -s ALLOW_MEMORY_GROWTH=1
CXX=emcc
CC=emcc
GIACOBJS = sym2poly.o gausspol.o threaded.o maple.o ti89.o lpsolve.o optimization.o mathml.o moyal.o misc.o permu.o quater.o desolve.o input_parser.o symbolic.o index.o modpoly.o modfactor.o ezgcd.o derive.o solve.o intg.o intgab.o risch.o lin.o series.o subst.o vecteur.o sparse.o csturm.o tex.o global.o ifactor.o alg_ext.o gauss.o isom.o help.o plot.o plot3d.o rpn.o prog.o pari.o cocoa.o lpsolve.o TmpLESystemSolver.o TmpFGLM.o unary.o usual.o identificateur.o gen.o input_lexer.o tinymt32.o opengl.o
LIBS = libpari.a libmpfi.a libmpfr.a libgmp.a libglpk.a --js-library time.js
giac.js:	$(GIACOBJS) 
	$(CXX) $(LDJSFLAGS) $(GIACOBJS) $(LIBS) -o giac.js #--preload-file doc/fr/keywords #  --closure 1 >& log
#	/bin/cp giacggb.js giac.js
#	/bin/rm giac.js.lz
#	lzma/lzma-native < giacggb.js > giacggb.js.compress
#	cp giacggb.js /shared/ggbjs/giac.js
#	scp giac.js malherbe.ujf-grenoble.fr:public_html/giac
#	scp giac.js.mem malherbe.ujf-grenoble.fr:public_html/giac
#	lzip giac.js
#	scp giac.js.lz malherbe.ujf-grenoble.fr:public_html/giac
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
libgiac.a:
	emar cru libgiac.a $(GIACOBJS)
icas.html:	$(GIACOBJS) icas.o
	$(CXX) $(LDJSFLAGS) $(GIACOBJS) $(LIBS) icas.o -o icas.html
test.html:	libgiac.a test.cc
	$(CXX) $(PREC) $(GIACOBJS) $(LIBS) test.cc -o test.html
clean:
	rm -f $(GIACOBJS) test.o icas.o
input_lexer.o: input_lexer.cc
	$(CXX) $(CXXFLAGSNOASM) -c input_lexer.cc
input_parser.o: input_parser.cc
	$(CXX) $(CXXFLAGSNOASM) -c input_parser.cc
.cc.o:
	$(CXX) $(CXXFLAGS) -c $<
.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

