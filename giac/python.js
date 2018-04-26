// CodeMirror, copyright (c) by Marijn Haverbeke and others
// Distributed under an MIT license: http://codemirror.net/LICENSE

(function(mod) {
  if (typeof exports == "object" && typeof module == "object") // CommonJS
    mod(require("../../lib/codemirror"));
  else if (typeof define == "function" && define.amd) // AMD
    define(["../../lib/codemirror"], mod);
  else // Plain browser env
    mod(CodeMirror);
})(function(CodeMirror) {
  "use strict";

  function wordRegexp(words) {
    return new RegExp("^((" + words.join(")|(") + "))\\b");
  }

  var wordOperators = wordRegexp(["and", "or", "not", "is"]);
  var commonKeywords = ["as", "assert", "break", "class", "continue",
                        "def", "del", "elif", "else", "except", "finally",
                        "for", "from", "global", "if", "import",
                        "lambda", "pass", "raise", "return",
                        "try", "while", "with", "yield", "in"];
  var commonBuiltins = ["abs", "all", "any", "bin", "bool", "bytearray", "callable", "chr",
                        "classmethod", "compile", "complex", "delattr", "dict", "dir", "divmod",
                        "enumerate", "eval", "filter", "float", "format", "frozenset",
                        "getattr", "globals", "hasattr", "hash", "help", "hex", "id",
                        "input", "int", "isinstance", "issubclass", "iter", "len",
                        "list", "locals", "map", "max", "memoryview", "min", "next",
                        "object", "oct", "open", "ord", "pow", "property", "range",
                        "repr", "reversed", "round", "set", "setattr", "slice",
                        "sorted", "staticmethod", "str", "sum", "super", "tuple",
                        "type", "vars", "zip", "__import__", "NotImplemented",
                        "Ellipsis", "__debug__"];

  function top(state) {
    return state.scopes[state.scopes.length - 1];
  }

  CodeMirror.defineMode("python", function(conf, parserConf) {
    var ERRORCLASS = "error";

    var delimiters = parserConf.delimiters || parserConf.singleDelimiters || /^[\(\)\[\]\{\}@,:`=;\.]/;
    //               (Backwards-compatiblity with old, cumbersome config system)
    var operators = [parserConf.singleOperators, parserConf.doubleOperators, parserConf.doubleDelimiters, parserConf.tripleDelimiters,
                     parserConf.operators || /^([-+*/%\/&|^]=?|[<>=]+|\/\/=?|\*\*=?|!=|[~!@])/]
    for (var i = 0; i < operators.length; i++) if (!operators[i]) operators.splice(i--, 1)

    var hangingIndent = parserConf.hangingIndent || conf.indentUnit;

    var myKeywords = commonKeywords, myBuiltins = commonBuiltins;
    if (parserConf.extra_keywords != undefined)
      myKeywords = myKeywords.concat(parserConf.extra_keywords);

    if (parserConf.extra_builtins != undefined)
      myBuiltins = myBuiltins.concat(parserConf.extra_builtins);

    var py3 = !(parserConf.version && Number(parserConf.version) < 3)
    if (py3) {
      // since http://legacy.python.org/dev/peps/pep-0465/ @ is also an operator
      var identifiers = parserConf.identifiers|| /^[_A-Za-z\u00A1-\uFFFF][_A-Za-z0-9\u00A1-\uFFFF]*/;
      myKeywords = myKeywords.concat(["nonlocal", "False", "True", "None", "async", "await"]);
      myBuiltins = myBuiltins.concat(["ascii", "bytes", "exec", "print"]);
      var stringPrefixes = new RegExp("^(([rbuf]|(br))?('{3}|\"{3}|['\"]))", "i");
    } else {
      var identifiers = parserConf.identifiers|| /^[_A-Za-z][_A-Za-z0-9]*/;
      myKeywords = myKeywords.concat(["exec", "print"]);
      myBuiltins = myBuiltins.concat(["apply", "basestring", "buffer", "cmp", "coerce", "execfile",
                                      "file", "intern", "long", "raw_input", "reduce", "reload",
                                      "unichr", "unicode", "xrange", "False", "True", "None"]);
      var stringPrefixes = new RegExp("^(([rubf]|(ur)|(br))?('{3}|\"{3}|['\"]))", "i");
    }
    myBuiltins = myBuiltins.concat(["ACOSH","ACOT","ACSC","ASEC","ASIN","ASINH","ATAN","ATANH","Airy_Ai","Airy_Bi","Archive","BesselJ","BesselY","Beta","BlockDiagonal","COND","COS","COSH","COT","CSC","CST","Celsius2Fahrenheit","Ci","Circle","ClrDraw","ClrGraph","ClrIO","Col","CopyVar","CyclePic","DIGITS","DOM_COMPLEX","DOM_FLOAT","DOM_FUNC","DOM_IDENT","DOM_INT","DOM_LIST","DOM_RAT","DOM_STRING","DOM_SYMBOLIC","DOM_int","DelFold","DelVar","Det","Dialog","Digits","Dirac","Disp","DispG","DispHome","DrawFunc","DrawInv","DrawParm","DrawPol","DrawSlp","DropDown","DrwCtour","ERROR","EXP","Ei","EndDlog","FALSE","Factor","Fahrenheit2Celsius","False","Fill","GF","Gamma","Gcd","GetFold","Graph","Heaviside","IFTE","Input","InputStr","Int","Inverse","JordanBlock","LN","LQ","LSQ","LU","Li","Line","LineHorz","LineTan","LineVert","NORMALD","NewFold","NewPic","Nullspace","Output","Ox_2d_unit_vector","Ox_3d_unit_vector","Oy_2d_unit_vector","Oy_3d_unit_vector","Oz_3d_unit_vector","Pause","Phi","Pi","PopUp","Psi","QR","Quo","REDIM","REPLACE","RandSeed","RclPic","Rem","Request","Resultant","Row","RplcPic","Rref","SCALE","SCALEADD","SCHUR","SIN","SVD","SVL","SWAPCOL","SWAPROW","SetFold","Si","SortA","SortD","StoPic","Store","TAN","TRUE","TeX","Text","Title","True","UTPC","UTPF","UTPN","UTPT","Unarchiv","VARS","VAS","VAS_positive","WAIT","Zeta","a2q","abcuv","about","abs","abscissa","accumulate_head_tail","acos","acos2asin","acos2atan","acosh","acot","acsc","add","additionally","adjoint_matrix","affix","algsubs","algvar","alog10","alors","altitude","and","angle","angle_radian","angleat","angleatraw","animate","animate3d","animation","ans","append","apply","approx","approx_mode","arc","arcLen","arccos","arccosh","archive","arclen","arcsin","arcsinh","arctan","arctanh","area","areaat","areaatraw","areaplot","arg","args","array","as_function_of","asc","asec","asin","asin2acos","asin2atan","asinh","assert","assign","assume","at","atan","atan2acos","atan2asin","atanh","atrig2ln","augment","auto_correlation","autosimplify","avance","avgRC","axes","back","backquote","backward","baisse_crayon","bar_plot","bartlett_hann_window","barycenter","base","basis","batons","begin","bernoulli","besselJ","besselY","betad","betad_cdf","betad_icdf","bezier","bezout_entiers","binomial","binomial_cdf","binomial_icdf","bisection_solver","bisector","bitand","bitor","bitxor","black","blackman_harris_window","blackman_window","bloc","blockmatrix","blue","bohman_window","border","boxwhisker","break","breakpoint","brent_solver","by","c1oc2","c1op2","cFactor","cSolve","cZeros","cache_tortue","camembert","canonical_form","cap","cap_flat_line","cap_round_line","cap_square_line","cas_setup","case","cat","catch","cauchy","cauchy_cdf","cauchy_icdf","cauchyd","cauchyd_cdf","cauchyd_icdf","cd","cdf","ceil","ceiling","center","center2interval","centered_cube","centered_tetrahedron","cfactor","cfsolve","changebase","char","charpoly","chinrem","chisquare","chisquare_cdf","chisquare_icdf","chisquared","chisquared_cdf","chisquared_icdf","chisquaret","choice","cholesky","choosebox","chr","chrem","circle","circumcircle","classes","click","close","coeff","coeffs","col","colDim","colNorm","colSwap","coldim","collect","colnorm","color","colspace","colswap","comDenom","comb","combine","comment","common_perpendicular","companion","compare","complex","complex_mode","complex_variables","complexroot","concat","cond","cone","confrac","conic","conj","conjugate_gradient","cont","contains","content","continue","contourplot","convert","convertir","convexhull","convolution","coordinates","copy","correlation","cos","cos2sintan","cosh","cosine_window","cot","cote","count","count_eq","count_inf","count_sup","courbe_parametrique","courbe_polaire","covariance","covariance_correlation","cpartfrac","crationalroot","crayon","cross","crossP","cross_correlation","cross_point","cross_ratio","crossproduct","csc","csolve","csv2gen","cube","cumSum","cumsum","cumulated_frequencies","curl","current_sheet","curvature","curve","cyan","cycle2perm","cycleinv","cycles2permu","cyclotomic","cylinder","dash_line","dashdot_line","dashdotdot_line","dayofweek","de","deSolve","debug","debut_enregistrement","default","degree","del","delcols","delrows","deltalist","denom","densityplot","derive","deriver","desolve","dessine_tortue","det","det_minor","developper","developper_transcendant","dfc","dfc2f","diag","diff","dim","display","disque","disque_centre","distance","distance2","distanceat","distanceatraw","div","divergence","divide","divis","division_point","divisors","divmod","divpc","dnewton_solver","do","dodecahedron","domain","dot","dotP","dot_paper","dotprod","double","droit","droite_tangente","dsolve","e","e2r","ecart_type","ecart_type_population","ecris","efface","egcd","egv","egvl","eigVc","eigVl","eigenvals","eigenvalues","eigenvectors","eigenvects","element","elif","eliminate","ellipse","else","end","end_for","end_if","end_while","entry","envelope","epaisseur","epaisseur_ligne_1","epaisseur_ligne_2","epaisseur_ligne_3","epaisseur_ligne_4","epaisseur_ligne_5","epaisseur_ligne_6","epaisseur_ligne_7","epaisseur_point_1","epaisseur_point_2","epaisseur_point_3","epaisseur_point_4","epaisseur_point_5","epaisseur_point_6","epaisseur_point_7","epsilon","epsilon2zero","equal","equal2diff","equal2list","equation","equilateral_triangle","erase","erase3d","erf","erfc","error","est_permu","et","euler","euler_gamma","eval","eval_level","evala","evalb","evalc","evalf","evalm","even","evolute","exact","exbisector","excircle","execute","exp","exp2list","exp2pow","exp2trig","expand","expexpand","expln","exponential","exponential_cdf","exponential_icdf","exponential_regression","exponential_regression_plot","exponentiald","exponentiald_cdf","exponentiald_icdf","expr","expression","extend","extract_measure","extrema","ezgcd","f2nd","fMax","fMin","fPart","faces","facteurs_premiers","factor","factor_xn","factorial","factoriser","factoriser_entier","factoriser_sur_C","factors","fadeev","faire","false","falsepos_solver","fclose","fcoeff","fdistrib","feuille","ffaire","ffonction","fft","fi","fieldplot","filled","fin_enregistrement","find","findhelp","fisher","fisher_cdf","fisher_icdf","fisherd","fisherd_cdf","fisherd_icdf","flatten","float","float2rational","floor","fonction	","fonction_derivee","fopen","for","format","forward","fourier_an","fourier_bn","fourier_cn","fpour","fprint","frac","fracmod","frame_2d","frame_3d","frames","frequencies","frobenius_norm","from","froot","fsi","fsi","fsolve","ftantque","fullparfrac","func","funcplot","function","function_diff","fxnd","gammad","gammad_cdf","gammad_icdf","gauche","gauss","gauss15","gauss_seidel_linsolve","gaussian_window","gaussjord","gaussquad","gbasis","gcd","gcdex","genpoly","geometric","geometric_cdf","geometric_icdf","getDenom","getKey","getNum","getType","gl_ortho","gl_quaternion","gl_rotation","gl_showaxes","gl_shownames","gl_texture","gl_x","gl_x_axis_color","gl_x_axis_name","gl_x_axis_unit","gl_xtick","gl_y","gl_y_axis_color","gl_y_axis_name","gl_y_axis_unit","gl_ytick","gl_z","gl_z_axis_color","gl_z_axis_name","gl_z_axis_unit","gl_ztick","gnuplot","goto","grad","gramschmidt","graph2tex","graph3d2tex","graphe","graphe3d","graphe_suite","greduce","green","grid_paper","groupermu","hadamard","half_cone","half_line","halftan","halftan_hyp2exp","halt","hamdist","hamming_window","hann_poisson_window","hann_window","harmonic_conjugate","harmonic_division","has","hasard","head","heading","heapify","heappop","heappush","hermite","hessenberg","hessian","heugcd","hexagon","hidden_name","highpass","hilbert","histogram","hold","homothety","horner","hybrid_solver","hybridj_solver","hybrids_solver","hybridsj_solver","hyp2exp","hyperbola","i","iPart","iabcuv","ibasis","ibpdv","ibpu","icdf","ichinrem","ichrem","icontent","icosahedron","id","identifier","identity","idivis","idn","iegcd","if","ifactor","ifactors","ifft","ifte","igamma","igcd","igcdex","ihermite","ilaplace","im","imag","image","implicitdiff","implicitplot","in","inString","in_ideal","incircle","indets","index","inequationplot","inf","infinity","input","inputform","insert","insmod","int","intDiv","integer","integrate","integrer","inter","interactive_odeplot","interactive_plotode","interp","intersect","interval","interval2center","inv","inverse","inversion","invisible_point","invlaplace","invztrans","iquo","iquorem","iratrecon","irem","isPrime","is_collinear","is_concyclic","is_conjugate","is_coplanar","is_cospheric","is_cycle","is_element","is_equilateral","is_harmonic","is_harmonic_circle_bundle","is_harmonic_line_bundle","is_included","is_inside","is_isosceles","is_orthogonal","is_parallel","is_parallelogram","is_permu","is_perpendicular","is_prime","is_pseudoprime","is_rectangle","is_rhombus","is_square","ismith","isobarycenter","isom","isopolygon","isosceles_triangle","isprime","ithprime","jacobi_linsolve","jacobi_symbol","jordan","jusqu_a","jusqua","jusque","keep_algext","keep_pivot","ker","kernel","kill","kolmogorovd","kolmogorovt","l1norm","l2norm","label","labels","lagrange","laguerre","laplace","laplacian","latex","lcm","lcoeff","ldegree","left","left_rectangle","legend","legendre","legendre_symbol","len","length","leve_crayon","lgcd","lhs","ligne_chapeau_carre","ligne_chapeau_plat","ligne_chapeau_rond","ligne_polygonale","ligne_polygonale_pointee","ligne_tiret","ligne_tiret_point","ligne_tiret_pointpoint","ligne_trait_plein","limit","limite","lin","line","line_inter","line_paper","line_segments","line_width_1","line_width_2","line_width_3","line_width_4","line_width_5","line_width_6","line_width_7","linear_interpolate","linear_regression","linear_regression_plot","lineariser","lineariser_trigo","linfnorm","linsolve","linspace","lis","lis_phrase","list","list2exp","list2mat","listplot","lll","ln","lname","lncollect","lnexpand","local","locus","log","log10","logarithmic_regression","logarithmic_regression_plot","logb","logistic_regression","logistic_regression_plot","lower","lowpass","lp_assume","lp_bestprojection","lp_binary","lp_binaryvariables","lp_breadthfirst","lp_depthfirst","lp_depthlimit","lp_firstfractional","lp_gaptolerance","lp_hybrid","lp_initialpoint","lp_integer","lp_integertolerance","lp_integervariables","lp_interiorpoint","lp_iterationlimit","lp_lastfractional","lp_maxcuts","lp_maximize","lp_method","lp_mostfractional","lp_nodelimit","lp_nodeselect","lp_nonnegative","lp_nonnegint","lp_pseudocost","lp_simplex","lp_timelimit","lp_variables","lp_varselect","lp_verbose","lpsolve","lsmod","lsq","lu","lvar","mRow","mRowAdd","magenta","makelist","makemat","makesuite","makevector","map","maple2mupad","maple2xcas","maple_ifactors","maple_mode","markov","mat2list","mathml","matpow","matrix","matrix_norm","max","maximize","maxnorm","mean","median","median_line","member","mgf","mid","middle_point","midpoint","min","minimax","minimize","minus","mkisom","mksa","mod","modgcd","mods","montre_tortue","moustache","moyal","moyenne","mul","mult_c_conjugate","mult_conjugate","multinomial","multiplier_conjugue","multiplier_conjugue_complexe","multiply","mupad2maple","mupad2xcas","nCr","nDeriv","nInt","nPr","nSolve","ncols","negbinomial","negbinomial_cdf","negbinomial_icdf","newList","newMat","newton","newton_solver","newtonj_solver","nextperm","nextprime","nlpsolve","nodisp","nom_cache","non","non_recursive_normal","nop","nops","norm","normal","normal_cdf","normal_icdf","normald","normald_cdf","normald_icdf","normalize","normalt","not","nprimes","nrows","nuage_points","nullspace","numer","octahedron","od","odd","odeplot","odesolve","of","op","open","open_polygon","option","or","ord","order_size","ordinate","orthocenter","orthogonal","osculating_circle","otherwise","ou","output","p1oc2","p1op2","pa2b2","pade","parabola","parallel","parallelepiped","parallelogram","parameq","parameter","paramplot","parfrac","pari","part","partfrac","parzen_window","pas","pas_de_cote","pcar","pcar_hessenberg","pcoef","pcoeff","pencolor","pendown","penup","perimeter","perimeterat","perimeteratraw","periodic","perm","perminv","permu2cycles","permu2mat","permuorder","perpen_bisector","perpendicular","peval","pi","piecewise","pivot","pixoff","pixon","plane","playsnd","plex","plot","plot3d","plotarea","plotcdf","plotcontour","plotdensity","plotfield","plotfunc","plotimplicit","plotinequation","plotlist","plotode","plotparam","plotpolar","plotproba","plotseq","plus_point","pmin","point","point2d","point3d","point_carre","point_croix","point_etoile","point_invisible","point_losange","point_milieu","point_plus","point_point","point_triangle","point_width_1","point_width_2","point_width_3","point_width_4","point_width_5","point_width_6","point_width_7","poisson","poisson_cdf","poisson_icdf","poisson_window","polar","polar_coordinates","polar_point","polarplot","pole","poly2symb","polyEval","polygon","polygone_rempli","polygonplot","polygonscatterplot","polyhedron","polynom","polynomial_regression","polynomial_regression_plot","position","poslbdLMQ","posubLMQ","potential","pour","pow","pow2exp","power_regression","power_regression_plot","powermod","powerpc","powexpand","powmod","prepend","preval","prevperm","prevprime","primpart","print","printf","prism","proc","product","program","projection","proot","propFrac","propfrac","psrgcd","ptayl","purge","pwd","pyramid","python","python_compat","q2a","qr","quadrant1","quadrant2","quadrant3","quadrant4","quadric","quadrilateral","quantile","quartile1","quartile3","quartiles","quest","quo","quorem","quote","r2e","radical_axis","radius","ramene","rand","randMat","randNorm","randPoly","randbinomial","randchisquare","randexp","randfisher","randgeometric","randint","randmarkov","randmatrix","randmultinomial","randnorm","random","randperm","randpoisson","randpoly","randseed","randstudent","randvector","range","rank","ranm","ranv","rassembler_trigo","rat_jordan","rational","rationalroot","ratnormal","rcl","rdiv","re","read","readrgb","readwav","real","realroot","reciprocation","rectangle","rectangle_droit","rectangle_gauche","rectangle_plein","rectangular_coordinates","recule","red","redim","reduced_conic","reduced_quadric","ref","reflection","regroup","rem","remain","remove","reorder","repeat","repete","repeter","replace","residue","resoudre","resoudre_dans_C","resoudre_systeme_lineaire","restart","resultant","return","reverse","reverse_rsolve","revert","revlex","revlist","rhombus","rhombus_point","rhs","riemann_window","right","right_rectangle","right_triangle","risch","rm_a_z","rm_all_vars","rmbreakpoint","rmmod","rmwatch","romberg","rombergm","rombergt","rond","root","rootof","roots","rotate","rotation","round","row","rowAdd","rowDim","rowNorm","rowSwap","rowdim","rownorm","rowspace","rowswap","rref","rsolve","same","sample","sans_factoriser","saute","sauve","save_history","scalarProduct","scalar_product","scale","scaleadd","scatterplot","schur","sec","secant_solver","segment","select","semi_augment","seq","seqplot","seqsolve","series","shift","shift_phase","shuffle","si","sign","signature","signe","similarity","simp2","simplex_reduce","simplifier","simplify","simpson","simult","sin","sin2costan","sincos","single_inter","sinh","sinon","size","sizes","slope","slopeat","slopeatraw","smith","smith","smod","snedecor","snedecor_cdf","snedecor_icdf","snedecord","snedecord_cdf","snedecord_icdf","solid_line","solve","somme","sommet","sort","sorta","sortd","soundsec","sphere","spline","split","sq","sqrfree","sqrt","square","square_point","srand","sst","sst_in","stack","star_point","start","stdDev","stddev","stddevp","steffenson_solver","step","sto","str","string","string","student","student_cdf","student_icdf","studentd","studentt","sturm","sturmab","sturmseq","style","subMat","subs","subsop","subst","substituer","subtype","sum","sum_riemann","suppress","surd","svd","swapcol","swaprow","switch","switch_axes","sylvester","symb2poly","symbol","syst2mat","tCollect","tExpand","table","tablefunc","tableseq","tabvar","tail","tan","tan2cossin2","tan2sincos","tan2sincos2","tangent","tangente","tanh","tantque","taux_accroissement","taylor","tchebyshev1","tchebyshev2","tcoeff","tcollect","tdeg","test","tetrahedron","texpand","textinput","then","thickness","thiele","threshold","throw","time","title","titre","tlin","to","tourne_droite","tourne_gauche","tpsolve","trace","trames","tran","translation","transpose","trapeze","trapezoid","triangle","triangle_paper","triangle_plein","triangle_point","triangle_window","trig2exp","trigcos","trigexpand","triginterp","trigsimplify","trigsin","trigtan","trn","true","trunc","truncate","try","tsimplify","tukey_window","type","ufactor","ugamma","unapply","unarchive","unfactored","uniform","uniform_cdf","uniform_icdf","uniformd","uniformd_cdf","uniformd_icdf","union","unitV","unquote","until","upper","user_operator","usimplify","valuation","vandermonde","var","variables_are_files","variance","vector","vector","vers","version","vertices","vertices_abc","vertices_abca","vpotential","watch","weibull","weibull_cdf","weibull_icdf","weibulld","weibulld_cdf","weibulld_icdf","welch_window","when","while","white","widget_size","wilcoxonp","wilcoxons","wilcoxont","write","writergb","writewav","wz_certificate","xcas","xcas_mode","xor","xyztrange","yellow","zeros","zip","ztrans",
				    "ΔLIST","ΠLIST","Σ","ΣLIST","∂","∫"]);
    var keywords = wordRegexp(myKeywords);
    var builtins = wordRegexp(myBuiltins);
    CodeMirror.registerHelper("hintWords", "python", commonKeywords.concat(commonBuiltins).concat(myBuiltins));

    // tokenizers
    function tokenBase(stream, state) {
      if (stream.sol()) state.indent = stream.indentation()
      // Handle scope changes
      if (stream.sol() && top(state).type == "py") {
        var scopeOffset = top(state).offset;
        if (stream.eatSpace()) {
          var lineOffset = stream.indentation();
          if (lineOffset > scopeOffset)
            pushPyScope(state);
          else if (lineOffset < scopeOffset && dedent(stream, state) && stream.peek() != "#")
            state.errorToken = true;
          return null;
        } else {
          var style = tokenBaseInner(stream, state);
          if (scopeOffset > 0 && dedent(stream, state))
            style += " " + ERRORCLASS;
          return style;
        }
      }
      return tokenBaseInner(stream, state);
    }

    function tokenBaseInner(stream, state) {
      if (stream.eatSpace()) return null;

      var ch = stream.peek();

      // Handle Comments
      if (ch == "#") {
        stream.skipToEnd();
        return "comment";
      }

      // Handle Number Literals
      if (stream.match(/^[0-9\.]/, false)) {
        var floatLiteral = false;
        // Floats
        if (stream.match(/^[\d_]*\.\d+(e[\+\-]?\d+)?/i)) { floatLiteral = true; }
        if (stream.match(/^[\d_]+\.\d*/)) { floatLiteral = true; }
        if (stream.match(/^\.\d+/)) { floatLiteral = true; }
        if (floatLiteral) {
          // Float literals may be "imaginary"
          stream.eat(/J/i);
          return "number";
        }
        // Integers
        var intLiteral = false;
        // Hex
        if (stream.match(/^0x[0-9a-f_]+/i)) intLiteral = true;
        // Binary
        if (stream.match(/^0b[01_]+/i)) intLiteral = true;
        // Octal
        if (stream.match(/^0o[0-7_]+/i)) intLiteral = true;
        // Decimal
        if (stream.match(/^[1-9][\d_]*(e[\+\-]?[\d_]+)?/)) {
          // Decimal literals may be "imaginary"
          stream.eat(/J/i);
          // TODO - Can you have imaginary longs?
          intLiteral = true;
        }
        // Zero by itself with no other piece of number.
        if (stream.match(/^0(?![\dx])/i)) intLiteral = true;
        if (intLiteral) {
          // Integer literals may be "long"
          stream.eat(/L/i);
          return "number";
        }
      }

      // Handle Strings
      if (stream.match(stringPrefixes)) {
        state.tokenize = tokenStringFactory(stream.current());
        return state.tokenize(stream, state);
      }

      for (var i = 0; i < operators.length; i++)
        if (stream.match(operators[i])) return "operator"

      if (stream.match(delimiters)) return "punctuation";

      if (state.lastToken == "." && stream.match(identifiers))
        return "property";

      if (stream.match(keywords) || stream.match(wordOperators))
        return "keyword";

      if (stream.match(builtins))
        return "builtin";

      if (stream.match(/^(self|cls)\b/))
        return "variable-2";

      if (stream.match(identifiers)) {
        if (state.lastToken == "def" || state.lastToken == "class")
          return "def";
        return "variable";
      }

      // Handle non-detected items
      stream.next();
      return ERRORCLASS;
    }

    function tokenStringFactory(delimiter) {
      while ("rubf".indexOf(delimiter.charAt(0).toLowerCase()) >= 0)
        delimiter = delimiter.substr(1);

      var singleline = delimiter.length == 1;
      var OUTCLASS = "string";

      function tokenString(stream, state) {
        while (!stream.eol()) {
          stream.eatWhile(/[^'"\\]/);
          if (stream.eat("\\")) {
            stream.next();
            if (singleline && stream.eol())
              return OUTCLASS;
          } else if (stream.match(delimiter)) {
            state.tokenize = tokenBase;
            return OUTCLASS;
          } else {
            stream.eat(/['"]/);
          }
        }
        if (singleline) {
          if (parserConf.singleLineStringErrors)
            return ERRORCLASS;
          else
            state.tokenize = tokenBase;
        }
        return OUTCLASS;
      }
      tokenString.isString = true;
      return tokenString;
    }

    function pushPyScope(state) {
      while (top(state).type != "py") state.scopes.pop()
      state.scopes.push({offset: top(state).offset + conf.indentUnit,
                         type: "py",
                         align: null})
    }

    function pushBracketScope(stream, state, type) {
      var align = stream.match(/^([\s\[\{\(]|#.*)*$/, false) ? null : stream.column() + 1
      state.scopes.push({offset: state.indent + hangingIndent,
                         type: type,
                         align: align})
    }

    function dedent(stream, state) {
      var indented = stream.indentation();
      while (state.scopes.length > 1 && top(state).offset > indented) {
        if (top(state).type != "py") return true;
        state.scopes.pop();
      }
      return top(state).offset != indented;
    }

    function tokenLexer(stream, state) {
      if (stream.sol()) state.beginningOfLine = true;

      var style = state.tokenize(stream, state);
      var current = stream.current();

      // Handle decorators
      if (state.beginningOfLine && current == "@")
        return stream.match(identifiers, false) ? "meta" : py3 ? "operator" : ERRORCLASS;

      if (/\S/.test(current)) state.beginningOfLine = false;

      if ((style == "variable" || style == "builtin")
          && state.lastToken == "meta")
        style = "meta";

      // Handle scope changes.
      if (current == "pass" || current == "return")
        state.dedent += 1;

      if (current == "lambda") state.lambda = true;
      if (current == ":" && !state.lambda && top(state).type == "py")
        pushPyScope(state);

      var delimiter_index = current.length == 1 ? "[({".indexOf(current) : -1;
      if (delimiter_index != -1)
        pushBracketScope(stream, state, "])}".slice(delimiter_index, delimiter_index+1));

      delimiter_index = "])}".indexOf(current);
      if (delimiter_index != -1) {
        if (top(state).type == current) state.indent = state.scopes.pop().offset - hangingIndent
        else return ERRORCLASS;
      }
      if (state.dedent > 0 && stream.eol() && top(state).type == "py") {
        if (state.scopes.length > 1) state.scopes.pop();
        state.dedent -= 1;
      }

      return style;
    }

    var external = {
      startState: function(basecolumn) {
        return {
          tokenize: tokenBase,
          scopes: [{offset: basecolumn || 0, type: "py", align: null}],
          indent: basecolumn || 0,
          lastToken: null,
          lambda: false,
          dedent: 0
        };
      },

      token: function(stream, state) {
        var addErr = state.errorToken;
        if (addErr) state.errorToken = false;
        var style = tokenLexer(stream, state);

        if (style && style != "comment")
          state.lastToken = (style == "keyword" || style == "punctuation") ? stream.current() : style;
        if (style == "punctuation") style = null;

        if (stream.eol() && state.lambda)
          state.lambda = false;
        return addErr ? style + " " + ERRORCLASS : style;
      },

      indent: function(state, textAfter) {
        if (state.tokenize != tokenBase)
          return state.tokenize.isString ? CodeMirror.Pass : 0;

        var scope = top(state), closing = scope.type == textAfter.charAt(0)
        if (scope.align != null)
          return scope.align - (closing ? 1 : 0)
        else
          return scope.offset - (closing ? hangingIndent : 0)
      },

      electricInput: /^\s*[\}\]\)]$/,
      closeBrackets: {triples: "'\""},
      lineComment: "#",
      fold: "indent"
    };
    return external;
  });

  CodeMirror.defineMIME("text/x-python", "python");

  var words = function(str) { return str.split(" "); };

  CodeMirror.defineMIME("text/x-cython", {
    name: "python",
    extra_keywords: words("by cdef cimport cpdef ctypedef enum except "+
                          "extern gil include nogil property public "+
                          "readonly struct union DEF IF ELIF ELSE")
  });

});
