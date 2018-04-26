#include "config.h"
#include "giac.h"

using namespace std;
using namespace giac;

int main(){
  //debug_infolevel=20;
  context ct;
  gen g("x^4-1",&ct);
  cout << g << endl;
  return 0;
  cout << "Enter expressions to evaluate" << endl;
  cout << "Example: factor(x^4-1); simplify(sin(3x)/sin(x))" << endl;
  cout << "int(1/(x^4-1)); int(1/(x^4+1)^4,x,0,inf)" << endl;
  cout << "f(x):=sin(x^2); f'(2); f'(y)" << endl;
  cout << "Enter 0 to stop" << endl;
  for (int i=0;i<10;++i){
    gen g;
    g=i; // cin >> g;
    if (is_zero(g))
      break;
    cout << g << "=" << caseval(g.print(&ct).c_str()) << endl;
    //cout << g << "=" << eval(g,1,&ct) << endl;
  }
}
