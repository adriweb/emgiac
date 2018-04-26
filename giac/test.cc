#include <iostream>

using namespace std;
typedef unsigned long long ulonglong;

int _Random(){
  static int r=0;
  r = unsigned ((1664525*ulonglong(r)+1013904223)%(ulonglong(1)<<31));
  return r;
}

unsigned Random(){
  static unsigned r=0;
  return r = 1664525*r+1013904223;
}

double R=2*2147483648.0;
double r(){
  return Random()/R;
}

int main(int argc,char **argv){
  for (int i=0;i<100;++i){
    cerr << r() << endl;
  }
}
