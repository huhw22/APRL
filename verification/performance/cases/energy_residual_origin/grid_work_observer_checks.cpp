// Analytic quadrature / pre-inlet endpoint tests. No PIC simulation.
#include "grid_work_audit.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

int main(int argc,char** argv) {
  try {
    if(argc!=2) throw std::runtime_error("Expected an empty output directory");
    const aprl::EBGridGeometry g(10,12,3,1.,1.,1.,1.e-10);
    aprl::EBFieldGrid f(g);
    f.ex().fill(2.); f.ey().fill(2.); f.ez().fill(2.);
    f.jx().fill(3.); f.jy().fill(3.); f.jz().fill(3.);
    aprl::GridWorkAudit audit(argv[1],0,2,2);
    audit.begin(f,0);
    f.ex().fill(4.); f.ey().fill(4.); f.ez().fill(4.);
    aprl::GridWorkAudit::beforeInlet(&audit,f);
    for(std::size_t j=0;j<f.ex().ny();++j)
      for(std::size_t i=0;i<f.ex().nx();++i) f.ex()(i,j,0)=1000.;
    for(std::size_t j=0;j<f.ey().ny();++j)
      for(std::size_t i=0;i<f.ey().nx();++i) f.ey()(i,j,0)=1000.;
    audit.complete(f); audit.close();
    std::ifstream input(std::string(argv[1])+"/mesh-work-rank-00000.csv");
    std::string line; std::getline(input,line); std::getline(input,line);
    std::istringstream row(line); std::vector<long double> values;
    while(std::getline(row,line,',')) values.push_back(std::stold(line));
    if(values.size()!=9 || values[0]!=1.) throw std::runtime_error("Invalid CSV");
    const long double cells=6*8*3,scale=g.dt;
    const long double inlet=6*8*2*.5L*3*(1000-4)*.5L*scale;
    const long double expected[8]={cells*18*scale,cells*27*scale,cells*36*scale,
      cells*18*scale,cells*27*scale,cells*36*scale,inlet,inlet};
    for(unsigned n=0;n<8;++n)
      if(std::abs(values[n+1]-expected[n])>1.e-14L*std::abs(expected[n]))
        throw std::runtime_error("Analytic J dot E quadrature mismatch");
    if(f.ex()(0,0,0)!=1000. || f.ey()(0,0,0)!=1000. || f.ez()(0,0,0)!=4. || f.jx()(0,0,0)!=3.)
      throw std::runtime_error("Observer mutated a field/current");
    std::cout << "PASS analytic uniform native/collocated old/mid/new work; inlet overwrite excluded and separately measured; fields/current unchanged\n";
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
