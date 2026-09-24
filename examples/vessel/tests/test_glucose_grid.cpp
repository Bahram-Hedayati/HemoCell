#include "glucoseGrid.h"
#include "surfaceRay.h"
#include "voxelWallContact.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <numeric>
#include <stdexcept>
using namespace vessel;
namespace {
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
Mesh cube(Vec low,Vec high,bool tx=false) {
    Mesh m;m.transmitter=tx;
    m.vertices={{{low[0],low[1],low[2]}},{{high[0],low[1],low[2]}},{{high[0],high[1],low[2]}},{{low[0],high[1],low[2]}},
                {{low[0],low[1],high[2]}},{{high[0],low[1],high[2]}},{{high[0],high[1],high[2]}},{{low[0],high[1],high[2]}}};
    m.faces={{{0,2,1}},{{0,3,2}},{{4,5,6}},{{4,6,7}},{{0,1,5}},{{0,5,4}},
             {{1,2,6}},{{1,6,5}},{{2,3,7}},{{2,7,6}},{{3,0,4}},{{3,4,7}}};return m;
}
double variance(const GlucoseGrid& g,double mean) {
    double sum=0;for(int i=0;i<int(g.amount.size());++i){double d=g.position(i)[0]-mean;sum+=g.amount[i]*d*d;}return sum/g.mass();
}
double center(const GlucoseGrid& g){double sum=0;for(int i=0;i<int(g.amount.size());++i)sum+=g.amount[i]*g.position(i)[0];return sum/g.mass();}
void gaussian(GlucoseGrid& g,double mean,double sigma) {
    for(int i=0;i<int(g.amount.size());++i){double x=(g.position(i)[0]-mean)/sigma;g.amount[i]=std::exp(-x*x/2);}
    double m=g.mass();for(double& q:g.amount)q/=m;
}
void diffusion() {
    GlucoseGrid g(128,3,3,5e-7,1e-7,9.2e-10);g.updateGeometry({});
    gaussian(g,64,4);const double initial=variance(g,64);
    std::vector<Vec> velocity(g.amount.size(),{{0,0,0}});
    for(int n=0;n<1000;++n)g.transport(velocity);
    require(std::abs(g.mass()-1)<1e-12,"Diffusion lost mass");
    require(std::abs((variance(g,64)-initial)-2*g.diffusionLattice*1000)<1e-9,"Incorrect physical diffusivity");
    std::cout<<"diffusion: variance increase="<<variance(g,64)-initial<<" (expected 0.736 LU^2)\n";
}
double advectionError(int refinement) {
    const int n=64*refinement;GlucoseGrid g(n,3,3,1./refinement,.2/refinement,.01);g.updateGeometry({});
    const double mean=20.*refinement,sigma=3.*refinement;gaussian(g,mean,sigma);
    std::vector<Vec> velocity(g.amount.size(),{{.15,0,0}});
    const int steps=40*refinement;for(int t=0;t<steps;++t)g.transport(velocity);
    const double targetMean=mean+.15*steps,targetVariance=sigma*sigma+2*g.diffusionLattice*steps;
    double norm=0,error=0;std::vector<double> exact(g.amount.size());
    for(int i=0;i<int(exact.size());++i){double d=g.position(i)[0]-targetMean;exact[i]=std::exp(-d*d/(2*targetVariance));norm+=exact[i];}
    for(int i=0;i<int(exact.size());++i)error+=std::abs(g.amount[i]-exact[i]/norm);
    require(std::abs(g.mass()-1)<1e-12,"Advection lost mass");
    std::cout<<"advection mean error [LU]: "<<center(g)-targetMean<<"\n";
    require(std::abs(center(g)-targetMean)<0.01,"Wrong advection displacement");
    return error;
}
void membranesAndMotion() {
    GlucoseGrid g(24,12,12,5e-7,1e-7,9.2e-10);
    Mesh body=cube({{8.2,3.2,3.2}},{{12.2,7.2,7.2}},true);g.updateGeometry({body});
    require(!g.plasma[g.index(10,5,5)],"Interior node not excluded");
    require(std::abs(g.sourceArea-96)<1e-10,"Surface area incorrect");
    require(std::abs(std::accumulate(g.sourceWeights.begin(),g.sourceWeights.end(),0.)-1)<1e-12,"Source not normalized");
    for(int n=0;n<20;++n){g.release(100*g.dt);g.transport(std::vector<Vec>(g.amount.size(),{{.02,0,0}}));}
    require(std::abs(g.mass()-2e-4)<1e-15,"Release dose incorrect");require(g.excludedMass()==0,"Glucose entered membrane");
    // An advancing membrane swallows a populated exterior node.
    std::fill(g.amount.begin(),g.amount.end(),0);g.amount[g.index(8,5,5)]=1;
    for(Vec& v:body.vertices)v[0]-=1.;
    g.updateGeometry({body});
    require(std::abs(g.mass()-1)<1e-12,"Moving membrane lost mass");
    require(g.excludedMass()==0 && g.amount[g.index(8,5,5)]==0,"Moving membrane retained interior glucose");
    require(g.remappedMolecules>0,"Remap not exercised");
    // Exercise the geometric fallback when the old six-link graph is isolated.
    // Visibility must keep evacuated mass on the exterior side of the old cube.
    GlucoseGrid isolated(24,12,12,1,.1,.01);
    Mesh moving=cube({{8.2,3.2,3.2}},{{8.8,7.2,7.2}});
    isolated.updateGeometry({moving});
    const int trapped=int(isolated.index(8,5,5));isolated.amount[trapped]=1;
    for(int a=0;a<3;++a) {
        isolated.open[trapped][a]=0;
        isolated.open[isolated.neighbor(trapped,a,-1)][a]=0;
    }
    for(Vec& vertex:moving.vertices)vertex[0]-=.4;
    isolated.updateGeometry({moving});isolated.validate();
    require(std::abs(isolated.mass()-1)<1e-12,"Exterior fallback lost mass");
    for(int i=0;i<int(isolated.amount.size());++i)
        if(isolated.position(i)[0]>=9)require(isolated.amount[i]==0,"Remap crossed an old membrane");
    // Thin membrane crossing a link whose BOTH endpoints are exterior.
    GlucoseGrid thin(24,8,8,1,.1,.1);std::vector<unsigned char> walls(thin.amount.size(),1);
    for(int i=0;i<int(walls.size());++i)if(thin.position(i)[0]==0||thin.position(i)[0]==23)walls[i]=0;
    thin.setWalls(walls);thin.updateGeometry({cube({{10.2,-2,-2}},{{10.3,10,10}})});
    require(thin.plasma[thin.index(10,4,4)] && thin.plasma[thin.index(11,4,4)],"Thin barrier test does not straddle exterior nodes");
    require(!thin.open[thin.index(10,4,4)][0],"Subgrid membrane crossing left open");
    thin.amount[thin.index(9,4,4)]=1;
    for(int n=0;n<50;++n)thin.transport(std::vector<Vec>(thin.amount.size(),{{.1,0,0}}));
    for(int i=0;i<int(thin.amount.size());++i)if(thin.position(i)[0]>=11)require(thin.amount[i]==0,"Glucose crossed thin membrane");
    require(std::abs(thin.mass()-1)<1e-12,"No-flux membrane lost mass");
    // Complete cell crossing the periodic seam: still closed, source counted once.
    GlucoseGrid periodic(24,12,12,1,.1,.01);
    periodic.updateGeometry({cube({{-2.2,3.2,3.2}},{{2.2,7.2,7.2}},true)});
    require(!periodic.plasma[periodic.index(0,5,5)]&&!periodic.plasma[periodic.index(23,5,5)],"Periodic cell exclusion failed");
    periodic.release(1);require(std::abs(periodic.mass()-1)<1e-12,"Periodic source duplicated dose");
}
void substepsAndZero() {
    GlucoseGrid g(24,3,3,1,1,1);g.updateGeometry({});
    std::vector<Vec> v(g.amount.size(),{{2,0,0}});g.transport(v);require(g.mass()==0,"Zero state changed");
    g.amount[g.index(8,1,1)]=1;g.transport(v);
    require(g.lastSubsteps>1,"CFL substeps not used");require(std::abs(g.mass()-1)<1e-12,"Substeps lost mass");
}
void wallContact() {
    // Actual failing vertex, immediately before entering solid node (47,38,48).
    const Vec start{{47.383146577908008,37.499993903119496,47.524176434652446}};
    const Vec proposed{{start[0]+.003,start[1]+.00002,start[2]}};
    const auto wall=[](int,int y,int){return y>=38;};
    bool touched=false;
    const auto end=hemo::constrainVoxelWallMotion(start,proposed,wall,touched);
    require(touched && end[1]<37.5,"Wall contact did not preserve fluid-side vertex");
    require(std::abs(end[0]-proposed[0])<1e-12,"Wall contact lost tangential motion");
    const auto reverse=hemo::constrainVoxelWallMotion(end,start,wall,touched);
    require(!touched && std::abs(reverse[1]-start[1])<1e-12,"Wall contact prevents detachment");
    const auto slab=[](int x,int,int){return x==2;};
    const auto thin=hemo::constrainVoxelWallMotion(Vec{{0,0,0}},Vec{{4,0,0}},slab,touched);
    require(touched && thin[0]<1.5,"Swept wall test allowed tunneling");
    const auto corner=[](int x,int y,int){return x==1 || y==1;};
    const auto diagonal=hemo::constrainVoxelWallMotion(Vec{{0,0,0}},Vec{{2,2,0}},corner,touched);
    require(touched && diagonal[0]<.5 && diagonal[1]<.5,"Corner contact allowed penetration");
    const auto negative=hemo::constrainVoxelWallMotion(Vec{{0,0,0}},Vec{{-2,0,0}},
        [](int x,int,int){return x<0;},touched);
    require(touched && negative[0]>-.5,"Negative-direction wall contact failed");
    const auto periodic=hemo::constrainVoxelWallMotion(Vec{{102.49,0,0}},Vec{{102.51,0,0}},
        [](int,int,int){return false;},touched);
    require(!touched && periodic[0]==102.51,"Open periodic seam was treated as a wall");
    bool rejected=false;
    try{hemo::constrainVoxelWallMotion(Vec{{2,0,0}},Vec{{0,0,0}},slab,touched);}
    catch(const std::runtime_error&){rejected=true;}
    require(rejected,"Invalid initial membrane silently repaired");
}
void grazingRays() {
    // Two triangles sharing a projected diagonal must contribute one crossing,
    // regardless of winding. Two *different* nearby surfaces must remain two.
    const Vec a{{8,2,2}},b{{8,6,2}},c{{8,6,6}},d{{8,2,6}};
    double hit;
    const int count=surfaceRayIntersection(a,b,c,0,4,4,hit)
                   +surfaceRayIntersection(a,c,d,0,4,4,hit);
    const int reversed=surfaceRayIntersection(c,b,a,0,4,4,hit)
                      +surfaceRayIntersection(d,c,a,0,4,4,hit);
    require(count==1 && reversed==1,"Shared triangle edge counted twice or missed");
    GlucoseGrid thin(24,12,12,1,.1,.01);
    thin.updateGeometry({cube({{8.2,2.2,2.2}},{{8.2+3.56037617e-8,7.2,7.2}})});
    require(!thin.open[thin.index(8,4,4)][0],"Close paired crossings were merged");

    // Real RBC 179 snapshot from the step-6018 failure. No simulator or MPI is
    // required: test every mesh ray plus the exact grazing ray in isolation.
    std::ifstream input(std::string(VESSEL_TEST_FIXTURE_DIR)+"/rbc179_step6018.mesh");
    require(bool(input),"Missing RBC grazing-ray regression fixture");
    Mesh body;body.id=179;std::size_t nv=0,nf=0;
    input>>nv>>nf;require(nv==642 && nf==1280,"Unexpected RBC fixture dimensions");
    body.vertices.resize(nv);body.faces.resize(nf);
    for(auto& vertex:body.vertices)input>>vertex[0]>>vertex[1]>>vertex[2];
    for(auto& face:body.faces)input>>face[0]>>face[1]>>face[2];
    require(bool(input),"Incomplete RBC regression fixture");
    std::vector<double> hits;
    for(const auto& f:body.faces)
        if(surfaceRayIntersection(body.vertices[f[0]],body.vertices[f[1]],body.vertices[f[2]],
                                  0,19+1.234e-8L,36+2.345e-8L,hit))hits.push_back(hit);
    require(hits.size()==2 && std::abs(hits[0]-hits[1])<1e-7,"Grazing ray lost its paired crossings");
    GlucoseGrid real(103,53,53,5e-7,1e-7,9.2e-10);real.updateGeometry({body});real.validate();
}
}
int main(){try{diffusion();double coarse=advectionError(1),fine=advectionError(2);
    std::cout<<"advection-diffusion L1 errors: "<<coarse<<", "<<fine<<"\n";
    require(fine<coarse*.7,"Transport does not converge with refinement");
    membranesAndMotion();substepsAndZero();wallContact();grazingRays();
    std::cout<<"All glucose grid tests passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}
