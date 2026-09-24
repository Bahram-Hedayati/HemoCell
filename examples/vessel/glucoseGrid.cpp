#include "glucoseGrid.h"
#include "surfaceRay.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <sstream>
#include <unordered_map>

namespace vessel {
namespace {
Vec sub(const Vec& a, const Vec& b) { return {{a[0]-b[0], a[1]-b[1], a[2]-b[2]}}; }
Vec cross(const Vec& a, const Vec& b) {
    return {{a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]}};
}
double dot(const Vec& a, const Vec& b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
double minmod(double a, double b) {
    return a*b <= 0 ? 0 : std::copysign(std::min(std::abs(a), std::abs(b)), a);
}
struct Triangle { Vec a, b, c; };
bool intersects(const Vec& start, const Vec& end, const Triangle& tri) {
    const Vec d=sub(end,start), e1=sub(tri.b,tri.a), e2=sub(tri.c,tri.a);
    const Vec p=cross(d,e2); const double det=dot(e1,p);
    if (std::abs(det)<1e-12) return false;
    const Vec s=sub(start,tri.a); const double u=dot(s,p)/det;
    if (u < -1e-10 || u > 1+1e-10) return false;
    const Vec q=cross(s,e1); const double v=dot(d,q)/det;
    if (v < -1e-10 || u+v > 1+1e-10) return false;
    const double t=dot(e2,q)/det;
    return t > 1e-7 && t < 1-1e-7;
}
}

GlucoseGrid::GlucoseGrid(int x, int y, int z, double spacing, double timestep, double d)
    : nx(x), ny(y), nz(z), dx(spacing), dt(timestep), diffusion(d),
      diffusionLattice(d*timestep/(spacing*spacing)) {
    if (x<3 || y<3 || z<3 || !std::isfinite(dx) || !std::isfinite(dt)
        || !std::isfinite(d) || dx<=0 || dt<=0 || d<0)
        throw std::runtime_error("Invalid glucose grid or physical parameters");
    const std::size_t n=std::size_t(nx)*ny*nz;
    if (n>std::size_t(std::numeric_limits<int>::max()/3))
        throw std::runtime_error("Glucose reference solver grid is too large");
    amount.assign(n,0); walls.assign(n,1); plasma=walls;
    open.resize(n); sourceWeights.assign(n,0); stage.resize(n); result.resize(n); slope.resize(n);
}
std::size_t GlucoseGrid::index(int x, int y, int z) const {
    x=(x%nx+nx)%nx;
    return x+std::size_t(nx)*(y+ny*z);
}
Vec GlucoseGrid::position(int i) const { return {{double(i%nx),double((i/nx)%ny),double(i/(nx*ny))}}; }
int GlucoseGrid::neighbor(int i,int a,int direction) const {
    Vec p=position(i); int c=int(p[a])+direction;
    if (a!=0 && (c<0 || c>=(a==1?ny:nz))) return -1;
    p[a]=c; return int(index(int(p[0]),int(p[1]),int(p[2])));
}
void GlucoseGrid::setWalls(const std::vector<unsigned char>& fluid) {
    if (fluid.size()!=amount.size()) throw std::runtime_error("Glucose wall mask size mismatch");
    walls=fluid;
}
bool GlucoseGrid::connected(int i,int a,int direction,
    const std::vector<std::array<unsigned char,3>>& links) const {
    const int j=neighbor(i,a,direction);
    return j>=0 && (direction>0?links[i][a]:links[j][a]);
}

void GlucoseGrid::updateGeometry(const std::vector<Mesh>& meshes) {
    const auto oldOpen=open;
    const auto oldPlasma=plasma;
    plasma=walls;
    for (auto& edge:open) edge={{1,1,1}};
    sourceWeights.assign(amount.size(),0); sourceArea=0;
    // Spatial bins for short source-to-plasma visibility tests.
    std::vector<Triangle> triangles;
    std::unordered_map<long long,std::vector<int>> bins;
    const int margin=5;
    const auto key=[&](int x,int y,int z) {
        return (x+margin)+static_cast<long long>(nx+2*margin)*
            ((y+margin)+static_cast<long long>(ny+2*margin)*(z+margin));
    };
    struct Patch { Vec center, normal; double area; };
    std::vector<Patch> patches;
    const int dims[3]={nx,ny,nz};
    for (const Mesh& mesh:meshes) {
        if (mesh.vertices.empty()) throw std::runtime_error("Empty membrane mesh");
        double volume=0;
        const Vec origin=mesh.vertices.front();
        for (const Face& f:mesh.faces) {
            for (int v:f) if (v<0 || std::size_t(v)>=mesh.vertices.size())
                throw std::runtime_error("Invalid membrane triangle index");
            volume+=dot(sub(mesh.vertices[f[0]],origin),
                        cross(sub(mesh.vertices[f[1]],origin),sub(mesh.vertices[f[2]],origin)))/6;
        }
        if (std::abs(volume)<1e-10) throw std::runtime_error("Degenerate membrane volume");
        if (mesh.transmitter) for (const Face& f:mesh.faces) {
            const Vec a=mesh.vertices[f[0]], b=mesh.vertices[f[1]], c=mesh.vertices[f[2]];
            Vec n=cross(sub(b,a),sub(c,a)); const double len=std::sqrt(dot(n,n));
            if (len<1e-12) throw std::runtime_error("Degenerate TX surface triangle");
            Vec center; for (int d=0;d<3;++d) {center[d]=(a[d]+b[d]+c[d])/3; n[d]*=(volume>0?1:-1)/len;}
            center[0]-=std::floor((center[0]+0.5)/nx)*nx;
            patches.push_back({center,n,len/2}); sourceArea+=len/2;
        }
        // Reconstruct complete periodic images before ray casting. Each cell
        // independently contributes its interior; overlapping cells form a union.
        for (int image=-1;image<=1;++image) {
            double minX=mesh.vertices.front()[0],maxX=minX;
            for(const Vec& v:mesh.vertices){minX=std::min(minX,v[0]);maxX=std::max(maxX,v[0]);}
            if(maxX+image*nx < -0.5 || minX+image*nx > nx-0.5)continue;
            std::vector<Triangle> imageTriangles;
            for (const Face& f:mesh.faces) {
                Triangle tri{mesh.vertices[f[0]],mesh.vertices[f[1]],mesh.vertices[f[2]]};
                tri.a[0]+=image*nx;tri.b[0]+=image*nx;tri.c[0]+=image*nx;
                const int loX=int(std::floor(std::min({tri.a[0],tri.b[0],tri.c[0]})));
                const int hiX=int(std::floor(std::max({tri.a[0],tri.b[0],tri.c[0]})));
                imageTriangles.push_back(tri);
                if (hiX < -margin || loX>=nx+margin) continue;
                const int id=int(triangles.size()); triangles.push_back(tri);
                int lo[3],hi[3];
                for(int d=0;d<3;++d) {
                    lo[d]=std::max(-margin,int(std::floor(std::min({tri.a[d],tri.b[d],tri.c[d]}))));
                    hi[d]=std::min(dims[d]+margin-1,int(std::floor(std::max({tri.a[d],tri.b[d],tri.c[d]}))));
                }
                for(int z=lo[2];z<=hi[2];++z) for(int y=lo[1];y<=hi[1];++y) for(int x=lo[0];x<=hi[0];++x)
                    bins[key(x,y,z)].push_back(id);
            }
            if (imageTriangles.empty()) continue;
            for (int axis=0;axis<3;++axis) {
                const int u=(axis+1)%3,v=(axis+2)%3;
                std::vector<std::vector<double>> rays(dims[u]*dims[v]);
                for (const Triangle& tri:imageTriangles) {
                    const int ulo=std::max(0,int(std::ceil(std::min({tri.a[u],tri.b[u],tri.c[u]})-1e-7)));
                    const int uhi=std::min(dims[u]-1,int(std::floor(std::max({tri.a[u],tri.b[u],tri.c[u]})+1e-7)));
                    const int vlo=std::max(0,int(std::ceil(std::min({tri.a[v],tri.b[v],tri.c[v]})-1e-7)));
                    const int vhi=std::min(dims[v]-1,int(std::floor(std::max({tri.a[v],tri.b[v],tri.c[v]})+1e-7)));
                    for(int j=vlo;j<=vhi;++j) for(int i=ulo;i<=uhi;++i) {
                        double hit;
                        if(surfaceRayIntersection(tri.a,tri.b,tri.c,axis,
                                                  i+1.234e-8L,j+2.345e-8L,hit))
                            rays[i+dims[u]*j].push_back(hit);
                    }
                }
                for(int j=0;j<dims[v];++j) for(int i=0;i<dims[u];++i) {
                    auto& hits=rays[i+dims[u]*j]; if(hits.empty()) continue;
                    std::sort(hits.begin(),hits.end());
                    // Keep paired close crossings. Shared-edge duplicates are
                    // excluded by the half-open triangle predicate above.
                    if(hits.size()%2) {
                        std::ostringstream error;
                        error << "Unclosed membrane ray: cell=" << mesh.id << " image=" << image
                              << " axis=" << axis << " transverse_indices=" << i << "," << j << " hits=";
                        for(double hit:hits)error << hit << ",";
                        throw std::runtime_error(error.str());
                    }
                    int p[3];p[u]=i;p[v]=j;
                    for(double hit:hits) {
                        const int k=int(std::floor(hit));
                        if(k>=0 && k<dims[axis]) {p[axis]=k;open[index(p[0],p[1],p[2])][axis]=0;}
                        // A periodic x face may intersect the image just below zero.
                        if(axis==0 && k==-1) {p[axis]=nx-1;open[index(p[0],p[1],p[2])][axis]=0;}
                    }
                    if(axis==0) for(std::size_t h=0;h<hits.size();h+=2) {
                        for(int k=std::max(0,int(std::ceil(hits[h])));k<=std::min(nx-1,int(std::floor(hits[h+1])));++k) {
                            p[axis]=k;plasma[index(p[0],p[1],p[2])]=0;
                        }
                    }
                }
            }
        }
    }
    for(int i=0;i<int(amount.size());++i) for(int a=0;a<3;++a) {
        int j=neighbor(i,a,1);
        if(!plasma[i] || j<0 || !plasma[j])open[i][a]=0;
    }
    // Moving-boundary remap: evacuate swallowed control volumes along the OLD
    // plasma graph. This cannot shortcut through an old impermeable membrane.
    // Newly uncovered control volumes initially contain no molecules.
    std::vector<double> transfers(amount.size(),0);
    std::vector<int> visited(amount.size(),0); int stamp=0;
    if(initialized) for(int i=0;i<int(amount.size());++i) if(!plasma[i] && amount[i]>0) {
        ++stamp; std::vector<int> frontier(1,i), destinations;visited[i]=stamp;
        while(destinations.empty() && !frontier.empty()) {
            std::vector<int> next;
            for(int from:frontier) for(int a=0;a<3;++a) for(int sign:{-1,1}) {
                if(!connected(from,a,sign,oldOpen))continue;
                const int to=neighbor(from,a,sign);if(visited[to]==stamp)continue;
                visited[to]=stamp;
                if(plasma[to])destinations.push_back(to);else next.push_back(to);
            }
            frontier.swap(next);
        }
        if(destinations.empty()) {
            // Six-link voxel connectivity can disappear in a narrow oblique
            // gap even though a continuous exterior path still exists. Search
            // nearby diagonal destinations, accepting only segments that do
            // not cross ANY of the old membrane triangles.
            const Vec start=position(i);
            for(int radius=1;radius<=3 && destinations.empty();++radius) {
                for(int z=int(start[2])-radius;z<=int(start[2])+radius;++z)
                for(int y=int(start[1])-radius;y<=int(start[1])+radius;++y)
                for(int x=int(start[0])-radius;x<=int(start[0])+radius;++x) {
                    if(y<0||y>=ny||z<0||z>=nz)continue;
                    const int to=int(index(x,y,z));if(!plasma[to]||!oldPlasma[to])continue;
                    const Vec end{{double(x),double(y),double(z)}};
                    bool clear=true;
                    for(const Mesh& old:previousMeshes) {
                        for(int image=-1;image<=1 && clear;++image)for(const Face& f:old.faces) {
                            Triangle tri{old.vertices[f[0]],old.vertices[f[1]],old.vertices[f[2]]};
                            tri.a[0]+=image*nx;tri.b[0]+=image*nx;tri.c[0]+=image*nx;
                            bool nearby=true;
                            for(int a=0;a<3;++a) {
                                if(std::max({tri.a[a],tri.b[a],tri.c[a]})<std::min(start[a],end[a])
                                   ||std::min({tri.a[a],tri.b[a],tri.c[a]})>std::max(start[a],end[a])){nearby=false;break;}
                            }
                            if(nearby&&intersects(start,end,tri)){clear=false;break;}
                        }
                        if(!clear)break;
                    }
                    if(clear)destinations.push_back(to);
                }
            }
        }
        if(destinations.empty())throw std::runtime_error("Moving membrane trapped glucose with no resolved exterior path; refine plasma grid/cell clearance");
        const double portion=amount[i]/destinations.size();
        for(int j:destinations)transfers[j]+=portion;
        remappedMolecules+=amount[i];amount[i]=0;
    }
    for(std::size_t i=0;i<amount.size();++i)amount[i]+=transfers[i];
    initialized=true;
    // One quadrature point per triangle, weighted by instantaneous triangle
    // area. Each patch has its own normalized, exterior-only deposition stencil.
    for(const Patch& patch:patches) {
        Vec start=patch.center;for(int a=0;a<3;++a)start[a]+=patch.normal[a]*1e-6;
        std::vector<std::pair<int,double>> candidates; double total=0;
        for(int radius=1;radius<=3 && candidates.empty();++radius) {
            for(int z=int(std::floor(start[2]))-radius;z<=int(std::ceil(start[2]))+radius;++z)
            for(int y=int(std::floor(start[1]))-radius;y<=int(std::ceil(start[1]))+radius;++y)
            for(int x=int(std::floor(start[0]))-radius;x<=int(std::ceil(start[0]))+radius;++x) {
                if(y<0||y>=ny||z<0||z>=nz)continue;
                const int id=int(index(x,y,z));if(!plasma[id])continue;
                Vec end{{double(x),double(y),double(z)}},d=sub(end,patch.center);
                const double d2=dot(d,d);if(dot(d,patch.normal)<1e-8 || d2>radius*radius*3.)continue;
                bool visible=true;
                for(int bz=int(std::floor(std::min(start[2],end[2])));bz<=int(std::floor(std::max(start[2],end[2])))&&visible;++bz)
                for(int by=int(std::floor(std::min(start[1],end[1])));by<=int(std::floor(std::max(start[1],end[1])))&&visible;++by)
                for(int bx=int(std::floor(std::min(start[0],end[0])));bx<=int(std::floor(std::max(start[0],end[0])))&&visible;++bx) {
                    auto bin=bins.find(key(bx,by,bz));if(bin==bins.end())continue;
                    for(int t:bin->second)if(intersects(start,end,triangles[t])){visible=false;break;}
                }
                if(visible){const double w=std::exp(-2*d2);candidates.push_back({id,w});total+=w;}
            }
        }
        if(candidates.empty())throw std::runtime_error("TX surface patch has no resolved exterior plasma support; refine grid or increase cell clearance");
        for(const auto& c:candidates)sourceWeights[c.first]+=patch.area/sourceArea*c.second/total;
    }
    previousMeshes=meshes;
    if(!patches.empty()) {
        const double sum=std::accumulate(sourceWeights.begin(),sourceWeights.end(),0.);
        for(double& w:sourceWeights)w/=sum;
    }
}
void GlucoseGrid::release(double molecules) {
    if(!std::isfinite(molecules)||molecules<0)throw std::runtime_error("Invalid glucose release amount");
    if(molecules==0)return;
    if(sourceArea<=0)throw std::runtime_error("No TX surface for glucose release");
    for(std::size_t i=0;i<amount.size();++i)amount[i]+=molecules*sourceWeights[i];
    emittedMolecules+=molecules;
}
void GlucoseGrid::euler(const std::vector<double>& q,std::vector<double>& out,
                       const std::vector<Vec>& velocity,double step) {
    out=q;
    for(int i=0;i<int(q.size());++i) for(int a=0;a<3;++a) {
        const int lo=neighbor(i,a,-1),hi=neighbor(i,a,1);
        slope[i][a]=(plasma[i] && connected(i,a,-1,open) && connected(i,a,1,open))?
            minmod(q[i]-q[lo],q[hi]-q[i]):0;
    }
    for(int i=0;i<int(q.size());++i)for(int a=0;a<3;++a)if(open[i][a]) {
        const int j=neighbor(i,a,1);const double u=(velocity[i][a]+velocity[j][a])/2;
        const double advected=u>=0?q[i]+slope[i][a]/2:q[j]-slope[j][a]/2;
        const double flux=step*(u*advected+diffusionLattice*(q[i]-q[j]));
        out[i]-=flux;out[j]+=flux;
    }
}
void GlucoseGrid::transport(const std::vector<Vec>& velocity) {
    if(velocity.size()!=amount.size())throw std::runtime_error("Glucose velocity size mismatch");
    double bound=0;
    for(int i=0;i<int(amount.size());++i)if(plasma[i]) {
        double outgoing=0;
        for(int a=0;a<3;++a)for(int sign:{-1,1})if(connected(i,a,sign,open)) {
            int j=neighbor(i,a,sign);double u=(velocity[i][a]+velocity[j][a])*sign/2;
            if(!std::isfinite(u))throw std::runtime_error("Nonfinite plasma velocity");
            outgoing+=2*std::max(0.,u)+diffusionLattice;
        }
        bound=std::max(bound,outgoing);
    }
    if(!std::isfinite(bound)||bound>8000)throw std::runtime_error("Invalid/excessive glucose CFL bound");
    lastSubsteps=std::max(1,int(std::ceil(bound/0.8)));
    if(lastSubsteps>10000)throw std::runtime_error("Glucose transport requires too many substeps");
    const double step=1./lastSubsteps;
    for(int k=0;k<lastSubsteps;++k) {
        euler(amount,stage,velocity,step);euler(stage,result,velocity,step);
        for(std::size_t i=0;i<amount.size();++i)amount[i]=(amount[i]+result[i])/2;
    }
    validate();
}
double GlucoseGrid::mass() const {return std::accumulate(amount.begin(),amount.end(),0.);}
double GlucoseGrid::excludedMass() const {
    double m=0;for(std::size_t i=0;i<amount.size();++i)if(!plasma[i])m+=std::abs(amount[i]);return m;
}
void GlucoseGrid::validate() const {
    for(std::size_t i=0;i<amount.size();++i)
        if(!std::isfinite(amount[i]) || amount[i]<0 || (!plasma[i] && amount[i]!=0))
            throw std::runtime_error("Glucose positivity/interior exclusion invariant failed");
}
}
