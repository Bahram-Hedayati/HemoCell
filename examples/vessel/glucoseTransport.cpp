#include "glucoseTransport.h"
#include "cellMechanics.h"
#include "palabos3D.h"
#include "palabos3D.hh"
#include <hdf5.h>
#include <hdf5_hl.h>
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

namespace vessel {
namespace {
void checked(herr_t status) { if(status<0)throw std::runtime_error("Glucose HDF5 operation failed"); }
struct H5File {
    hid_t id;
    explicit H5File(const std::string& path,bool read=false):id(read?H5Fopen(path.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT):H5Fcreate(path.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT)) {
        if(id<0)throw std::runtime_error("Cannot open glucose HDF5 file: "+path);
    }
    ~H5File(){H5Fclose(id);}
};
std::string iterationName(unsigned int iter) {std::ostringstream s;s<<std::setw(12)<<std::setfill('0')<<iter;return s.str();}
struct VertexRecord { int cell,type,vertex; double p[3]; };
}
GlucoseTransport::GlucoseTransport(hemo::HemoCell& simulation):sim(simulation),box(simulation.lattice->getBoundingBox()) {
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    nx=box.getNx();ny=box.getNy();nz=box.getNz();size=std::size_t(nx)*ny*nz;
    dx=hemo::param::dx;dt=hemo::param::dt;
    diffusion=(*sim.cfg)["glucose"]["diffusionCoefficient"].read<double>();
    rate=(*sim.cfg)["glucose"]["releaseRate"].read<double>();
    if(!std::isfinite(rate)||rate<0||!std::isfinite(diffusion)||diffusion<=0)
        throw std::runtime_error("Glucose rate must be finite/nonnegative and diffusion finite/positive");
    if(size>std::size_t(std::numeric_limits<int>::max()/3))throw std::runtime_error("Glucose grid exceeds MPI count limit");
    for(unsigned int t=0;t<sim.cellfields->size();++t) {
        std::vector<Face> faces;
        for(const auto& tri:(*sim.cellfields)[t]->triangle_list)
            faces.push_back({{int(tri[0]),int(tri[1]),int(tri[2])}});
        std::map<std::pair<int,int>,std::pair<int,int>> edges;
        for(const Face& f:faces)for(int a=0;a<3;++a) {
            const int i=f[a],j=f[(a+1)%3];
            auto& edge=edges[std::minmax(i,j)];++edge.first;edge.second+=i<j?1:-1;
        }
        for(const auto& edge:edges)if(edge.second.first!=2||edge.second.second!=0)
            throw std::runtime_error("Glucose membranes must be closed and consistently oriented");
        topology.push_back(faces);vertexCounts.push_back((*sim.cellfields)[t]->numVertex);
    }
    outputDirectory=plb::global::directories().getOutputDir()+"/glucose";
    rootAction([&](){
        grid.reset(new GlucoseGrid(nx,ny,nz,dx,dt,diffusion));velocity.resize(size);
        hemo::mkpath(outputDirectory.c_str(),0777);
    });
}
void GlucoseTransport::rootAction(const std::function<void()>& action) {
    std::string error;
    if(rank==0)try{action();}catch(const std::exception& e){error=e.what();}
    int n=int(error.size());MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
    if(n){error.resize(n);MPI_Bcast(&error[0],n,MPI_CHAR,0,MPI_COMM_WORLD);throw std::runtime_error(error);}
}
void GlucoseTransport::gatherFluid(bool includeWalls) {
    std::vector<double> local(size*3,0),all(rank==0?size*3:0);
    std::vector<unsigned char> mask(includeWalls?size:0,0),combined(rank==0&&includeWalls?size:0);
    for(plb::plint id:sim.lattice->getLocalInfo().getBlocks()) {
        auto& block=sim.lattice->getComponent(id);const auto loc=block.getLocation();
        const auto bulk=plb::SmartBulk3D(sim.lattice->getMultiBlockManagement(),id).getBulk();
        for(int z=bulk.z0;z<=bulk.z1;++z)for(int y=bulk.y0;y<=bulk.y1;++y)for(int x=bulk.x0;x<=bulk.x1;++x) {
            const std::size_t i=(x-box.x0)+std::size_t(nx)*((y-box.y0)+ny*(z-box.z0));
            auto& cell=block.get(x-loc.x,y-loc.y,z-loc.z);
            if(cell.getDynamics().isBoundary())continue;
            if(includeWalls)mask[i]=1;
            plb::Array<double,3> u;cell.computeVelocity(u);
            for(int a=0;a<3;++a)local[3*i+a]=u[a];
        }
    }
    MPI_Reduce(local.data(),rank==0?all.data():nullptr,int(size*3),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    if(includeWalls)MPI_Reduce(mask.data(),rank==0?combined.data():nullptr,int(size),MPI_UNSIGNED_CHAR,MPI_MAX,0,MPI_COMM_WORLD);
    if(rank==0){
        for(std::size_t i=0;i<size;++i)for(int a=0;a<3;++a)velocity[i][a]=all[3*i+a];
        if(includeWalls)grid->setWalls(combined);
    }
}
std::vector<Mesh> GlucoseTransport::gatherMeshes() {
    std::vector<VertexRecord> local;
    for(plb::plint id:sim.cellfields->immersedParticles->getLocalInfo().getBlocks()) {
        auto& pf=sim.cellfields->immersedParticles->getComponent(id);
        for(const auto& particle:pf.particles)if(pf.isContainedABS(particle.sv.position,pf.localDomain)) {
            VertexRecord r{};r.cell=sim.cellfields->base_cell_id(particle.sv.cellId);
            r.type=particle.sv.celltype;r.vertex=particle.sv.vertexId;
            r.p[0]=particle.sv.position[0]-box.x0;r.p[1]=particle.sv.position[1]-box.y0;r.p[2]=particle.sv.position[2]-box.z0;
            local.push_back(r);
        }
    }
    int ranks;MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    const int bytes=int(local.size()*sizeof(VertexRecord));
    std::vector<int> counts(ranks),offsets(ranks);
    MPI_Gather(&bytes,1,MPI_INT,counts.data(),1,MPI_INT,0,MPI_COMM_WORLD);
    int total=0;if(rank==0)for(int r=0;r<ranks;++r){offsets[r]=total;total+=counts[r];}
    std::vector<VertexRecord> all(total/sizeof(VertexRecord));
    MPI_Gatherv(local.data(),bytes,MPI_BYTE,all.data(),counts.data(),offsets.data(),MPI_BYTE,0,MPI_COMM_WORLD);
    std::vector<Mesh> meshes;
    rootAction([&](){
        std::map<int,std::map<int,VertexRecord>> cells;
        for(const auto& r:all) {
            if(!cells[r.cell].emplace(r.vertex,r).second)throw std::runtime_error("Duplicate owned glucose membrane vertex");
        }
        for(const auto& cell:cells) {
            const int type=cell.second.begin()->second.type;
            if(int(cell.second.size())!=vertexCounts.at(type))throw std::runtime_error("Incomplete glucose membrane; cannot allow transport through a missing surface");
            Mesh mesh;mesh.id=cell.first;mesh.transmitter=(*sim.cellfields)[type]->name=="TX";
            mesh.faces=topology[type];mesh.vertices.resize(vertexCounts[type]);
            for(int v=0;v<vertexCounts[type];++v) {
                const auto& r=cell.second.at(v);for(int a=0;a<3;++a)mesh.vertices[v][a]=r.p[a];
                if(v>0)mesh.vertices[v][0]-=std::round((mesh.vertices[v][0]-mesh.vertices[0][0])/nx)*nx;
            }
            meshes.push_back(std::move(mesh));
        }
        int tx=0;for(const auto& m:meshes)tx+=m.transmitter;
        if(tx!=1)throw std::runtime_error("Glucose release requires exactly one complete TX mesh");
    });
    return meshes;
}
void GlucoseTransport::initialize() {
    gatherFluid(true);const auto meshes=gatherMeshes();
    rootAction([&](){grid->updateGeometry(meshes);});
}
void GlucoseTransport::advance() {
    rootAction([&](){grid->release(rate*dt);grid->transport(velocity);});
}
void GlucoseTransport::updateMovingGeometry() {
    gatherFluid(false);const auto meshes=gatherMeshes();
    rootAction([&](){grid->updateGeometry(meshes);grid->validate();});
}
void GlucoseTransport::writeOutput() {
    rootAction([&](){
        grid->validate();const double expected=rate*sim.iter*dt,mass=grid->mass();
        if(std::abs(mass-grid->emittedMolecules)>1e-10*std::max(1.,grid->emittedMolecules)
           ||std::abs(expected-grid->emittedMolecules)>1e-10*std::max(1.,expected))
            throw std::runtime_error("Glucose molecule balance failed");
        const std::string base="Glucose."+iterationName(sim.iter);
        const double volume=dx*dx*dx,time=sim.iter*dt;
        std::vector<double> concentration(size),molar(size);
        for(std::size_t i=0;i<size;++i){concentration[i]=grid->amount[i]/volume;molar[i]=concentration[i]/6.02214076e23;}
        {
            H5File file(outputDirectory+"/"+base+".h5");const hsize_t dims[3]={hsize_t(nz),hsize_t(ny),hsize_t(nx)};
            checked(H5LTmake_dataset_double(file.id,"Concentration_molecules_per_m3",3,dims,concentration.data()));
            checked(H5LTmake_dataset_double(file.id,"Concentration_mol_per_m3",3,dims,molar.data()));
            checked(H5LTmake_dataset(file.id,"PlasmaMask",3,dims,H5T_NATIVE_UCHAR,grid->plasma.data()));
            checked(H5LTset_attribute_double(file.id,"/","time_s",&time,1));
            checked(H5LTset_attribute_double(file.id,"/","dx_m",&dx,1));
            checked(H5LTset_attribute_double(file.id,"/","diffusion_m2_per_s",&diffusion,1));
            checked(H5LTset_attribute_double(file.id,"/","released_molecules",&grid->emittedMolecules,1));
        }
        std::ofstream xmf(outputDirectory+"/"+base+".xmf");if(!xmf)throw std::runtime_error("Cannot write glucose XMF");
        xmf<<std::setprecision(17)<<"<?xml version=\"1.0\"?>\n<Xdmf Version=\"2.2\"><Domain><Grid Name=\"Glucose\" GridType=\"Uniform\">\n"
           <<"<Topology TopologyType=\"3DCoRectMesh\" Dimensions=\""<<nz<<" "<<ny<<" "<<nx<<"\"/>\n"
           <<"<Geometry GeometryType=\"ORIGIN_DXDYDZ\"><DataItem Dimensions=\"3\" Format=\"XML\">"
           <<box.z0*dx<<" "<<box.y0*dx<<" "<<box.x0*dx<<"</DataItem>"
           <<"<DataItem Dimensions=\"3\" Format=\"XML\">"<<dx<<" "<<dx<<" "<<dx<<"</DataItem></Geometry>\n";
        for(const std::string name:{"Concentration_molecules_per_m3","Concentration_mol_per_m3","PlasmaMask"})
            xmf<<"<Attribute Name=\""<<name<<"\" AttributeType=\"Scalar\" Center=\"Node\"><DataItem Dimensions=\""<<nz<<" "<<ny<<" "<<nx
               <<"\" NumberType=\""<<(name=="PlasmaMask"?"UInt":"Float")<<"\" Precision=\""<<(name=="PlasmaMask"?1:8)
               <<"\" Format=\"HDF\">"<<base<<".h5:/"<<name<<"</DataItem></Attribute>\n";
        xmf<<"</Grid></Domain></Xdmf>\n";
        const auto limits=std::minmax_element(concentration.begin(),concentration.end());
        if(firstOutput) {
            std::vector<std::string> history;
            if(sim.cfg->checkpointed) {
                std::ifstream old(outputDirectory+"/mass_balance.csv");std::string line;
                std::getline(old,line); // header
                while(std::getline(old,line)) {
                    std::istringstream row(line);unsigned int iteration;
                    if(row>>iteration && iteration<sim.iter)history.push_back(line);
                }
            }
        std::ofstream out(outputDirectory+"/mass_balance.csv");
        if(!out)throw std::runtime_error("Cannot write glucose diagnostics");
        out<<"iteration,time_s,released_molecules,expected_molecules,total_molecules,mass_error_molecules,excluded_molecules,min_concentration_molecules_per_m3,max_concentration_molecules_per_m3,remapped_molecules,tx_area_m2,substeps\n";
            for(const auto& row:history)out<<row<<"\n";
            firstOutput=false;
        }
        std::ofstream csv(outputDirectory+"/mass_balance.csv",std::ios::app);
        if(!csv)throw std::runtime_error("Cannot append glucose mass balance");
        csv<<std::setprecision(17)<<sim.iter<<","<<time<<","<<grid->emittedMolecules<<","<<expected<<","<<mass<<","<<mass-grid->emittedMolecules
           <<","<<grid->excludedMass()<<","<<*limits.first<<","<<*limits.second<<","<<grid->remappedMolecules<<","<<grid->sourceArea*dx*dx<<","<<grid->lastSubsteps<<"\n";
        hemo::hlog<<"(Glucose) t="<<time<<" s, released="<<grid->emittedMolecules<<", plasma="<<mass<<", excluded="<<grid->excludedMass()<<" molecules"<<std::endl;
    });
}
void GlucoseTransport::saveCheckpoint() {
    rootAction([&](){
        const std::string file=hemo::global.checkpointDirectory+"glucose.h5";
        {
            H5File out(file+".tmp");hsize_t n=size;
            checked(H5LTmake_dataset_double(out.id,"amount",1,&n,grid->amount.data()));
            const double metadata[]={double(sim.iter),double(nx),double(ny),double(nz),dx,dt,diffusion,rate,grid->emittedMolecules,grid->remappedMolecules};
            checked(H5LTset_attribute_double(out.id,"/","state_v1",metadata,10));
        }
        std::ifstream previous(file);
        if(previous.good()) {
            previous.close();
            if(std::rename(file.c_str(),(file+".old").c_str())!=0)
                throw std::runtime_error("Cannot preserve previous glucose checkpoint");
        }
        if(std::rename((file+".tmp").c_str(),file.c_str())!=0)throw std::runtime_error("Cannot install glucose checkpoint");
    });
}
void GlucoseTransport::loadCheckpoint() {
    rootAction([&](){
        H5File file(hemo::global.checkpointDirectory+"glucose.h5",true);
        hsize_t shape=0;H5T_class_t type;size_t width;
        int dimensions=0;checked(H5LTget_dataset_ndims(file.id,"amount",&dimensions));
        if(dimensions!=1)throw std::runtime_error("Invalid glucose checkpoint shape");
        checked(H5LTget_dataset_info(file.id,"amount",&shape,&type,&width));
        if(shape!=size)throw std::runtime_error("Glucose checkpoint grid mismatch");
        int attributeDimensions=0;checked(H5LTget_attribute_ndims(file.id,"/","state_v1",&attributeDimensions));
        if(attributeDimensions!=1)throw std::runtime_error("Invalid glucose checkpoint metadata");
        hsize_t attributeSize=0;checked(H5LTget_attribute_info(file.id,"/","state_v1",&attributeSize,&type,&width));
        if(attributeSize!=10)throw std::runtime_error("Unknown glucose checkpoint version");
        double state[10];checked(H5LTget_attribute_double(file.id,"/","state_v1",state));
        const double expected[]={double(sim.iter),double(nx),double(ny),double(nz),dx,dt,diffusion,rate};
        for(int i=0;i<8;++i)if(state[i]!=expected[i])throw std::runtime_error("Glucose checkpoint iteration/parameters mismatch");
        checked(H5LTread_dataset_double(file.id,"amount",grid->amount.data()));
        grid->emittedMolecules=state[8];grid->remappedMolecules=state[9];grid->validate();
    });
}
}
