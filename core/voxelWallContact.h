#ifndef HEMOCELL_VOXEL_WALL_CONTACT_H
#define HEMOCELL_VOXEL_WALL_CONTACT_H

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace hemo {
// A kinematic contact constraint on the same nearest-node voxel geometry used
// by HemoCell. This is not an elastic wall force or a sub-voxel wall model.
// isSolid receives integer lattice coordinates, including halo coordinates.
template<class Solid>
std::array<double,3> constrainVoxelWallMotion(
    const std::array<double,3>& start, const std::array<double,3>& proposed,
    Solid isSolid, bool& contacted)
{
    using Point = std::array<double,3>;
    Point position=start, remaining;
    contacted=false;
    for(int a=0;a<3;++a) {
        if(!std::isfinite(start[a]) || !std::isfinite(proposed[a]))
            throw std::runtime_error("Nonfinite membrane wall-contact position");
        if(std::abs(start[a])>std::numeric_limits<int>::max()/2.
           || std::abs(proposed[a])>std::numeric_limits<int>::max()/2.)
            throw std::runtime_error("Membrane position exceeds wall-grid index range");
        remaining[a]=proposed[a]-start[a];
    }
    const auto node=[](const Point& p) {
        return std::array<int,3>{{int(std::floor(p[0]+.5)),
                                  int(std::floor(p[1]+.5)),int(std::floor(p[2]+.5))}};
    };
    auto cell=node(position);
    if(isSolid(cell[0],cell[1],cell[2]))
        throw std::runtime_error("Membrane starts inside a solid voxel; use a valid initial state/checkpoint");

    // Each contact suppresses at least one displacement component. Up to three
    // contacts, followed by a final unrestricted tangential segment, suffice.
    for(int contact=0;contact<4;++contact) {
        cell=node(position);
        std::array<int,3> direction;
        Point next, delta;
        for(int a=0;a<3;++a) {
            direction[a]=(remaining[a]>0)-(remaining[a]<0);
            if(direction[a]) {
                next[a]=(cell[a]+.5*direction[a]-position[a])/remaining[a];
                delta[a]=1/std::abs(remaining[a]);
            } else next[a]=delta[a]=std::numeric_limits<double>::infinity();
        }
        bool hit=false;
        for(int traversed=0;traversed<100000;++traversed) {
            const double t=std::min({next[0],next[1],next[2]});
            if(t>1) {
                for(int a=0;a<3;++a)position[a]+=remaining[a];
                const auto end=node(position);
                if(isSolid(end[0],end[1],end[2]))
                    throw std::runtime_error("Wall-contact rounding placed a vertex inside solid");
                return position;
            }
            std::array<bool,3> crossed{{false,false,false}};
            for(int a=0;a<3;++a)crossed[a]=next[a]==t;
            // Check every voxel touched at an edge/corner, preventing a
            // diagonal step from slipping between solid voxels.
            std::array<bool,3> blocked{{false,false,false}};
            for(int subset=1;subset<8;++subset) {
                auto candidate=cell;bool relevant=true;
                for(int a=0;a<3;++a)if(subset&(1<<a)) {
                    if(!crossed[a]){relevant=false;break;}
                    candidate[a]+=direction[a];
                }
                if(relevant && isSolid(candidate[0],candidate[1],candidate[2])) {
                    hit=true;
                    for(int a=0;a<3;++a)if(subset&(1<<a))blocked[a]=true;
                }
            }
            if(hit) {
                contacted=true;
                for(int a=0;a<3;++a) {
                    position[a]+=t*remaining[a];
                    if(blocked[a]) {
                        // A finite lattice-space clearance also survives
                        // conversion to absolute coordinates at a block offset.
                        position[a]=cell[a]+.5*direction[a]-direction[a]*1e-9;
                        remaining[a]=0;
                    } else remaining[a]*=1-t;
                }
                break;
            }
            for(int a=0;a<3;++a)if(crossed[a]) {
                cell[a]+=direction[a];next[a]+=delta[a];
            }
            if(traversed==99999)
                throw std::runtime_error("Excessive membrane displacement in wall-contact traversal");
        }
        if(!hit)break;
    }
    throw std::runtime_error("Membrane wall-contact traversal did not converge");
}
}
#endif
