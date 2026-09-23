#ifndef VESSEL_GLUCOSE_GRID_H
#define VESSEL_GLUCOSE_GRID_H

#include <array>
#include <cstddef>
#include <vector>

namespace vessel {
using Vec = std::array<double, 3>;
using Face = std::array<int, 3>;
struct Mesh {
    int id = 0;
    bool transmitter = false;
    std::vector<Vec> vertices;
    std::vector<Face> faces;
};

// Coordinates and velocities are in fluid lattice units. State is expected
// molecules per full grid control volume; concentration = amount / dx^3.
// The scalar is advanced independently of HemoCell's optional CEPAC lattice.
class GlucoseGrid {
public:
    GlucoseGrid(int nx, int ny, int nz, double dx, double dt, double diffusion);
    std::size_t index(int x, int y, int z) const;
    int neighbor(int i, int axis, int direction) const;
    Vec position(int i) const;
    void setWalls(const std::vector<unsigned char>& fluid);
    void updateGeometry(const std::vector<Mesh>& meshes);
    void release(double molecules);
    void transport(const std::vector<Vec>& velocity);
    double mass() const;
    double excludedMass() const;
    void validate() const;

    int nx, ny, nz;
    double dx, dt, diffusion, diffusionLattice;
    std::vector<double> amount;
    std::vector<unsigned char> plasma;
    std::vector<std::array<unsigned char, 3>> open;
    std::vector<double> sourceWeights;
    double sourceArea = 0; // lattice area
    double remappedMolecules = 0;
    double emittedMolecules = 0;
    int lastSubsteps = 0;
private:
    std::vector<unsigned char> walls;
    std::vector<double> stage, result;
    std::vector<Vec> slope;
    bool initialized = false;
    std::vector<Mesh> previousMeshes;
    bool connected(int i, int axis, int direction,
                   const std::vector<std::array<unsigned char, 3>>& links) const;
    void euler(const std::vector<double>& input, std::vector<double>& output,
               const std::vector<Vec>& velocity, double step);
};
}
#endif
