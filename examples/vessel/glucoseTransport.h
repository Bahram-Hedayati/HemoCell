#ifndef VESSEL_GLUCOSE_TRANSPORT_H
#define VESSEL_GLUCOSE_TRANSPORT_H
#include "glucoseGrid.h"
#include "hemocell.h"
#include <functional>
#include <memory>
#include <string>

namespace vessel {
// Reference implementation: HemoCell remains MPI distributed; glucose geometry
// and finite-volume updates are gathered to rank zero. No halo is counted twice.
class GlucoseTransport {
public:
    explicit GlucoseTransport(hemo::HemoCell& simulation);
    void initialize();
    void advance();
    void updateMovingGeometry();
    void writeOutput();
    void saveCheckpoint();
    void loadCheckpoint();
private:
    hemo::HemoCell& sim;
    plb::Box3D box;
    int rank, nx, ny, nz;
    double rate, diffusion, dx, dt;
    std::size_t size;
    std::unique_ptr<GlucoseGrid> grid;
    std::vector<Vec> velocity;
    std::vector<std::vector<Face>> topology;
    std::vector<int> vertexCounts;
    std::string outputDirectory;
    bool firstOutput = true;
    void rootAction(const std::function<void()>& action);
    void gatherFluid(bool includeWalls);
    std::vector<Mesh> gatherMeshes();
};
}
#endif
