/*
This file is part of the HemoCell library

HemoCell is developed and maintained by the Computational Science Lab
in the University of Amsterdam. Any questions or remarks regarding this library
can be sent to: info@hemocell.eu

When using the HemoCell library in scientific work please cite the
corresponding paper: https://doi.org/10.3389/fphys.2017.00563

The HemoCell library is free software: you can redistribute it and/or
modify it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <hemocell.h>
#include <helper/voxelizeDomain.h>
#include "rbcHighOrderModel.h"
#include "pltSimpleModel.h"
#include "cellInfo.h"
#include "writeCellInfoCSV.h"
#include "glucoseTransport.h"
#include "palabos3D.h"
#include "palabos3D.hh"
#include <fstream>
#include <stdexcept>
#include <mpi.h>

using namespace hemo;

namespace {
// Run on every rank: the HemoCell count helpers perform MPI collectives.
void checkMachines(HemoCell& simulation) {
    const pluint tx = CellInformationFunctionals::getNumberOfCellsFromType(&simulation, "TX");
    const pluint rx = CellInformationFunctionals::getNumberOfCellsFromType(&simulation, "RX");
    if (tx != 1 || rx != 1) {
        throw std::runtime_error("Vessel requires exactly one complete TX and one complete RX. "
                                 "Check position files, wall clearance and checkpoint compatibility.");
    }
    hlog << "(Vessel) iteration=" << simulation.iter
         << " time=" << simulation.iter * param::dt << " s; TX=" << tx << " RX=" << rx
         << " RBC=" << CellInformationFunctionals::getNumberOfCellsFromType(&simulation, "RBC")
         << " PLT=" << CellInformationFunctionals::getNumberOfCellsFromType(&simulation, "PLT") << endl;
    unsigned long long wallContacts=0;
    MPI_Allreduce(&simulation.cellfields->wallContactCorrections,&wallContacts,1,
                  MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
    hlog << "(Vessel) wall-contact corrections since startup=" << wallContacts << endl;
}

void checkMachineInput(const std::string& role) {
    std::ifstream input(role + ".pos");
    int count = 0;
    T x, y, z, a, b, c;
    std::string extra;
    if (!(input >> count) || count != 1 || !(input >> x >> y >> z >> a >> b >> c)
        || (input >> extra)) {
        throw std::runtime_error(role + ".pos must contain a count of 1 and exactly one x y z rotationX rotationY rotationZ row.");
    }
}
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <configuration.xml>" << std::endl;
        return 1;
    }
    const bool membranesOnly = argc == 3 && std::string(argv[2]) == "--diagnose-membranes-only";
    if (argc > 2 && !membranesOnly) {
        std::cerr << "Unknown diagnostic option" << std::endl;
        return 1;
    }
    HemoCell hemocell(argv[1], argc, argv);
    Config& cfg = *hemocell.cfg;
    hemocell.outputInSiUnits = true;
    if (hemo::global.enableCEPACfield) {
        throw std::runtime_error("Vessel plasma-only glucose uses its own transport solver; keep enableCEPACfield=0.");
    }
    if (!cfg.checkpointed) {
        checkMachineInput("TX");
        checkMachineInput("RX");
    }
    const unsigned int tmax = cfg["sim"]["tmax"].read<unsigned int>();
    const unsigned int tmeas = cfg["sim"]["tmeas"].read<unsigned int>();
    const unsigned int tcsv = cfg["sim"]["tcsv"].read<unsigned int>();
    const unsigned int tcheckpoint = cfg["sim"]["tcheckpoint"].read<unsigned int>();
    const int materialEvery = cfg["ibm"]["stepMaterialEvery"].read<int>();
    const int particleEvery = cfg["ibm"]["stepParticleEvery"].read<int>();
    if (!tmeas || !tcsv || !tcheckpoint || materialEvery <= 0 || particleEvery <= 0
        || materialEvery % particleEvery != 0) {
        throw std::runtime_error("Output/update intervals must be positive; material interval must be divisible by particle interval.");
    }

    std::auto_ptr<MultiScalarField3D<int>> flagMatrix;
    std::auto_ptr<VoxelizedDomain3D<T>> voxelizedDomain;
    getFlagMatrixFromSTL(cfg["domain"]["geometry"].read<string>(),
        cfg["domain"]["fluidEnvelope"].read<int>(), cfg["domain"]["refDirN"].read<int>(),
        cfg["domain"]["refDir"].read<int>(), voxelizedDomain, flagMatrix,
        cfg["domain"]["blockSize"].read<int>(), cfg["domain"]["particleEnvelope"].read<int>());
    param::lbm_pipe_parameters(cfg, flagMatrix.get());
    param::printParameters();
    hemocell.lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(
        voxelizedDomain->getMultiBlockManagement(),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
        new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0 / param::tau));
    defineDynamics(*hemocell.lattice, *flagMatrix, hemocell.lattice->getBoundingBox(),
                   new BounceBack<T, DESCRIPTOR>(1.0), 0);
    hemocell.lattice->toggleInternalStatistics(false);
    hemocell.lattice->periodicity().toggleAll(false);
    hemocell.latticeEquilibrium(1., plb::Array<T,3>(0., 0., 0.));
    hemocell.lattice->initialize();
    hemocell.initializeCellfield();
    hemocell.cellfields->preserveMembranesAtWalls = true;

    // Independent material files permit later geometry/mechanics customization.
    for (const std::string name : {"RBC", "TX", "RX"}) {
        hemocell.addCellType<RbcHighOrderModel>(name, RBC_FROM_SPHERE);
        hemocell.setMaterialTimeScaleSeparation(name, materialEvery);
        hemocell.setInitialMinimumDistanceFromSolid(name, 0.5); // micrometres
    }
    hemocell.addCellType<PltSimpleModel>("PLT", ELLIPSOID_FROM_SPHERE);
    hemocell.setMaterialTimeScaleSeparation("PLT", materialEvery);
    hemocell.setParticleVelocityUpdateTimeScaleSeparation(particleEvery);
    hemocell.setRepulsion(cfg["domain"]["kRep"].read<T>(), cfg["domain"]["RepCutoff"].read<T>());
    hemocell.setRepulsionTimeScaleSeperation(materialEvery);
    const vector<int> outputs = {OUTPUT_POSITION, OUTPUT_TRIANGLES, OUTPUT_VELOCITY,
                                 OUTPUT_CELL_ID, OUTPUT_FORCE};
    for (const std::string name : {"RBC", "TX", "RX", "PLT"}) {
        hemocell.setOutputs(name, outputs);
    }
    hemocell.setFluidOutputs({OUTPUT_VELOCITY, OUTPUT_DENSITY, OUTPUT_BOUNDARY});
    hemocell.setSystemPeriodicity(0, true);
    if (cfg.checkpointed) hemocell.loadCheckPoint();
    else hemocell.loadParticles();
    checkMachines(hemocell);

    const T force = 8. * param::nu_lbm * (param::u_lbm_max * 0.5)
                    / param::pipe_radius / param::pipe_radius;
    const auto drive = [&]() {
        setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
            DESCRIPTOR<T>::ExternalField::forceBeginsAt, plb::Array<T,3>(force, 0., 0.));
    };
    drive();
    if (!cfg.checkpointed) {
        // Warm up only plasma, then output the initial machine positions at t=0.
        for (plint i = 0; i < cfg["parameters"]["warmup"].read<plint>(); ++i)
            hemocell.lattice->collideAndStream();
    }
    vessel::GlucoseTransport glucose(hemocell);
    if (membranesOnly) {
        hlog << "(Diagnostic) Membrane-only run: glucose evolution/output disabled; not a molecular simulation." << endl;
        glucose.checkMembranes();
    } else {
        glucose.initialize();
        if (cfg.checkpointed) glucose.loadCheckpoint();
    }
    hemocell.writeOutput();
    if (!membranesOnly) glucose.writeOutput();
    while (hemocell.iter < tmax) {
        if (!membranesOnly) glucose.advance(); // release and transport on the current plasma geometry
        hemocell.iterate(); // all four cell types use the same mobile IBM coupling
        drive();
        if (membranesOnly) glucose.checkMembranes();
        else glucose.updateMovingGeometry(); // conservatively evacuate newly covered nodes
        if (hemocell.iter % tmeas == 0) {
            checkMachines(hemocell);
            hemocell.writeOutput();
            if (!membranesOnly) glucose.writeOutput();
        }
        if (hemocell.iter % tcsv == 0) writeCellInfo_CSV(hemocell);
        if (!membranesOnly && hemocell.iter % tcheckpoint == 0) {
            hemocell.saveCheckPoint();
            glucose.saveCheckpoint();
        }
    }
    checkMachines(hemocell);
    if (hemocell.iter % tmeas != 0) {
        hemocell.writeOutput();
        if (!membranesOnly) glucose.writeOutput();
    }
    hlog << "(Vessel) Simulation finished." << endl;
    return 0;
}
