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
#include "fluidInfo.h"
#include "particleInfo.h"
#include "writeCellInfoCSV.h"
#include <fenv.h>

#include "palabos3D.h"
#include "palabos3D.hh"

#include <cmath>

typedef double T;

using namespace hemo;

int main(int argc, char *argv[]) {
  if(argc < 2) {
      std::cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  HemoCell hemocell(argv[1], argc, argv);

  hemocell.outputInSiUnits = false; // Temporary: output raw lattice values for CEPAC validation

  Config * cfg = hemocell.cfg;


  hlogfile << "(PipeFlow) (Geometry) reading and voxelizing STL file " << (*cfg)["domain"]["geometry"].read<string>() << endl;

  std::auto_ptr<MultiScalarField3D<int>> flagMatrix;
  std::auto_ptr<VoxelizedDomain3D<T>> voxelizedDomain;
  getFlagMatrixFromSTL((*cfg)["domain"]["geometry"].read<string>(),
                       (*cfg)["domain"]["fluidEnvelope"].read<int>(),
                       (*cfg)["domain"]["refDirN"].read<int>(),
                       (*cfg)["domain"]["refDir"].read<int>(),
                       voxelizedDomain, flagMatrix,
                       (*cfg)["domain"]["blockSize"].read<int>(),
                       (*cfg)["domain"]["particleEnvelope"].read<int>());

  param::lbm_pipe_parameters((*cfg),flagMatrix.get());
  param::tau_CEPAC = 1.0;
  param::printParameters();

  hemocell.lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(
            voxelizedDomain.get()->getMultiBlockManagement(),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
            new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0/param::tau));

  defineDynamics(*hemocell.lattice, *flagMatrix.get(), (*hemocell.lattice).getBoundingBox(), new BounceBack<T, DESCRIPTOR>(1.), 0);

  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.lattice->periodicity().toggleAll(false);
  hemocell.latticeEquilibrium(1.,plb::Array<T, 3>(0.,0.,0.));

  //Driving Force
  T poiseuilleForce =  8 * param::nu_lbm * (param::u_lbm_max * 0.5) / param::pipe_radius / param::pipe_radius;

  hemocell.lattice->initialize();

  //Adding all the cells
  hemocell.initializeCellfield();

  // Impermeable vessel wall for CEPAC scalar.
  // The same flagMatrix used for the fluid identifies solid nodes with flag 0.
  defineDynamics(*hemocell.cellfields->CEPACfield, *flagMatrix.get(), hemocell.cellfields->CEPACfield->getBoundingBox(), new BounceBack<T, CEPAC_DESCRIPTOR>(1.0), 0);

  hemocell.addCellType<RbcHighOrderModel>("RBC", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("RBC", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setInitialMinimumDistanceFromSolid("RBC", 0.5); //Micrometer! not LU

  // Synthetic cell: mechanically identical to RBC
  hemocell.addCellType<RbcHighOrderModel>("SYN", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("SYN", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setInitialMinimumDistanceFromSolid("SYN", 0.5);

  hemocell.addCellType<PltSimpleModel>("PLT", ELLIPSOID_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("PLT", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());

  // plint cx = hemocell.cellfields->CEPACfield->getNx() / 4;
  // plint cy = hemocell.cellfields->CEPACfield->getNy() / 2;
  // plint cz = hemocell.cellfields->CEPACfield->getNz() / 2;

  // Small initial concentration blob near tube centreline
  // Box3D CEPACsource(cx, cx + 1, cy - 1, cy + 1, cz - 1, cz + 1);

  // Entire scalar field starts at zero
  plb::initializeAtEquilibrium(*hemocell.cellfields->CEPACfield, hemocell.cellfields->CEPACfield->getBoundingBox(), 1.0, plb::Array<T,3>(0.0, 0.0, 0.0));

  hemocell.cellfields->CEPACfield->initialize();

  hemocell.setParticleVelocityUpdateTimeScaleSeparation((*cfg)["ibm"]["stepParticleEvery"].read<int>());

  //hemocell.setRepulsion((*cfg)["domain"]["kRep"].read<T>(), (*cfg)["domain"]["RepCutoff"].read<T>());
  //hemocell.setRepulsionTimeScaleSeperation((*cfg)["ibm"]["stepMaterialEvery"].read<int>());

  vector<int> outputs = {OUTPUT_POSITION,OUTPUT_TRIANGLES,OUTPUT_FORCE,OUTPUT_FORCE_VOLUME,OUTPUT_FORCE_BENDING,OUTPUT_FORCE_LINK,OUTPUT_FORCE_AREA,OUTPUT_FORCE_VISC};
  hemocell.setOutputs("RBC", outputs);
  hemocell.setOutputs("SYN", outputs);
  hemocell.setOutputs("PLT", outputs);

  outputs = {OUTPUT_VELOCITY,OUTPUT_DENSITY,OUTPUT_FORCE};
  hemocell.setFluidOutputs(outputs);

  outputs = {OUTPUT_DENSITY};
  hemocell.setCEPACOutputs(outputs);

  // Turn on periodicity in the X direction
  hemocell.setSystemPeriodicity(0, true);
  hemocell.cellfields->CEPACfield->periodicity().toggle(0, true);

  // Enable boundary particles
  //hemocell.enableBoundaryParticles((*cfg)["domain"]["kRep"].read<T>(), (*cfg)["domain"]["BRepCutoff"].read<T>(),(*cfg)["ibm"]["stepMaterialEvery"].read<int>());

  //loading the cellfield
  if (not cfg->checkpointed) {
    hemocell.loadParticles();
    hemocell.writeOutput();
  } else {
    hemocell.loadCheckPoint();
  }

  //Restructure atomic blocks on processors when possible
  //hemocell.doRestructure(false); // cause errors

  if (hemocell.iter == 0) {
    hlog << "(PipeFlow) fresh start: warming up cell-free fluid domain for "  << (*cfg)["parameters"]["warmup"].read<plint>() << " iterations..." << endl;
	setExternalVector(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
                    DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                    plb::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));
    for (plint itrt = 0; itrt < (*cfg)["parameters"]["warmup"].read<plint>(); ++itrt) {
      hemocell.lattice->collideAndStream();
    }
  }

  unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();
  unsigned int tcheckpoint = (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();
  unsigned int tcsv = (*cfg)["sim"]["tcsv"].read<unsigned int>();

  // Temporary moving-source test
  const unsigned int releaseEvery = 10;
  const T sourceDensity = 1.05;

  hlog << "(PipeFlow) Starting simulation..." << endl;

  while (hemocell.iter < tmax ) {
    hemocell.iterate();

	// Temporary CEPAC source following the SYN cell
	if (hemocell.iter % releaseEvery == 0) {

		std::map<int, CellInformation> cellInfo;
		CellInformationFunctionals::calculateCellInformation(
		    &hemocell,
		    cellInfo
		);

		pluint synType = (*hemocell.cellfields)["SYN"]->ctype;

		for (const auto& pair : cellInfo) {

		    const CellInformation& cinfo = pair.second;

		    if (cinfo.centerLocal &&
		        cinfo.cellType == synType) {

		        // SYN center is already expressed in lattice units
		        plint sx = (plint)std::round(cinfo.position[0]);
		        plint sy = (plint)std::round(cinfo.position[1]);
		        plint sz = (plint)std::round(cinfo.position[2]);

		        // Small 3x3x3 source surrounding the SYN center
		        Box3D movingSource(
		            sx - 1, sx + 1,
		            sy - 1, sy + 1,
		            sz - 1, sz + 1
		        );

		        plb::initializeAtEquilibrium(
		            *hemocell.cellfields->CEPACfield,
		            movingSource,
		            sourceDensity,
		            plb::Array<T,3>(0.0, 0.0, 0.0)
		        );
		    }
		}
	}


    //Set driving force as required after each iteration
    setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                plb::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));

    if (hemocell.iter % tmeas == 0) {

        // Get current information for all cells
	std::map<int, CellInformation> cellInfo;
	CellInformationFunctionals::calculateCellInformation(&hemocell, cellInfo);

	// Numeric type-ID corresponding to SYN
	pluint synType = (*hemocell.cellfields)["SYN"]->ctype;

	// Find the SYN cell and print its current center
	for (const auto& pair : cellInfo) {

	    const CellInformation& cinfo = pair.second;

	    if (cinfo.centerLocal && cinfo.cellType == synType) {

		hlog << "SYN center = ("
		     << cinfo.position[0] << ", "
		     << cinfo.position[1] << ", "
		     << cinfo.position[2] << ")"
		     << endl;
	    }
	}


        hlog << "(main) Stats. @ " <<  hemocell.iter << " (" << hemocell.iter * param::dt << " s):" << endl;
        hlog << "\t # of cells: " << CellInformationFunctionals::getTotalNumberOfCells(&hemocell);
        hlog << " | # of RBC: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "RBC");
        hlog << ", SYN: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "SYN") << endl;
        hlog << ", PLT: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "PLT") << endl;
        FluidStatistics finfo = FluidInfo::calculateVelocityStatistics(&hemocell); T toMpS = param::dx / param::dt;
        hlog << "\t Velocity  -  max.: " << finfo.max * toMpS << " m/s, mean: " << finfo.avg * toMpS<< " m/s, rel. app. viscosity: " << (param::u_lbm_max*0.5) / finfo.avg << endl;
        ParticleStatistics pinfo = ParticleInfo::calculateForceStatistics(&hemocell); T topN = param::df * 1.0e12;
        hlog << "\t Force  -  min.: " << pinfo.min * topN << " pN, max.: " << pinfo.max * topN << " pN (" << pinfo.max << " lf), mean: " << pinfo.avg * topN << " pN" << endl;

        // Additional useful stats, if needed
        //finfo = FluidInfo::calculateForceStatistics(&hemocell);
        //Set force as required after this function;
        // setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
        //           DESCRIPTOR<T>::ExternalField::forceBeginsAt,
        //           hemo::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));
        // pcout << "Fluid force, Minimum: " << finfo.min << " Maximum: " << finfo.max << " Average: " << finfo.avg << endl;
        // ParticleStatistics pinfo = ParticleInfo::calculateVelocityStatistics(&hemocell);
        // pcout << "Particle velocity, Minimum: " << pinfo.min << " Maximum: " << pinfo.max << " Average: " << pinfo.avg << endl;
        hemocell.writeOutput();
    }

    if (hemocell.iter % tcsv == 0) {
      hlog << "Saving simple mean cell values to CSV at timestep " << hemocell.iter << endl;
      writeCellInfo_CSV(hemocell);
    }

    if (hemocell.iter % tcheckpoint == 0) {
       hemocell.saveCheckPoint();
    }
  }

  hlog << "(main) Simulation finished :) " << endl;

  return 0;
}
