#include "ParticleExporter_NPY.h"
#include <Utilities/Logger.h>
#include <Utilities/FileSystem.h>
#include "SPlisHSPlasH/Simulation.h"
#include "Simulator/SceneConfiguration.h"
#include <regex>

using namespace SPH;
using namespace Utilities;

ParticleExporter_NPY::ParticleExporter_NPY(SimulatorBase *base) :
	ExporterBase(base)
{
	m_outfile = nullptr;
    Simulation* sim = Simulation::getCurrent();
    m_sim2D = sim->is2DSimulation();
    m_numFields = 3 + 3 + 3 + 3 + 6;
}

ParticleExporter_NPY::~ParticleExporter_NPY(void)
{
}

void ParticleExporter_NPY::init(const std::string& outputPath)
{
    std::string sceneFile = SceneConfiguration::getCurrent()->getSceneFile();
    Real pr = SceneConfiguration::getCurrent()->getScene().particleRadius;
    sceneFile = FileSystem::getFileName(sceneFile);
	m_exportPath = FileSystem::normalizePath(outputPath + "/output/" + sceneFile + "/npy_" + std::to_string(pr));
}

void ParticleExporter_NPY::step(const unsigned int frame)
{
	if (!m_active)
		return;

    std::cout << "Writing frame " << frame << std::endl;

    Simulation* sim = Simulation::getCurrent();
	m_sim2D = sim->is2DSimulation();
    Real pr = sim->getParticleRadius();

    std::string frame_string = std::to_string(frame);
    int num_zeros = 4 - frame_string.length();
    frame_string.insert(0, num_zeros, '0');

	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		FluidModel* model = sim->getFluidModel(i);
        const unsigned int numParticles = model->numActiveParticles();
        
        Real particle_radius = sim->getParticleRadius();
		std::string fileName = "Particle";
		if (!m_base->getValue<bool>(SimulatorBase::EXPORT_OBJECT_SPLITTING))
		{
            
			fileName = fileName + "_" + model->getId() + "_" + frame_string;
			std::string exportFileName = FileSystem::normalizePath(m_exportPath + "/" + fileName);
			writeParticles(exportFileName + ".npy", model);
		}
		else
		{
			std::cout << "Error! Not Implemented error" << std::endl;
            std::exit(0);
		}
	}

    for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++) 
    {
        BoundaryModel* model = sim->getBoundaryModel(i);

        std::string fileName = "Particle_Boundary_" + std::to_string(i);
        fileName += "_" + frame_string;
        std::string exportFileName = FileSystem::normalizePath(m_exportPath + "/" + fileName);
		writeBoundaryParticles(exportFileName + ".npy", model);
    }

}

void ParticleExporter_NPY::reset()
{
}

void ParticleExporter_NPY::setActive(const bool active)
{
	ExporterBase::setActive(active);
	if (m_active)
		FileSystem::makeDirs(m_exportPath);
}

void ParticleExporter_NPY::createParticleFile(const std::string& fileName, FluidModel* model)
{
}

void ParticleExporter_NPY::writeParticles(const std::string& fileName, FluidModel* model)
{
	const unsigned int numParticles = model->numActiveParticles();
    const unsigned int offset = 12;
    // pos, vel, acc_v, acc_p, mass, volume, density, density advected, aii, pressure

    Real* data = (Real *) std::malloc(sizeof(Real) * numParticles * m_numFields);

    const FieldDescription& aiiField                    = model->getField("a_ii");
    const FieldDescription& pressureField               = model->getField("pressure");
    const FieldDescription& pressureAccelerationField   = model->getField("pressure_acceleration");
    const FieldDescription& densityAdvectedField        = model->getField("advected density");

    for (int i = 0; i < numParticles; i++) {
        const auto& pos_i       = model->getPosition(i);
        const auto& vel_i       = model->getVelocity(i);
        const auto& acc_v_i     = model->getAcceleration(i);
        Vector3r acc_p_i((Real*)pressureAccelerationField.getFct(i));
        
        const auto& mass_i      = model->getMass(i);
        const auto& volume_i    = model->getVolume(i);
        const auto& rho_i       = model->getDensity(i);
        Real density_Advected_i = *((Real*)densityAdvectedField.getFct(i));
        Real aii_i              = *((Real*)aiiField.getFct(i));
        Real p_i                = *((Real*)pressureField.getFct(i));

        for (int j = 0; j < 3; j++) {
            data[i * m_numFields + j + 0] = pos_i[j];
            data[i * m_numFields + j + 3] = vel_i[j];
            data[i * m_numFields + j + 6] = acc_v_i[j];
            data[i * m_numFields + j + 9] = acc_p_i[j];
        }
        
        data[i * m_numFields + offset + 0] = mass_i;
        data[i * m_numFields + offset + 1] = volume_i;
        data[i * m_numFields + offset + 2] = rho_i;
        data[i * m_numFields + offset + 3] = density_Advected_i;
        data[i * m_numFields + offset + 4] = aii_i;
        data[i * m_numFields + offset + 5] = p_i;
    }

#ifdef USE_DOUBLE
        npy::npy_data_ptr<double> d;
        d.data_ptr = data;
        d.shape = {numParticles, m_numFields};
        d.fortran_order = false; // optional
#else
        npy::npy_data_ptr<float> d;
        d.data_ptr = data;
        d.shape = {numParticles, m_numFields};
        d.fortran_order = false; // optional
#endif

        // const std::string path{filename};
        npy::write_npy(fileName, d);
    
    free(data);
    data = NULL;

}

void ParticleExporter_NPY::writeBoundaryParticles(const std::string& fileName, BoundaryModel* model)
{
    Simulation* sim = Simulation::getCurrent();
    if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012) {

        
        BoundaryModel_Akinci2012 *bm_neighbor = static_cast<BoundaryModel_Akinci2012*>(model);
        const unsigned int numParticles = bm_neighbor->numberOfParticles(); 
        unsigned int numBoundaryFields = 7; // pos, vel, volume
        Real* data = (Real *) std::malloc(sizeof(Real) * numParticles * numBoundaryFields);

        for (int i = 0; i < numParticles; i++) {
            const auto& pos_i       = bm_neighbor->getPosition(i);
            const auto& vel_i       = bm_neighbor->getVelocity(i);
            const auto& vol_i       = bm_neighbor->getVolume(i);

            for (int j = 0; j < 3; j++) {
                data[i * numBoundaryFields + j + 0] = pos_i[j];
                data[i * numBoundaryFields + j + 3] = vel_i[j];
            }
            data[i * numBoundaryFields + 6] = vol_i;   
        }

    #ifdef USE_DOUBLE
            npy::npy_data_ptr<double> d;
            d.data_ptr = data;
            d.shape = {numParticles, numBoundaryFields};
            d.fortran_order = false; // optional
    #else
            npy::npy_data_ptr<float> d;
            d.data_ptr = data;
            d.shape = {numParticles, numBoundaryFields};
            d.fortran_order = false; // optional
    #endif

            // const std::string path{filename};
            npy::write_npy(fileName, d);
        
        free(data);
        data = NULL;
    }
    else {
        std::cout << "Unsupported boundary discretization for exporter" << std::endl;
        std::exit(0);
    }
}
