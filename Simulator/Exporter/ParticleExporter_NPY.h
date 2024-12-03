#ifndef __ParticleExporter_NPY_h__
#define __ParticleExporter_NPY_h__

#include "ExporterBase.h"
#include "npy.hpp"
#include <fstream>

namespace SPH
{
	/** \brief Particle exporter for the NPY format.
	*/
	class ParticleExporter_NPY : public ExporterBase
	{
	protected: 
		std::string m_exportPath;
		std::ofstream *m_outfile;
		std::vector<std::string> m_attributes;
        std::vector<std::string> m_fieldNames;
        bool m_sim2D;
        unsigned int m_numFields;


		void createParticleFile(const std::string& fileName, FluidModel* model);
		void writeParticles(const std::string& fileName, FluidModel* model);
		void writeBoundaryParticles(const std::string& fileName, BoundaryModel* model);

	public:
		ParticleExporter_NPY(SimulatorBase *base);
		ParticleExporter_NPY(const ParticleExporter_NPY&) = delete;
		ParticleExporter_NPY& operator=(const ParticleExporter_NPY&) = delete;
		virtual ~ParticleExporter_NPY(void);

		virtual void init(const std::string& outputPath);
		virtual void step(const unsigned int frame);
		virtual void reset();
		virtual void setActive(const bool active); 
	};
}

#endif
