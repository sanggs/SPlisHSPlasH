import pysplishsplash as sph
import numpy as np
import meshio
import sys
import os
import json

class CustomExporter(sph.ExporterBase):
    def __init__(self, base):
        sph.ExporterBase.__init__(self, base)
        self.exportPath = ""
        self.firstFrame = True
        
    def init(self, outputPath): 
        self.exportPath = os.path.join(outputPath, "rigid_bodies")
        
    def step(self, frame): 
        if (not self.getActive()):
            return;
        
        fileName = "rb_data_"
        exportFileName = os.path.join(self.exportPath, fileName)
        self.writeRigidBodies(exportFileName, frame)
        
    def reset(self): 
        self.firstFrame = True
        return
        
    def setActive(self, active): 
        sph.ExporterBase.setActive(self,active)
        if (self.getActive()):
            os.makedirs(self.exportPath, exist_ok = True)
            
    def writeRigidBodies(self, fileName, frame):
        sim = sph.Simulation.getCurrent()
        nBoundaryModels = sim.numberOfBoundaryModels()

        # check if we have a static model
        isStatic = True
        for i in range(0, nBoundaryModels):
            bm = sim.getBoundaryModel(i)
            if (bm.getRigidBodyObject().isDynamic()):
                isStatic = False
                break
                
           
        if (self.firstFrame or not isStatic):
            for i in range(0, nBoundaryModels):
                bm = sim.getBoundaryModel(i)
                rbo = bm.getRigidBodyObject()
                vertices = np.array(rbo.getVertexBuffer(), copy=False)
                faces = np.array(rbo.getFaceBuffer(), copy=False)
                tris = faces.reshape((-1,3))
                cells = [("triangle", tris)]
                
                meshio.write_points_cells(fileName + str(i) + "_" + str(frame) + ".obj", vertices, cells)
        self.firstFrame = False

    
def main():
    base = sph.Exec.SimulatorBase()
    
    base.init(sys.argv, "[Python] SPlisHSPlasH")

    gui = sph.GUI.Simulator_GUI_imgui(base)
    base.setGui(gui)
    
    base.activateExporter("NPY Exporter", True)	
        
    # init the simulation
    base.initSimulation()

    base.runSimulation()

    export_info = {}
    sim = sph.Simulation.getCurrent()
    scene_file_path = sph.Exec.SceneConfiguration.getSceneFile(sph.Exec.SceneConfiguration.getCurrent())
    scene_file_name = sph.Utilities.FileSystem.getFileName(scene_file_path)
    output_path = os.path.join(base.getOutputPath(), 'output')
    nBoundaryModels = sim.numberOfBoundaryModels()

    # check if we have a static model
    isStatic = True
    for i in range(0, nBoundaryModels):
        bm = sim.getBoundaryModel(i)
        if (bm.getRigidBodyObject().isDynamic()):
            isStatic = False
            break
    assert(isStatic)
    
    global_min = None
    global_max = None

    for i in range(0, nBoundaryModels):
        bm = sim.getBoundaryModel(i)
        bp_num = bm.numberOfParticles()

        particles = []
        for j in range(bp_num):
            particles.append(bm.getPosition(j))
        
        particles = np.array(particles)

        if (global_min is None):
            global_min = np.min(particles, axis=0)
            global_max = np.max(particles, axis=0)
        else:
            local_min = np.min(particles, axis=0)
            local_max = np.max(particles, axis=0)
            for k in range(3):
                global_min[k] = min(global_min[k], local_min[k])
                global_max[k] = max(global_max[k], local_max[k])
    
    pr = sim.getParticleRadius()
    pr_str = f'{pr:8.6f}'

    assert(global_min is not None)
    assert(global_max is not None)

    export_info['folder_name'] = f'npy_{pr_str}'
    export_info['particle_radius'] = float(f'{pr_str}')
    export_info['bounds'] = [global_min.tolist(), global_max.tolist()]
    
    info_file = os.path.join(output_path, scene_file_name)
    info_file = os.path.join(info_file, f'info_{pr_str}.json')

    with open(info_file, 'w') as f:
        json.dump(export_info, f, indent=True)
    
    base.cleanup()

if __name__ == "__main__":
	main()
