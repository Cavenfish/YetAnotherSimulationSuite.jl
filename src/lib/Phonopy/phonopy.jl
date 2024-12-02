
# Communication between phonopy and JMD goes through the phonopy yaml file.
# Python should be used for all non-force calculations.

function phonopy_addForces(yamlFile, forces)
  phonopy    = pyimport("phonopy")
  obj        = phonopy.load(yamlFile)
  obj.forces = forces
  obj.save()
end