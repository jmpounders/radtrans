cmake_minimum_required (VERSION 2.8)
project ( Transport CXX )

include ( LocalConfig.cmake )

set ( transport_SRC main.cpp
                    associatedlegendre.cpp
                    dataset.cpp
                    dofobj.cpp
		    element.cpp
                    fixedsource.cpp
                    global.cpp
                    hdf5interface.cpp
		    inputparser.cpp
                    logpolicy.cpp
                    log.cpp
		    legendre.cpp
		    material.cpp
                    mathematics.cpp
		    mesh.cpp
                    meshfactory.cpp
                    moabmesh.cpp
		    node.cpp
		    output.cpp
		    outputgenerator.cpp
                    perfstats.cpp
                    solutionmanager.cpp
		    solverbase.cpp
                    solverlocalmoc.cpp
                    solverregmoc.cpp
                    sweeper.cpp
                    timing.cpp
                    transient.cpp
                    transient_ndadaptive.cpp
                    transient_uts.cpp
		    transportproblem.cpp)

include_directories ( ../include ${hdfPath}/include ${moabPath}/include ${mpiPath}/include)
link_directories( ${hdfPath}/lib ${moabPath}/lib ${mpiPath}/lib ${perfTools}/lib )

add_executable ( Transport ${transport_SRC} )
target_link_libraries( Transport MOAB mpi mpicxx hdf5_cpp hdf5)
