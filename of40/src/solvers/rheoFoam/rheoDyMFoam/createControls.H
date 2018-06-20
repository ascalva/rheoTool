/**
 * author: Alberto Serrano - ascalva@gmail.com
 */

#include "createControl.H"
#include "createTimeControls.H"
#include "readTimeControls.H"

bool correctPhi
(
    mesh
    .solutionDict()
    .subDict("SIMPLE")
    .lookupOrDefault("correctPhi", false)
);

bool checkMeshCourantNo
(
    mesh
    .solutionDict()
    .subDict("SIMPLE")
    .lookupOrDefault("checkMeshCourantNo", false)
);
