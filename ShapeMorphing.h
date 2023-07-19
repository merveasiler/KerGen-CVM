// @author Merve Asiler

#pragma once

#include "Mesh.h"

void morphByKernel(Mesh* sourceMesh, Mesh* targetMesh, vector<Mesh*>& interMeshes, double centerSource[3], double centerTarget[3]);

void morphByLerp(Mesh* sourceMesh, Mesh* targetMesh, vector<Mesh*>& interMeshes);


