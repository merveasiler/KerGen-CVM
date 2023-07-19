// @author Merve Asiler

#pragma once

#include "Painter.h"
#include "Scene.h"
#include "Mesh.h"

// SCENE DRAWING MANAGER METHODs

void drawMeshToScene(Mesh* mesh);

void drawDoubleMeshToScene(Mesh* mesh1, Mesh* mesh2, vector<double> colorSource);

void drawMeshToScene(string meshName);

void drawRotatedMeshToScene(string meshName);

void drawMultipleMeshToScene(vector<tuple<Mesh*, MaterialSetting*>> mesh_mat_set);

void drawNextToNextSketchToScene(vector<tuple<tuple<Mesh*, MaterialSetting*>, tuple<Mesh*, MaterialSetting*>>> mesh_mat_tuple_set);

void drawMultipleScenes(vector<tuple<tuple<Mesh*, MaterialSetting*>, tuple<Mesh*, MaterialSetting*>>> outputs, vector<double*> positions, double sceneSize);

void drawSegmentationToScene(Mesh* mesh, vector<double*>& coveringNodes, vector<int*>& trisIds, vector<int>& numOfTrisPerNode, float radius);

void drawMeshOnSphere(Mesh* mesh, Mesh* originalMesh, double center[3], float radius);
