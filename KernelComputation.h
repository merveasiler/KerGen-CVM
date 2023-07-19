// @author Merve Asiler

#pragma once
#pragma comment(lib, "boost_filesystem-vc140-mt.lib")

#include <string>
#include <vector>
using namespace std;

class Mesh;
class KernelExpansion;

void doExperimentForPaper(string meshName);

void ComputeKernel(string meshName, string algoType);

void ComputeBatchKernel(string inputFolderName, string outputFolderName, string algoType);

Mesh Run(vector<KernelExpansion*>& kernelExpansions, Mesh& mesh, std::ofstream& outputFile, double& elapsedTime, string algoType);

void CompareKernelQuality(Mesh groundTruth, Mesh kernel, string algoType, std::ofstream& outputFile, double& volDiffPercentage, double* hausdorffDistances);

void FindKernelPoint_SDLP(string meshName);

void SphericalParametrize(string meshName);

void ShapeMorphByKernel(string sourceMeshName, string targetMeshName);

void ShapeMorphByLerp(string sourceMeshName, string targetMeshName);

void computeCenterOfKernel(Mesh& mesh, double center[3]);

vector<double> produceColorSource(Mesh* ground_truth, Mesh* exp_mesh);



