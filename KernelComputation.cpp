// @author Merve Asiler

#include "KernelComputation.h"
#include "SphericalParametrization.h"
#include "ShapeMorphing.h"
#include "KernelExpansion.h"
#include "KerGen_Ideal.h"
#include "KernelByCGAL.h"
#include "sdlp.h"
#include "SceneManager.h"
#include "CommonUtils.h"
#include "CGALUtils.h"

#include <boost/filesystem.hpp>
#include <ctime>
#include <fstream>

using namespace boost::filesystem;

/* ***************************** GLOBAL VARIABLES ***************************** */
/* */ double executionCount = 10;
/* */ bool qualityComparison = true;
/* **************************************************************************** */

void doExperimentForPaper(string meshName) {

	/************************************************ READ MESH ************************************************/
	Mesh mesh;
	if (meshName.substr(meshName.length() - 3, 3) == "off")
		mesh.loadOff(meshName.c_str());
	else
		mesh.loadObj(meshName.c_str());

	/*************************************** INGREDIENTS / PREPARATIONS  ***************************************/
		// ... for kernel computation
	vector<string> algoTypes{ "kernel_by_cgal", "kergen" };	// initial one is assumed to be groundTruth
	vector<Mesh> kernels;
	kernels.resize(algoTypes.size());	// number of experiments: 1 + 4 + 4
		// ... for positions of the shapes on the scene
	vector<double*> positions;
	Grid boundingBox(&mesh, NULL, NULL, NULL, -1);	// compute the bounding box of the shape
	double* minborder = boundingBox.getMinCoords();
	double* maxborder = boundingBox.getMaxCoords();
	double objectWidth[3], breakWidth[3];
	for (int j = 0; j < 3; j++) {
		objectWidth[j] = maxborder[j] - minborder[j];
		breakWidth[j] = objectWidth[j] / 2.0;
	}
	double sceneSizeUnit = objectWidth[0] > objectWidth[1] ? objectWidth[0] : objectWidth[1];
	double totalSceneSize = (sceneSizeUnit + sceneSizeUnit / 10.0) * (kernels.size() / 2);
		// ... for shape materials
	vector<tuple<tuple<Mesh*, MaterialSetting*>, tuple<Mesh*, MaterialSetting*>>> mesh_mat_sets;
	MaterialSetting* kernelMatSetting = new MaterialSetting(0, 0, 1, 0);
	MaterialSetting* meshMatSetting = new MaterialSetting(1, 1, 1, 0.5);
		// ... for printing results&statistics
	string tempFileName = "temp_file.txt";
	std::ofstream outputFile(tempFileName);

	/************************************************ ALGORITHMS ************************************************/
	for (int i = 0; i < algoTypes.size(); i++) {
		vector<KernelExpansion*> kernelExpansions;
		string algoType = algoTypes[i];
		double elapsedTime = 0;
		double volDiffPercentage = 0;
		double hausdorffDistances[3] = { 0, 0, 0 };

		// Choose algoType
		if (algoType == "kergen_ideal")
			kernelExpansions.push_back(new KerGen_Ideal(mesh));
		else if (algoType == "kernel_by_cgal")
			kernelExpansions.push_back(new KernelByCGAL(mesh));
		else;

		// Run
		kernels[i] = Run(kernelExpansions, mesh, outputFile, elapsedTime, algoType);
		//double angles[3] = { 0, PI/4.0, 0 };
		//kernels[i].rotate(angles);
		//mesh.rotate(angles);
		mesh_mat_sets.push_back(make_tuple(make_tuple(&kernels[i], kernelMatSetting), make_tuple(&mesh, meshMatSetting)));

		if (i == 0)
			positions.push_back(new double[3]{ 0, 0, 0 });	// on the corner scene
		else
			positions.push_back(new double[3]{ 1.5 * i * sceneSizeUnit, 0, 0 });

		/*
		if (i == 0)
			positions.push_back(new double[3]{ 0, -1.5 * (sceneSizeUnit / 2), 0 });	// on the corner scene
		else if (i < 3)
			positions.push_back(new double[3]{ 1.5 * i * sceneSizeUnit, 0, 0 });
		else
			positions.push_back(new double[3]{ 1.5 * (i - 2) * sceneSizeUnit, -1.5 * sceneSizeUnit, 0 });
		*/

		// Quality
		if (qualityComparison && kernels[i].getNumOfVerts() > 0 && i > 0)
			CompareKernelQuality(kernels[0], kernels[i], algoType, outputFile, volDiffPercentage, hausdorffDistances);

		kernelExpansions.clear();
	}

	// Print
	outputFile.close();
	std::ifstream inputFile(tempFileName);
	std::stringstream ss;
	ss << inputFile.rdbuf();
	cout << ss.str();
	inputFile.close();
	remove(tempFileName);

	/*************************************************** DRAW ***************************************************/
	if (kernels[0].getNumOfVerts() > 0)
		drawMultipleScenes(mesh_mat_sets, positions, totalSceneSize);

	/************************************************* CLEAN-UP *************************************************/
	for (int i = 0; i < algoTypes.size(); i++)
		delete[] positions[i];
	delete kernelMatSetting;
	delete meshMatSetting;

}

void ComputeKernel(string meshName, string algoType) {

	// Read mesh
	Mesh mesh;
	if (meshName.substr(meshName.length() - 3, 3) == "off")
		mesh.loadOff(meshName.c_str());
	else
		mesh.loadObj(meshName.c_str());

	// Ingredients:
	double elapsedTime = 0;
	double volDiffPercentage = 0;
	double hausdorffDistances[3] = { 0, 0, 0 };
	Mesh kernel;
	string tempFileName = "temp_file.txt";
	std::ofstream outputFile(tempFileName);

	// Compute kernel
	vector<KernelExpansion*> kernelExpansions;
	for (int i = 0; i < executionCount; i++) {
		if (algoType == "kergen_ideal")
			kernelExpansions.push_back(new KerGen_Ideal(mesh));
		else if (algoType == "kernel_by_cgal")
			kernelExpansions.push_back(new KernelByCGAL(mesh));
		else;
	}

	// Run
	kernel = Run(kernelExpansions, mesh, outputFile, elapsedTime, algoType);

	// Quality
	if (qualityComparison && kernel.getNumOfVerts() > 0) {
		KernelByCGAL kernelByCGAL(mesh);
		kernelByCGAL.expandKernel();
		Mesh groundTruth = kernelByCGAL.getKernel();
		CompareKernelQuality(groundTruth, kernel, algoType, outputFile, volDiffPercentage, hausdorffDistances);
	}

	// Print
	outputFile.close();
	std::ifstream inputFile(tempFileName);
	std::stringstream ss;
	ss << inputFile.rdbuf();
	cout << ss.str();
	inputFile.close();
	remove(tempFileName);

	// Draw
	if (kernel.getNumOfVerts() > 0) {
		MaterialSetting* kernelMatSetting = new MaterialSetting(0, 0, 1, 0);
		MaterialSetting* meshMatSetting = new MaterialSetting(1, 1, 1, 0.5);
		vector<tuple<Mesh*, MaterialSetting*>> mesh_mat_set = { make_tuple(&kernel, kernelMatSetting), make_tuple(&mesh, meshMatSetting) };
		drawMultipleMeshToScene(mesh_mat_set);
	}

}

void ComputeBatchKernel(string inputFolderName, string outputFolderName, string algoType) {

	// Detect the algo type
	string extension, statistics_file_name;
	statistics_file_name = "KernelResults_" + algoType + ".txt";
	extension = "_" + algoType + "_kernel.off";

	double avgTime_star = 0, avgTime_nonstar = 0;
	int numOfStarShapes = 0, numOfNonStarShapes = 0, numOfNonManifoldShapes = 0;
	double avgVolDiffPercentage = 0;				// for star-shapes only
	double avgHausdorffDistances[3] = { 0, 0, 0 };	// for star-shapes only

	// Read folder, fecth mesh names, compute kernels
	int mesh_rank_id = 0;
	path p(inputFolderName);
	for (auto i = directory_iterator(p); i != directory_iterator(); i++, mesh_rank_id++)
	{
		// Fetch the mesh name
		string meshName;
		if (!is_directory(i->path())) //we eliminate directories
			meshName = i->path().filename().string();
		else
			continue;

		// if this is a previously created kernel file
		if (meshName.length() > 11 && meshName.substr(meshName.length() - 11, 11) == "_kernel.off")
			continue;

		// Open the file to append statistics & results
		std::ofstream outputFile(outputFolderName + "/" + statistics_file_name, ios_base::app);
		cout << meshName << endl;

		// Read mesh
		Mesh mesh;
		if (meshName.substr(meshName.length() - 3, 3) == "off")
			mesh.loadOff((inputFolderName + "/" + meshName).c_str());
		else
			mesh.loadObj((inputFolderName + "/" + meshName).c_str());

		if (mesh.isManifold() == false) {
			outputFile << "NOT MANIFOLD: " << mesh_rank_id << ": " << meshName << " !" << endl;
			cout << "Not manifold: " << mesh_rank_id << ": " << meshName << " !" << endl;
			numOfNonManifoldShapes++;
			continue;
		}

		cout << "Processing: " << mesh_rank_id << ": " << meshName << endl;
		outputFile << meshName << endl;

		// Choose the method
		double elapsedTime = 0;
		Mesh kernel;
		vector<KernelExpansion*> kernelExpansions;
		for (int i = 0; i < executionCount; i++) {
			if (algoType == "batch_kergen_ideal")
				kernelExpansions.push_back(new KerGen_Ideal(mesh));
			else if (algoType == "batch_kernel_cgal")
				kernelExpansions.push_back(new KernelByCGAL(mesh));
			else;
		}

		// Run
		kernel = Run(kernelExpansions, mesh, outputFile, elapsedTime, algoType);

		// Quality
		if (qualityComparison && kernel.getNumOfVerts() > 0) {
			KernelByCGAL kernelByCGAL(mesh);
			kernelByCGAL.expandKernel();
			Mesh groundTruth = kernelByCGAL.getKernel();
			double volDiffPercentage;
			double hausdorffDistances[3];
			CompareKernelQuality(groundTruth, kernel, algoType, outputFile, volDiffPercentage, hausdorffDistances);
			avgVolDiffPercentage += volDiffPercentage;
			for (int k = 0; k < 3; k++)
				avgHausdorffDistances[k] += hausdorffDistances[k];
		}

		// Output notes
		if (kernel.getNumOfVerts() > 0) {
			kernel.writeOff(outputFolderName + "/" + meshName.substr(0, meshName.length() - 4) + extension);
			avgTime_star += elapsedTime;
			numOfStarShapes++;
		}
		else {
			avgTime_nonstar += elapsedTime;
			numOfNonStarShapes++;
		}

		outputFile.close();
	}

	std::ofstream outputFile(outputFolderName + "/" + statistics_file_name, ios_base::app);
	outputFile << endl << endl;
	outputFile << "AVERAGE kernel computation has been completed in " << avgTime_star / numOfStarShapes << " milisecond(s) for " << numOfStarShapes << " star-shapes." << endl;
	outputFile << "AVERAGE kernel computation has been completed in " << avgTime_nonstar / numOfNonStarShapes << " milisecond(s) for " << numOfNonStarShapes << " nonstar-shapes." << endl;
	outputFile << "AVERAGE percentage of volume differences in computed kernels is " << avgVolDiffPercentage / numOfStarShapes << "% for " << numOfStarShapes << " star-shapes." << endl;
	outputFile << "AVERAGE hausdorff distance as 'grountTruth -> computed' is " << avgHausdorffDistances[0] / numOfStarShapes << " for " << numOfStarShapes << " star-shapes." << endl;
	outputFile << "AVERAGE hausdorff distance as 'computed -> grountTruth' is " << avgHausdorffDistances[1] / numOfStarShapes << " for " << numOfStarShapes << " star-shapes." << endl;
	outputFile << "AVERAGE hausdorff distance as 'symmetric' is " << avgHausdorffDistances[2] / numOfStarShapes << " for " << numOfStarShapes << " star-shapes." << endl;
	outputFile << "There has not been done kernel computation operation for " << numOfNonManifoldShapes << " non-manifold shapes." << endl;
	outputFile.close();

}

Mesh Run(vector<KernelExpansion*>& kernelExpansions, Mesh& mesh, std::ofstream& outputFile, double& elapsedTime, string algoType) {

	// Execute by <executionCount>-many times
	Mesh kernel;
	double totalTime = 0;
	for (int i = 0; i < kernelExpansions.size(); i++) {
		KernelExpansion* kernelExpansion = kernelExpansions[i];
		clock_t begin = clock();
		kernelExpansion->expandKernel();
		kernel = kernelExpansion->getKernel();
		clock_t end = clock();
		totalTime += double(end - begin) * 1000.0 / CLOCKS_PER_SEC;
		delete kernelExpansion;
	}

	// Print the statistics
	outputFile << "Mesh: [faces: " << mesh.getNumOfTris() << "], [edges: " << mesh.getNumOfEdges() << "], [vertices: " << mesh.getNumOfVerts() << "]" << endl;
	if (kernel.getNumOfVerts() > 0)
		outputFile << "Kernel: [faces: " << kernel.getNumOfTris() << "], [edges: " << kernel.getNumOfEdges() << "], [vertices: " << kernel.getNumOfVerts() << "]" << endl;
	else
		outputFile << "Kernel is empty!" << endl;

	// Print the average time
	elapsedTime = totalTime / executionCount;
	outputFile << "Kernel computation has been completed in " << elapsedTime << " milisecond(s) by KerGen." << endl;

	return kernel;

}

void CompareKernelQuality(Mesh groundTruth, Mesh kernel, string algoType, std::ofstream& outputFile, double& volDiffPercentage, double* hausdorffDistances) {

	double groundTruthVolume = groundTruth.computeVolume();
	double volume = 0;
	if (kernel.getNumOfVerts() > 0)
		volume = kernel.computeVolume();
	//if (volume < EPSILON)
	//	outputFile << "comparison is invalid due to ZERO volume!" << endl;

	outputFile << "Volume of the kernel computed by <" << algoType << "> : " << volume << " out of " << groundTruthVolume << "." << endl;
	double volumeDiff = groundTruthVolume - volume;
	if (abs(volumeDiff) < EPSILON)
		volumeDiff = 0;
	volDiffPercentage = (100.0 * abs(volumeDiff)) / groundTruthVolume;
	outputFile << ">>>> Volume difference : " << volumeDiff << endl;
	outputFile << ">>>> Percentage of difference : " << volDiffPercentage << "%" << endl;

	double* hd;
	if (volume < EPSILON) {
		hd = new double[3];
		for (int k = 0; k < 3; k++)
			hd[k] = 0;
	}
	else
		hd = computeHausdorffDistance(groundTruth, kernel);

	for (int k = 0; k < 3; k++) {
		if (hd[k] < EPSILON)
			hd[k] = 0;
		hausdorffDistances[k] = hd[k];
	}
	outputFile << "Hausdorff distances : " << endl;
	outputFile << "\t\t    grountTruth -> computed : " << hd[0] << endl;
	outputFile << "\t\t    computed -> groundTruth : " << hd[1] << endl;
	outputFile << "\t\t    symmetric : " << hd[2] << " ." << endl << endl;
	delete[] hd;

}

void FindKernelPoint_SDLP(string meshName) {

	Mesh* mesh = new Mesh();
	if (meshName.substr(meshName.length() - 3, 3) == "off")
		mesh->loadOff(meshName.c_str());
	else
		mesh->loadObj(meshName.c_str());

	double extremeDirection[3] = { 0, 0, 1 };

	double* kernel_point = sdlpMain(*mesh, extremeDirection);
	if (kernel_point != NULL) {
		cout << "Final kernel point: " << kernel_point[0] << " " << kernel_point[1] << " " << kernel_point[2] << endl;
		delete[] kernel_point;
	}
	delete mesh;
}

void SphericalParametrize(string meshName) {

	Mesh mesh;
	if (meshName.substr(meshName.length() - 3, 3) == "off")
		mesh.loadOff(meshName.c_str());
	else
		mesh.loadObj(meshName.c_str());

	int resolution = 10;
	double radius[1] = { 1.0 };
	
	double center[3];
	computeCenterOfKernel(mesh, center);

	/*
	double extremeDirection[3] = { 0, 0, 1 };
	double* kernelPoint = sdlpMain(mesh, extremeDirection);
	double* center = kernelPoint;
	*/

	Mesh sphericalMesh;
	parametrizeByKernel(mesh, sphericalMesh, center, radius, resolution);
	drawMeshOnSphere(&sphericalMesh, &mesh, center, radius[0]);

}

void ShapeMorphByKernel(string sourceMeshName, string targetMeshName) {

	Mesh* sourceMesh = new Mesh();
	Mesh* targetMesh = new Mesh();

	if (sourceMeshName.substr(sourceMeshName.length() - 3, 3) == "off")
		sourceMesh->loadOff(sourceMeshName.c_str());
	else
		sourceMesh->loadObj(sourceMeshName.c_str());
	
	if (targetMeshName.substr(targetMeshName.length() - 3, 3) == "off")
		targetMesh->loadOff(targetMeshName.c_str());
	else
		targetMesh->loadObj(targetMeshName.c_str());

	int numOfInterMeshes = 5;
	vector<Mesh*> interMeshes;
	for (int i = 0; i < numOfInterMeshes; i++)
		interMeshes.push_back(new Mesh());

	double centerSource[3] = { 0, 0, 0 };
	double centerTarget[3] = { 0, 0, 0 };
	//double centerSource[3], centerTarget[3];
	//computeCenterOfKernel(sourceMesh, centerSource);
	//computeCenterOfKernel(targetMesh, centerTarget);
	
	morphByKernel(sourceMesh, targetMesh, interMeshes, centerSource, centerTarget);

	vector<tuple<Mesh*, MaterialSetting*>> outputs;
	MaterialSetting* meshMatSetting = new MaterialSetting(1, 1, 1, 0);
	outputs.push_back(make_tuple(sourceMesh, meshMatSetting));
	for (int i = 0; i < interMeshes.size(); i++)
		outputs.push_back(make_tuple(interMeshes[i], meshMatSetting));

	drawMultipleMeshToScene(outputs);

	delete sourceMesh;
	delete targetMesh;	
}

void ShapeMorphByLerp(string sourceMeshName, string targetMeshName) {

	Mesh* sourceMesh = new Mesh();
	Mesh* targetMesh = new Mesh();

	if (sourceMeshName.substr(sourceMeshName.length() - 3, 3) == "off")
		sourceMesh->loadOff(sourceMeshName.c_str());
	else
		sourceMesh->loadObj(sourceMeshName.c_str());

	if (targetMeshName.substr(targetMeshName.length() - 3, 3) == "off")
		targetMesh->loadOff(targetMeshName.c_str());
	else
		targetMesh->loadObj(targetMeshName.c_str());

	int numOfInterMeshes = 5;
	vector<Mesh*> interMeshes;
	for (int i = 0; i < numOfInterMeshes; i++)
		interMeshes.push_back(new Mesh());

	morphByLerp(sourceMesh, targetMesh, interMeshes);

	vector<tuple<Mesh*, MaterialSetting*>> outputs;
	MaterialSetting* meshMatSetting = new MaterialSetting(1, 1, 1, 0);
	outputs.push_back(make_tuple(sourceMesh, meshMatSetting));
	for (int i = 0; i < interMeshes.size(); i++)
		outputs.push_back(make_tuple(interMeshes[i], meshMatSetting));

	drawMultipleMeshToScene(outputs);

	delete sourceMesh;
	delete targetMesh;
}

void computeCenterOfKernel(Mesh& mesh, double center[3]) {

	KernelExpansion* kernelExpansion(new KerGen_Ideal(mesh));
	kernelExpansion->expandKernel();
	Mesh kernel = kernelExpansion->getKernel();

	for (int k = 0; k < 3; k++)
		center[k] = 0;
	for (int i = 0; i < kernel.getNumOfVerts(); i++) {
		for (int k = 0; k < 3; k++)
			center[k] += kernel.getVertex(i).coords[k];
	}
	for (int k = 0; k < 3; k++)
		center[k] /= kernel.getNumOfVerts();

}

vector<double> produceColorSource(Mesh* ground_truth, Mesh* exp_mesh) {

	vector<double> distances;
	double maxDistance = -numeric_limits<double>::infinity();

	for (int i = 0; i < ground_truth->getNumOfVerts(); i++) {

		Vertex v1 = ground_truth->getVertex(i);		
		double minDistance = numeric_limits<double>::infinity();

		for (int j = 0; j < exp_mesh->getNumOfVerts(); j++) {
	
			Vertex v2 = exp_mesh->getVertex(j);
			double* diff = diffVects(v1.coords, v2.coords);
			double distance = computeLength(diff);
			delete[] diff;

			if (distance < minDistance)
				minDistance = distance;
		}

		if (minDistance < EPSILON)
			minDistance = 0;
		distances.push_back(minDistance);
		if (minDistance > maxDistance)
			maxDistance = minDistance;
	}

	if (!(maxDistance < EPSILON)) {
		// normalize
		for (int i = 0; i < distances.size(); i++)
			distances[i] /= maxDistance;
	}

	return distances;
}

