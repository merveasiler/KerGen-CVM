// @author Merve Asiler

#include "SceneManager.h"
#include "CommonUtils.h"

#include <string>

void drawMeshToScene(string meshName) {

	Mesh* mesh = new Mesh;
	if (meshName.substr(meshName.length() - 3, 3) == "off")
		mesh->loadOff(meshName.c_str());
	else
		mesh->loadObj(meshName.c_str());

	drawMeshToScene(mesh);
	delete mesh;

}

void drawMeshToScene(Mesh* mesh) {

	Scene* scene = new Scene();
	Painter* painter = new Painter();
	SoSeparator* res = new SoSeparator();
	painter->getShapeSep(mesh, res);
	painter->drawTriangulation(mesh, res);
	scene->makeScene(res);
	delete scene;
	delete painter;

}

void drawRotatedMeshToScene(string meshName) {

	Mesh* mesh = new Mesh;
	if (meshName.substr(meshName.length() - 3, 3) == "off")
		mesh->loadOff(meshName.c_str());
	else
		mesh->loadObj(meshName.c_str());

	double angle_X = 0.0 * (PI / 180.0);
	double angle_Y = 0.0 * (PI / 180.0);
	double angle_Z = 195.0 * (PI / 180.0);

	// Rotate mesh
	Mesh* rotatedMesh = new Mesh();;
	double rotation_X[3][3] = { {1, 0, 0}, {0, cos(angle_X), -sin(angle_X)}, {0, sin(angle_X), cos(angle_X)} };
	double rotation_Y[3][3] = { {cos(angle_Y), 0, sin(angle_Y)}, {0, 1, 0}, {-sin(angle_Y), 0, cos(angle_Y)} };
	double rotation_Z[3][3] = { {cos(angle_Z), -sin(angle_Z), 0}, {sin(angle_Z), cos(angle_Z), 0}, {0, 0, 1} };
	double R[3][3] = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				R[i][k] += rotation_Y[i][j] * rotation_Z[j][k];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				R[i][k] += rotation_X[i][j] * R[j][k];

	for (int i = 0; i < mesh->getNumOfVerts(); i++) {
		double r[3] = { 0, 0, 0 };
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				r[j] += R[j][k] * mesh->getVertex(i).coords[k];
		rotatedMesh->addVertex(r[0], r[1], r[2]);
	}
	for (int i = 0; i < mesh->getNumOfTris(); i++)
		rotatedMesh->addTriangle(mesh->getTriangle(i).corners[0], mesh->getTriangle(i).corners[1], mesh->getTriangle(i).corners[2]);

	// Write to file
	ofstream outFile;
	outFile.open(meshName.substr(0, meshName.length() - 4) + "_rotated.off");
	outFile << "OFF\n";
	outFile << mesh->getNumOfVerts() << " " << mesh->getNumOfTris() << " 0\n";
	for (int i = 0; i < mesh->getNumOfVerts(); i++)
		outFile << rotatedMesh->getVertex(i).coords[0] << " " << rotatedMesh->getVertex(i).coords[1] << " " << rotatedMesh->getVertex(i).coords[2] << "\n";
	for (int i = 0; i < mesh->getNumOfTris(); i++)
		outFile << "3 " << rotatedMesh->getTriangle(i).corners[0] << " " << rotatedMesh->getTriangle(i).corners[1] << " " << rotatedMesh->getTriangle(i).corners[2] << "\n";
	outFile.close();


	Scene* scene = new Scene();
	Painter* painter = new Painter();
	SoSeparator* res = new SoSeparator();
	painter->getShapeSep(rotatedMesh, res);
	painter->drawTriangulation(rotatedMesh, res);
	scene->makeScene(res);


	delete scene;
	delete painter;
	delete mesh;

}

void drawDoubleMeshToScene(Mesh* mesh1, Mesh* mesh2, vector<double> colorSource) {

	Scene* scene = new Scene();
	Painter* painter = new Painter();

	SoSeparator* res1 = new SoSeparator();
	painter->getColorfulShapeSep(mesh1, res1, colorSource);
	SoSeparator* res2 = new SoSeparator();
	painter->getShapeSep(mesh2, res2);
	//painter->drawTriangulation(mesh2, res2);
	
	vector<SoSeparator*> resSet;
	resSet.push_back(res1);
	resSet.push_back(res2);
	scene->makeScene(resSet);

	delete scene;
	delete painter;

}

void drawMultipleMeshToScene(vector<tuple<Mesh*, MaterialSetting*>> mesh_mat_set) {

	Scene* scene = new Scene();
	Painter* painter = new Painter();

	vector<SoSeparator*> resSet;
	for (int i = 0; i < mesh_mat_set.size(); i++) {
		Mesh* mesh;
		MaterialSetting* mat;
		tie(mesh, mat) = mesh_mat_set[i];
		SoSeparator* res = new SoSeparator();
		/*
		if (i == 0)
			painter->getShapeSepByEdges(mesh, mat, res);
		else {
		*/
			painter->getShapeSep(mesh, mat, res);
			//if (i > 0)
			painter->drawTriangulation(mesh, res);
		//}
		resSet.push_back(res);
	}

	scene->makeScene(resSet);

	resSet.clear();
	delete scene;
	delete painter;

}

void drawNextToNextSketchToScene(vector<tuple<tuple<Mesh*, MaterialSetting*>, tuple<Mesh*, MaterialSetting*>>> mesh_mat_tuple_set) {
	
	Scene* scene = new Scene();
	Painter* painter = new Painter();

	vector<SoSeparator*> resSet;
	for (int i = 0; i < mesh_mat_tuple_set.size(); i++) {
		tuple<Mesh*, MaterialSetting*> tuple1, tuple2;
		tie(tuple1, tuple2) = mesh_mat_tuple_set[i];

		Mesh* mesh1, * mesh2;
		MaterialSetting* mat1, * mat2;
		tie(mesh1, mat1) = tuple1;
		tie(mesh2, mat2) = tuple2;

		SoSeparator* res1 = new SoSeparator();
		painter->getShapeSep(mesh1, mat1, res1);
		SoSeparator* res2 = new SoSeparator();
		painter->getShapeSep(mesh2, mat2, res2);

		resSet.push_back(res1);
		resSet.push_back(res2);
	}
	scene->makeScene(resSet);

	resSet.clear();
	delete scene;
	delete painter;

}

void drawMultipleScenes(vector<tuple<tuple<Mesh*, MaterialSetting*>, tuple<Mesh*, MaterialSetting*>>> outputs, vector<double*> positions, double sceneSize) {

	Scene* scene = new Scene();
	vector<SoSeparator*> resSets;
	SoSeparator* cornerResSet = NULL;
	Painter* painter = new Painter();
	
	// SCENES
	for (int i = 0; i < outputs.size(); i++) {
		tuple<tuple<Mesh*, MaterialSetting*>, tuple<Mesh*, MaterialSetting*>> kernel_mesh_tuple = outputs[i];
		tuple<Mesh*, MaterialSetting*> kernel_mat_set, mesh_mat_set;
		tie(kernel_mat_set, mesh_mat_set) = kernel_mesh_tuple;

		Mesh* mesh[2];
		MaterialSetting* mesh_mat[2];
		tie(mesh[0], mesh_mat[0]) = kernel_mat_set;
		tie(mesh[1], mesh_mat[1]) = mesh_mat_set;

		SoSeparator* resSet = new SoSeparator;
		if (i == 0)
			cornerResSet = resSet;
		else
			resSets.push_back(resSet);
		
		// compute position
		SoTransform* transform = new SoTransform();
		transform->translation.setValue(positions[i][0], positions[i][1], positions[i][2]);
		resSet->addChild(transform);
		
		for (int j = 0; j < 2; j++) {
			SoSeparator* res = new SoSeparator();
			//if (i==1 && j==0)
			//	painter->getShapeSepByEdges(mesh[j], mesh_mat[j], res);
			//else
				painter->getShapeSep(mesh[j], mesh_mat[j], res);
				if (j > 0)
					painter->drawTriangulation(mesh[j], res);
			resSet->addChild(res);
		}
	}

	scene->makeMultipleScene(resSets, cornerResSet, sceneSize);

	// clean-up
	for (int i = 0; i < outputs.size()-1; i++)
		resSets[i]->removeAllChildren();
	resSets.clear();
	delete painter;
	delete scene;

}

// Each coveringNodes consists of a 3D coordinate and represents a segment 
// Each element of trisIds consist of a set of triangles belonging to corresponding segment in coveringNodes
// Each element of numOfTrisPerNodes corresponds to the size of the corresponding triangle list in trisIds
void drawSegmentationToScene(Mesh* mesh, vector<double*>& coveringNodes, vector<int*>& trisIds, vector<int>& numOfTrisPerNode, float radius) {
	Scene* scene = new Scene();
	Painter* painter = new Painter();
	SoSeparator* res = new SoSeparator();

	painter->getShapeSep(mesh, coveringNodes, trisIds, numOfTrisPerNode, radius, res);
	scene->makeScene(res);

	delete scene;
	delete painter;
}

void drawMeshOnSphere(Mesh* mesh, Mesh* originalMesh, double center[3], float radius) {
	Scene* scene = new Scene();
	Painter* painter = new Painter();
	SoSeparator* res1 = new SoSeparator();
	SoSeparator* res2 = new SoSeparator();
	SoSeparator* res3 = new SoSeparator();
	
	vector<double> colorSource;
	
	painter->getShapeSep(mesh, res2);
	painter->getShapeSep(originalMesh, res1);
	painter->drawTriangulation(mesh, res2);
	painter->drawTriangulation(originalMesh, res1);
	
	//painter->getColorfulShapeSep(mesh, res2, colorSource);
	//painter->getColorfulShapeSep(originalMesh, res1, colorSource);
	//painter->drawColorfulTriangulation(mesh, res2);
	//painter->drawColorfulTriangulation(originalMesh, res1);
	//painter->getShapeSep(center, radius, res2);

	vector<SoSeparator*> resSet;
	resSet.push_back(res2);
	resSet.push_back(res1);
	resSet.push_back(res3);
	scene->makeScene(resSet);

	resSet.clear();
	delete scene;
	delete painter;
}
