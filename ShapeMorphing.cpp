// @author Merve Asiler

#include "ShapeMorphing.h"

double computeWidth(Mesh* mesh) {

	double min_x = numeric_limits<double>::infinity();
	double max_x = -numeric_limits<double>::infinity();
	for (int j = 0; j < mesh->getNumOfVerts(); j++) {
		if (mesh->getVertex(j).coords[0] < min_x)
			min_x = mesh->getVertex(j).coords[0];
		if (mesh->getVertex(j).coords[0] > max_x)
			max_x = mesh->getVertex(j).coords[0];
	}
	
	return max_x - min_x;

}

void morphByKernel(Mesh* sourceMesh, Mesh* targetMesh, vector<Mesh*>& interMeshes, double centerSource[3], double centerTarget[3]) {

	double meshWidth = computeWidth(sourceMesh);
	double breakWidth = meshWidth;

	double time_step = 1.0 / (interMeshes.size() + 1);
	// compute translation first
	double translation[3];
	for (int i = 0; i < 3; i++)
		translation[i] = centerTarget[i] - centerSource[i];
	interMeshes.push_back(new Mesh());	// recompute targetMesh

	// define source and target mesh vertices w.r.t. center of their kernels
	vector <double*> distance_vectors_source, distance_vectors_target;
	vector<double> angles;
	for (int i = 0; i < sourceMesh->getNumOfVerts(); i++) {
		Vertex vs = sourceMesh->getVertex(i);
		Vertex vt = targetMesh->getVertex(i);

		double* ray_dir_source = new double[4];
		double* ray_dir_target = new double[4];
		for (int j = 0; j < 3; j++) {
			ray_dir_source[j] = vs.coords[j] - centerSource[j];
			ray_dir_target[j] = vt.coords[j] - centerTarget[j];
		}

		ray_dir_source[3] = computeLength(ray_dir_source);	// only for the first 3 elements
		ray_dir_target[3] = computeLength(ray_dir_target);	// only for the first 3 elements
		normalize(ray_dir_source);							// only for the first 3 elements
		normalize(ray_dir_target);							// only for the first 3 elements

		distance_vectors_source.push_back(ray_dir_source);
		distance_vectors_target.push_back(ray_dir_target);
		
		double cosAngle = dotProduct(ray_dir_source, ray_dir_target);
		angles.push_back(acos(cosAngle));
		//angles.push_back(0);
	}

	// interpolate
	for (int m = 0; m < interMeshes.size(); m++) {
		Mesh* mesh = interMeshes[m];
		double time = time_step * (m + 1);

		for (int i = 0; i < sourceMesh->getNumOfVerts(); i++) {
			// angle
			double* sourceDir = distance_vectors_source[i];
			double* targetDir = distance_vectors_target[i];
			double angle = angles[i];

			double sourceCoeff = sin((1.0 - time) * angle) / sin(angle);
			double targetCoeff = sin(time * angle) / sin(angle);

			double point_on_sphere[3], direction[3];
			for (int k = 0; k < 3; k++) {
				point_on_sphere[k] = (centerSource[k] + sourceDir[k]) * sourceCoeff + (centerTarget[k] + targetDir[k]) * targetCoeff;
				direction[k] = point_on_sphere[k] - centerTarget[k];
			}
			normalize(direction);

			double totalDistanceDiff = distance_vectors_target[i][3] - distance_vectors_source[i][3];
			double distance = time * totalDistanceDiff;

			double center[3], point[3];
			for (int k = 0; k < 3; k++) {
				center[k] = centerSource[k] + time * translation[k];
				point[k] = center[k] + direction[k] * (distance + distance_vectors_source[i][3]);
			}
			mesh->addVertex(point[0] + ((meshWidth + breakWidth) * (m + 1)), point[1], point[2]);
		}
	
		for (int t = 0; t < sourceMesh->getNumOfTris(); t++) {
			Triangle triangle = sourceMesh->getTriangle(t);
			mesh->addTriangle(triangle.corners[0], triangle.corners[1], triangle.corners[2]);
		}
	}

}

void morphByLerp(Mesh* sourceMesh, Mesh* targetMesh, vector<Mesh*>& interMeshes) {

	double meshWidth = computeWidth(sourceMesh);
	double breakWidth = meshWidth;
	double time_step = 1.0 / (interMeshes.size() + 1);
	interMeshes.push_back(new Mesh());	// recompute targetMesh

	for (int m = 0; m < interMeshes.size(); m++) {
		Mesh* mesh = interMeshes[m];
		double time = time_step * (m + 1);

		for (int i = 0; i < sourceMesh->getNumOfVerts(); i++) {
			Vertex vs = sourceMesh->getVertex(i);
			Vertex vt = targetMesh->getVertex(i);

			double distance[3];
			for (int k = 0; k < 3; k++)
				distance[k] = vt.coords[k] - vs.coords[k];

			double point[3];
			for (int k = 0; k < 3; k++)
				point[k] = vs.coords[k] + distance[k] * time;

			mesh->addVertex(point[0] + ((meshWidth + breakWidth) * (m + 1)), point[1], point[2]);
		}

		for (int t = 0; t < sourceMesh->getNumOfTris(); t++) {
			Triangle triangle = sourceMesh->getTriangle(t);
			mesh->addTriangle(triangle.corners[0], triangle.corners[1], triangle.corners[2]);
		}
	}

}