// @author Merve Asiler

#include "SphericalParametrization.h"

void findBoundingBox(Mesh& mesh, double min[3], double max[3]) {

	for (int k = 0; k < 3; k++) {
		min[k] = numeric_limits<double>::infinity();
		max[k] = -numeric_limits<double>::infinity();
	}

	for (int i = 0; i < mesh.getNumOfVerts(); i++) {
		Vertex v = mesh.getVertex(i);
		for (int k = 0; k < 3; k++) {
			if (v.coords[k] < min[k])
				min[k] = v.coords[k];
			if (v.coords[k] > max[k])
				max[k] = v.coords[k];
		}
	}

}

void computeColor(double coords[3], double min[3], double max[3], double color[3]) {

	double rgb_length[3];
	for (int k = 0; k < 3; k++) {
		rgb_length[k] = max[k] - min[k];
		color[k] = (coords[k] - min[k]) / rgb_length[k];
	}

}

void parametrizeByKernel(Mesh& mesh, Mesh& sphericalMesh, double center[3], double radius[1], int resolution) {

	double min[3], max[3];
	findBoundingBox(mesh, min, max);
	radius[0] = (max[0] - min[0]) / 2;

	// project mesh vertices onto the sphere
	for (int i = 0; i < mesh.getNumOfVerts(); i++) {
		Vertex v = mesh.getVertex(i);
		double rayDirection[3];
		for (int j = 0; j < 3; j++)
			rayDirection[j] = v.coords[j] - center[j];
		normalize(rayDirection);
		// new coordinates:
		double newcoords[3];
		for (int j = 0; j < 3; j++)
			newcoords[j] = center[j] + rayDirection[j] * radius[0];
		sphericalMesh.addVertex(newcoords[0], newcoords[1], newcoords[2]);
		double color[3];
		computeColor(v.coords, min, max, color);
		sphericalMesh.addVertexColor(i, color);
		mesh.addVertexColor(i, color);
	}

	
	for (int i = 0; i < mesh.getNumOfTris(); i++) {
		Triangle t = mesh.getTriangle(i);
		sphericalMesh.addTriangle(t.corners[0], t.corners[1], t.corners[2]);
	}
	
	/*
	for (int i = 0; i < mesh->getNumOfEdges(); i++) {
		Edge e = mesh->getEdge(i);
		Vertex v1 = mesh->getVertex(e.endVerts[0]);
		Vertex v2 = mesh->getVertex(e.endVerts[1]);
		double* diff = diffVects(v2.coords, v1.coords);
		double length = computeLength(diff);
		normalize(diff);
		double unitLength = length / resolution;
		double startCoords[3] = { v1.coords[0], v1.coords[1], v1.coords[2] };
		int startVertexId = v1.idx, endVertexId;
		for (int j = 0; j < resolution; j++) {

			double endCoords[3], rayDirection[3];
			for (int k = 0; k < 3; k++) {
				endCoords[k] = startCoords[k] + diff[k] * unitLength;	// original coords
				rayDirection[k] = endCoords[k] - center[k];
			}
			normalize(rayDirection);
			for (int k = 0; k < 3; k++)
				endCoords[k] = center[k] + rayDirection[k] * radius;	// spherical coords

			if (j == resolution - 1)
				endVertexId = v2.idx;
			else {
				sphericalMesh->addVertex(endCoords[0], endCoords[1], endCoords[2]);
				endVertexId = sphericalMesh->getNumOfVerts() - 1;
				double color[3];
				computeColor(endCoords, min, max, color);
				sphericalMesh->addVertexColor(endVertexId, color);
			}
			sphericalMesh->addEdge(startVertexId, endVertexId);

			startVertexId = endVertexId;
			for (int k = 0; k < 3; k++)
				startCoords[k] = endCoords[k];
		}
	}
	*/
	
	double width = max[0] - min[0] + 10.0;
	double height = 10.0;
	for (int i = 0; i < mesh.getNumOfVerts(); i++) {
		Vertex v = mesh.getVertex(i);
		double coords[3] = { v.coords[0] - width, v.coords[1], v.coords[2] + height};
		mesh.changeVertexCoords(i, coords);
	}
}