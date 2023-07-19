// @author Merve Asiler

#include "KernelExpansion.h"
#include "BaseGeoOpUtils.h"
#include "sdlp.h"

KernelExpansion::KernelExpansion(const Mesh& hostMesh, bool defineInitialPoint) {

	this->hostMeshptr = &hostMesh;

	computeHalfSpacesFromTriangles(hostMeshptr->getAllTris(), hostMeshptr->getAllVerts(), this->halfSpaceSet);

	// FIND EXTREME KERNEL BOUNDING CORNERS BY USING ALL OF THE 6 EXTREME DIRECTIONS, TAKE THEIR AVERAGE
	extremeCorners[0] = new double[3];	extremeCorners[1] = new double[3];
	double extremeDirections[6][3] = { {-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1} };
	for (int i = 0; i < 6; i++) {
		this->extremePoints[i] = nullptr;
		this->extremePoints[i] = sdlpMain(extremeDirections[i], halfSpaceSet);	// compute initial kernel point at the given extreme direction
		if (!this->extremePoints[i])
			return;
		this->extremeCorners[i % 2][int(i / 2)] = this->extremePoints[i][int(i / 2)];
	}

	if (defineInitialPoint) {
		this->initialPoint = new double[3]{ 0, 0, 0 };
		for (int i = 0; i < 3; i++) {
			for (int p = 0; p < 6; p++)
				this->initialPoint[i] += this->extremePoints[p][i];
			this->initialPoint[i] /= 6;
		}
			//this->initialPoint[i] = (this->extremeCorners[0][i] + this->extremeCorners[1][i]) / 2.0;	// compute initial kernel point	
	}
	else
		this->initialPoint = nullptr;

}

KernelExpansion::KernelExpansion(const Mesh& hostMesh) {

	this->hostMeshptr = &hostMesh;
	this->initialPoint = nullptr;
	this->extremeCorners[0] = nullptr;	this->extremeCorners[1] = nullptr;
	computeHalfSpacesFromTriangles(hostMesh.getAllTris(), hostMesh.getAllVerts(), this->halfSpaceSet);

}

KernelExpansion::~KernelExpansion() {

	if (initialPoint) {
		delete[] initialPoint;
		initialPoint = nullptr;
	}

	if (extremeCorners[0]) {
		delete[] extremeCorners[0];
		extremeCorners[0] = nullptr;
		delete[] extremeCorners[1];
		extremeCorners[1] = nullptr;
	}

	for (int i = 0; i < 6; i++) {
		if (extremePoints[i]) {
			delete[] extremePoints[i];
			extremePoints[i] = nullptr;
		}
	}

	for (int i = 0; i < handledCells.size(); i++)
		delete[] handledCells[i];
	handledCells.clear();
	handledGrids.clear();

	halfSpaceSet.clear();

}

Mesh& KernelExpansion::getKernel() {
	return kernel;
}

double* KernelExpansion::getInitialKernelPoint() {
	return initialPoint;
}

vector<HalfSpace>& KernelExpansion::getHalfSpaceSet() {
	return halfSpaceSet;
}

Grid& KernelExpansion::getGrid() {
	return grid;
}

void KernelExpansion::decideOnNeighbor2(int neigh_i, int neigh_j, int neigh_k, double neighVert1_scalar, double neighVert2_scalar, double neighVert3_scalar, double neighVert4_scalar) {
	/*
	if (isCellValidForQueue(neigh_i, neigh_j, neigh_k)) {

		// if any vertex of the given face is a kernel point, then we should check this neighbor cell sharing the face
		if (neighVert1_scalar <= 0 ||
			neighVert2_scalar <= 0 ||
			neighVert3_scalar <= 0 ||
			neighVert4_scalar <= 0) {	// scalarValue is located.

			int* neigh_ijk = new int[3]{ neigh_i, neigh_j, neigh_k };
			possibleKernelCells.push(neigh_ijk);
			handledCells.push_back(neigh_ijk);
		}
	}
	*/
}

void KernelExpansion::decideOnNeighbor(int neigh_i, int neigh_j, int neigh_k,
	double coords[8][4], double* scalars[8],
	int ind1, int ind2, int ind3, int ind4,
	Grid* grid) {

	if (isCellValidForQueue(neigh_i, neigh_j, neigh_k, grid)) {

		// if any vertex of the given face is a kernel point, then we should check this neighbor cell sharing the face
		if (coords[ind1][3] <= 0 || coords[ind2][3] <= 0 || coords[ind3][3] <= 0 || coords[ind4][3] <= 0) {

			int ind[4] = { ind1, ind2, ind3, ind4 };
			double kernelPoint[3];
			for (int i = 0; i < 4; i++)
				if (coords[ind[i]][3] <= 0) {
					for (int j = 0; j < 3; j++)
						kernelPoint[j] = coords[ind[i]][j];
					break;
				}

			int* neigh_ijk = new int[3]{ neigh_i, neigh_j, neigh_k };
			fillContainers(new Grid(*grid), neigh_ijk, kernelPoint);
		}
		/*
		else {

			int ind[4] = { ind1, ind2, ind3, ind4};
			vector<int> excluders[4];
			for (int i = 0; i < 4; i++)
				excluders[i] = detectExcluderHalfSpaces(scalars[ind[i]]);

			for (int a = 0; a < excluders[0].size(); a++)
				for (int b = 0; b < excluders[1].size(); b++)
					if (excluders[0][a] == excluders[1][b])
						for (int c = 0; c < excluders[2].size(); c++)
							if (excluders[0][a] == excluders[2][c])
								for (int d = 0; d < excluders[3].size(); d++)
									if (excluders[0][a] == excluders[3][d])
										return;

			double kernelPoint[3] = { coords[ind1][0], coords[ind1][1], coords[ind1][2] };
			int* neigh_ijk = new int[3]{ neigh_i, neigh_j, neigh_k };
			fillContainers(new Grid(*grid), neigh_ijk, kernelPoint);
		}
		*/
	}

}

bool KernelExpansion::isCellValidForQueue(int  neighbor_i, int neighbor_j, int neighbor_k, Grid* grid) {

	int neighbor_ijk[3] = { neighbor_i, neighbor_j, neighbor_k };
	/*
		for (int i = 0; i < 3; i++)
			if (neighbor_ijk[i] >= grid->numOfCells[i] || neighbor_ijk[i] < 0)	// there is no such a neighbor because of invalid index
				return false;
	*/
	// did we previously check this cell?

	for (int i = 0; i < handledCells.size(); i++) {
		if (handledGrids[i] == grid->gridId &&
			handledCells[i][0] == neighbor_ijk[0] && handledCells[i][1] == neighbor_ijk[1] && handledCells[i][2] == neighbor_ijk[2])
			return false;
	}

	return true;
}

void KernelExpansion::fillContainers(Grid* grid, int* initialCell, double* startPoint) {

	if (!initialCell)
		initialCell = grid->findHomeCell(startPoint);	// find the cell which the initial point belongs to and append it into the queue
	possibleKernelCells.push(initialCell);					// this is a queue consisting of cells whose corners will be checked for being inside/outside the kernel
	ownerGridForCells.push(grid);
	kernelPointsForCells.push(new double[3]{ startPoint[0], startPoint[1], startPoint[2] });
	handledCells.push_back(initialCell);
	handledGrids.push_back(grid->gridId);

}

bool KernelExpansion::isInKernel(double* scalarVector) {
	for (int i = 0; i < halfSpaceSet.size(); i++)
		if (scalarVector[i] > EPSILON_3)
			return false;	// out
	return true;		// in
}

vector<int> KernelExpansion::detectExcluderHalfSpaces(double* scalarVector) {
	vector<int> excluderIds;
	for (int i = 0; i < halfSpaceSet.size(); i++)
		if (scalarVector[i] > EPSILON_3)
			excluderIds.push_back(i);
	return excluderIds;
}

double KernelExpansion::findIdealCelSize_OuterToInner() {

	double leastCoordinates[3] = { extremeCorners[0][0], extremeCorners[0][1], extremeCorners[0][2] };
	double mostCoordinates[3] = { extremeCorners[1][0], extremeCorners[1][1], extremeCorners[1][2] };

	double coords[8][3] = { {leastCoordinates[0], leastCoordinates[1], leastCoordinates[2]},
							{mostCoordinates[0], leastCoordinates[1], leastCoordinates[2]},
							{mostCoordinates[0], leastCoordinates[1], mostCoordinates[2]},
							{leastCoordinates[0], leastCoordinates[1], mostCoordinates[2]},
							{leastCoordinates[0], mostCoordinates[1], leastCoordinates[2]},
							{mostCoordinates[0], mostCoordinates[1], leastCoordinates[2]},
							{mostCoordinates[0], mostCoordinates[1], mostCoordinates[2]},
							{leastCoordinates[0], mostCoordinates[1], mostCoordinates[2]}
	};

	double closestPoints[8][3];

	for (int k = 0; k < 3; k++)
		closestPoints[0][k] = (extremePoints[0][k] + extremePoints[2][k] + extremePoints[4][k]) / 3.0;
	for (int k = 0; k < 3; k++)
		closestPoints[1][k] = (extremePoints[1][k] + extremePoints[2][k] + extremePoints[4][k]) / 3.0;
	for (int k = 0; k < 3; k++)
		closestPoints[2][k] = (extremePoints[1][k] + extremePoints[2][k] + extremePoints[5][k]) / 3.0;
	for (int k = 0; k < 3; k++)
		closestPoints[3][k] = (extremePoints[0][k] + extremePoints[2][k] + extremePoints[5][k]) / 3.0;
	for (int k = 0; k < 3; k++)
		closestPoints[4][k] = (extremePoints[0][k] + extremePoints[3][k] + extremePoints[4][k]) / 3.0;
	for (int k = 0; k < 3; k++)
		closestPoints[5][k] = (extremePoints[1][k] + extremePoints[3][k] + extremePoints[4][k]) / 3.0;
	for (int k = 0; k < 3; k++)
		closestPoints[6][k] = (extremePoints[1][k] + extremePoints[3][k] + extremePoints[5][k]) / 3.0;
	for (int k = 0; k < 3; k++)
		closestPoints[7][k] = (extremePoints[0][k] + extremePoints[3][k] + extremePoints[5][k]) / 3.0;


	double min_cell_size = numeric_limits<double>::infinity();
	for (int i = 0; i < 8; i++) {
		double diffVect[3];
		for (int k = 0; k < 3; k++)
			diffVect[k] = closestPoints[i][k] - coords[i][k];

		normalize(diffVect);
		Line ray(coords[i], diffVect);
		double maxDistance = 0;
		bool valid = false;
		double* scalars = findValueVectorSatisfiedByPoint(coords[i], halfSpaceSet);
		for (int j = 0; j < halfSpaceSet.size(); j++) {
			if (scalars[j] > 0) {
				double t = findLinePlaneIntersection(ray, halfSpaceSet[j]);
				if (t > maxDistance) {
					maxDistance = t;
					valid = true;
				}
			}
		}
		delete[] scalars;

		if (valid && maxDistance < min_cell_size)
			min_cell_size = maxDistance;
	}

	return (min_cell_size / sqrt(3.0));
}

double KernelExpansion::findIdealCelSize_InnerToOuter() {

	double leastCoordinates[3] = { extremeCorners[0][0], extremeCorners[0][1], extremeCorners[0][2] };
	double mostCoordinates[3] = { extremeCorners[1][0], extremeCorners[1][1], extremeCorners[1][2] };

	double coords[8][3] = { {leastCoordinates[0], leastCoordinates[1], leastCoordinates[2]},
							{mostCoordinates[0], leastCoordinates[1], leastCoordinates[2]},
							{mostCoordinates[0], leastCoordinates[1], mostCoordinates[2]},
							{leastCoordinates[0], leastCoordinates[1], mostCoordinates[2]},
							{leastCoordinates[0], mostCoordinates[1], leastCoordinates[2]},
							{mostCoordinates[0], mostCoordinates[1], leastCoordinates[2]},
							{mostCoordinates[0], mostCoordinates[1], mostCoordinates[2]},
							{leastCoordinates[0], mostCoordinates[1], mostCoordinates[2]}
	};

	//double temp[8][3] = { {1, 1, 1}, {1, 1, -1}, {1, -1, 1}, {-1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}, {-1, -1, -1} };

	double min_cell_size = numeric_limits<double>::infinity();
	for (int i = 0; i < 8; i++) {
		double diffVect[3];
		for (int k = 0; k < 3; k++)
			diffVect[k] = coords[i][k] - initialPoint[k];
		//diffVect[k] = temp[i][k];

		normalize(diffVect);
		Line ray(initialPoint, diffVect);
		double minDistance = numeric_limits<double>::infinity();
		for (int j = 0; j < halfSpaceSet.size(); j++) {
			double t = findLinePlaneIntersection(ray, halfSpaceSet[j]);
			if (t == 0) {
				minDistance = numeric_limits<double>::infinity();
				break;
			}
			if (t > 0 && t < minDistance)
				minDistance = t;
		}

		if (minDistance < min_cell_size)
			min_cell_size = minDistance;
	}

	return cbrt(min_cell_size);

}

void KernelExpansion::findSensitiveIdealCelSizePerCoordinate(double* minIntervals) {

	double xyzVals[3][7];
	for (int d = 0; d < 3; d++) {
		for (int i = 0; i < 7; i++) {
			xyzVals[d][i] = numeric_limits<double>::infinity();
		}
	}

	double* basePoints[7];
	for (int i = 0; i < 6; i++)
		basePoints[i] = extremePoints[i];
	basePoints[6] = this->initialPoint;

	// increasingly order extreme point coordinates 
	for (int d = 0; d < 3; d++) {
		for (int i = 0; i < 7; i++) {
			for (int j = 0; j < 7; j++) {
				if (basePoints[i][d] < xyzVals[d][j]) {
					for (int k = 6; k > j; k--)	// shift
						xyzVals[d][k] = xyzVals[d][k - 1];
					xyzVals[d][j] = basePoints[i][d];
					break;
				}
			}
		}
	}

	double maxIntervals[3], avgIntervals[3] = { 0, 0, 0 };;
	for (int d = 0; d < 3; d++) {
		minIntervals[d] = numeric_limits<double>::infinity();
		maxIntervals[d] = -numeric_limits<double>::infinity();
	}

	for (int x = 0; x < 6; x++) {
		double interval = abs(xyzVals[0][x + 1] - xyzVals[0][x]);
		avgIntervals[0] += interval;
		if (interval > EPSILON) {
			if (interval < minIntervals[0])
				minIntervals[0] = interval;
			else if (interval > maxIntervals[0])
				maxIntervals[0] = interval;
		}
	}
	for (int y = 0; y < 6; y++) {
		double interval = abs(xyzVals[1][y + 1] - xyzVals[1][y]);
		avgIntervals[1] += interval;
		if (interval > EPSILON) {
			if (interval < minIntervals[1])
				minIntervals[1] = interval;
			else if (interval > maxIntervals[1])
				maxIntervals[1] = interval;
		}
	}
	for (int z = 0; z < 6; z++) {
		double interval = abs(xyzVals[2][z + 1] - xyzVals[2][z]);
		avgIntervals[2] += interval;
		if (interval > EPSILON) {
			if (interval < minIntervals[2])
				minIntervals[2] = interval;
			else if (interval > maxIntervals[2])
				maxIntervals[2] = interval;
		}
	}


	for (int d = 0; d < 3; d++)
		avgIntervals[d] /= 6;

	minIntervals[0] = avgIntervals[0] / 2;
	minIntervals[1] = avgIntervals[1] / 2;
	minIntervals[2] = avgIntervals[2] / 2;

	/*
	minIntervals[0] = maxIntervals[0];
	minIntervals[1] = maxIntervals[1];
	minIntervals[2] = maxIntervals[2];
	*/

}

void KernelExpansion::findIdealCelSizePerCoordinate(double* minIntervals) {

	double xyzVals[3][6];
	for (int d = 0; d < 3; d++) {
		for (int i = 0; i < 6; i++) {
			xyzVals[d][i] = numeric_limits<double>::infinity();
		}
	}

	// increasingly order extreme point coordinates 
	for (int d = 0; d < 3; d++) {
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				if (extremePoints[i][d] < xyzVals[d][j]) {
					for (int k = 5; k > j; k--)	// shift
						xyzVals[d][k] = xyzVals[d][k - 1];
					xyzVals[d][j] = extremePoints[i][d];
					break;
				}
			}
		}
	}

	double maxIntervals[3], avgIntervals[3] = { 0, 0, 0 };;
	for (int d = 0; d < 3; d++) {
		minIntervals[d] = numeric_limits<double>::infinity();
		maxIntervals[d] = -numeric_limits<double>::infinity();
	}

	for (int x = 0; x < 5; x++) {
		double interval = abs(xyzVals[0][x + 1] - xyzVals[0][x]);
		avgIntervals[0] += interval;
		if (interval > EPSILON) {
			if (interval < minIntervals[0])
				minIntervals[0] = interval;
			else if (interval > maxIntervals[0])
				maxIntervals[0] = interval;
		}
	}
	for (int y = 0; y < 5; y++) {
		double interval = abs(xyzVals[1][y + 1] - xyzVals[1][y]);
		avgIntervals[1] += interval;
		if (interval > EPSILON) {
			if (interval < minIntervals[1])
				minIntervals[1] = interval;
			else if (interval > maxIntervals[1])
				maxIntervals[1] = interval;
		}
	}
	for (int z = 0; z < 5; z++) {
		double interval = abs(xyzVals[2][z + 1] - xyzVals[2][z]);
		avgIntervals[2] += interval;
		if (interval > EPSILON) {
			if (interval < minIntervals[2])
				minIntervals[2] = interval;
			else if (interval > maxIntervals[2])
				maxIntervals[2] = interval;
		}
	}


	for (int d = 0; d < 3; d++)
		avgIntervals[d] /= 5;

	minIntervals[0] = avgIntervals[0];
	minIntervals[1] = avgIntervals[1];
	minIntervals[2] = avgIntervals[2];

	/*
	minIntervals[0] = maxIntervals[0];
	minIntervals[1] = maxIntervals[1];
	minIntervals[2] = maxIntervals[2];
	*/

}

double KernelExpansion::findIdealCelSize() {

	double xyzVals[3][6];
	for (int d = 0; d < 3; d++) {
		for (int i = 0; i < 6; i++) {
			xyzVals[d][i] = numeric_limits<double>::infinity();
		}
	}

	// increasingly order extreme point coordinates 
	for (int d = 0; d < 3; d++) {
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				if (extremePoints[i][d] < xyzVals[d][j]) {
					for (int k = 5; k > j; k--)	// shift
						xyzVals[d][k] = xyzVals[d][k - 1];
					xyzVals[d][j] = extremePoints[i][d];
					break;
				}
			}
		}
	}

	double minInterval = numeric_limits<double>::infinity();
	double maxInterval = -numeric_limits<double>::infinity();
	double avgInterval = 0;
	for (int x = 0; x < 5; x++) {
		double interval = abs(xyzVals[0][x + 1] - xyzVals[0][x]);
		avgInterval += interval;
		if (interval > 0 && interval < minInterval)
			minInterval = interval;
	}
	for (int y = 0; y < 5; y++) {
		double interval = abs(xyzVals[0][y + 1] - xyzVals[0][y]);
		avgInterval += interval;
		if (interval > 0 && interval < minInterval)
			minInterval = interval;
	}
	for (int z = 0; z < 5; z++) {
		double interval = abs(xyzVals[0][z + 1] - xyzVals[0][z]);
		avgInterval += interval;
		if (interval > 0 && interval < minInterval)
			minInterval = interval;
	}

	avgInterval /= 15;
	return avgInterval;

}

void KernelExpansion::checkKernelForNonKernelVertices() {

	int num_of_non_kernel_points = 0;
	for (int i = 0; i < kernel.getNumOfVerts(); i++) {
		if (findClosestValueSatisfiedByPoint(kernel.getVertex(i).coords, halfSpaceSet) > 3 * EPSILON ||
			findClosestValueSatisfiedByPoint(kernel.getVertex(i).coords, halfSpaceSet) > 3 * EPSILON ||
			findClosestValueSatisfiedByPoint(kernel.getVertex(i).coords, halfSpaceSet) > 3 * EPSILON) {
			cout << "non-kernel point" << endl;
			num_of_non_kernel_points++;
		}
	}
	cout << "num of non-kernel points is: " << num_of_non_kernel_points << endl;

}

void KernelExpansion::checkKernelForIrregularTriangles() {

	int num_of_irregular_triangles = 0;
	for (int i = 0; i < kernel.getNumOfTris(); i++) {
		if (kernel.getTriangle(i).corners[0] == kernel.getTriangle(i).corners[1] ||
			kernel.getTriangle(i).corners[1] == kernel.getTriangle(i).corners[2] ||
			kernel.getTriangle(i).corners[2] == kernel.getTriangle(i).corners[0]) {
			cout << "errorrrrr: triangle no: " << i << endl;
			num_of_irregular_triangles++;
		}
	}
	cout << "num of irregular triangles is: " << num_of_irregular_triangles << endl;

}



Grid::Grid(const Mesh* hostMeshPtr, double* leastCoordinates, double* mostCoordinates, int* gridDimension, double cellSizeRatio) {

	if (hostMeshPtr)
		this->defineGridByMeshCoordinates(*hostMeshPtr);							// USE MESH BOUNDARY BOX TO DEFINE GRID
	else
		this->defineGridByExternalCoordinates(leastCoordinates, mostCoordinates);	// USE THE GIVEN BOUNDARY BOX TO DEFINE GRID

	if (gridDimension)
		this->partitionGridByNumOfCells(gridDimension);
	else if (cellSizeRatio > 0)
		this->partitionGridByDiagonalRatio(cellSizeRatio);
	else
		this->noPartitionGrid();
}

void Grid::defineGridByMeshCoordinates(const Mesh& hostMesh) {

	for (int j = 0; j < 3; j++) {
		this->minGridCoords[j] = numeric_limits<double>::infinity();
		this->maxGridCoords[j] = -numeric_limits<double>::infinity();
	}

	// find the "most" and "least" coordinates of the mesh
	for (int i = 0; i < hostMesh.getNumOfVerts(); i++) {
		Vertex vertex = hostMesh.getVertex(i);
		for (int j = 0; j < 3; j++) {
			if (this->minGridCoords[j] > vertex.coords[j])
				this->minGridCoords[j] = vertex.coords[j];
			if (this->maxGridCoords[j] < vertex.coords[j])
				this->maxGridCoords[j] = vertex.coords[j];
		}
	}

}

void Grid::defineGridByExternalCoordinates(double* leastCoordinates, double* mostCoordinates) {

	for (int j = 0; j < 3; j++) {
		this->minGridCoords[j] = leastCoordinates[j];
		this->maxGridCoords[j] = mostCoordinates[j];
	}

}

void Grid::partitionGridByCellSize(double cellSize) {

	for (int j = 0; j < 3; j++)
		this->cellSize[j] = cellSize;

	for (int j = 0; j < 3; j++) {
		this->minGridCoords[j] -= this->cellSize[j] / 2;
		this->maxGridCoords[j] += this->cellSize[j] / 2;
	}

	// compute number of cells of grid at each direction
	for (int j = 0; j < 3; j++)
		this->numOfCells[j] = (this->maxGridCoords[j] - this->minGridCoords[j]) / this->cellSize[j];

}

void Grid::partitionGridByCellSize(double cellSize[3], double* centerPoint, bool update_borders) {

	for (int j = 0; j < 3; j++)
		this->cellSize[j] = cellSize[j];

	double minCornerOfCenter[3], maxCornerOfCenter[3];

	while (minCornerOfCenter[0] > this->minGridCoords[0] || minCornerOfCenter[1] > this->minGridCoords[1] || minCornerOfCenter[2] > this->minGridCoords[2]) {
		for (int j = 0; j < 3; j++)
			minCornerOfCenter[j] -= cellSize[j];
	}

	while (maxCornerOfCenter[0] < this->maxGridCoords[0] || maxCornerOfCenter[1] < this->maxGridCoords[1] || maxCornerOfCenter[2] < this->maxGridCoords[2]) {
		for (int j = 0; j < 3; j++)
			maxCornerOfCenter[j] += cellSize[j];
	}

	if (update_borders) {
		for (int j = 0; j < 3; j++) {
			minCornerOfCenter[j] -= cellSize[j];
			maxCornerOfCenter[j] += cellSize[j];

			this->minGridCoords[j] = minCornerOfCenter[j];
			this->maxGridCoords[j] = maxCornerOfCenter[j];
		}
	}

	// compute number of cells of grid at each direction
	for (int j = 0; j < 3; j++)
		this->numOfCells[j] = (this->maxGridCoords[j] - this->minGridCoords[j]) / this->cellSize[j];

}

void Grid::partitionGridByDiagonalRatio(double ratio) {

	// compute each edge length of the grid
	double edgeLength[3];
	for (int j = 0; j < 3; j++)
		edgeLength[j] = this->maxGridCoords[j] - this->minGridCoords[j];

	// compute diagonal's length & cellSize
	double diagonal = sqrt((edgeLength[0] * edgeLength[0]) + (edgeLength[1] * edgeLength[1]) + (edgeLength[2] * edgeLength[2]));
	double cellSize = diagonal * ratio;

	partitionGridByCellSize(cellSize);

}

void Grid::partitionGridByNumOfCells(int gridDimension[3]) {

	if (gridDimension[0] == gridDimension[1] && gridDimension[1] == gridDimension[2]) {
		partitionGridIntoEqualPieces(gridDimension[0]);
		return;
	}

	for (int d = 0; d < 3; d++) {
		this->numOfCells[d] = gridDimension[d];
		this->cellSize[d] = (this->maxGridCoords[d] - this->minGridCoords[d]) / this->numOfCells[d];
	}

}

void Grid::partitionGridIntoEqualPieces(int numOfPieces) {

	// assumes the grid is a cube!!!

	// find the minimum cell size ( depends on floating point error)
	double cellSize = numeric_limits<double>::infinity();
	for (int j = 0; j < 3; j++) {
		double temp = (this->maxGridCoords[j] - this->minGridCoords[j]) / numOfPieces;
		if (temp < cellSize)
			cellSize = temp;
	}

	partitionGridByCellSize(cellSize);

}

void Grid::noPartitionGrid() {

	for (int d = 0; d < 3; d++) {
		this->numOfCells[d] = 1;
		this->cellSize[d] = (this->maxGridCoords[d] - this->minGridCoords[d]) / this->numOfCells[d];
	}

}

Grid::Grid(const Grid& grid) {

	for (int i = 0; i < 3; i++) {
		numOfCells[i] = grid.numOfCells[i];
		cellSize[i] = grid.cellSize[i];
		minGridCoords[i] = grid.minGridCoords[i];
		maxGridCoords[i] = grid.maxGridCoords[i];
		MIN_NUM_OF_CELLS = grid.MIN_NUM_OF_CELLS;
	}

	gridId = grid.gridId;
	generationId = grid.generationId;
}

int Grid::getNumOfCells() {

	return numOfCells[0] * numOfCells[1] * numOfCells[2];
}

double Grid::getCellSize() {

	return cellSize[0];
}

double* Grid::getMinCoords() {

	return minGridCoords;
}

double* Grid::getMaxCoords() {

	return maxGridCoords;
}

int* Grid::findHomeCell(double* point) {

	int* cell_indices = new int[3];

	for (int d = 0; d < 3; d++) {	// d: dimension
		cell_indices[d] = numOfCells[d];
		double coordValue = minGridCoords[d];
		for (int i = 0; i < numOfCells[d] + 1; i++) {
			if (point[d] < coordValue) {
				cell_indices[d] = i - 1;
				break;
			}
			coordValue += cellSize[d];
		}
	}

	return cell_indices;
}

double* Grid::findCellCoords(int i, int j, int k) {

	int cell_indices[3] = { i, j, k };
	double* cornerCoords = new double[3];
	for (int c = 0; c < 3; c++)
		cornerCoords[c] = minGridCoords[c] + (cell_indices[c] * cellSize[c]);

	return cornerCoords;
}

void Grid::findCellCoords(int i, int j, int k, double* coords) {

	int cell_indices[3] = { i, j, k };
	for (int c = 0; c < 3; c++)
		coords[c] = minGridCoords[c] + (cell_indices[c] * cellSize[c]);

}


