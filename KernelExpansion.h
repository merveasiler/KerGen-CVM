// @author Merve Asiler

#pragma once

#include "Mesh.h"
#include "BasicGeometricElements.h"
#include <queue>

struct Grid {
	int numOfCells[3];	// number of cells for each of the three dimension
	double cellSize[3];	// length of the edges for each of the three dimension, normally equal because of cube
	double minGridCoords[3], maxGridCoords[3];
	int MIN_NUM_OF_CELLS = 4;
	int gridId = 0;
	int generationId = 0;

	Grid() {};
	Grid(const Mesh* hostMeshPtr, double* leastCoordinates, double* mostCoordinates, int* gridDimension, double cellSizeRatio);
	~Grid() {};
	Grid(const Grid& grid);
	void defineGridByMeshCoordinates(const Mesh& hostMesh);
	void defineGridByExternalCoordinates(double* leastCoordinates, double* mostCoordinates);
	void partitionGridByCellSize(double cellSize);
	void partitionGridByCellSize(double cellSize[3], double* centerPoint, bool update_borders);
	void partitionGridByDiagonalRatio(double ratio);
	void partitionGridByNumOfCells(int gridDimension[3]);
	void partitionGridIntoEqualPieces(int numOfPieces);
	void noPartitionGrid();
	int getNumOfCells();
	double getCellSize();
	double* getMinCoords();
	double* getMaxCoords();
	int* findHomeCell(double* point);
	double* findCellCoords(int i, int j, int k);
	void findCellCoords(int i, int j, int k, double* coords);

};

class KernelExpansion {

protected:
	double* initialPoint = nullptr;
	double* extremeCorners[2] = { nullptr, nullptr };
	double* extremePoints[6] = { nullptr, nullptr, nullptr, nullptr, nullptr, nullptr };
	const Mesh* hostMeshptr;
	Mesh kernel;

	Grid grid;
	double EPSILON_3 = 3 * EPSILON;

	queue<int*> possibleKernelCells;		// the cells which may include some pieces of the kernel
	vector<int*> handledCells;				// the cells which were previously checked or in queue to be checked for being inside kernel 
	vector<int> handledGrids;				// the grids ids showing the ownership of the corresponding cell in handledCells
	queue<double*> kernelPointsForCells;	// holds a kernel point for each corresponding cell in possibleKernelCells
	queue<Grid*> ownerGridForCells;			// holds grids in which cells are defined
	vector<HalfSpace> halfSpaceSet;

	void decideOnNeighbor2(int neigh_i, int neigh_j, int neigh_k, double neighVert1_scalar, double neighVert2_scalar, double neighVert3_scalar, double neighVert4_scalar);
	void decideOnNeighbor(int neigh_i, int neigh_j, int neigh_k, double coords[8][4], double* scalars[8], int ind1, int ind2, int ind3, int ind4, Grid* grid);
	bool isCellValidForQueue(int  neighbor_i, int neighbor_j, int neighbor_k, Grid* grid);
	void fillContainers(Grid* grid, int* initialCell, double* initialPoint);
	bool isInKernel(double* scalarVector);
	vector<int> detectExcluderHalfSpaces(double* scalarVector);
	double findIdealCelSize_OuterToInner();
	double findIdealCelSize_InnerToOuter();
	double findIdealCelSize();
	void findIdealCelSizePerCoordinate(double* minIntervals);
	void findSensitiveIdealCelSizePerCoordinate(double* minIntervals);

public:
	KernelExpansion(const Mesh& hostMesh, bool defineInitialPoint);
	KernelExpansion(const Mesh& hostMesh);
	~KernelExpansion();
	Mesh& getKernel();
	double* getInitialKernelPoint();
	vector<HalfSpace>& getHalfSpaceSet();
	Grid& getGrid();
	virtual void expandKernel() = 0;
	void checkKernelForNonKernelVertices();
	void checkKernelForIrregularTriangles();
};
