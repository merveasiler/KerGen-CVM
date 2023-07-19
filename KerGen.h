// @author Merve Asiler

#pragma once

#include "Mesh.h"
#include "BasicGeometricElements.h"
#include "FilteredPredicates.h"
#include "KernelExpansion.h"
#include <queue>

enum class ProcessColor {
	RED = 2,
	WHITE = 4,
	GREEN = 6
};

struct EdgePartnerTriple {

	int partnerId;
	double edgeDirection[3];
	double startPoint[3];
	int startPointId;

	EdgePartnerTriple(int partnerId, double* edgeDirection, double* startPoint, int startPointId) {
		this->partnerId = partnerId;
		for (int i = 0; i < 3; i++) {
			this->edgeDirection[i] = edgeDirection[i];
			this->startPoint[i] = startPoint[i];
		}
		this->startPointId = startPointId;
	};

	~EdgePartnerTriple() {
	}

};

class KerGen : public KernelExpansion {

	vector<queue<EdgePartnerTriple>> edgePartners;	// "id" of the other plane to construct and edge line with this & "direction" to walk on the line
	vector<vector<int>> edgePartnerIds;
	vector<vector<int>> vertexParentIds;
	vector<ProcessColor> isKernelFace;
	double _EPSILON = 3 * EPSILON;
	double BIG_EPSILON = 1e-8;

	vector<double> initialize(double* point);
	int findTheClosestHalfSpace(double* point, vector<double>& scalarsVector);
	int findTheClosestHalfSpace(double* point, int id);
	vector<int> findTheClosestHalfSpace(int vertexId, double* lineDirection, int lineParent1Id, int lineParent2Id, double* newpoint);
	void orderTheFaces(int base_id, int partner_id, vector<int> next_partner_ids, double* startPoint, double* currentEdgeDirection);
	bool isWalkedEdge(int hs1_id, int hs2_id);
	void findEdgeDirection(int hs1_id, int hs2_id, bool should_revert, double* edgeDirection, double* directioner);
	bool isValidEdge(double* startPoint, double* edgeDirection);
	void filterRepetitions(double* distances, vector<double>& scalarsVector);
	void filterRepetitions2(double* distances, vector<double>& scalarsVector);


	bool isLineDirectionValid(double* startPoint, double* lineDirection);
	bool isLineValid(double* startPoint, double* lineDirection);
	bool shiftOnAnotherPlane(double* point, double* rootpoint, int& shiftId);
	bool shiftOnAnotherPlane(double* point, double* rootpoint, int& shiftId, int baseId);
	bool isPointOnPlane(double* point, int planeId);
	bool isPointInFrontOfThePlane(double* point, int planeId);
	bool areTheSamePlanes(int plane1Id, int plane2Id);
	void findInitialVertexAndLines(double* point);
	bool identifyInitialPlanes(double* vertex);
	vector<int> findFirstPlane(double* point, double& distanceToFirst);
	vector<int> findSecondPlane(double* point, int id, vector<double>& theClosestDirections, vector<double>& intersectionLineDirections, double& theClosestDistance);
	vector<int> findThirdPlane(double* startPoint, double* lineDirection, int lineParent1Id, int lineParent2Id, double* newpoint);
	bool isPreviousLine(int hs1_id, int hs2_id);
	void recordLineInfo(int hs1_id, int hs2_id, double* lineDir, double* point);
	bool recordWalkedLineInfo(int hs1_id, int hs2_id);
	int recordVertexInfo(double* vertex, vector<int> thirdParentIds, vector<int> lineParentIds);
	vector<int> findClosestPlaneOnLine(int startVertexId, double* lineDirection, int lineParent1Id, int lineParent2Id, double* vertex);
	vector<int> findClosestPlaneOnLine2(int startVertexId, double* lineDirection, int lineParent1Id, int lineParent2Id, double* vertex);

public:
	KerGen(const Mesh& hostMesh);
	~KerGen();
	void expandKernel();

	void findInitialPoint_1(double* point);
	void findInitialPoint_2(double* point);
	void findInitialPoint_3(double* point);
	void findInitialPoint_4(double* point);
	void findInitialPoint_5(double* point);
};