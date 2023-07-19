// @author Merve Asiler

#pragma once

#include "Mesh.h"
#include "BasicGeometricElements.h"
#include "FilteredPredicates.h"
#include "KernelExpansion.h"
#include <queue>

enum class ProcessColor_Modified {
	RED = 2,
	WHITE = 4,
	GREEN = 6
};

struct EdgePartnerTriple_Modified {

	int partnerPlaneId;
	int backPlaneId;
	double edgeDirection[3];
	double startPoint[3];
	int startPointId;

	EdgePartnerTriple_Modified(int partnerPlaneId, int backPlaneId, double* edgeDirection, double* startPoint, int startPointId) {
		this->partnerPlaneId = partnerPlaneId;
		this->backPlaneId = backPlaneId;
		for (int i = 0; i < 3; i++) {
			this->edgeDirection[i] = edgeDirection[i];
			this->startPoint[i] = startPoint[i];
		}
		this->startPointId = startPointId;
	};

	~EdgePartnerTriple_Modified() {
	}

};

class KerGen_Modified : public KernelExpansion {

	vector<queue<EdgePartnerTriple_Modified>> edgePartners;	// "id" of the other plane to construct and edge line with this & "direction" to walk on the line
	vector<vector<int>> edgePartnerIds;
	vector<vector<int>> vertexParentIds;
	vector<ProcessColor_Modified> isKernelFace;
	double BIG_EPSILON = 1e-8;
	double _EPSILON = 3 * EPSILON;

	void initialize();
	bool isLineDirectionValid(double* startPoint, double* lineDirection);
	bool isLineValid(double* startPoint, double* lineDirection);
	bool correctLineDirection(double* startPoint, double* lineDirection, int referencePlaneId);
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
	void recordLineInfo(int hs1_id, int hs2_id, int hs3_id, double* lineDir, double* vertex);
	int recordVertexInfo(double* vertex, vector<int> thirdParentIds, vector<int> lineParentIds);

	vector<int> findTheClosestHalfSpace(int vertexId, double* lineDirection, int lineParent1Id, int lineParent2Id, int backPlaneId, double* newpoint);
	vector<int> orderFromFrontToBack(int base_id, int partner_id, vector<int> next_partner_ids, double* startPoint);
	void identifyNextEdges(int base_plane_id, int partner_plane_id, vector<int> ordered_next_partner_ids, double* startpoint);

public:
	KerGen_Modified(const Mesh& hostMesh);
	~KerGen_Modified();
	void expandKernel();
};