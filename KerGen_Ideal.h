// @author Merve Asiler

#pragma once

#include "Mesh.h"
#include "BasicGeometricElements.h"
#include "FilteredPredicates.h"
#include "KernelExpansion.h"
#include <queue>

enum class ProcessColor_Ideal {
	RED = 2,
	WHITE = 4,
	GREEN = 6
};

struct EdgePartnerTriple_Ideal {

	int partnerId;
	double edgeDirection[3];
	double startPoint[3];
	int startPointId;
	int backPlaneId;

	EdgePartnerTriple_Ideal(int partnerId, double* edgeDirection, double* startPoint, int startPointId, int backPlaneId) {
		this->partnerId = partnerId;
		for (int i = 0; i < 3; i++) {
			this->edgeDirection[i] = edgeDirection[i];
			this->startPoint[i] = startPoint[i];
		}
		this->startPointId = startPointId;
		this->backPlaneId = backPlaneId;
	};

	~EdgePartnerTriple_Ideal() {
	}

};

class KerGen_Ideal : public KernelExpansion {

	vector<queue<EdgePartnerTriple_Ideal>> edgePartners;	// "id" of the other plane to construct and edge line with this & "direction" to walk on the line
	vector<vector<int>> edgePartnerIds;
	vector<vector<int>> vertexParentIds;
	vector<ProcessColor_Ideal> isKernelFace;
	double BIG_EPSILON = 1e-12;
	
	vector<double> initialize(double* point);
	int findTheClosestHalfSpace2(double* point, vector<double>& scalarsVector);
	int findTheClosestHalfSpace2(double* point, int id);
	int findTheClosestHalfSpace(double* point, vector<double>& scalarsVector);
	int findTheClosestHalfSpace(double* point, int id);
	vector<int> findTheClosestHalfSpace(int vertexId, double* lineDirection, int lineParent1Id, int lineParent2Id, int backPlaneId, double* newpoint);
	void orderTheFaces(int base_id, int partner_id, vector<int> next_partner_ids, double* startPoint, double* currentEdgeDirection);
	bool isRecordedEdge(int hs1_id, int hs2_id);
	void filterRepetitions2(double* distances, vector<double>& scalarsVector);

public:
	KerGen_Ideal(const Mesh& hostMesh);
	~KerGen_Ideal();
	void expandKernel();

	void findInitialPoint_1(double* point);

};