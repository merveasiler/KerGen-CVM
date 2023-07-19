// @author Merve Asiler

#pragma once

#include "Mesh.h"
#include "BasicGeometricElements.h"
#include "KernelExpansion.h"
#include <queue>

enum class ProcessColor_Robust {
	RED = 2,
	WHITE = 4,
	GREEN = 6
};

struct VertexPartners_Robust {

	int partner1_id;
	int partner2_id;
	int vertex_id;

	VertexPartners_Robust(int partner1_id, int partner2_id, int vertex_id) {
		this->partner1_id = partner1_id;
		this->partner2_id = partner2_id;
		this->vertex_id = vertex_id;
	};
};

class KerGen_Robust : public KernelExpansion {

	vector<vector<VertexPartners_Robust>> vertexPartnerships;
	vector<ProcessColor_Robust> faceProcess;
	queue<int> tobeProcessedFaceIds;
	vector<vector<int>> generatorsOfVertices;

	void initialize();
	int findClosestHalfSpace(double* point);
	int findClosestHalfSpace(double* point, int id);
	vector<int> findClosestHalfSpace(double* point, int plane1_id, int plane2_id);
	vector<int> findClosestHalfSpace(int base_plane_id, int partner_plane_id, int back_plane_id, int vertex_id);
	int computeKernelVertex(int plane1_id, int plane2_id, vector<int> others);
	vector<int> findFrontestPlane(int base_plane_id, int partner_plane_id, vector<int> latest_plane_ids);
	bool addIntoProcessQueueAfterCheck(int plane1_id, int plane2_id, int plane3_id, int plane4_id, int vertex_id);

public:
	KerGen_Robust(const Mesh& hostMesh);
	~KerGen_Robust();
	void expandKernel();

	void findInitialPoint_1(double* point);
};