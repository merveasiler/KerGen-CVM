// @author Merve Asiler

#pragma once

#include "Mesh.h"
#include "BasicGeometricElements.h"
#include "KernelExpansion.h"
#include <queue>

enum class ProcessColor_Slow {
	RED = 2,
	WHITE = 4,
	GREEN = 6
};

struct VertexPartners_Slow {

	int partner1_id;
	int partner2_id;
	Line common_line_with_partner1;

	VertexPartners_Slow(int partner1_id, int partner2_id, Line common_line_with_partner1) {
		this->partner1_id = partner1_id;
		this->partner2_id = partner2_id;
		this->common_line_with_partner1 = common_line_with_partner1;
	};
};

class KerGen_Slow : public KernelExpansion {

	vector<vector<VertexPartners_Slow>> vertexPartnerships;
	vector<ProcessColor_Slow> faceProcess;
	queue<int> tobeProcessedFaceIds;
	vector<vector<int>> generatorsOfVertices;

	void initialize();
	int findClosestHalfSpace(double* point);
	int findClosestHalfSpace(double* point, int id);
	vector<int> findClosestHalfSpace(double* point, int base_plane_id, int partner_plane_id);
	vector<int> findClosestHalfSpace(double* point, int base_plane_id, int partner_plane_id, int back_plane_id, Line line);
	int recordFoundVertex(double* vertex);
	vector<int> orderFromFrontToBack(double* point, int base_plane_id, int partner_plane_id, vector<int> latest_plane_ids);
	bool addIntoProcessQueueAfterCheck(int plane1_id, int plane2_id, int plane3_id, int base_plane_id, double* point);

	bool ClosestPlaneToPoint_3D(int candidate_plane_id, double* point, double* candidate_point);
	bool ClosestPlaneToPoint_2D(int candidate_plane_id, int base_plane_id, double* point, double* candidate_point);
	bool ClosestPlaneToPoint_1D(int candidate_plane_id, int base_plane_id, int partner_plane_id, Line line, double* candidate_point);
	bool ClosestPlaneToPoint_0D(int candidate_plane_id, int base_plane_id, int partner_plane_id, int back_plane_id, Line line, double* candidate_point);

public:
	KerGen_Slow(const Mesh& hostMesh);
	~KerGen_Slow();
	void expandKernel();

	void findInitialPoint_1(double* point);
};