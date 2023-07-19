// @author Merve Asiler

#pragma once

#include "Mesh.h"
#include "BasicGeometricElements.h"
#include "KernelExpansion.h"
#include <queue>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/enum.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Quotient.h>

typedef CGAL::Simple_cartesian<double>							C;
typedef CGAL::Exact_predicates_exact_constructions_kernel		EK;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel	IK;
typedef CGAL::Exact_kernel_selector<C>::Exact_kernel			XK;


enum class ProcessColor_Optimal {
	RED = 2,
	WHITE = 4,
	GREEN = 6
};

struct VertexPartners_Optimal {

	int partner1_id;
	int partner2_id;
	Line common_line_with_partner1;

	VertexPartners_Optimal(int partner1_id, int partner2_id, Line common_line_with_partner1) {
		this->partner1_id = partner1_id;
		this->partner2_id = partner2_id;
		this->common_line_with_partner1 = common_line_with_partner1;
	};
};

class KerGen_Optimal : public KernelExpansion  {

	vector<vector<VertexPartners_Optimal>> vertexPartnerships;
	vector<ProcessColor_Optimal> faceProcess;
	queue<int> tobeProcessedFaceIds;
	vector<vector<int>> generatorsOfVertices;
	double BIG_EPSILON = 1e-12;

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
	bool ClosestPlaneToPoint_0D(int candidate_plane_id, int base_plane_id, int partner_plane_id, int back_plane_id, Line line, double* candidate_point, double& smallest_t);

public:
	KerGen_Optimal(const Mesh& hostMesh);
	~KerGen_Optimal();
	void expandKernel();

	void findInitialPoint_1(double* point);
};