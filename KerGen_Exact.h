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
typedef EK::Point_3												Point_3;
typedef EK::Plane_3												Plane_3;
typedef EK::Line_3												Line_3;
typedef EK::Vector_3											Vector_3;
typedef EK::RT													RT;


enum class ProcessColor_Exact {
	RED = 2,
	WHITE = 4,
	GREEN = 6
};

struct VertexPartners_Exact {

	int partner1_id;
	int partner2_id;
	int vertex_id;

	VertexPartners_Exact(int partner1_id, int partner2_id, int vertex_id) {
		this->partner1_id = partner1_id;
		this->partner2_id = partner2_id;
		this->vertex_id = vertex_id;
	};
};

class KerGen_Exact : public KernelExpansion {

	vector<vector<VertexPartners_Exact>> vertexPartnerships;
	vector<ProcessColor_Exact> faceProcess;
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

	EK::Point_3 CPoint(double* point);
	EK::Plane_3 CPlane(HalfSpace halfspace);
	bool ParallelPlanes(Point_3 plane1_p1, Point_3 plane1_p2, Point_3 plane1_p3,
		Point_3 plane2_p1, Point_3 plane2_p2, Point_3 plane2_p3);
	bool SamePlanes(Point_3 plane1_p1, Point_3 plane1_p2, Point_3 plane1_p3,
		Point_3 plane2_p1, Point_3 plane2_p2, Point_3 plane2_p3);
	CGAL::Sign IntersectablePlane(Point_3 plane1_p1, Point_3 plane1_p2, Point_3 plane1_p3,
		Point_3 plane2_p1, Point_3 plane2_p2, Point_3 plane2_p3,
		Point_3 candidate_plane_p1, Point_3 candidate_plane_p2, Point_3 candidate_plane_p3);
	CGAL::Sign IntersectableFrontPlane(Point_3 plane1_p1, Point_3 plane1_p2, Point_3 plane1_p3,
		Point_3 plane2_p1, Point_3 plane2_p2, Point_3 plane2_p3,
		Point_3 candidate_plane_p1, Point_3 candidate_plane_p2, Point_3 candidate_plane_p3,
		Point_3 back_plane_p1, Point_3 back_plane_p2, Point_3 back_plane_p3);
	CGAL::Sign ClosestPlaneToPoint_3D(Point_3 plane1_p1, Point_3 plane1_p2, Point_3 plane1_p3,
		Point_3 plane2_p1, Point_3 plane2_p2, Point_3 plane2_p3,
		Point_3 point);
	CGAL::Sign ClosestPlaneToPoint_2D(Point_3 base_plane_p1, Point_3 base_plane_p2, Point_3 base_plane_p3,
		Point_3 plane1_p1, Point_3 plane1_p2, Point_3 plane1_p3,
		Point_3 plane2_p1, Point_3 plane2_p2, Point_3 plane2_p3,
		Point_3 initialpoint);
	CGAL::Comparison_result ClosestPlaneToPoint_1D(Point_3 plane1_p1, Point_3 plane1_p2, Point_3 plane1_p3,
		Point_3 plane2_p1, Point_3 plane2_p2, Point_3 plane2_p3,
		Point_3 candidate_plane_p1, Point_3 candidate_plane_p2, Point_3 candidate_plane_p3,
		Point_3 front_plane_p1, Point_3 front_plane_p2, Point_3 front_plane_p3,
		Point_3 back_plane_p1, Point_3 back_plane_p2, Point_3 back_plane_p3);
	CGAL::Comparison_result ClosestPlaneToPoint_UOD(Point_3 plane1_p1, Point_3 plane1_p2, Point_3 plane1_p3,
		Point_3 plane2_p1, Point_3 plane2_p2, Point_3 plane2_p3,
		Point_3 candidate_plane_p1, Point_3 candidate_plane_p2, Point_3 candidate_plane_p3,
		Point_3 front_plane_p1, Point_3 front_plane_p2, Point_3 front_plane_p3,
		Point_3 initialpoint);
	CGAL::Sign FrontierPlane(Point_3 base_plane_p1, Point_3 base_plane_p2, Point_3 base_plane_p3,
		Point_3 partner_plane_p1, Point_3 partner_plane_p2, Point_3 partner_plane_p3,
		Point_3 candidate_plane_p1, Point_3 candidate_plane_p2, Point_3 candidate_plane_p3,
		Point_3 reference_plane_p1, Point_3 reference_plane_p2, Point_3 reference_plane_p3);

public:
	KerGen_Exact(const Mesh& hostMesh);
	~KerGen_Exact();
	void expandKernel();

	void findInitialPoint_1(double* point);
};