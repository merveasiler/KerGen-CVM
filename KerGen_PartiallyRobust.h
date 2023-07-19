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


enum class ProcessColor_PartiallyRobust {
	RED = 2,
	WHITE = 4,
	GREEN = 6
};

struct VertexPartners_PartiallyRobust {

	int partner1_id;
	int partner2_id;

	VertexPartners_PartiallyRobust(int partner1_id, int partner2_id) {
		this->partner1_id = partner1_id;
		this->partner2_id = partner2_id;
	};
};

class KerGen_PartiallyRobust : public KernelExpansion {

	vector<XK::Plane_3> cgalPlaneSet;
	vector<vector<VertexPartners_PartiallyRobust>> vertexPartnerships;
	vector<ProcessColor_PartiallyRobust> faceProcess;
	queue<int> tobeProcessedFaceIds;
	vector<vector<int>> generatorsOfVertices;

	void initialize();
	int findClosestHalfSpace(double* point);
	int findClosestHalfSpace(double* point, int id);
	vector<int> findClosestHalfSpace(int base_plane_id, int partner_plane_id);
	vector<int> findClosestHalfSpace(int base_plane_id, int partner_plane_id, int back_plane_id);
	int computeKernelVertex(int plane1_id, int plane2_id, int plane3_id);
	vector<int> orderFromFrontToBack(int base_plane_id, int partner_plane_id, vector<int> latest_plane_ids);
	bool addIntoProcessQueueAfterCheck(int plane1_id, int plane2_id, int plane3_id, int base_plane_id);

	bool ClosestPlaneToPoint_3D(int candidate_plane_id, XK::Point_3 point);
	bool ClosestPlaneToPoint_2D(int candidate_plane_id, int base_plane_id, XK::Point_3 point);
	bool ClosestPlaneToPoint_1D(int candidate_plane_id, int base_plane_id, int partner_plane_id, XK::Line_3 line);
	bool ClosestPlaneToPoint_0D(int candidate_plane_id, int base_plane_id, int partner_plane_id, int back_plane_id, XK::Line_3 line);
	XK::Point_3 CPoint(double* point);
	XK::Plane_3 CPlane(HalfSpace halfspace);

public:
	KerGen_PartiallyRobust(const Mesh& hostMesh);
	~KerGen_PartiallyRobust();
	void expandKernel();

	void findInitialPoint_1(double* point);
};