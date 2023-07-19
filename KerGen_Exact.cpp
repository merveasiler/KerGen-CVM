// @author Merve Asiler

#include "KerGen_Exact.h"
#include "BaseGeoOpUtils.h"
#include "sdlp.h"
#include "CGALUtils.h"

KerGen_Exact::KerGen_Exact(const Mesh& hostMesh) :
	KernelExpansion(hostMesh, true) {

}

KerGen_Exact::~KerGen_Exact() {

}

void KerGen_Exact::expandKernel() {

	//cout << "Number of bounding planes: " << halfSpaceSet.size() << " out of " << hostMeshptr->getNumOfTris() << " triangles." << endl;

	initialize();

	double point[3];
	if (this->initialPoint == NULL)
		return;		// NOT STAR-SHAPED!
	for (int i = 0; i < 3; i++)
		point[i] = this->initialPoint[i];

	clock_t begin1 = clock();
	int closest1_id = findClosestHalfSpace(point);
	clock_t end1 = clock();
	cout << "PREDICATE TIME-1: " << double(end1 - begin1) * 1000 / CLOCKS_PER_SEC << "   " << closest1_id << endl;

	clock_t begin3 = clock();
	int closest2_id = findClosestHalfSpace(point, closest1_id);
	clock_t end3 = clock();
	cout << "PREDICATE TIME-2: " << double(end3 - begin3) * 1000 / CLOCKS_PER_SEC << "   " << closest2_id << endl;

	clock_t begin4 = clock();
	vector<int> closest3_ids = findClosestHalfSpace(point, closest1_id, closest2_id);
	clock_t end4 = clock();
	cout << "PREDICATE TIME-3: " << double(end4 - begin4) * 1000 / CLOCKS_PER_SEC << "   " << closest3_ids[0] << endl;

	vector<int> ordered_closest3_ids = findFrontestPlane(closest1_id, closest2_id, closest3_ids);
	int closest3_id = ordered_closest3_ids[0];
	int vertex_id = computeKernelVertex(closest1_id, closest2_id, closest3_ids);
	addIntoProcessQueueAfterCheck(closest1_id, closest2_id, closest3_id, closest3_id, vertex_id);

	if (ordered_closest3_ids.size() == 1)
		addIntoProcessQueueAfterCheck(closest2_id, closest3_id, closest1_id, closest1_id, vertex_id);
	else {
		ordered_closest3_ids.push_back(closest2_id);
		for (int i = ordered_closest3_ids.size() - 1, j = i - 1; j >= 1; j--) {
			if (addIntoProcessQueueAfterCheck(ordered_closest3_ids[i], ordered_closest3_ids[j], ordered_closest3_ids[j - 1], closest1_id, vertex_id))
				i = j;
			if (j == 1)
				addIntoProcessQueueAfterCheck(ordered_closest3_ids[i], ordered_closest3_ids[j - 1], closest1_id, closest1_id, vertex_id);
		}
	}

	int number_of_kernel_faces = 0;
	int number_of_kernel_edges = 0;

	while (!tobeProcessedFaceIds.empty()) {
		int face_id = tobeProcessedFaceIds.front();
		//cout << face_id << endl;
		for (int partner = 0; ; partner++) {
			// initializations
			VertexPartners_Exact& vp = vertexPartnerships[face_id][partner];
			int partner1_id = vp.partner1_id;
			int partner2_id = vp.partner2_id;
			int vertex_id = vp.vertex_id;

			//cout << face_id << " " << partner1_id << " " << partner2_id << " " << vertex_id << endl;
			number_of_kernel_edges++;

			vector<int> closest_plane_ids = findClosestHalfSpace(face_id, partner1_id, partner2_id, vertex_id);
			vector<int> ordered_closest_plane_ids = findFrontestPlane(face_id, partner1_id, closest_plane_ids);

			int next_partner_id = ordered_closest_plane_ids[0];
			int initial_partner_id = vertexPartnerships[face_id][0].partner1_id;
			if (abs(dotProduct(halfSpaceSet[next_partner_id].ABCD, halfSpaceSet[initial_partner_id].ABCD) - 1.0) < EPSILON)
				break;	// loop completed

			int found_vertex_id = computeKernelVertex(face_id, partner1_id, closest_plane_ids);
			vertexPartnerships[face_id].push_back(VertexPartners_Exact(next_partner_id, partner1_id, found_vertex_id));

			if (ordered_closest_plane_ids.size() == 1)
				addIntoProcessQueueAfterCheck(partner1_id, next_partner_id, face_id, face_id, vertex_id);
			else {
				ordered_closest_plane_ids.push_back(partner1_id);
				for (int i = ordered_closest_plane_ids.size() - 1, j = i - 1; j >= 1; j--) {
					if (addIntoProcessQueueAfterCheck(ordered_closest_plane_ids[i], ordered_closest_plane_ids[j], ordered_closest_plane_ids[j - 1], face_id, vertex_id))
						i = j;
					if (j == 1)
						addIntoProcessQueueAfterCheck(ordered_closest_plane_ids[i], ordered_closest_plane_ids[j - 1], face_id, face_id, vertex_id);
				}
			}

		}

		faceProcess[face_id] = ProcessColor_Exact::GREEN;
		tobeProcessedFaceIds.pop();
		number_of_kernel_faces++;

	}


	//	cout << "[number of kernel edges: " << number_of_kernel_edges << "], number of kernel faces: " << number_of_kernel_faces << "]" << endl;

	if (kernel.getNumOfVerts() > 0)
		kernel = computeConvexHull(kernel.getAllVerts());

}

void KerGen_Exact::initialize() {

	for (int i = 0; i < halfSpaceSet.size(); i++) {
		vertexPartnerships.push_back(vector<VertexPartners_Exact>());
		faceProcess.push_back(ProcessColor_Exact::RED);
	}

}

int KerGen_Exact::findClosestHalfSpace(double* point) {

	// initialize the closest plane temporarily
	int closest_plane_id = 0;

	// update the closest plane if necessary
	for (int i = 1; i < halfSpaceSet.size(); i++) {
		if (ClosestPlaneToPoint_3D(CPoint(halfSpaceSet[closest_plane_id].point1), CPoint(halfSpaceSet[closest_plane_id].point2), CPoint(halfSpaceSet[closest_plane_id].point3),
			CPoint(halfSpaceSet[i].point1), CPoint(halfSpaceSet[i].point2), CPoint(halfSpaceSet[i].point3),
			CPoint(point)) == CGAL::LARGER)
			closest_plane_id = i;
	}

	return closest_plane_id;

}

int KerGen_Exact::findClosestHalfSpace(double* point, int base_plane_id) {

	// initialize the closest plane temporarily
	int closest_plane_id = 0;
	for (int i = 0; i < halfSpaceSet.size(); i++) {
		if (i == base_plane_id)
			continue;
		if (ParallelPlanes(CPoint(halfSpaceSet[base_plane_id].point1), CPoint(halfSpaceSet[base_plane_id].point2), CPoint(halfSpaceSet[base_plane_id].point3),
			CPoint(halfSpaceSet[i].point1), CPoint(halfSpaceSet[i].point2), CPoint(halfSpaceSet[i].point3)))
			continue;
		closest_plane_id = i;
		break;
	}

	// update the closest plane if necessary
	for (int i = closest_plane_id + 1; i < halfSpaceSet.size(); i++) {
		if (i == base_plane_id)
			continue;
		if (ParallelPlanes(CPoint(halfSpaceSet[base_plane_id].point1), CPoint(halfSpaceSet[base_plane_id].point2), CPoint(halfSpaceSet[base_plane_id].point3),
			CPoint(halfSpaceSet[i].point1), CPoint(halfSpaceSet[i].point2), CPoint(halfSpaceSet[i].point3)))
			continue;
		if (ClosestPlaneToPoint_2D(CPoint(halfSpaceSet[base_plane_id].point1), CPoint(halfSpaceSet[base_plane_id].point2), CPoint(halfSpaceSet[base_plane_id].point3),
			CPoint(halfSpaceSet[closest_plane_id].point1), CPoint(halfSpaceSet[closest_plane_id].point2), CPoint(halfSpaceSet[closest_plane_id].point3),
			CPoint(halfSpaceSet[i].point1), CPoint(halfSpaceSet[i].point2), CPoint(halfSpaceSet[i].point3),
			CPoint(point)) == CGAL::LARGER)
			closest_plane_id = i;
	}

	return closest_plane_id;

}

vector<int> KerGen_Exact::findClosestHalfSpace(double* point, int base_plane_id, int partner_plane_id) {

	vector<int> closest_plane_ids;
	// initialize the closest plane temporarily
	int closest_plane_id = 0;
	for (int i = 0; i < halfSpaceSet.size(); i++) {
		if (i == base_plane_id || i == partner_plane_id)
			continue;
		if (IntersectablePlane(CPoint(halfSpaceSet[base_plane_id].point1), CPoint(halfSpaceSet[base_plane_id].point2), CPoint(halfSpaceSet[base_plane_id].point3),
			CPoint(halfSpaceSet[partner_plane_id].point1), CPoint(halfSpaceSet[partner_plane_id].point2), CPoint(halfSpaceSet[partner_plane_id].point3),
			CPoint(halfSpaceSet[i].point1), CPoint(halfSpaceSet[i].point2), CPoint(halfSpaceSet[i].point3))
			== CGAL::NEGATIVE)
			continue;
		closest_plane_id = i;
		closest_plane_ids.push_back(i);
		break;
	}

	// update the closest plane if necessary
	for (int i = closest_plane_id + 1; i < halfSpaceSet.size(); i++) {
		if (i == base_plane_id || i == partner_plane_id)
			continue;
		CGAL::Comparison_result cr = ClosestPlaneToPoint_UOD(
			CPoint(halfSpaceSet[base_plane_id].point1), CPoint(halfSpaceSet[base_plane_id].point2), CPoint(halfSpaceSet[base_plane_id].point3),
			CPoint(halfSpaceSet[partner_plane_id].point1), CPoint(halfSpaceSet[partner_plane_id].point2), CPoint(halfSpaceSet[partner_plane_id].point3),
			CPoint(halfSpaceSet[i].point1), CPoint(halfSpaceSet[i].point2), CPoint(halfSpaceSet[i].point3),
			CPoint(halfSpaceSet[closest_plane_id].point1), CPoint(halfSpaceSet[closest_plane_id].point2), CPoint(halfSpaceSet[closest_plane_id].point3),
			CPoint(point));
		if (cr == CGAL::LARGER) {
			closest_plane_id = i;
			closest_plane_ids.clear();
			closest_plane_ids.push_back(i);
		}
		if (cr == CGAL::EQUAL)
			closest_plane_ids.push_back(i);
	}

	return closest_plane_ids;

}

vector<int> KerGen_Exact::findClosestHalfSpace(int base_plane_id, int partner_plane_id, int back_plane_id, int vertex_id) {

	vector<int> closest_plane_ids;
	// initialize the closest plane temporarily
	int closest_plane_id = 0;
	for (int i = 0; i < halfSpaceSet.size(); i++) {
		bool is_generator = false;
		for (int g = 0; g < generatorsOfVertices[vertex_id].size(); g++)
			if (i == generatorsOfVertices[vertex_id][g]) {
				is_generator = true;
				break;
			}
		if (is_generator)
			continue;

		if (IntersectableFrontPlane(CPoint(halfSpaceSet[base_plane_id].point1), CPoint(halfSpaceSet[base_plane_id].point2), CPoint(halfSpaceSet[base_plane_id].point3),
			CPoint(halfSpaceSet[partner_plane_id].point1), CPoint(halfSpaceSet[partner_plane_id].point2), CPoint(halfSpaceSet[partner_plane_id].point3),
			CPoint(halfSpaceSet[i].point1), CPoint(halfSpaceSet[i].point2), CPoint(halfSpaceSet[i].point3),
			CPoint(halfSpaceSet[back_plane_id].point1), CPoint(halfSpaceSet[back_plane_id].point2), CPoint(halfSpaceSet[back_plane_id].point3))
			!= CGAL::NEGATIVE)
			continue;
		closest_plane_id = i;
		closest_plane_ids.push_back(i);
		break;
	}
	//cout << "closest-1: " << closest_plane_id << endl;
	// update the closest plane if necessary
	for (int i = closest_plane_id + 1; i < halfSpaceSet.size(); i++) {
		bool is_generator = false;
		for (int g = 0; g < generatorsOfVertices[vertex_id].size(); g++)
			if (i == generatorsOfVertices[vertex_id][g]) {
				is_generator = true;
				break;
			}
		if (is_generator)
			continue;

		CGAL::Comparison_result cr = ClosestPlaneToPoint_1D(
			CPoint(halfSpaceSet[base_plane_id].point1), CPoint(halfSpaceSet[base_plane_id].point2), CPoint(halfSpaceSet[base_plane_id].point3),
			CPoint(halfSpaceSet[partner_plane_id].point1), CPoint(halfSpaceSet[partner_plane_id].point2), CPoint(halfSpaceSet[partner_plane_id].point3),
			CPoint(halfSpaceSet[i].point1), CPoint(halfSpaceSet[i].point2), CPoint(halfSpaceSet[i].point3),
			CPoint(halfSpaceSet[closest_plane_id].point1), CPoint(halfSpaceSet[closest_plane_id].point2), CPoint(halfSpaceSet[closest_plane_id].point3),
			CPoint(halfSpaceSet[back_plane_id].point1), CPoint(halfSpaceSet[back_plane_id].point2), CPoint(halfSpaceSet[back_plane_id].point3));
		if (cr == CGAL::LARGER) {
			closest_plane_id = i;
			closest_plane_ids.clear();
			closest_plane_ids.push_back(i);
		}
		else if (cr == CGAL::EQUAL)
			closest_plane_ids.push_back(i);
	}
	//cout << "closest-2: " << closest_plane_id << endl;
	return closest_plane_ids;
}

vector<int> KerGen_Exact::findFrontestPlane(int base_plane_id, int partner_plane_id, vector<int> latest_plane_ids) {

	vector<int> planesFromFrontToBack;
	for (int i = 0; i < latest_plane_ids.size(); i++)
		planesFromFrontToBack.push_back(-1);

	// update the frontest plane if necessary
	for (int i = 0; i < latest_plane_ids.size(); i++) {
		bool is_frontier = true;
		int rank_id = 0;
		for (int j = 0; j < latest_plane_ids.size(); j++) {
			if (CGAL::parallel(CPlane(halfSpaceSet[latest_plane_ids[i]]), CPlane(halfSpaceSet[latest_plane_ids[j]])))
				continue;
			int candidate_plane_id = latest_plane_ids[i];
			int reference_plane_id = latest_plane_ids[j];
			if (FrontierPlane(CPoint(halfSpaceSet[base_plane_id].point1), CPoint(halfSpaceSet[base_plane_id].point2), CPoint(halfSpaceSet[base_plane_id].point3),
				CPoint(halfSpaceSet[partner_plane_id].point1), CPoint(halfSpaceSet[partner_plane_id].point2), CPoint(halfSpaceSet[partner_plane_id].point3),
				CPoint(halfSpaceSet[candidate_plane_id].point1), CPoint(halfSpaceSet[candidate_plane_id].point2), CPoint(halfSpaceSet[candidate_plane_id].point3),
				CPoint(halfSpaceSet[reference_plane_id].point1), CPoint(halfSpaceSet[reference_plane_id].point2), CPoint(halfSpaceSet[reference_plane_id].point3))
				!= CGAL::NEGATIVE) {
				is_frontier = false;
				rank_id++;
			}
		}

		planesFromFrontToBack[rank_id] = latest_plane_ids[i];
	}

	vector<int> trimmed;
	for (int i = 0; i < latest_plane_ids.size(); i++) {
		if (planesFromFrontToBack[i] == -1)
			break;
		trimmed.push_back(planesFromFrontToBack[i]);
	}
	return trimmed;

}

int KerGen_Exact::computeKernelVertex(int plane1_id, int plane2_id, vector<int> others) {

	auto auto_line = CGAL::intersection(CPlane(halfSpaceSet[plane1_id]), CPlane(halfSpaceSet[plane2_id]));
	EK::Line_3* line = boost::get<EK::Line_3>(&*auto_line);

	auto auto_point = CGAL::intersection(CPlane(halfSpaceSet[others[0]]), *line);
	EK::Point_3* point = boost::get<EK::Point_3>(&*auto_point);
	double found_point[3] = { CGAL::to_double(point->x()), CGAL::to_double(point->y()), CGAL::to_double(point->z()) };
/*
	double* scalar = findValueVectorSatisfiedByPoint(found_point, halfSpaceSet);
	for (int i = 0; i < halfSpaceSet.size(); i++)
		if (scalar[i] > EPSILON)
			cout << "NON-KERNEL\n";
	delete[] scalar;
*/
	int vertex_id = -1;
	for (int v = 0; v < kernel.getNumOfVerts(); v++)
		if (isTripleSame(kernel.getVertex(v).coords, found_point)) {
			vertex_id = v;
			break;
		}

	if (vertex_id == -1) {
		vertex_id = kernel.getNumOfVerts();
		kernel.addVertex(CGAL::to_double(point->x()), CGAL::to_double(point->y()), CGAL::to_double(point->z()));
		vector<int> vertex_generator_ids;
		vertex_generator_ids.push_back(plane1_id);
		vertex_generator_ids.push_back(plane2_id);
		for (int i = 0; i < others.size(); i++)
			vertex_generator_ids.push_back(others[i]);
		generatorsOfVertices.push_back(vertex_generator_ids);
	}

	return vertex_id;
}

bool KerGen_Exact::addIntoProcessQueueAfterCheck(int plane1_id, int plane2_id, int plane3_id, int plane4_id, int vertex_id) {

	bool is_valid = true;
	Plane plane_set1[2] = {Plane(halfSpaceSet[plane1_id]), Plane(halfSpaceSet[plane2_id]) };
	Plane plane_set2[2] = {Plane(halfSpaceSet[plane1_id]), Plane(halfSpaceSet[plane3_id]) };
	Line* line1 = find2PlaneIntersection(plane_set1);
	Line* line2 = find2PlaneIntersection(plane_set2);
	normalize(line1->directionVector);
	normalize(line2->directionVector);
	if (abs(dotProduct(line1->directionVector, line2->directionVector) - 1.0) < EPSILON)
		is_valid = false;
	else {
		if (faceProcess[plane1_id] == ProcessColor_Exact::RED) {
			faceProcess[plane1_id] = ProcessColor_Exact::WHITE;
			tobeProcessedFaceIds.push(plane1_id);
			vertexPartnerships[plane1_id].push_back(VertexPartners_Exact(plane2_id, plane4_id, vertex_id));
		}
	}

	delete line1;
	delete line2;

	return is_valid;

}




void KerGen_Exact::findInitialPoint_1(double* point) {

	//computeHalfSpacesFromTriangles(hostMeshptr->getAllTris(), hostMeshptr->getAllVerts(), halfSpaceSet);

	this->initialPoint = new double[18];
	double extremeDirections[6][3] = { {-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1} };
	for (int i = 0; i < 6; i++) {
		double* extremePoint = sdlpMain(extremeDirections[i], halfSpaceSet);	// compute initial kernel point at the given extreme direction
		if (extremePoint) {
			for (int j = 0; j < 3; j++)
				this->initialPoint[i * 3 + j] = extremePoint[j];
			delete[] extremePoint;
		}
		else {
			delete[] this->initialPoint;
			this->initialPoint = nullptr;
			break;
		}

	}

	if (this->initialPoint) {
		for (int i = 0; i < 3; i++)
			point[i] = 0;
		// find the center of the kernel's bonding box
		for (int i = 0; i < 6; i++)
			point[i / 2] += this->initialPoint[i * 3 + i / 2];
		for (int i = 0; i < 3; i++)
			point[i] /= 2.0;
	}
	else {
		for (int i = 0; i < 3; i++)
			point[i] = numeric_limits<double>::infinity();
	}

}




bool KerGen_Exact::ParallelPlanes(Point_3 plane1_p1, Point_3 plane1_p2, Point_3 plane1_p3,
					Point_3 plane2_p1, Point_3 plane2_p2, Point_3 plane2_p3) {
	
	Plane_3 plane1(plane1_p1, plane1_p2, plane1_p3);
	Plane_3 plane2(plane2_p1, plane2_p2, plane2_p3);

	if (CGAL::parallel(plane1, plane2))
		return true;
	return false;

}

bool KerGen_Exact::SamePlanes(Point_3 plane1_p1, Point_3 plane1_p2, Point_3 plane1_p3,
				Point_3 plane2_p1, Point_3 plane2_p2, Point_3 plane2_p3) {

	Plane_3 plane1(plane1_p1, plane1_p2, plane1_p3);
	Plane_3 plane2(plane2_p1, plane2_p2, plane2_p3);

	if (CGAL::parallel(plane1, plane2)) {
		if (CGAL::scalar_product(plane1.orthogonal_vector(), plane2.orthogonal_vector()) > 0)
			return true;
	}
	return false;

}

CGAL::Sign KerGen_Exact::IntersectablePlane(	Point_3 plane1_p1, Point_3 plane1_p2, Point_3 plane1_p3,
								Point_3 plane2_p1, Point_3 plane2_p2, Point_3 plane2_p3,
								Point_3 candidate_plane_p1, Point_3 candidate_plane_p2, Point_3 candidate_plane_p3) {

	Plane_3 plane1(plane1_p1, plane1_p2, plane1_p3);
	Plane_3 plane2(plane2_p1, plane2_p2, plane2_p3);

	Plane_3 candidate_plane(candidate_plane_p1, candidate_plane_p2, candidate_plane_p3);
	if (CGAL::parallel(candidate_plane, plane1) || CGAL::parallel(candidate_plane, plane2))
		return CGAL::NEGATIVE;	// not valid

	auto auto_line = CGAL::intersection(plane1, plane2);	// these two planes are assummed to be already intersectable
	Line_3* line = boost::get<Line_3>(&*auto_line);

	auto auto_candidate_point = CGAL::intersection(candidate_plane, *line);
	Point_3* candidate_point = boost::get<Point_3>(&*auto_candidate_point);

	if (candidate_point == NULL)
		return CGAL::NEGATIVE;	// not valid

	return CGAL::POSITIVE;

}

CGAL::Sign KerGen_Exact::IntersectableFrontPlane(	Point_3 plane1_p1, Point_3 plane1_p2, Point_3 plane1_p3,
									Point_3 plane2_p1, Point_3 plane2_p2, Point_3 plane2_p3,
									Point_3 candidate_plane_p1, Point_3 candidate_plane_p2, Point_3 candidate_plane_p3,
									Point_3 back_plane_p1, Point_3 back_plane_p2, Point_3 back_plane_p3) {

	Plane_3 plane1(plane1_p1, plane1_p2, plane1_p3);
	Plane_3 plane2(plane2_p1, plane2_p2, plane2_p3);

	Plane_3 candidate_plane(candidate_plane_p1, candidate_plane_p2, candidate_plane_p3);
	if (CGAL::parallel(candidate_plane, plane1) || CGAL::parallel(candidate_plane, plane2))
		return CGAL::POSITIVE;	// not valid

	auto auto_line = CGAL::intersection(plane1, plane2);	// these two planes are assummed to be already intersectable
	Line_3* line = boost::get<Line_3>(&*auto_line);

	auto auto_candidate_point = CGAL::intersection(candidate_plane, *line);
	Point_3* candidate_point = boost::get<Point_3>(&*auto_candidate_point);

	if (candidate_point == NULL)
		return CGAL::POSITIVE;	// not valid

	CGAL::Sign sign = CGAL::orientation(back_plane_p1, back_plane_p2, back_plane_p3, *candidate_point);
	if (sign == CGAL::POSITIVE)
		return CGAL::POSITIVE;
	if (sign == CGAL::COPLANAR) {
		Plane_3 back_plane(back_plane_p1, back_plane_p2, back_plane_p3);
		if (CGAL::parallel(back_plane, candidate_plane) &&
			CGAL::scalar_product(back_plane.orthogonal_vector(), candidate_plane.orthogonal_vector()) > 0)
			return CGAL::POSITIVE;
	}
	return CGAL::NEGATIVE;

}

CGAL::Sign KerGen_Exact::ClosestPlaneToPoint_3D(	Point_3 plane1_p1, Point_3 plane1_p2, Point_3 plane1_p3,
									Point_3 plane2_p1, Point_3 plane2_p2, Point_3 plane2_p3,
									Point_3 point) {

	Plane_3 plane1(plane1_p1, plane1_p2, plane1_p3);
	Plane_3 plane2(plane2_p1, plane2_p2, plane2_p3);

	RT distance1 = CGAL::squared_distance(plane1, point);
	RT distance2 = CGAL::squared_distance(plane2, point);
	return CGAL::compare(distance1, distance2);

}

CGAL::Sign KerGen_Exact::ClosestPlaneToPoint_2D(	Point_3 base_plane_p1, Point_3 base_plane_p2, Point_3 base_plane_p3,
									Point_3 plane1_p1, Point_3 plane1_p2, Point_3 plane1_p3,
									Point_3 plane2_p1, Point_3 plane2_p2, Point_3 plane2_p3,
									Point_3 initialpoint) {

	Plane_3 base_plane(base_plane_p1, base_plane_p2, base_plane_p3);
	Plane_3 plane1(plane1_p1, plane1_p2, plane1_p3);
	Plane_3 plane2(plane2_p1, plane2_p2, plane2_p3);		// it is assummed that none of these planes are parallel

	Point_3 point = base_plane.projection(initialpoint);

	auto auto_line1 = CGAL::intersection(base_plane, plane1);
	Line_3* line1 = boost::get<Line_3>(&*auto_line1);

	auto auto_line2 = CGAL::intersection(base_plane, plane2);
	Line_3* line2 = boost::get<Line_3>(&*auto_line2);

	RT distance1 = CGAL::squared_distance(*line1, point);
	RT distance2 = CGAL::squared_distance(*line2, point);
	return CGAL::compare(distance1, distance2);

}

CGAL::Comparison_result KerGen_Exact::ClosestPlaneToPoint_1D(	Point_3 plane1_p1, Point_3 plane1_p2, Point_3 plane1_p3,
												Point_3 plane2_p1, Point_3 plane2_p2, Point_3 plane2_p3,
												Point_3 candidate_plane_p1, Point_3 candidate_plane_p2, Point_3 candidate_plane_p3,
												Point_3 front_plane_p1, Point_3 front_plane_p2, Point_3 front_plane_p3,
												Point_3 back_plane_p1, Point_3 back_plane_p2, Point_3 back_plane_p3) {

	Plane_3 plane1(plane1_p1, plane1_p2, plane1_p3);
	Plane_3 plane2(plane2_p1, plane2_p2, plane2_p3);

	Plane_3 candidate_plane(candidate_plane_p1, candidate_plane_p2, candidate_plane_p3);
	if (CGAL::parallel(candidate_plane, plane1) || CGAL::parallel(candidate_plane, plane2))
		return CGAL::SMALLER;	// not valid

	auto auto_line = CGAL::intersection(plane1, plane2);
	Line_3* line = boost::get<Line_3>(&*auto_line);

	auto auto_candidate_point = CGAL::intersection(candidate_plane, *line);
	Point_3* candidate_point = boost::get<Point_3>(&*auto_candidate_point);

	if (candidate_point == NULL)
		return CGAL::SMALLER;

	CGAL::Sign sign = CGAL::orientation(back_plane_p1, back_plane_p2, back_plane_p3, *candidate_point);
	//cout << "SIGN: " << sign << endl;
	//Plane_3 back_plane2(back_plane_p1, back_plane_p2, back_plane_p3);
	//cout << "PLACEMENT: " << back_plane2.a() * candidate_point->x() + back_plane2.b() * candidate_point->y() + back_plane2.c() * candidate_point->z() + back_plane2.d() << endl;
	if (sign == CGAL::POSITIVE)
		return CGAL::SMALLER;
	if (sign == CGAL::COPLANAR) {
		Plane_3 back_plane(back_plane_p1, back_plane_p2, back_plane_p3);
		if (CGAL::parallel(back_plane, candidate_plane) &&
			CGAL::scalar_product(back_plane.orthogonal_vector(), candidate_plane.orthogonal_vector()) > 0)
			return CGAL::SMALLER;
	}

	Plane_3 front_plane(front_plane_p1, front_plane_p2, front_plane_p3);
	auto auto_front_point = CGAL::intersection(front_plane, *line);
	Point_3* front_point = boost::get<Point_3>(&*auto_front_point);

	Plane_3 back_plane(back_plane_p1, back_plane_p2, back_plane_p3);
	auto auto_back_point = CGAL::intersection(back_plane, *line);
	Point_3* back_point = boost::get<Point_3>(&*auto_back_point);

	RT distance1 = CGAL::squared_distance(*back_point, *front_point);
	RT distance2 = CGAL::squared_distance(*back_point, *candidate_point);
	return CGAL::compare(distance1, distance2);

}

CGAL::Comparison_result KerGen_Exact::ClosestPlaneToPoint_UOD(Point_3 plane1_p1, Point_3 plane1_p2, Point_3 plane1_p3,
												Point_3 plane2_p1, Point_3 plane2_p2, Point_3 plane2_p3,
												Point_3 candidate_plane_p1, Point_3 candidate_plane_p2, Point_3 candidate_plane_p3,
												Point_3 front_plane_p1, Point_3 front_plane_p2, Point_3 front_plane_p3,
												Point_3 initialpoint) {

	Plane_3 plane1(plane1_p1, plane1_p2, plane1_p3);
	Plane_3 plane2(plane2_p1, plane2_p2, plane2_p3);

	Plane_3 candidate_plane(candidate_plane_p1, candidate_plane_p2, candidate_plane_p3);
	if (CGAL::parallel(candidate_plane, plane1) || CGAL::parallel(candidate_plane, plane2))
		return CGAL::SMALLER;	// not valid

	auto auto_line = CGAL::intersection(plane1, plane2);
	Line_3* line = boost::get<Line_3>(&*auto_line);

	auto auto_candidate_point = CGAL::intersection(candidate_plane, *line);
	Point_3* candidate_point = boost::get<Point_3>(&*auto_candidate_point);

	if (candidate_point == NULL)
		return CGAL::SMALLER;

	Plane_3 front_plane(front_plane_p1, front_plane_p2, front_plane_p3);
	auto auto_front_point = CGAL::intersection(front_plane, *line);
	Point_3* front_point = boost::get<Point_3>(&*auto_front_point);

	Point_3 point = plane1.projection(initialpoint);
	point = line->projection(point);

	RT distance1 = CGAL::squared_distance(point, *front_point);
	RT distance2 = CGAL::squared_distance(point, *candidate_point);
	return CGAL::compare(distance1, distance2);

}

CGAL::Sign KerGen_Exact::FrontierPlane(	Point_3 base_plane_p1, Point_3 base_plane_p2, Point_3 base_plane_p3,
							Point_3 partner_plane_p1, Point_3 partner_plane_p2, Point_3 partner_plane_p3,
							Point_3 candidate_plane_p1, Point_3 candidate_plane_p2, Point_3 candidate_plane_p3,
							Point_3 reference_plane_p1, Point_3 reference_plane_p2, Point_3 reference_plane_p3) {

	Plane_3 base_plane(base_plane_p1, base_plane_p2, base_plane_p3);
	Plane_3 partner_plane(partner_plane_p1, partner_plane_p2, partner_plane_p3);
	Plane_3 candidate_plane(candidate_plane_p1, candidate_plane_p2, candidate_plane_p3);
	Plane_3 reference_plane(reference_plane_p1, reference_plane_p2, reference_plane_p3);

	if (CGAL::parallel(candidate_plane, reference_plane))	// then they are the same planes
		return CGAL::POSITIVE;

	auto auto_line = CGAL::intersection(base_plane, candidate_plane);
	Line_3* line = boost::get<Line_3>(&*auto_line);
	auto auto_found_point = CGAL::intersection(*line, partner_plane);
	Point_3* found_point = boost::get<Point_3>(&*auto_found_point);

	Vector_3 direction = line->to_vector();
	Point_3 test_point = *found_point + direction;

	if (CGAL::orientation(partner_plane_p1, partner_plane_p2, partner_plane_p3, test_point) != CGAL::NEGATIVE) {
		direction = line->opposite().to_vector();
		test_point = *found_point + direction;
	}
	return CGAL::orientation(reference_plane_p1, reference_plane_p2, reference_plane_p3, test_point);

}

EK::Point_3 KerGen_Exact::CPoint(double* point) {

	return EK::Point_3(point[0], point[1], point[2]);
}

EK::Plane_3 KerGen_Exact::CPlane(HalfSpace halfspace) {

	return EK::Plane_3(CPoint(halfspace.point1), CPoint(halfspace.point2), CPoint(halfspace.point3));
}






