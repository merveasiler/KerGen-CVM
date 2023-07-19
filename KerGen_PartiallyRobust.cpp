// @author Merve Asiler

#include "KerGen_PartiallyRobust.h"
#include "BaseGeoOpUtils.h"
#include "sdlp.h"
#include "CGALUtils.h"

KerGen_PartiallyRobust::KerGen_PartiallyRobust(const Mesh& hostMesh) :
	KernelExpansion(hostMesh, true) {

}

KerGen_PartiallyRobust::~KerGen_PartiallyRobust() {

}

void KerGen_PartiallyRobust::expandKernel() {

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
	vector<int> closest3_ids = findClosestHalfSpace(closest1_id, closest2_id);
	clock_t end4 = clock();
	cout << "PREDICATE TIME-3: " << double(end4 - begin4) * 1000 / CLOCKS_PER_SEC << "   " << closest3_ids[0] << endl;

	vector<int> ordered_closest3_ids = orderFromFrontToBack(closest1_id, closest2_id, closest3_ids);
	int closest3_id = ordered_closest3_ids[0];
	int vertex_id = computeKernelVertex(closest1_id, closest2_id, closest3_id);
	addIntoProcessQueueAfterCheck(closest1_id, closest2_id, closest3_id, closest3_id);

	if (ordered_closest3_ids.size() == 1)
		addIntoProcessQueueAfterCheck(closest2_id, closest3_id, closest1_id, closest1_id);
	else {
		ordered_closest3_ids.push_back(closest2_id);
		int i = ordered_closest3_ids.size() - 1;
		for (int j = i - 1; j >= 1; j--) {
			if (addIntoProcessQueueAfterCheck(ordered_closest3_ids[i], ordered_closest3_ids[j], ordered_closest3_ids[j - 1], closest1_id))
				i = j;
			if (j == 1)
				addIntoProcessQueueAfterCheck(ordered_closest3_ids[i], ordered_closest3_ids[j - 1], closest1_id, closest1_id);
		}
	}

	int number_of_kernel_faces = 0;
	int number_of_kernel_edges = 0;

	while (!tobeProcessedFaceIds.empty()) {
		int face_id = tobeProcessedFaceIds.front();
		for (int partner = 0; ; partner++) {
			// initializations
			VertexPartners_PartiallyRobust& vp = vertexPartnerships[face_id][partner];
			int partner1_id = vp.partner1_id;
			int partner2_id = vp.partner2_id;

			cout << face_id << " " << partner1_id << " " << partner2_id << endl;
			number_of_kernel_edges++;

			vector<int> closest_plane_ids = findClosestHalfSpace(face_id, partner1_id, partner2_id);
			vector<int> ordered_closest_plane_ids = orderFromFrontToBack(face_id, partner1_id, closest_plane_ids);

			int next_partner_id = ordered_closest_plane_ids[0];
			int initial_partner_id = vertexPartnerships[face_id][0].partner1_id;
			if (abs(dotProduct(halfSpaceSet[next_partner_id].ABCD, halfSpaceSet[initial_partner_id].ABCD) - 1.0) < EPSILON)
				break;	// loop completed

			int found_vertex_id = computeKernelVertex(face_id, partner1_id, next_partner_id);
			vertexPartnerships[face_id].push_back(VertexPartners_PartiallyRobust(next_partner_id, partner1_id));

			if (ordered_closest_plane_ids.size() == 1)
				addIntoProcessQueueAfterCheck(partner1_id, next_partner_id, face_id, face_id);
			else {
				ordered_closest_plane_ids.push_back(partner1_id);
				int i = ordered_closest_plane_ids.size() - 1;
				for (int j = i - 1; j >= 1; j--) {
					if (addIntoProcessQueueAfterCheck(ordered_closest_plane_ids[i], ordered_closest_plane_ids[j], ordered_closest_plane_ids[j - 1], face_id))
						i = j;
					if (j == 1)
						addIntoProcessQueueAfterCheck(ordered_closest_plane_ids[i], ordered_closest_plane_ids[j - 1], face_id, face_id);
				}
			}

		}

		faceProcess[face_id] = ProcessColor_PartiallyRobust::GREEN;
		tobeProcessedFaceIds.pop();
		number_of_kernel_faces++;

	}


	//	cout << "[number of kernel edges: " << number_of_kernel_edges << "], number of kernel faces: " << number_of_kernel_faces << "]" << endl;

	if (kernel.getNumOfVerts() > 0)
		kernel = computeConvexHull(kernel.getAllVerts());

}

void KerGen_PartiallyRobust::initialize() {

	for (int i = 0; i < halfSpaceSet.size(); i++) {
		vertexPartnerships.push_back(vector<VertexPartners_PartiallyRobust>());
		faceProcess.push_back(ProcessColor_PartiallyRobust::RED);
		cgalPlaneSet.push_back(XK::Plane_3(CPoint(halfSpaceSet[i].point1), CPoint(halfSpaceSet[i].point2), CPoint(halfSpaceSet[i].point3)));
	}

}

int KerGen_PartiallyRobust::findClosestHalfSpace(double* point) {

	// initialize the closest plane temporarily
	int closest_plane_id = 0;
	XK::Point_3 cpoint = CPoint(point);

	// update the closest plane if necessary
	for (int i = 0; i < halfSpaceSet.size(); i++) {
		if (ClosestPlaneToPoint_3D(i, cpoint)) {
			closest_plane_id = i;
			break;
		}
	}

	return closest_plane_id;

}

int KerGen_PartiallyRobust::findClosestHalfSpace(double* point, int base_plane_id) {

	XK::Point_3 cpoint = CPoint(point);
	XK::Plane_3 base_plane = cgalPlaneSet[base_plane_id];
	cpoint = base_plane.projection(cpoint);

	// initialize the closest plane temporarily
	int closest_plane_id = 0;
	for (int i = 0; i < halfSpaceSet.size(); i++) {
		if (ClosestPlaneToPoint_2D(i, base_plane_id, cpoint)) {
			closest_plane_id = i;
			break;
		}
	}

	return closest_plane_id;

}

vector<int> KerGen_PartiallyRobust::findClosestHalfSpace(int base_plane_id, int partner_plane_id) {

	XK::Plane_3 base_plane = cgalPlaneSet[base_plane_id];
	XK::Plane_3 partner_plane = cgalPlaneSet[partner_plane_id];
	auto auto_line = CGAL::intersection(base_plane, partner_plane);
	XK::Line_3* line = boost::get<XK::Line_3>(&*auto_line);

	// initialize the closest plane temporarily
	vector<int> closest_plane_ids;
	for (int i = 0; i < halfSpaceSet.size(); i++) {
		if (ClosestPlaneToPoint_1D(i, base_plane_id, partner_plane_id, *line))
			closest_plane_ids.push_back(i);
	}
	/*
		auto auto_candidate_point = CGAL::intersection(*line, cgalPlaneSet[closest_plane_ids[0]]);
		XK::Point_3* candidate_point = boost::get<XK::Point_3>(&*auto_candidate_point);
		double found_point[3] = {CGAL::to_double(candidate_point->x()), CGAL::to_double(candidate_point->y()), CGAL::to_double(candidate_point->z()) };
		double* scalar = findValueVectorSatisfiedByPoint(found_point, halfSpaceSet);
		for (int i = 0; i < halfSpaceSet.size(); i++)
			if (scalar[i] > EPSILON)
				cout << "NON-KERNEL\n";
		delete[] scalar;
	*/
	return closest_plane_ids;

}

vector<int> KerGen_PartiallyRobust::findClosestHalfSpace(int base_plane_id, int partner_plane_id, int back_plane_id) {

	XK::Plane_3 base_plane = cgalPlaneSet[base_plane_id];
	XK::Plane_3 partner_plane = cgalPlaneSet[partner_plane_id];
	auto auto_line = CGAL::intersection(base_plane, partner_plane);
	XK::Line_3* line = boost::get<XK::Line_3>(&*auto_line);

	// initialize the closest plane temporarily
	vector<int> closest_plane_ids;
	for (int i = 0; i < halfSpaceSet.size(); i++) {
		if (ClosestPlaneToPoint_0D(i, base_plane_id, partner_plane_id, back_plane_id, *line))
			closest_plane_ids.push_back(i);
	}
	/*
		auto auto_candidate_point = CGAL::intersection(*line, cgalPlaneSet[closest_plane_ids[0]]);
		XK::Point_3* candidate_point = boost::get<XK::Point_3>(&*auto_candidate_point);
		double found_point[3] = { CGAL::to_double(candidate_point->x()), CGAL::to_double(candidate_point->y()), CGAL::to_double(candidate_point->z()) };
		double* scalar = findValueVectorSatisfiedByPoint(found_point, halfSpaceSet);
		for (int i = 0; i < halfSpaceSet.size(); i++)
			if (scalar[i] > EPSILON)
				cout << "NON-KERNEL\n";
		delete[] scalar;
	*/
	return closest_plane_ids;
}

vector<int> KerGen_PartiallyRobust::orderFromFrontToBack(int base_plane_id, int partner_plane_id, vector<int> latest_plane_ids) {

	vector<int> planesFromFrontToBack;
	for (int i = 0; i < latest_plane_ids.size(); i++)
		planesFromFrontToBack.push_back(-1);

	XK::Point_3 subject_point;
	bool is_subject_point_picked = false;
	vector<XK::Line_3> lines;

	// update the frontest plane if necessary
	for (int i = 0; i < latest_plane_ids.size(); i++) {
		int candidate_plane_id = latest_plane_ids[i];

		XK::Plane_3 base_plane = cgalPlaneSet[base_plane_id];
		XK::Plane_3 candidate_plane = cgalPlaneSet[candidate_plane_id];
		XK::Plane_3 partner_plane = cgalPlaneSet[partner_plane_id];

		auto auto_line = CGAL::intersection(base_plane, candidate_plane);
		XK::Line_3* line = boost::get<XK::Line_3>(&*auto_line);
		auto auto_found_point = CGAL::intersection(*line, partner_plane);
		XK::Point_3* found_point = boost::get<XK::Point_3>(&*auto_found_point);

		if (!is_subject_point_picked) {
			subject_point = *found_point;
			is_subject_point_picked = true;
		}
		else {
			if (abs(CGAL::to_double(subject_point.x() - found_point->x())) < EPSILON &&
				abs(CGAL::to_double(subject_point.y() - found_point->y())) < EPSILON &&
				abs(CGAL::to_double(subject_point.z() - found_point->z())) < EPSILON)
				;	// OK
			else
				latest_plane_ids[i] = -1;	// REMOVE
		}

		lines.push_back(*line);
	}

	for (int i = 0; i < latest_plane_ids.size(); i++) {
		if (latest_plane_ids[i] == -1)
			continue;
		int rank_id = 0;
		int candidate_plane_id = latest_plane_ids[i];

		XK::Vector_3 direction = lines[i].to_vector();
		XK::Point_3 test_point = subject_point + direction;
		if (CGAL::orientation(CPoint(halfSpaceSet[partner_plane_id].point1),
			CPoint(halfSpaceSet[partner_plane_id].point2),
			CPoint(halfSpaceSet[partner_plane_id].point3),
			test_point) == CGAL::POSITIVE) {
			direction = lines[i].opposite().to_vector();
			test_point = subject_point + direction;
		}

		for (int j = 0; j < latest_plane_ids.size(); j++) {
			int reference_plane_id = latest_plane_ids[j];
			if (candidate_plane_id == reference_plane_id)
				continue;
			if (reference_plane_id == -1)
				continue;

			if (CGAL::parallel(cgalPlaneSet[candidate_plane_id], cgalPlaneSet[reference_plane_id])) {
				if (CGAL::scalar_product(cgalPlaneSet[candidate_plane_id].orthogonal_vector(),
					cgalPlaneSet[reference_plane_id].orthogonal_vector()) > 0)
					// same planes;
					continue;
				else
					cout << "DUE TO WRONG SELECTION OF PARTNER PLANE\n";
			}

			CGAL::Sign sign = CGAL::orientation(CPoint(halfSpaceSet[reference_plane_id].point1),
				CPoint(halfSpaceSet[reference_plane_id].point2),
				CPoint(halfSpaceSet[reference_plane_id].point3),
				test_point);
			if (sign == CGAL::POSITIVE)
				rank_id++;
		}

		planesFromFrontToBack[rank_id] = candidate_plane_id;
	}

	vector<int> trimmed;
	for (int i = 0; i < latest_plane_ids.size(); i++) {
		if (planesFromFrontToBack[i] == -1)
			break;
		trimmed.push_back(planesFromFrontToBack[i]);
	}
	return trimmed;

}

int KerGen_PartiallyRobust::computeKernelVertex(int plane1_id, int plane2_id, int plane3_id) {

	auto auto_line = CGAL::intersection(cgalPlaneSet[plane1_id], cgalPlaneSet[plane2_id]);
	XK::Line_3* line = boost::get<XK::Line_3>(&*auto_line);
	auto auto_point = CGAL::intersection(cgalPlaneSet[plane3_id], *line);
	XK::Point_3* point = boost::get<XK::Point_3>(&*auto_point);
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
	}

	return vertex_id;
}

bool KerGen_PartiallyRobust::addIntoProcessQueueAfterCheck(int plane1_id, int plane2_id, int plane3_id, int base_plane_id) {

	auto auto_line1 = CGAL::intersection(cgalPlaneSet[plane1_id], cgalPlaneSet[plane2_id]);
	XK::Line_3* line1 = boost::get<XK::Line_3>(&*auto_line1);

	auto auto_line2 = CGAL::intersection(cgalPlaneSet[plane1_id], cgalPlaneSet[plane3_id]);
	XK::Line_3* line2 = boost::get<XK::Line_3>(&*auto_line2);

	bool is_valid = false;

	if (!CGAL::parallel(*line1, *line2)) {
		if (faceProcess[plane1_id] == ProcessColor_PartiallyRobust::RED) {
			faceProcess[plane1_id] = ProcessColor_PartiallyRobust::WHITE;
			tobeProcessedFaceIds.push(plane1_id);
			vertexPartnerships[plane1_id].push_back(VertexPartners_PartiallyRobust(plane2_id, base_plane_id));
		}
		is_valid = true;
	}

	return is_valid;

}

void KerGen_PartiallyRobust::findInitialPoint_1(double* point) {

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



bool KerGen_PartiallyRobust::ClosestPlaneToPoint_3D(int candidate_plane_id, XK::Point_3 point) {

	XK::Plane_3 candidate_plane = cgalPlaneSet[candidate_plane_id];
	XK::Point_3 candidate_point = candidate_plane.projection(point);

	for (int i = 0; i < halfSpaceSet.size(); i++) {
		CGAL::Sign sign = CGAL::orientation(CPoint(halfSpaceSet[i].point1), CPoint(halfSpaceSet[i].point2), CPoint(halfSpaceSet[i].point3), candidate_point);
		if (sign == CGAL::POSITIVE)
			return false;
	}

	return true;

}

bool KerGen_PartiallyRobust::ClosestPlaneToPoint_2D(int candidate_plane_id, int base_plane_id, XK::Point_3 point) {

	XK::Plane_3 candidate_plane = cgalPlaneSet[candidate_plane_id];
	XK::Plane_3 base_plane = cgalPlaneSet[base_plane_id];

	if (CGAL::parallel(candidate_plane, base_plane))
		return false;

	auto auto_line = CGAL::intersection(base_plane, candidate_plane);
	XK::Line_3* line = boost::get<XK::Line_3>(&*auto_line);

	if (line == NULL)
		return false;

	XK::Point_3 candidate_point = line->projection(point);

	for (int i = 0; i < halfSpaceSet.size(); i++) {
		CGAL::Sign sign = CGAL::orientation(CPoint(halfSpaceSet[i].point1), CPoint(halfSpaceSet[i].point2), CPoint(halfSpaceSet[i].point3), candidate_point);
		if (sign == CGAL::POSITIVE)
			return false;
	}

	return true;

}

bool KerGen_PartiallyRobust::ClosestPlaneToPoint_1D(int candidate_plane_id, int base_plane_id, int partner_plane_id, XK::Line_3 line) {

	XK::Plane_3 candidate_plane = cgalPlaneSet[candidate_plane_id];
	XK::Plane_3 base_plane = cgalPlaneSet[base_plane_id];
	XK::Plane_3 partner_plane = cgalPlaneSet[partner_plane_id];

	if (CGAL::parallel(candidate_plane, base_plane) || CGAL::parallel(candidate_plane, partner_plane))
		return false;

	auto auto_candidate_point = CGAL::intersection(candidate_plane, line);
	XK::Point_3* candidate_point = boost::get<XK::Point_3>(&*auto_candidate_point);

	if (candidate_point == NULL)
		return false;

	for (int i = 0; i < halfSpaceSet.size(); i++) {
		CGAL::Sign sign = CGAL::orientation(CPoint(halfSpaceSet[i].point1), CPoint(halfSpaceSet[i].point2), CPoint(halfSpaceSet[i].point3), *candidate_point);
		if (sign == CGAL::POSITIVE)
			return false;
	}

	return true;

}

bool KerGen_PartiallyRobust::ClosestPlaneToPoint_0D(int candidate_plane_id, int base_plane_id, int partner_plane_id, int back_plane_id, XK::Line_3 line) {

	XK::Plane_3 candidate_plane = cgalPlaneSet[candidate_plane_id];
	XK::Plane_3 base_plane = cgalPlaneSet[base_plane_id];
	XK::Plane_3 partner_plane = cgalPlaneSet[partner_plane_id];
	XK::Plane_3 back_plane = cgalPlaneSet[back_plane_id];

	if (CGAL::parallel(candidate_plane, base_plane) || CGAL::parallel(candidate_plane, partner_plane))
		return false;

	if (abs(dotProduct(halfSpaceSet[candidate_plane_id].ABCD, halfSpaceSet[back_plane_id].ABCD) - 1.0) < EPSILON)
		return false;

	auto auto_candidate_point = CGAL::intersection(candidate_plane, line);
	XK::Point_3* candidate_point = boost::get<XK::Point_3>(&*auto_candidate_point);
	if (candidate_point == NULL)
		return false;

	CGAL::Sign sign = CGAL::orientation(CPoint(halfSpaceSet[back_plane_id].point1),
		CPoint(halfSpaceSet[back_plane_id].point2),
		CPoint(halfSpaceSet[back_plane_id].point3), *candidate_point);
	if (sign != CGAL::NEGATIVE)
		return false;

	for (int i = 0; i < halfSpaceSet.size(); i++) {
		CGAL::Sign sign = CGAL::orientation(CPoint(halfSpaceSet[i].point1), CPoint(halfSpaceSet[i].point2), CPoint(halfSpaceSet[i].point3), *candidate_point);
		if (sign == CGAL::POSITIVE)
			return false;
	}

	return true;

}

XK::Point_3 KerGen_PartiallyRobust::CPoint(double* point) {

	return XK::Point_3(point[0], point[1], point[2]);
}

XK::Plane_3 KerGen_PartiallyRobust::CPlane(HalfSpace halfspace) {

	return XK::Plane_3(CPoint(halfspace.point1), CPoint(halfspace.point2), CPoint(halfspace.point3));
}






