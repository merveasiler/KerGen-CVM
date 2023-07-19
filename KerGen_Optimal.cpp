// @author Merve Asiler

#include "KerGen_Optimal.h"
#include "BaseGeoOpUtils.h"
#include "sdlp.h"
#include "CGALUtils.h"

KerGen_Optimal::KerGen_Optimal(const Mesh& hostMesh) :
	KernelExpansion(hostMesh, true) {

}

KerGen_Optimal::~KerGen_Optimal() {

}

void KerGen_Optimal::expandKernel() {

	//cout << "Number of bounding planes: " << halfSpaceSet.size() << " out of " << hostMeshptr->getNumOfTris() << " triangles." << endl;

	initialize();

	double point[3];
	if (this->initialPoint == NULL)
		return;		// NOT STAR-SHAPED!
	for (int i = 0; i < 3; i++)
		point[i] = this->initialPoint[i];

	/*
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
	*/

	int closest1_id = findClosestHalfSpace(point);
	int closest2_id = findClosestHalfSpace(point, closest1_id);
	vector<int> closest3_ids = findClosestHalfSpace(point, closest1_id, closest2_id);

	vector<int> ordered_closest3_ids = orderFromFrontToBack(point, closest1_id, closest2_id, closest3_ids);
	int closest3_id = ordered_closest3_ids[0];
	int vertex_id = recordFoundVertex(point);
	addIntoProcessQueueAfterCheck(closest1_id, closest2_id, -1, closest3_id, point);

	if (ordered_closest3_ids.size() == 1)
		addIntoProcessQueueAfterCheck(closest2_id, closest3_id, -1, closest1_id, point);
	else {
		ordered_closest3_ids.push_back(closest2_id);
		int i = ordered_closest3_ids.size() - 1;
		for (int j = i - 1; j >= 1; j--) {
			if (addIntoProcessQueueAfterCheck(ordered_closest3_ids[i], ordered_closest3_ids[j], ordered_closest3_ids[j - 1], closest1_id, point))
				i = j;
			if (j == 1)
				addIntoProcessQueueAfterCheck(ordered_closest3_ids[i], ordered_closest3_ids[j - 1], -1, closest1_id, point);
		}
	}

	int number_of_kernel_faces = 0;
	int number_of_kernel_edges = 0;

	while (!tobeProcessedFaceIds.empty()) {
		int face_id = tobeProcessedFaceIds.front();
		for (int partner = 0; ; partner++) {
			// initializations
			VertexPartners_Optimal& vp = vertexPartnerships[face_id][partner];
			int partner1_id = vp.partner1_id;
			int partner2_id = vp.partner2_id;
			Line common_line_with_partner1 = vp.common_line_with_partner1;

			cout << face_id << " " << partner1_id << " " << partner2_id << endl;
			number_of_kernel_edges++;

			double point[3];
			vector<int> closest_plane_ids = findClosestHalfSpace(point, face_id, partner1_id, partner2_id, common_line_with_partner1);
			vector<int> ordered_closest_plane_ids = orderFromFrontToBack(point, face_id, partner1_id, closest_plane_ids);

			int next_partner_id = ordered_closest_plane_ids[0];
			int initial_partner_id = vertexPartnerships[face_id][0].partner1_id;
			if (abs(dotProduct(halfSpaceSet[next_partner_id].ABCD, halfSpaceSet[initial_partner_id].ABCD) - 1.0) < BIG_EPSILON)
				break;	// loop completed

			int latest_num_of_verts = kernel.getNumOfVerts();
			int found_vertex_id = recordFoundVertex(point);

			double next_line_direction[3];
			crossProduct(halfSpaceSet[face_id].ABCD, halfSpaceSet[next_partner_id].ABCD, next_line_direction);
			Line next_line(point, next_line_direction);
			vertexPartnerships[face_id].push_back(VertexPartners_Optimal(next_partner_id, partner1_id, next_line));

			if (latest_num_of_verts > found_vertex_id)
				continue;	// previously processed vertex

			if (ordered_closest_plane_ids.size() == 1)
				addIntoProcessQueueAfterCheck(partner1_id, next_partner_id, -1, face_id, point);
			else {
				ordered_closest_plane_ids.push_back(partner1_id);
				int i = ordered_closest_plane_ids.size() - 1;
				for (int j = i - 1; j >= 1; j--) {
					if (addIntoProcessQueueAfterCheck(ordered_closest_plane_ids[i], ordered_closest_plane_ids[j], ordered_closest_plane_ids[j - 1], face_id, point))
						i = j;
					if (j == 1)
						addIntoProcessQueueAfterCheck(ordered_closest_plane_ids[i], ordered_closest_plane_ids[j - 1], -1, face_id, point);
				}
			}

		}

		faceProcess[face_id] = ProcessColor_Optimal::GREEN;
		tobeProcessedFaceIds.pop();
		number_of_kernel_faces++;

	}


	//	cout << "[number of kernel edges: " << number_of_kernel_edges << "], number of kernel faces: " << number_of_kernel_faces << "]" << endl;

	if (kernel.getNumOfVerts() > 0)
		kernel = computeConvexHull(kernel.getAllVerts());

}

void KerGen_Optimal::initialize() {

	for (int i = 0; i < halfSpaceSet.size(); i++) {
		vertexPartnerships.push_back(vector<VertexPartners_Optimal>());
		faceProcess.push_back(ProcessColor_Optimal::RED);
	}

}

int KerGen_Optimal::findClosestHalfSpace(double* point) {

	int closest_plane_id;
	double closest_point[3];

	// test each plane one by one whether it can be base plane or not (base plane: the first plane )
	// when you found one, stop
	for (int i = 0; i < halfSpaceSet.size(); i++) {
		if (ClosestPlaneToPoint_3D(i, point, closest_point)) {
			closest_plane_id = i;
			break;
		}
	}

	// replace the initial point with its projection on the selected base plane
	for (int i = 0; i < 3; i++)
		point[i] = closest_point[i];

	return closest_plane_id;

}

int KerGen_Optimal::findClosestHalfSpace(double* point, int base_plane_id) {

	int closest_plane_id;
	double closest_point[3];

	// test each plane one by one whether it can be a partner to the given base plane, or not
	// when you found one, stop
	for (int i = 0; i < halfSpaceSet.size(); i++) {
		if (ClosestPlaneToPoint_2D(i, base_plane_id, point, closest_point)) {
			closest_plane_id = i;
			break;
		}
	}

	// replace the given base point with its projection on the common line of the base plane and the selected partner plane
	for (int i = 0; i < 3; i++)
		point[i] = closest_point[i];

	return closest_plane_id;

}

vector<int> KerGen_Optimal::findClosestHalfSpace(double* point, int base_plane_id, int partner_plane_id) {

	// find the common line of the base plane and partner plane
	double common_line_direction[3];
	crossProduct(halfSpaceSet[base_plane_id].ABCD, halfSpaceSet[partner_plane_id].ABCD, common_line_direction);
	Line line(point, common_line_direction);

	vector<int> closest_plane_ids;
	double closest_point[3], candidate_point[3];
	bool point_recorded = false;

	// test each plane one by one whether it can be a third plane to construct a vertex with the base plane and partner plane, or not
	// find each one that satisfy the condition
	for (int i = 0; i < halfSpaceSet.size(); i++) {
		if (ClosestPlaneToPoint_1D(i, base_plane_id, partner_plane_id, line, candidate_point)) {

			// we find the intersections (vertices) at both end of the common line
			// we record the ones locating at only one of those points
			// simply we record the group of planes that giving the point we first met

			if (point_recorded) {
				if (isTripleSame(closest_point, candidate_point))	// if it is not the first point we met, then eliminate it
					closest_plane_ids.push_back(i);
			}
			else {
				for (int j = 0; j < 3; j++)
					closest_point[j] = candidate_point[j];
				point_recorded = true;
				closest_plane_ids.push_back(i);
			}
		}
	}

	// replace the given point with the vertex at the intersection point of the three planes
	for (int i = 0; i < 3; i++)
		point[i] = closest_point[i];

	return closest_plane_ids;

}

vector<int> KerGen_Optimal::findClosestHalfSpace(double* point, int base_plane_id, int partner_plane_id, int back_plane_id, Line line) {

	vector<int> closest_plane_ids;
	double closest_point[3], candidate_point[3];
	double smallest_distance = numeric_limits<double>::infinity();

	// test each plane one by one whether it can be a third plane to construct a vertex with the base plane and partner plane
	//		and not locating at the point of back_plane, or not
	// find each one that satisfy the condition
	for (int i = 0; i < halfSpaceSet.size(); i++) {
		double previous_distance = smallest_distance;
		if (ClosestPlaneToPoint_0D(i, base_plane_id, partner_plane_id, back_plane_id, line, candidate_point, smallest_distance)) {

			// if previously found a plane at this distance, just add the new plane id
			if (abs(previous_distance - smallest_distance) < BIG_EPSILON)
				closest_plane_ids.push_back(i);

			// if the new distance is much smaller than the previous distance, then reset everything, initialize the new foundings
			else {
				closest_plane_ids.clear();
				closest_plane_ids.push_back(i);

				for (int j = 0; j < 3; j++)
					closest_point[j] = candidate_point[j];
			}
		}
	}

	// replace the given point with the vertex at the intersection point of the three planes
	for (int i = 0; i < 3; i++)
		point[i] = closest_point[i];

	return closest_plane_ids;
}

vector<int> KerGen_Optimal::orderFromFrontToBack(double* point, int base_plane_id, int partner_plane_id, vector<int> latest_plane_ids) {

	if (latest_plane_ids.size() == 1)
		return latest_plane_ids;

	// reserve an empty array 
	vector<int> planesFromFrontToBack;
	for (int i = 0; i < latest_plane_ids.size(); i++)
		planesFromFrontToBack.push_back(-1);

	// first prepare the common lines with the base plane, we will use these lines
	vector<Line> lines;	// save those lines into this vector
	for (int i = 0; i < latest_plane_ids.size(); i++) {
		int candidate_plane_id = latest_plane_ids[i];

		double candidate_line_direction[3];
		crossProduct(halfSpaceSet[base_plane_id].ABCD, halfSpaceSet[candidate_plane_id].ABCD, candidate_line_direction);
		Line candidate_line(point, candidate_line_direction);
		lines.push_back(candidate_line);
	}

	// we will use a test point on the common line of each candidate-base double
	// this test point will be just a bit ahead of the intersection point of candidate-base-partner triple
	for (int i = 0; i < latest_plane_ids.size(); i++) {
		int candidate_plane_id = latest_plane_ids[i];
		int rank_id = 0;

		// initialize the test point, note that this point may be wrong if the direction vector was taken reversed
		double test_point[3];
		for (int j = 0; j < 3; j++)
			test_point[j] = point[j] + lines[i].directionVector[j];

		// check if the test point is wrong or not by controlling whether it is behid the partner plane or not. It should fall in front of the partner plane
		double scalar = halfSpaceSet[partner_plane_id].ABCD[3];
		for (int j = 0; j < 3; j++)
			scalar += halfSpaceSet[partner_plane_id].ABCD[j] * test_point[j];
		if (scalar > BIG_EPSILON) {	// if it falls behind the partner plane, recompute the test point by taking the reverse direction of the common line direction
			for (int j = 0; j < 3; j++)
				test_point[j] = point[j] - lines[i].directionVector[j];
		}

		// now check whether the test point falls behid the other candidate planes or not
		// if falls behind, then it means the current candidate plane is at the backside of that other candidate plane
		//	in that case, increase its rank, it is at least one plane back (en az 1 plane arka tarafta)
		for (int j = 0; j < latest_plane_ids.size(); j++) {
			int reference_plane_id = latest_plane_ids[j];

			double scalar = halfSpaceSet[reference_plane_id].ABCD[3];
			for (int j = 0; j < 3; j++)
				scalar += halfSpaceSet[reference_plane_id].ABCD[j] * test_point[j];
			if (scalar > BIG_EPSILON)
				rank_id++;
		}

		planesFromFrontToBack[rank_id] = candidate_plane_id;
	}

	// if some of the planes are the same, we take only one of them, which means there occurs some empty seats towards the last seats, trim them
	vector<int> trimmed;
	for (int i = 0; i < latest_plane_ids.size(); i++) {
		if (planesFromFrontToBack[i] == -1)
			break;
		trimmed.push_back(planesFromFrontToBack[i]);
	}
	return trimmed;

}

int KerGen_Optimal::recordFoundVertex(double* vertex) {

	// if this vertex was found before
	int vertex_id = -1;
	for (int v = 0; v < kernel.getNumOfVerts(); v++)
		if (isTripleSame(kernel.getVertex(v).coords, vertex)) {
			vertex_id = v;
			break;
		}

	if (vertex_id == -1) {
		// check if this new vertex inside the kernel or not, theoretically it should be, but it may not be due to float arithmetic
		double* scalar = findValueVectorSatisfiedByPoint(vertex, halfSpaceSet);
		for (int i = 0; i < halfSpaceSet.size(); i++)
			if (scalar[i] > EPSILON) {	// not inside the kernel
				// satisfy robustness by updating partner plane with the i^th plane
				cout << "NON-KERNEL\n";
				delete[] scalar;
				return -i - 1;	// return the id of the i^th plane in negative not to make confusion with a regular kernel vertex id
			}
		delete[] scalar;

		// there is no problem, record the kernel vertex
		vertex_id = kernel.getNumOfVerts();
		kernel.addVertex(vertex[0], vertex[1], vertex[2]);
	}

	return vertex_id;
}

bool KerGen_Optimal::addIntoProcessQueueAfterCheck(int plane1_id, int plane2_id, int plane3_id, int base_plane_id, double* point) {

	// we will add plane1 and plane2 neighborhood into the queue according to the below tests:

	// first prepare the common line of plane 1 and plane2
	double line1_direction[3];
	crossProduct(halfSpaceSet[plane1_id].ABCD, halfSpaceSet[plane2_id].ABCD, line1_direction);
	normalize(line1_direction);

	bool is_valid = false;

	if (plane3_id == -1) {
		is_valid = true;	// no need to extra check, this is already decided to be a valid neighborgood to add into the queue
	}
	else {
		// check whether the common line of plane1 and plane2 is the same of the common line of plane1 and plane3,
		// if so, then it means plane2 is invalid because of book-type neighborhood with plane3. 
		// (note that we can say so since we already ordered them before and plane2 was behind of plane3)
		double line2_direction[3];
		crossProduct(halfSpaceSet[plane1_id].ABCD, halfSpaceSet[plane3_id].ABCD, line2_direction);
		normalize(line2_direction);

		// are the common lines the same (namely parallel)?
		if (abs(abs(dotProduct(line1_direction, line2_direction)) - 1.0) > EPSILON)
			is_valid = true;
	}

	if (is_valid) {	// if not parallel then we can safely add this neighborhood into queue if there is no record related with plane1 id
		if (faceProcess[plane1_id] == ProcessColor_Optimal::RED) {
			faceProcess[plane1_id] = ProcessColor_Optimal::WHITE;
			tobeProcessedFaceIds.push(plane1_id);

			Line line(point, line1_direction);
			vertexPartnerships[plane1_id].push_back(VertexPartners_Optimal(plane2_id, base_plane_id, line));
		}
		is_valid = true;
	}

	return is_valid;

}



bool KerGen_Optimal::ClosestPlaneToPoint_3D(int candidate_plane_id, double* point, double* candidate_point) {

	// find the projection of the given initial point on the candidate base plane
	projectVertexOnPlane(halfSpaceSet[candidate_plane_id], point, candidate_point);

	// check if this point inside the kernel or not
	double* scalar = findValueVectorSatisfiedByPoint(candidate_point, halfSpaceSet);
	for (int i = 0; i < halfSpaceSet.size(); i++)
		if (scalar[i] > EPSILON) {
			delete[] scalar;
			return false;
		}

	delete[] scalar;
	return true;

}

bool KerGen_Optimal::ClosestPlaneToPoint_2D(int candidate_plane_id, int base_plane_id, double* point, double* candidate_point) {

	// check if the candidate plane is parallel to the base plane or not
	if (abs(abs(dotProduct(halfSpaceSet[candidate_plane_id].ABCD, halfSpaceSet[base_plane_id].ABCD)) - 1.0) < EPSILON)
		return false;	// parallel planes

	// purpose: project the given point onto the common line of candidate plane and base plane
	// how: 
	//	1. Find the direction of the common of base plane and candidate plane
	//	2. Find the direction which is perpendicular to the common line and locating on the base plane
	//		issue 1
	double common_line_direction[3];
	crossProduct(halfSpaceSet[base_plane_id].ABCD, halfSpaceSet[candidate_plane_id].ABCD, common_line_direction);
	//		issue 2
	double ray_direction[3];
	crossProduct(halfSpaceSet[base_plane_id].ABCD, common_line_direction, ray_direction);
	Line ray(point, ray_direction);

	// find the intersection parameter (distance) from the given point to the target point (which is on the common line of planes)
	double t = findLinePlaneIntersection(ray, halfSpaceSet[candidate_plane_id]);

	// compute the target (projection) point by using the found parameter
	for (int i = 0; i < 3; i++)
		candidate_point[i] = point[i] + t * ray_direction[i];

	// check if this point inside the kernel or not
	double* scalar = findValueVectorSatisfiedByPoint(candidate_point, halfSpaceSet);
	for (int i = 0; i < halfSpaceSet.size(); i++)
		if (scalar[i] > EPSILON) {
			delete[] scalar;
			return false;
		}

	delete[] scalar;
	return true;

}

bool KerGen_Optimal::ClosestPlaneToPoint_1D(int candidate_plane_id, int base_plane_id, int partner_plane_id, Line line, double* candidate_point) {

	// check if the candidate plane is parallel to the base plane or not
	if (abs(abs(dotProduct(halfSpaceSet[candidate_plane_id].ABCD, halfSpaceSet[base_plane_id].ABCD)) - 1.0) < EPSILON)
		return false;	// parallel planes
	// check if the candidate plane is parallel to the partner plane or not
	if (abs(abs(dotProduct(halfSpaceSet[candidate_plane_id].ABCD, halfSpaceSet[partner_plane_id].ABCD)) - 1.0) < EPSILON)
		return false;	// parallel planes

	// find the intersection parameter (distance) from the given point to the target point (which is on the intersection of the three planes)
	double t = findLinePlaneIntersection(line, halfSpaceSet[candidate_plane_id]);

	// the intersection of the three planes may be an empty set, check this
	if (t == numeric_limits<double>::infinity())
		return false;	// non-intersecting planes

	// compute the target (intersection) point by using the found parameter
	for (int i = 0; i < 3; i++)
		candidate_point[i] = line.point[i] + line.directionVector[i] * t;

	// check if this point inside the kernel or not
	double* scalar = findValueVectorSatisfiedByPoint(candidate_point, halfSpaceSet);
	for (int i = 0; i < halfSpaceSet.size(); i++)
		if (scalar[i] > EPSILON) {
			delete[] scalar;
			return false;
		}

	delete[] scalar;
	return true;

}

bool KerGen_Optimal::ClosestPlaneToPoint_0D(int candidate_plane_id, int base_plane_id, int partner_plane_id, int back_plane_id, Line line, double* candidate_point, double& smallest_t) {

	// check if the candidate plane is parallel to the base plane or not
	if (abs(abs(dotProduct(halfSpaceSet[candidate_plane_id].ABCD, halfSpaceSet[base_plane_id].ABCD)) - 1.0) < BIG_EPSILON)
		return false;	// parallel planes

	// check if the candidate plane is parallel to the partner plane or not
	if (abs(abs(dotProduct(halfSpaceSet[candidate_plane_id].ABCD, halfSpaceSet[partner_plane_id].ABCD)) - 1.0) < BIG_EPSILON)
		return false;	// parallel planes

	// find the intersection parameter (distance) from the given point to the target point (which is on the intersection of the three planes)
	double t = findLinePlaneIntersection(line, halfSpaceSet[candidate_plane_id]);

	// the intersection of the three planes may be an empty set, check this
	if (t == numeric_limits<double>::infinity())
		return false;	// non-intersecting planes

	// compute the target (intersection) point by using the found parameter
	for (int i = 0; i < 3; i++)
		candidate_point[i] = line.point[i] + line.directionVector[i] * t;

	// is this point behind the back plane? we do not want so
	double scalar = halfSpaceSet[back_plane_id].ABCD[3];
	for (int i = 0; i < 3; i++)
		scalar += halfSpaceSet[back_plane_id].ABCD[i] * candidate_point[i];
	if (scalar > -EPSILON)
		return false;
	
	// check if this point inside the kernel or not
	double* scalars = findValueVectorSatisfiedByPoint(candidate_point, halfSpaceSet);
	for (int i = 0; i < halfSpaceSet.size(); i++)
		if (scalars[i] > EPSILON) {
			delete[] scalars;
			return false;
		}
	delete[] scalars;
	

	// update t value for the normalized line direction vector
	double line_length = computeLength(line.directionVector);
	t *= line_length;
	// if this distance is larger than the previous ones, directly eliminate this candidate plane
	if (abs(t) > smallest_t && abs(abs(t) - smallest_t) > BIG_EPSILON)
		return false;

	
		// we want the other end point which is not locating at the back plane position, check this
		if (isTripleSame(line.point, candidate_point)) {
			// we are at the back plane position since we are inside the if clause
			// still do not give up, the position may be constituting a V-form combination of planes (explained in the paper), check this

			// check if the candidate plane is the same of the back plane or not
			if (isTripleSame(halfSpaceSet[candidate_plane_id].ABCD, halfSpaceSet[back_plane_id].ABCD))
				return false;	// the same planes

			// check if this is a contradictory plane (behind the back plane), or suitable to constitute v-type neighborhood
			vector<int> left_planes = { partner_plane_id, candidate_plane_id };
			vector<int> ordered_planes = orderFromFrontToBack(candidate_point, base_plane_id, back_plane_id, left_planes);
			if (ordered_planes[0] == partner_plane_id)	// if there is v-type neighborhood, candidate plane should in front of the partner plane
				return false;

		}
	
	smallest_t = abs(t);	// update the smallest distance
	return true;

}


void KerGen_Optimal::findInitialPoint_1(double* point) {

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




