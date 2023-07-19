// @author Merve Asiler

#include "KerGen_Robust.h"
#include "BaseGeoOpUtils.h"
#include "sdlp.h"
#include "RobustOperations.h"
#include "CGALUtils.h"

KerGen_Robust::KerGen_Robust(const Mesh& hostMesh) :
	KernelExpansion(hostMesh, true) {
}

KerGen_Robust ::~KerGen_Robust() {

}

C::Point_3 CPoint(double* point) {

	return C::Point_3(point[0], point[1], point[2]);
}

C::Plane_3 CPlane(HalfSpace halfspace) {

	return C::Plane_3(CPoint(halfspace.point1), CPoint(halfspace.point2), CPoint(halfspace.point3));
}

void KerGen_Robust::expandKernel() {

	//cout << "Number of bounding planes: " << halfSpaceSet.size() << " out of " << hostMeshptr->getNumOfTris() << " triangles." << endl;

	initialize();
	
	double point[3];
	if (this->initialPoint == NULL)
		return;		// NOT STAR-SHAPED!
	for (int i = 0; i < 3; i++)
		point[i] = this->initialPoint[i];

	clock_t begin = clock();
	int theClosestId11 = findClosestPlaneToPointin3D(halfSpaceSet, point);
	clock_t end = clock();
	cout << "PREDICATE TIME-1-V1: " << double(end - begin) * 1000 / CLOCKS_PER_SEC << "   " << theClosestId11 << endl;

	clock_t begin2 = clock();
	int theClosestId21 = findClosestPlaneToPointin2D(halfSpaceSet, theClosestId11, point);
	clock_t end2 = clock();
	cout << "PREDICATE TIME-2-V1: " << double(end2 - begin2) * 1000 / CLOCKS_PER_SEC  << "   " << theClosestId21 << endl;

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

	while (!tobeProcessedFaceIds.empty()){
		int face_id = tobeProcessedFaceIds.front();
		//cout << face_id << endl;
		for (int partner = 0; ; partner++) {
			// initializations
			VertexPartners_Robust& vp = vertexPartnerships[face_id][partner];
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
			vertexPartnerships[face_id].push_back(VertexPartners_Robust(next_partner_id, partner1_id, found_vertex_id));

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
		
		faceProcess[face_id] = ProcessColor_Robust::GREEN;
		tobeProcessedFaceIds.pop();
		number_of_kernel_faces++;
			
	}
	

	//	cout << "[number of kernel edges: " << number_of_kernel_edges << "], number of kernel faces: " << number_of_kernel_faces << "]" << endl;

	if (kernel.getNumOfVerts() > 0)
		kernel = computeConvexHull(kernel.getAllVerts());

}

void KerGen_Robust::initialize() {
	
	for (int i = 0; i < halfSpaceSet.size(); i++) {
		vertexPartnerships.push_back(vector<VertexPartners_Robust>());
		faceProcess.push_back(ProcessColor_Robust::RED);
	}

}

int KerGen_Robust::findClosestHalfSpace(double* point) {
	
	// initialize the closest plane temporarily
	int closest_plane_id = 0;

	// update the closest plane if necessary
	for (int i = 1; i < halfSpaceSet.size(); i++) {
		CP_3D cp3d;
		if (cp3d(	CPoint(halfSpaceSet[closest_plane_id].point1), CPoint(halfSpaceSet[closest_plane_id].point2), CPoint(halfSpaceSet[closest_plane_id].point3),
					CPoint(halfSpaceSet[i].point1), CPoint(halfSpaceSet[i].point2), CPoint(halfSpaceSet[i].point3), 
					CPoint(point) ) == CGAL::LARGER)
			closest_plane_id = i;
	}

	return closest_plane_id;

}

int KerGen_Robust::findClosestHalfSpace(double* point, int base_plane_id) {

	// initialize the closest plane temporarily
	int closest_plane_id = 0;
	for (int i = 0; i < halfSpaceSet.size(); i++) {
		if (i == base_plane_id)
			continue;
		PP pp;
		if (pp(	CPoint(halfSpaceSet[base_plane_id].point1), CPoint(halfSpaceSet[base_plane_id].point2), CPoint(halfSpaceSet[base_plane_id].point3),
				CPoint(halfSpaceSet[i].point1), CPoint(halfSpaceSet[i].point2), CPoint(halfSpaceSet[i].point3)) )
			continue;
		closest_plane_id = i;
		break;
	}

	// update the closest plane if necessary
	for (int i = closest_plane_id + 1; i < halfSpaceSet.size(); i++) {
		if (i == base_plane_id)
			continue;
		PP pp;
		if (pp(	CPoint(halfSpaceSet[base_plane_id].point1), CPoint(halfSpaceSet[base_plane_id].point2), CPoint(halfSpaceSet[base_plane_id].point3),
				CPoint(halfSpaceSet[i].point1), CPoint(halfSpaceSet[i].point2), CPoint(halfSpaceSet[i].point3)) )
			continue;
		CP_2D cp2d;
		if (cp2d(	CPoint(halfSpaceSet[base_plane_id].point1), CPoint(halfSpaceSet[base_plane_id].point2), CPoint(halfSpaceSet[base_plane_id].point3),
					CPoint(halfSpaceSet[closest_plane_id].point1), CPoint(halfSpaceSet[closest_plane_id].point2), CPoint(halfSpaceSet[closest_plane_id].point3),
					CPoint(halfSpaceSet[i].point1), CPoint(halfSpaceSet[i].point2), CPoint(halfSpaceSet[i].point3),
					CPoint(point)) == CGAL::LARGER)
			closest_plane_id = i;
	}

	return closest_plane_id;

}

vector<int> KerGen_Robust::findClosestHalfSpace(double* point, int base_plane_id, int partner_plane_id) {

	vector<int> closest_plane_ids;
	// initialize the closest plane temporarily
	int closest_plane_id = 0;
	for (int i = 0; i < halfSpaceSet.size(); i++) {
		if (i == base_plane_id || i == partner_plane_id)
			continue;
		IP ip;
		if (ip(	CPoint(halfSpaceSet[base_plane_id].point1), CPoint(halfSpaceSet[base_plane_id].point2), CPoint(halfSpaceSet[base_plane_id].point3),
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
		CP_UOD cpuod;
		CGAL::Comparison_result cr = cpuod(
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

vector<int> KerGen_Robust::findClosestHalfSpace(int base_plane_id, int partner_plane_id, int back_plane_id, int vertex_id) {

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

		IFP ifp;
		if (ifp(	CPoint(halfSpaceSet[base_plane_id].point1), CPoint(halfSpaceSet[base_plane_id].point2), CPoint(halfSpaceSet[base_plane_id].point3),
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

		CP_1D cp1d;
		CGAL::Comparison_result cr = cp1d(
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

vector<int> KerGen_Robust::findFrontestPlane(int base_plane_id, int partner_plane_id, vector<int> latest_plane_ids) {
	
	vector<int> planesFromFrontToBack;
	for (int i = 0; i < latest_plane_ids.size(); i++)
		planesFromFrontToBack.push_back(-1);

	// update the frontest plane if necessary
	for (int i = 0; i < latest_plane_ids.size(); i++) {
		bool is_frontier = true;
		int rank_id = 0;
		for (int j = 0; j < latest_plane_ids.size(); j++) {
			if ( CGAL::parallel(CPlane(halfSpaceSet[latest_plane_ids[i]]), CPlane(halfSpaceSet[latest_plane_ids[j]])) )
				continue;
			int candidate_plane_id = latest_plane_ids[i];
			int reference_plane_id = latest_plane_ids[j];
			FP fp;
			if (fp(	CPoint(halfSpaceSet[base_plane_id].point1), CPoint(halfSpaceSet[base_plane_id].point2), CPoint(halfSpaceSet[base_plane_id].point3),
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

int KerGen_Robust::computeKernelVertex(int plane1_id, int plane2_id, vector<int> others) {

	auto auto_line = CGAL::intersection(CPlane(halfSpaceSet[plane1_id]), CPlane(halfSpaceSet[plane2_id]));
	C::Line_3* line = boost::get<C::Line_3>(&*auto_line);

	auto auto_point = CGAL::intersection(CPlane(halfSpaceSet[others[0]]), *line);
	C::Point_3* point = boost::get<C::Point_3>(&*auto_point);
	double found_point[3] = {point->x(), point->y(), point->z()};
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
		kernel.addVertex(point->x(), point->y(), point->z());
		vector<int> vertex_generator_ids;
		vertex_generator_ids.push_back(plane1_id);
		vertex_generator_ids.push_back(plane2_id);
		for (int i = 0; i < others.size(); i++)
			vertex_generator_ids.push_back(others[i]);
		generatorsOfVertices.push_back(vertex_generator_ids);
	}

	return vertex_id;
}

bool KerGen_Robust::addIntoProcessQueueAfterCheck(int plane1_id, int plane2_id, int plane3_id, int plane4_id, int vertex_id) {

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
		if (faceProcess[plane1_id] == ProcessColor_Robust::RED) {
			faceProcess[plane1_id] = ProcessColor_Robust::WHITE;
			tobeProcessedFaceIds.push(plane1_id);
			vertexPartnerships[plane1_id].push_back(VertexPartners_Robust(plane2_id, plane4_id, vertex_id));
		}
	}

	delete line1;
	delete line2;

	return is_valid;

}




void KerGen_Robust::findInitialPoint_1(double* point) {

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




