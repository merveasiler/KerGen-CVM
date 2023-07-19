// @author Merve Asiler

#include "KerGen_Modified.h"
#include "BaseGeoOpUtils.h"
#include "sdlp.h"
#include "CGALUtils.h"

#include <cinolib/predicates.h>

using namespace cinolib;
/*
	if (fabs(d) < TOLL)
		return INTERSECT;
	else if (d > 0)
		return ABOVE;
	else
		return BELOW;
*/

KerGen_Modified::KerGen_Modified(const Mesh& hostMesh) :
	KernelExpansion(hostMesh, true) {

}

KerGen_Modified::~KerGen_Modified() {

}

void KerGen_Modified::expandKernel() {

	//cout << "Number of bounding planes: " << halfSpaceSet.size() << " out of " << hostMeshptr->getNumOfTris() << " triangles." << endl;

	double point[3];
	if (this->initialPoint == NULL)
		return;		// NOT STAR-SHAPED!
	for (int i = 0; i < 3; i++)
		point[i] = this->initialPoint[i];

	//clock_t begin = clock();
	initialize();
	findInitialVertexAndLines(point);
	//clock_t end = clock();
	//cout << "TIME: " << double(end - begin) * 1000 / CLOCKS_PER_SEC << endl;

	int number_of_kernel_faces = 0;
	int number_of_kernel_edges = 0;

	for (bool any_white_left = false; ; any_white_left = false) {
		for (int i = 0; i < halfSpaceSet.size(); i++) {
			if (isKernelFace[i] == ProcessColor_Modified::WHITE) {
				while (!edgePartners[i].empty()) {
					// initializations
					EdgePartnerTriple_Modified& ept = edgePartners[i].front();
					int partnerPlaneId = ept.partnerPlaneId;
					int backPlaneId = ept.backPlaneId;
					double edgeDirection[3] = { ept.edgeDirection[0], ept.edgeDirection[1], ept.edgeDirection[2] };
					double startPoint[3] = { ept.startPoint[0], ept.startPoint[1], ept.startPoint[2] };
					int startPointId = ept.startPointId;
					edgePartners[i].pop();

					// find the corner on the given edge
					int previousNumOfVerts = kernel.getNumOfVerts();
					double kernelVertex[3];
					vector<int> theClosestIds = findTheClosestHalfSpace(startPointId, edgeDirection, i, partnerPlaneId, backPlaneId, kernelVertex);
					if (theClosestIds.size() == 0) {
						cout << "invalid edge\n";
						continue;
					}

					if (kernel.getNumOfVerts() > previousNumOfVerts /* || isTripleSame(kernelVertex, startPoint)*/) {
						vector<int> orderedClosestIds = orderFromFrontToBack(i, partnerPlaneId, theClosestIds, kernelVertex);
						identifyNextEdges(i, partnerPlaneId, orderedClosestIds, kernelVertex);
						//orderTheFaces(i, partnerId, theClosestIds, kernelVertex, edgeDirection);
					}

					number_of_kernel_edges++;
				}
				any_white_left = true;
				isKernelFace[i] = ProcessColor_Modified::GREEN;
				number_of_kernel_faces++;
				break;
			}
		}
		if (!any_white_left)
			break;
	}

	//	cout << "[number of kernel edges: " << number_of_kernel_edges << "], number of kernel faces: " << number_of_kernel_faces << "]" << endl;

	// checkKernelForNonKernelVertices();

	if (kernel.getNumOfVerts() > 0)
		kernel = computeConvexHull(kernel.getAllVerts());

}

void KerGen_Modified::initialize() {

	for (int i = 0; i < halfSpaceSet.size(); i++) {
		edgePartners.push_back(queue<EdgePartnerTriple_Modified>());
		isKernelFace.push_back(ProcessColor_Modified::RED);
		edgePartnerIds.push_back(vector<int>());
	}

}

bool KerGen_Modified::isLineDirectionValid(double* startPoint, double* lineDirection) {

	double testPoint[3];
	for (int k = 0; k < 3; k++)
		testPoint[k] = startPoint[k] + lineDirection[k] * BIG_EPSILON;

	double* scalars = findValueVectorSatisfiedByPoint(testPoint, halfSpaceSet);

	bool isValid = true;
	for (int k = 0; k < halfSpaceSet.size(); k++)
		if (scalars[k] > BIG_EPSILON) {
			isValid = false;
			break;
		}
	delete[] scalars;

	return isValid;
}

bool KerGen_Modified::isLineValid(double* startPoint, double* lineDirection) {

	double testPoint[3];
	for (int k = 0; k < 3; k++)
		testPoint[k] = startPoint[k] + lineDirection[k] * BIG_EPSILON;

	double* scalars = findValueVectorSatisfiedByPoint(testPoint, halfSpaceSet);

	bool isValid = true;
	for (int k = 0; k < halfSpaceSet.size(); k++)
		if (scalars[k] > 3 * BIG_EPSILON) {
			isValid = false;
			break;
		}
	delete[] scalars;

	if (!isValid) {	// revert the line direction
		for (int k = 0; k < 3; k++) {
			lineDirection[k] *= -1;
			testPoint[k] = startPoint[k] + lineDirection[k] * BIG_EPSILON;
		}

		double* scalars = findValueVectorSatisfiedByPoint(testPoint, halfSpaceSet);

		isValid = true;
		for (int k = 0; k < halfSpaceSet.size(); k++)
			if (scalars[k] > 3 * BIG_EPSILON) {
				isValid = false;
				break;
			}
		delete[] scalars;
	}

	return isValid;
}

bool KerGen_Modified::correctLineDirection(double* startPoint, double* lineDirection, int referencePlaneId) {

	double testPoint[3];
	for (int k = 0; k < 3; k++)
		testPoint[k] = startPoint[k] + lineDirection[k];
	
	HalfSpace referencePlane = halfSpaceSet[referencePlaneId];
	double d = orient3d(referencePlane.point1, referencePlane.point2, referencePlane.point3, testPoint);

	if (!(fabs(d) < BIG_EPSILON) && d < 0) {
		for (int k = 0; k < 3; k++)
			lineDirection[k] *= -1;
		return true;	// reverted
	}

	return false;	// not reverted

}

bool KerGen_Modified::shiftOnAnotherPlane(double* point, double* rootpoint, int& shiftId) {

	bool shouldShift = false;
	double* scalars = findValueVectorSatisfiedByPoint(point, halfSpaceSet);
	double theInnerMostDistance = BIG_EPSILON;
	int theInnerMostId = -1;
	for (int i = 0; i < halfSpaceSet.size(); i++)
		if (scalars[i] > theInnerMostDistance) {

			shouldShift = true;
			if (isPointOnPlane(rootpoint, i)) {
				theInnerMostDistance = scalars[i];
				theInnerMostId = i;
			}

		}
	delete[] scalars;

	shiftId = theInnerMostId;
	return shouldShift;

}

bool KerGen_Modified::shiftOnAnotherPlane(double* point, double* rootpoint, int& shiftId, int baseId) {

	bool shouldShift = false;
	double* scalars = findValueVectorSatisfiedByPoint(point, halfSpaceSet);
	double theInnerMostDistance = -1;
	int theInnerMostId = -1;
	for (int i = 0; i < halfSpaceSet.size(); i++) {
		if (scalars[i] > BIG_EPSILON) {
			shouldShift = true;
			if (isPointOnPlane(rootpoint, i)) {
				double lineDir[3];
				crossProduct(halfSpaceSet[baseId].ABCD, halfSpaceSet[i].ABCD, lineDir);
				if (isZeroVector(lineDir))
					continue;
				normalize(lineDir);
				double rayDir[3];
				crossProduct(halfSpaceSet[baseId].ABCD, lineDir, rayDir);
				normalize(rayDir);

				Line ray(point, rayDir);
				double t = findLinePlaneIntersection(ray, halfSpaceSet[i]);

				if (abs(t) > theInnerMostDistance) {
					theInnerMostDistance = abs(t);
					theInnerMostId = i;
				}
			}
		}
	}
	delete[] scalars;

	shiftId = theInnerMostId;
	return shouldShift;

}

bool KerGen_Modified::isPointOnPlane(double* point, int planeId) {

	double distance = halfSpaceSet[planeId].ABCD[3];
	for (int d = 0; d < 3; d++)
		distance += halfSpaceSet[planeId].ABCD[d] * point[d];

	if (abs(distance) < BIG_EPSILON)
		return true;
	return false;

}

bool KerGen_Modified::isPointInFrontOfThePlane(double* point, int planeId) {

	double distance = halfSpaceSet[planeId].ABCD[3];
	for (int d = 0; d < 3; d++)
		distance += halfSpaceSet[planeId].ABCD[d] * point[d];

	if (distance < BIG_EPSILON)
		return true;
	return false;

}

bool KerGen_Modified::areTheSamePlanes(int plane1Id, int plane2Id) {

	if (isTripleSame(halfSpaceSet[plane1Id].ABCD, halfSpaceSet[plane2Id].ABCD)) {
		double* point = halfSpaceSet[plane1Id].point1;
		if (isPointOnPlane(point, plane2Id))
			return true;
	}

	return false;

}

bool KerGen_Modified::isPreviousLine(int hs1_id, int hs2_id) {

	bool walkedEdge = false;
	for (int j = 0; j < edgePartnerIds[hs1_id].size(); j++)
		if (edgePartnerIds[hs1_id][j] == hs2_id) {
			walkedEdge = true;
			break;
		}
	return walkedEdge;
}

void KerGen_Modified::recordLineInfo(int hs1_id, int hs2_id, int hs3_id, double* lineDir, double* point) {

	edgePartners[hs1_id].push(EdgePartnerTriple_Modified(hs2_id, hs3_id, lineDir, point, kernel.getNumOfVerts() - 1));
	edgePartnerIds[hs1_id].push_back(hs2_id);
	edgePartnerIds[hs2_id].push_back(hs1_id);
	if (isKernelFace[hs1_id] != ProcessColor_Modified::GREEN)
		isKernelFace[hs1_id] = ProcessColor_Modified::WHITE;
	if (isKernelFace[hs2_id] != ProcessColor_Modified::GREEN)
		isKernelFace[hs2_id] = ProcessColor_Modified::WHITE;

}

int KerGen_Modified::recordVertexInfo(double* vertex, vector<int> thirdParentIds, vector<int> lineParentIds) {

	int vertexId = -1;
	for (int v = 0; v < kernel.getNumOfVerts(); v++)
		if (isTripleSame(kernel.getVertex(v).coords, vertex)) {
			vertexId = v;
			break;
		}

	if (vertexId == -1) {
		vertexId = kernel.getNumOfVerts();
		kernel.addVertex(vertex[0], vertex[1], vertex[2]);
		vector<int> parentIdsForNewVertex;
		vertexParentIds.push_back(parentIdsForNewVertex);

		for (int i = 0; i < thirdParentIds.size(); i++)
			vertexParentIds[vertexId].push_back(thirdParentIds[i]);
		for (int i = 0; i < lineParentIds.size(); i++)
			vertexParentIds[vertexId].push_back(lineParentIds[i]);
	}

	return vertexId;
}




void KerGen_Modified::findInitialVertexAndLines(double* point) {

	double distanceToFirst;
	vector<int> firstClosestIds = findFirstPlane(point, distanceToFirst);

	for (int i = 0; i < firstClosestIds.size(); i++) {
		int plane1Id = firstClosestIds[i];

		for (int d = 0; d < 3; d++)
			point[d] = point[d] + halfSpaceSet[plane1Id].ABCD[i] * distanceToFirst;

		vector<double> pointerToLineDirections, intersectionLineDirections;
		double distanceToSecond;
		vector<int> secondClosestIds = findSecondPlane(point, plane1Id, pointerToLineDirections, intersectionLineDirections, distanceToSecond);

		for (int j = 0; j < secondClosestIds.size(); j++) {
			int plane2Id = secondClosestIds[j];

			double pointerToLineDir[3], firstLineDir[3];
			for (int d = 0; d < 3; d++) {
				pointerToLineDir[d] = pointerToLineDirections[j * 3 + d];
				firstLineDir[d] = intersectionLineDirections[j * 3 + d];
			}

			for (int d = 0; d < 3; d++)
				point[d] = point[d] + pointerToLineDir[d] * distanceToSecond;

			double vertex[3];
			vector<int> thirdClosestIds = findThirdPlane(point, firstLineDir, plane1Id, plane2Id, vertex);

			if (identifyInitialPlanes(vertex))
				return;


		}
	}

}

bool KerGen_Modified::identifyInitialPlanes(double* vertex) {

	vector<int> planeIds = vertexParentIds[0];
	cout << planeIds.size() << endl;

	for (int i = 0; i < planeIds.size(); i++) {
		int planeId1 = planeIds[i];
		bool plane1_valid = false;

		for (int j = i; j < planeIds.size(); j++) {
			int planeId2 = planeIds[j];
			bool plane2_valid = false;

			if (isTripleSame(halfSpaceSet[planeId1].ABCD, halfSpaceSet[planeId2].ABCD))
				continue;

			double lineDir1[3];
			crossProduct(halfSpaceSet[planeId1].ABCD, halfSpaceSet[planeId2].ABCD, lineDir1);
			normalize(lineDir1);

			if (!plane1_valid) {
				if (isLineValid(vertex, lineDir1)) {
					double perpDir1[3];
					crossProduct(halfSpaceSet[planeId1].ABCD, lineDir1, perpDir1);
					normalize(perpDir1);
					double point[3];
					for (int d = 0; d < 3; d++)
						point[d] = vertex[d] + BIG_EPSILON * lineDir1[d];
					if (isLineValid(point, perpDir1))
						plane1_valid = true;
					else
						continue;
				}
				else
					continue;
			}

			double perpDir2[3];
			crossProduct(halfSpaceSet[planeId2].ABCD, lineDir1, perpDir2);
			normalize(perpDir2);
			double point[3];
			for (int d = 0; d < 3; d++)
				point[d] = vertex[d] + BIG_EPSILON * lineDir1[d];
			if (isLineValid(point, perpDir2))
				plane2_valid = true;
			else
				continue;

			for (int k = j; k < planeIds.size(); k++) {
				int planeId3 = planeIds[k];
				bool plane3_valid = false;

				if (isTripleSame(halfSpaceSet[planeId1].ABCD, halfSpaceSet[planeId3].ABCD))
					continue;
				if (isTripleSame(halfSpaceSet[planeId2].ABCD, halfSpaceSet[planeId3].ABCD))
					continue;

				double lineDir2[3];
				crossProduct(halfSpaceSet[planeId1].ABCD, halfSpaceSet[planeId3].ABCD, lineDir2);
				normalize(lineDir2);

				if (isLineValid(vertex, lineDir2)) {
					double perpDir3[3];
					crossProduct(halfSpaceSet[planeId1].ABCD, lineDir2, perpDir3);
					normalize(perpDir3);
					double point[3];
					for (int d = 0; d < 3; d++)
						point[d] = vertex[d] + BIG_EPSILON * lineDir2[d];
					if (isLineValid(point, perpDir3))
						plane3_valid = true;
				}

				double lineDir3[3];
				crossProduct(halfSpaceSet[planeId2].ABCD, halfSpaceSet[planeId3].ABCD, lineDir3);
				normalize(lineDir3);

				if (isLineValid(vertex, lineDir3)) {
					double perpDir3[3];
					crossProduct(halfSpaceSet[planeId2].ABCD, lineDir3, perpDir3);
					normalize(perpDir3);
					double point[3];
					for (int d = 0; d < 3; d++)
						point[d] = vertex[d] + BIG_EPSILON * lineDir3[d];
					if (isLineValid(point, perpDir3))
						plane3_valid = true;
					else
						continue;
				}
				else
					continue;

				// if I come here, I successfully found 3 planes
				correctLineDirection(vertex, lineDir1, planeId3);
				correctLineDirection(vertex, lineDir2, planeId2);
				correctLineDirection(vertex, lineDir3, planeId1);

				recordLineInfo(planeId1, planeId2, planeId3, lineDir1, vertex);
				recordLineInfo(planeId1, planeId3, planeId2, lineDir2, vertex);
				recordLineInfo(planeId2, planeId3, planeId1, lineDir3, vertex);

				return true;

			}

		}

	}

	return false;
}

vector<int> KerGen_Modified::findFirstPlane(double* point, double& distanceToFirst) {

	double* scalars = findValueVectorSatisfiedByPoint(point, halfSpaceSet);

	distanceToFirst = numeric_limits<double>::infinity();
	vector<int> theClosestIds;
	for (int i = 0; i < halfSpaceSet.size(); i++) {
		if (abs(scalars[i]) <= distanceToFirst + BIG_EPSILON) {
			if (abs(scalars[i]) < distanceToFirst - BIG_EPSILON)
				theClosestIds.clear();
			theClosestIds.push_back(i);
			distanceToFirst = abs(scalars[i]);
		}
	}
	delete[] scalars;

	return theClosestIds;

}

vector<int> KerGen_Modified::findSecondPlane(double* point, int id, vector<double>& pointerToLineDirections, vector<double>& intersectionLineDirections, double& distanceToSecond) {

	distanceToSecond = numeric_limits<double>::infinity();
	int theClosestId = -1;
	vector<int> theClosestIds;
	for (int s = 0, i = 0, bound = id; s < 2; s++) {
		for (; i < bound; i++) {
			if (isTripleSame(halfSpaceSet[id].ABCD, halfSpaceSet[i].ABCD))
				continue;
			double lineDirection[3], rayDirection[3];
			crossProduct(halfSpaceSet[i].ABCD, halfSpaceSet[id].ABCD, lineDirection);
			normalize(lineDirection);
			crossProduct(halfSpaceSet[id].ABCD, lineDirection, rayDirection);
			normalize(rayDirection);
			Line ray(point, rayDirection);
			double t = findLinePlaneIntersection(ray, halfSpaceSet[i]);
			if (abs(t) <= distanceToSecond + BIG_EPSILON) {
				if (abs(t) < distanceToSecond - BIG_EPSILON) {
					theClosestIds.clear();
					pointerToLineDirections.clear();
					intersectionLineDirections.clear();
				}
				theClosestIds.push_back(i);
				distanceToSecond = abs(t);
				for (int k = 0; k < 3; k++) {
					pointerToLineDirections.push_back(rayDirection[k]);
					intersectionLineDirections.push_back(lineDirection[k]);
				}
			}
		}
		i = id + 1;
		bound = halfSpaceSet.size();
	}

	return theClosestIds;

}

vector<int> KerGen_Modified::findThirdPlane(double* startpoint, double* lineDirection, int lineParent1Id, int lineParent2Id, double* newpoint) {

	double theClosestDistance = numeric_limits<double>::infinity();
	vector<int> theClosestIds;
	vector<int> theSameOfLineParents;

	Line line(startpoint, lineDirection);

	for (int i = 0; i < halfSpaceSet.size(); i++) {

		if (areTheSamePlanes(lineParent1Id, i)) {
			theSameOfLineParents.push_back(i);
			continue;
		}
		if (areTheSamePlanes(lineParent2Id, i)) {
			theSameOfLineParents.push_back(i);
			continue;
		}

		double t = findLinePlaneIntersection(line, halfSpaceSet[i]);

		if (t >= -BIG_EPSILON && t <= theClosestDistance + BIG_EPSILON) {
			if (t < theClosestDistance - BIG_EPSILON)
				theClosestIds.clear();
			theClosestIds.push_back(i);
			theClosestDistance = t;
		}

	}

	for (int k = 0; k < 3; k++)
		newpoint[k] = startpoint[k] + lineDirection[k] * theClosestDistance;

	recordVertexInfo(newpoint, theClosestIds, theSameOfLineParents);

	return theClosestIds;

}

vector<int> KerGen_Modified::findTheClosestHalfSpace(int startVertexId, double* lineDirection, int lineParent1Id, int lineParent2Id, int backPlaneId, double* vertex) {

	double theClosestDistance = numeric_limits<double>::infinity();
	vector<int> theClosestIds;
	double* startpoint = kernel.getVertex(startVertexId).coords;
	vector<int> theSameOfLineParents;
	double closestPoint[3];
	bool closestPlaneInitialized = false;

	//correctLineDirection(startpoint, lineDirection, backPlaneId);
	Line line(startpoint, lineDirection);

	for (int i = 0; i < halfSpaceSet.size(); i++) {

		if (areTheSamePlanes(lineParent1Id, i)) {
			theSameOfLineParents.push_back(i);
			continue;
		}
		if (areTheSamePlanes(lineParent2Id, i)) {
			theSameOfLineParents.push_back(i);
			continue;
		}

		double t = findLinePlaneIntersection(line, halfSpaceSet[i]);

		/*
		for (int d = 0; d < 3; d++)
			vertex[d] = startpoint[d] + lineDirection[d] * t;
		
		HalfSpace backPlane = halfSpaceSet[backPlaneId];
		double d = orient3d(backPlane.point1, backPlane.point2, backPlane.point3, vertex);

		if (fabs(d) < BIG_EPSILON || d < 0)
			continue;
		*/
		if (t >= -EPSILON && t <= theClosestDistance) {
			bool previous_parent = false;
			if (t < EPSILON) {
				int n = vertexParentIds[startVertexId].size();
				for (int j = 0; j < n; j++)
					if (i == vertexParentIds[startVertexId][j] || isTripleSame(halfSpaceSet[i].ABCD, halfSpaceSet[vertexParentIds[startVertexId][j]].ABCD)) {
						//vertexParentIds[vertexId].push_back(i);
						previous_parent = true;
						break;
					}
			}
			if (previous_parent)
				continue;

			if (t < theClosestDistance)
				theClosestIds.clear();
			theClosestIds.push_back(i);
			theClosestDistance = t;
			//for (int d = 0; d < 3; d++)
			//	closestPoint[d] = vertex[d];
		}
		
		/*
		if (closestPlaneInitialized) {
			HalfSpace referencePlane = halfSpaceSet[theClosestIds[0]];
			double d = orient3d(referencePlane.point1, referencePlane.point2, referencePlane.point3, vertex);
			if (d >= 0 && d <= BIG_EPSILON)
				theClosestIds.push_back(i);
			else if (d > BIG_EPSILON) {
				theClosestIds.clear();
				theClosestIds.push_back(i);
				for (int d = 0; d < 3; d++)
					closestPoint[d] = vertex[d];
			}
		}
		else {
			for (int d = 0; d < 3; d++)
				closestPoint[d] = vertex[d];
			closestPlaneInitialized = true;
			theClosestIds.push_back(i);
		}
		*/
	}

	//for (int d = 0; d < 3; d++)
	//	vertex[d] = closestPoint[d];

	for (int d = 0; d < 3; d++)
		vertex[d] = startpoint[d] + lineDirection[d] * theClosestDistance;

	int vertexId = recordVertexInfo(vertex, theClosestIds, theSameOfLineParents);
/*
	double* sca = findValueVectorSatisfiedByPoint(vertex, halfSpaceSet);
	for (int s = 0; s < halfSpaceSet.size(); s++)
		if (sca[s] > BIG_EPSILON)
			cout << theClosestDistance << "   NOOOOOOOOOOOOOOO\n";
	delete[] sca;
	cout << vertexId << endl;
*/
	kernel.addEdge(startVertexId, vertexId);

	return theClosestIds;

}

vector<int> KerGen_Modified::orderFromFrontToBack(int base_plane_id, int partner_plane_id, vector<int> next_partner_ids, double* startpoint) {

	if (next_partner_ids.size() == 1)
		return next_partner_ids;
	cout << next_partner_ids.size() << endl;
	// reserve an empty array 
	vector<int> planesFromFrontToBack;
	for (int i = 0; i < next_partner_ids.size(); i++)
		planesFromFrontToBack.push_back(-1);

	// first prepare the common lines with the base plane, we will use these lines
	vector<Line> lines;	// save those lines into this vector
	for (int i = 0; i < next_partner_ids.size(); i++) {
		int candidate_plane_id = next_partner_ids[i];

		double candidate_line_direction[3];
		crossProduct(halfSpaceSet[base_plane_id].ABCD, halfSpaceSet[candidate_plane_id].ABCD, candidate_line_direction);
		normalize(candidate_line_direction);
		correctLineDirection(startpoint, candidate_line_direction, partner_plane_id);
		Line candidate_line(startpoint, candidate_line_direction);
		lines.push_back(candidate_line);
	}

	// we will use a test point on the common line of each candidate-base double
	// this test point will be just a bit ahead of the intersection point of candidate-base-partner triple
	for (int i = 0; i < next_partner_ids.size(); i++) {
		int candidate_plane_id = next_partner_ids[i];
		int rank_id = 0;

		// initialize the test point, note that this point may be wrong if the direction vector was taken reversed
		double test_point[3];
		for (int j = 0; j < 3; j++)
			test_point[j] = startpoint[j] + lines[i].directionVector[j];

		// now check whether the test point falls behid the other candidate planes or not
		// if falls behind, then it means the current candidate plane is at the backside of that other candidate plane
		//	in that case, increase its rank, it is at least one plane back (en az 1 plane arka tarafta)
		for (int j = 0; j < next_partner_ids.size(); j++) {
			int reference_plane_id = next_partner_ids[j];
			HalfSpace reference_plane = halfSpaceSet[reference_plane_id];

			double d = orient3d(reference_plane.point1, reference_plane.point2, reference_plane.point3, test_point);
			if (!(fabs(d) < BIG_EPSILON) && d < 0)
				rank_id++;
		}

		planesFromFrontToBack[rank_id] = candidate_plane_id;
	}

	// if some of the planes are the same, we take only one of them, which means there occurs some empty seats towards the last seats, trim them
	vector<int> trimmed;
	for (int i = 0; i < next_partner_ids.size(); i++) {
		if (planesFromFrontToBack[i] == -1)
			break;
		trimmed.push_back(planesFromFrontToBack[i]);
	}
	return trimmed;

}

void KerGen_Modified::identifyNextEdges(int base_plane_id, int partner_plane_id, vector<int> ordered_next_partner_ids, double* startpoint) {

	int next_partner_id = ordered_next_partner_ids[0];
	double nextLineDir[3];
	crossProduct(halfSpaceSet[base_plane_id].ABCD, halfSpaceSet[next_partner_id].ABCD, nextLineDir);
	bool shouldRevert = correctLineDirection(startpoint, nextLineDir, partner_plane_id);
	recordLineInfo(base_plane_id, next_partner_id, partner_plane_id, nextLineDir, startpoint);

	vector<double> foundLineDirs;
	ordered_next_partner_ids.push_back(partner_plane_id);
	
	for (int i = 0, j = 1; j < ordered_next_partner_ids.size(); j++) {
		int plane1Id = ordered_next_partner_ids[i];
		int plane2Id = ordered_next_partner_ids[j];

		if (!isPreviousLine(plane1Id, plane2Id)) {
			double lineDir[3];
			//crossProduct(halfSpaceSet[plane1Id].ABCD, halfSpaceSet[plane2Id].ABCD, lineDir);
			//normalize(lineDir);
			//correctLineDirection(startpoint, lineDir, base_plane_id);
			
			if (shouldRevert)
				crossProduct(halfSpaceSet[plane2Id].ABCD, halfSpaceSet[plane1Id].ABCD, lineDir);
			else
				crossProduct(halfSpaceSet[plane1Id].ABCD, halfSpaceSet[plane2Id].ABCD, lineDir);
			if (isZeroVector(lineDir))
				continue;
			normalize(lineDir);
			

			bool already_found = false;
			for (int f = 0; f < foundLineDirs.size(); f+=3) {
				double foundLineDir[3] = { foundLineDirs[f], foundLineDirs[f + 1], foundLineDirs[f + 2] };
				if (isTripleSame(lineDir, foundLineDir) || isTripleSame(foundLineDir, lineDir)) {
					already_found = true;
					break;
				}
			}
			if (!already_found) {
				for (int d = 0; d < 3; d++)
					foundLineDirs.push_back(lineDir[d]);
				recordLineInfo(plane1Id, plane2Id, base_plane_id, lineDir, startpoint);
				i = j;
			}
			else
				;
		}
		else
			i = j;
	}

	ordered_next_partner_ids.clear();
	foundLineDirs.clear();

}


