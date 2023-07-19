// @author Merve Asiler

#include "KerGen.h"
#include "BaseGeoOpUtils.h"
#include "sdlp.h"
#include "CGALUtils.h"

KerGen::KerGen(const Mesh& hostMesh) :
	KernelExpansion(hostMesh, true) {

}

KerGen::~KerGen() {

}

void KerGen::expandKernel() {

	//cout << "Number of bounding planes: " << halfSpaceSet.size() << " out of " << hostMeshptr->getNumOfTris() << " triangles." << endl;

	double point[3];
	if (this->initialPoint == NULL)
		return;		// NOT STAR-SHAPED!
	for (int i = 0; i < 3; i++)
		point[i] = this->initialPoint[i];

	//clock_t begin = clock();
	vector<double> scalarsVector = initialize(point);
	int theClosestId1 = findTheClosestHalfSpace(point, scalarsVector);
	//clock_t end = clock();
	//cout << "NON-PREDICATE TIME: " << double(end - begin) * 1000 / CLOCKS_PER_SEC << endl;

	//clock_t begin2 = clock();
	int theClosestId2 = findTheClosestHalfSpace(point, theClosestId1);
	//clock_t end2 = clock();
	//cout << "NON-PREDICATE TIME-2: " << double(end2 - begin2) * 1000 / CLOCKS_PER_SEC << endl;
	//cout << theClosestId2 << endl;
	int number_of_kernel_faces = 0;
	int number_of_kernel_edges = 0;
	
	//cout << theClosestId1 << " " << theClosestId2 << " " << point[0] << " " << point[1] << " " << point[2] << endl;
	
	for (bool any_white_left = false; ; any_white_left = false) {
		for (int i = 0; i < halfSpaceSet.size(); i++) {
			if (isKernelFace[i] == ProcessColor::WHITE) {
				while (!edgePartners[i].empty()) {
					// initializations
					EdgePartnerTriple& ept = edgePartners[i].front();
					int partnerId = ept.partnerId;
					double edgeDirection[3] = { ept.edgeDirection[0], ept.edgeDirection[1], ept.edgeDirection[2] };
					double startPoint[3] = { ept.startPoint[0], ept.startPoint[1], ept.startPoint[2] };
					int startPointId = ept.startPointId;
					edgePartners[i].pop();

					// find the corner on the given edge
					int previousNumOfVerts = kernel.getNumOfVerts();
					double kernelVertex[3];
					vector<int> theClosestIds = findTheClosestHalfSpace(startPointId, edgeDirection, i, partnerId, kernelVertex);
					if (theClosestIds.size() == 0) {
						cout << "invalid edge\n";
						return;
					}

					if (kernel.getNumOfVerts() > previousNumOfVerts || isTripleSame(kernelVertex, startPoint))
						orderTheFaces(i, partnerId, theClosestIds, kernelVertex, edgeDirection);

					number_of_kernel_edges++;
				}
				any_white_left = true;
				isKernelFace[i] = ProcessColor::GREEN;
				number_of_kernel_faces++;
				break;
			}
		}
		if (!any_white_left)
			break;
	}

	//cout << kernel.getNumOfVerts() << " " << number_of_kernel_faces << " " << number_of_kernel_edges << endl;

	//	cout << "[number of kernel edges: " << number_of_kernel_edges << "], number of kernel faces: " << number_of_kernel_faces << "]" << endl;

	// checkKernelForNonKernelVertices();

	if (kernel.getNumOfVerts() > 0)
		kernel = computeConvexHull(kernel.getAllVerts());

}

vector<double> KerGen::initialize(double* point) {

	double* scalarsArray = findValueVectorSatisfiedByPoint(point, halfSpaceSet);
	vector<double> scalarsVector;
	filterRepetitions2(scalarsArray, scalarsVector);

	for (int i = 0; i < halfSpaceSet.size(); i++) {
		edgePartners.push_back(queue<EdgePartnerTriple>());
		isKernelFace.push_back(ProcessColor::RED);
		edgePartnerIds.push_back(vector<int>());
	}

	return scalarsVector;
}

int KerGen::findTheClosestHalfSpace(double* point, vector<double>& scalarsVector) {

	double theClosestDistance = numeric_limits<double>::infinity();
	vector<int> theClosestIds;
	for (int i = 0; i < halfSpaceSet.size(); i++) {
		if (abs(scalarsVector[i]) <= theClosestDistance) {
			if (abs(scalarsVector[i]) < theClosestDistance)
				theClosestIds.clear();
			theClosestIds.push_back(i);
			theClosestDistance = abs(scalarsVector[i]);
		}
	}

	int theClosestId = theClosestIds[0];
	scalarsVector.clear();
	for (int i = 0; i < 3; i++)
		point[i] = point[i] + halfSpaceSet[theClosestId].ABCD[i] * theClosestDistance;
	return theClosestId;
}

int KerGen::findTheClosestHalfSpace(double* point, int id) {

	double theClosestDistance = numeric_limits<double>::infinity();
	int theClosestId = -1;
	double theClosestDirection[3], intersectionLineDirection[3];
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
			if (t < theClosestDistance) {
				theClosestId = i;
				theClosestDistance = t;
				for (int k = 0; k < 3; k++) {
					theClosestDirection[k] = rayDirection[k];
					intersectionLineDirection[k] = lineDirection[k];
				}
			}
		}
		i = id + 1;
		bound = halfSpaceSet.size();
	}

	for (int i = 0; i < 3; i++)
		point[i] = point[i] + theClosestDirection[i] * theClosestDistance;

	edgePartners[id].push(EdgePartnerTriple(theClosestId, intersectionLineDirection, point, 0));
	edgePartnerIds[id].push_back(theClosestId);
	edgePartnerIds[theClosestId].push_back(id);
	isKernelFace[id] = ProcessColor::WHITE;
	isKernelFace[theClosestId] = ProcessColor::WHITE;

	vector<int> parentIdsForNewVertex;
	parentIdsForNewVertex.push_back(id);
	parentIdsForNewVertex.push_back(theClosestId);
	vertexParentIds.push_back(parentIdsForNewVertex);
	kernel.addVertex(point[0], point[1], point[2]);

	return theClosestId;

}

vector<int> KerGen::findTheClosestHalfSpace(int vertexId, double* lineDirection, int lineParent1Id, int lineParent2Id, double* newpoint) {

	double theClosestDistance = numeric_limits<double>::infinity();
	vector<int> theClosestIds;
	Line line(kernel.getVertex(vertexId).coords, lineDirection);

	for (int i = 0; i < halfSpaceSet.size(); i++) {
		bool previous_parent = false;
		for (int j = 0; j < vertexParentIds[vertexId].size(); j++)
			if (i == vertexParentIds[vertexId][j]) {
				previous_parent = true;
				break;
			}
		if (previous_parent)
			continue;

		double t = findLinePlaneIntersection(line, halfSpaceSet[i]);
		if (t >= 0 && t <= theClosestDistance) {
			if (t < EPSILON) {
				int n = vertexParentIds[vertexId].size();
				for (int j = 0; j < n; j++)
					if (i == vertexParentIds[vertexId][j] || isTripleSame(halfSpaceSet[i].ABCD, halfSpaceSet[vertexParentIds[vertexId][j]].ABCD)) {
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
		}
	}

	for (int i = 0; i < 3; i++)
		newpoint[i] = kernel.getVertex(vertexId).coords[i] + lineDirection[i] * theClosestDistance;
	double t = findClosestValueSatisfiedByPoint(newpoint, halfSpaceSet);
	if (t > 3 * EPSILON)
		return theClosestIds;

	int startId = vertexId;
	if (theClosestDistance > 0) {
		bool previousVertex = false;
		for (int v = 0; v < kernel.getNumOfVerts(); v++)
			if (isTripleSame(kernel.getVertex(v).coords, newpoint)) {
				previousVertex = true;
				vertexId = v;
				break;
			}
		if (!previousVertex) {
			kernel.addVertex(newpoint[0], newpoint[1], newpoint[2]);
			vector<int> parentIdsForNewVertex;
			parentIdsForNewVertex.push_back(lineParent1Id);
			parentIdsForNewVertex.push_back(lineParent2Id);
			vertexParentIds.push_back(parentIdsForNewVertex);
			vertexId = kernel.getNumOfVerts() - 1;
		}
	}

	for (int i = 0; i < theClosestIds.size(); i++)
		vertexParentIds[vertexId].push_back(theClosestIds[i]);

	kernel.addEdge(startId, vertexId);

	return theClosestIds;

}

void KerGen::orderTheFaces(int base_id, int partner_id, vector<int> next_partner_ids, double* startPoint, double* currentEdgeDirection) {

	vector<int> ordered_next_partner_ids;
	vector<double> cosAngles;
	bool should_revert = false;

	double currentEdgeDirection2[3];
	crossProduct(halfSpaceSet[base_id].ABCD, halfSpaceSet[partner_id].ABCD, currentEdgeDirection2);
	if (dotProduct(currentEdgeDirection, currentEdgeDirection2) < 0) {
		multVect(currentEdgeDirection, -1.0, currentEdgeDirection);
		should_revert = true;
	}

	// ordering
	for (int i = 0; i < next_partner_ids.size(); i++) {

		double nextEdgeDirection[3];
		crossProduct(halfSpaceSet[base_id].ABCD, halfSpaceSet[next_partner_ids[i]].ABCD, nextEdgeDirection);
		normalize(nextEdgeDirection);
		double cosAngle = dotProduct(currentEdgeDirection, nextEdgeDirection);
		cosAngles.push_back(cosAngle);
		ordered_next_partner_ids.push_back(next_partner_ids[i]);

		for (int j = 0; j < cosAngles.size() - 1; j++) {
			if (abs(cosAngle - cosAngles[j]) < EPSILON) {
				if (isTripleSame(halfSpaceSet[next_partner_ids[i]].ABCD, halfSpaceSet[ordered_next_partner_ids[j]].ABCD)) {
					cosAngles.pop_back();
					ordered_next_partner_ids.pop_back();
					break;
				}
			}
			if (cosAngle < cosAngles[j]) {
				for (int k = cosAngles.size() - 1; k >= j; k--) {
					cosAngles[k] = cosAngles[k - 1];
					ordered_next_partner_ids[k] = ordered_next_partner_ids[k - 1];
				}
				cosAngles[j] = cosAngle;
				ordered_next_partner_ids[j] = next_partner_ids[i];
				break;
			}
		}
	}

	ordered_next_partner_ids.insert(ordered_next_partner_ids.begin(), base_id);
	ordered_next_partner_ids.push_back(partner_id);

	for (int i = 0, j = 1; j < ordered_next_partner_ids.size(); j++) {
		if (isWalkedEdge(ordered_next_partner_ids[i], ordered_next_partner_ids[j])) {
			double edgeDirection[3];
			if (i == 0)
				findEdgeDirection(ordered_next_partner_ids[i], ordered_next_partner_ids[j], should_revert, edgeDirection, NULL);
			else
				findEdgeDirection(ordered_next_partner_ids[i], ordered_next_partner_ids[j], should_revert, edgeDirection, halfSpaceSet[base_id].ABCD);
			if (isValidEdge(startPoint, edgeDirection))
				edgePartners[ordered_next_partner_ids[i]].push(EdgePartnerTriple(ordered_next_partner_ids[j], edgeDirection, startPoint, kernel.getNumOfVerts() - 1));
			else
				continue;
		}
		i = j;
	}

	ordered_next_partner_ids.clear();
	cosAngles.clear();
}

bool KerGen::isWalkedEdge(int hs1_id, int hs2_id) {

	bool walkedEdge = false;
	for (int j = 0; j < edgePartnerIds[hs1_id].size(); j++)
		if (edgePartnerIds[hs1_id][j] == hs2_id) {
			walkedEdge = true;
			break;
		}
	if (!walkedEdge) {
		edgePartnerIds[hs1_id].push_back(hs2_id);
		edgePartnerIds[hs2_id].push_back(hs1_id);
		if (isKernelFace[hs1_id] != ProcessColor::GREEN)
			isKernelFace[hs1_id] = ProcessColor::WHITE;
		if (isKernelFace[hs2_id] != ProcessColor::GREEN)
			isKernelFace[hs2_id] = ProcessColor::WHITE;
		return true;
	}
	return false;

}

void KerGen::findEdgeDirection(int hs1_id, int hs2_id, bool should_revert, double* edgeDirection, double* directioner) {

	crossProduct(halfSpaceSet[hs1_id].ABCD, halfSpaceSet[hs2_id].ABCD, edgeDirection);
	normalize(edgeDirection);
	if (directioner) {
		if (dotProduct(edgeDirection, directioner) > 0)
			multVect(edgeDirection, -1.0, edgeDirection);
	}
	else {
		if (should_revert)
			multVect(edgeDirection, -1.0, edgeDirection);
	}
}

bool KerGen::isValidEdge(double* startPoint, double* edgeDirection) {
	return true;
	double testPoint[3];
	for (int k = 0; k < 3; k++)
		testPoint[k] = startPoint[k] + EPSILON * edgeDirection[k];
	double t = findClosestValueSatisfiedByPoint(testPoint, halfSpaceSet);
	if (t > 0.5 * EPSILON)
		return false;
	return true;

}

void KerGen::findInitialPoint_1(double* point) {

	computeHalfSpacesFromTriangles(hostMeshptr->getAllTris(), hostMeshptr->getAllVerts(), halfSpaceSet);

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

void KerGen::findInitialPoint_2(double* point) {

	Vertex v = hostMeshptr->getVertex(0);
	for (int k = 0; k < 3; k++)
		point[k] = v.coords[k];

	for (int i = 0; i < hostMeshptr->getNumOfTris(); i++) {
		Triangle tri = hostMeshptr->getTriangle(i);
		HalfSpace halfSpace(hostMeshptr->getVertex(tri.corners[0]).coords, tri.normal, false);
		halfSpaceSet.push_back(halfSpace);

		double distance = halfSpaceSet[i].ABCD[3];
		for (int k = 0; k < 3; k++)
			distance += halfSpaceSet[i].ABCD[k] * point[k];

		if (distance > EPSILON) {
			double targetpoint[3];
			for (int k = 0; k < 3; k++)
				targetpoint[k] = -halfSpaceSet[i].ABCD[k] * distance + point[k];

			while (true) {
				bool satisfied = true;
				for (int j = 0; j < i; j++) {
					distance = halfSpaceSet[j].ABCD[3];
					for (int k = 0; k < 3; k++)
						distance += halfSpaceSet[j].ABCD[k] * targetpoint[k];
					if (distance > EPSILON) {
						satisfied = false;
						break;
					}
				}

				if (satisfied) {
					for (int k = 0; k < 3; k++)
						point[k] = targetpoint[k];
					break;
				}
				cout << "p: " << point[0] << " " << point[1] << " " << point[2] << endl;
				cout << "t: " << targetpoint[0] << " " << targetpoint[1] << " " << targetpoint[2] << endl;

				double goingDir[3];
				for (int k = 0; k < 3; k++)
					goingDir[k] = targetpoint[k] - point[k];
				double min_t = computeLength(goingDir);

				int theClosestId = -1;
				normalize(goingDir);
				Line goingLine(point, goingDir);
				for (int j = 0; j < i; j++) {
					distance = halfSpaceSet[j].ABCD[3];
					for (int k = 0; k < 3; k++)
						distance += halfSpaceSet[j].ABCD[k] * targetpoint[k];
					if (distance > EPSILON) {
						double t = findLinePlaneIntersection(goingLine, halfSpaceSet[j]);
						if (abs(t) < EPSILON) {
							min_t = t;
							theClosestId = j;
						}
						else if (t < min_t + EPSILON) {
							if (t < min_t)
								min_t = t;
							theClosestId = j;
						}
					}
				}
				for (int k = 0; k < 3; k++)
					point[k] += goingDir[k] * min_t;

				satisfied = true;
				for (int j = 0; j <= i; j++) {
					distance = halfSpaceSet[j].ABCD[3];
					for (int k = 0; k < 3; k++)
						distance += halfSpaceSet[j].ABCD[k] * point[k];
					if (distance > EPSILON) {
						satisfied = false;
						break;
					}
				}

				if (satisfied)
					break;

				double lineDir[3], goingDir2[3], controlpoint[3], controlDistance = -1.0;
				crossProduct(halfSpaceSet[theClosestId].ABCD, halfSpaceSet[i].ABCD, lineDir);
				normalize(lineDir);
				crossProduct(halfSpaceSet[theClosestId].ABCD, lineDir, goingDir2);
				normalize(goingDir2);
				Line goingLine2(point, goingDir2);
				double t = findLinePlaneIntersection(goingLine2, halfSpaceSet[i]);
				for (int k = 0; k < 3; k++) {
					targetpoint[k] = goingDir2[k] * t + point[k];
					controlpoint[k] = point[k] + goingDir2[k] * EPSILON;
				}

				/*
				if (prevId >= 0) {
					controlDistance = halfSpaceSet[prevId].ABCD[3];
					for (int k = 0; k < 3; k++)
						controlDistance += controlpoint[k] * halfSpaceSet[prevId].ABCD[k];
				}

				double zeroVector[3] = { 0, 0, 0 };
				if (isTripleSame(lineDir, zeroVector) ||
					isTripleSame(goingDir2, zeroVector) ||
					t == numeric_limits<double>::infinity() ||
					controlDistance > EPSILON ||
					t < EPSILON
					) {
					for (int k = 0; k < 3; k++)
						point[k] = numeric_limits<double>::infinity();
					return;
				}

				prevId = theClosestId;
				*/
			}
		}
	}
}

void KerGen::findInitialPoint_3(double* point) {

	Vertex v = hostMeshptr->getVertex(0);
	for (int k = 0; k < 3; k++)
		point[k] = v.coords[k];

	for (int i = 0; i < hostMeshptr->getNumOfTris(); i++) {

		// construct the halfspace
		Triangle tri = hostMeshptr->getTriangle(i);
		HalfSpace halfSpace(hostMeshptr->getVertex(tri.corners[0]).coords, tri.normal, false);
		halfSpaceSet.push_back(halfSpace);

		// calculate the distance of current point to i'th plane
		double distance = halfSpaceSet[i].ABCD[3];
		for (int k = 0; k < 3; k++)
			distance += halfSpaceSet[i].ABCD[k] * point[k];
		// if distance > 0, then it means point does not satisfy the i'th halfspace
		if (distance > 0) {
			// project the point onto i'th plane, now it is the new point
			for (int k = 0; k < 3; k++)
				point[k] = -halfSpaceSet[i].ABCD[k] * distance + point[k];

			// we made the previous point satisfy all the halfpaces upto i'th, check them for the new point
			for (int j = 0; j < i; j++) {
				// calculate the distance of the new point to j'th plane
				double distance = halfSpaceSet[j].ABCD[3];
				for (int k = 0; k < 3; k++)
					distance += halfSpaceSet[j].ABCD[k] * point[k];
				// if distance > 0, then it means that new point does not satisfy the j'th halfspace
				// but the previous point was satisfying j'th half-space
				// then we understand that j'th halfspace passes through between the previous point and the new point (current point)
				if (distance > EPSILON) {
					double lineDir[3], goingDir[3];
					crossProduct(halfSpaceSet[j].ABCD, halfSpaceSet[i].ABCD, lineDir);
					crossProduct(lineDir, halfSpaceSet[i].ABCD, goingDir);
					Line goingLine(point, goingDir);
					double t = findLinePlaneIntersection(goingLine, halfSpaceSet[j]);

					double zeroVector[3] = { 0, 0, 0 };
					if (isTripleSame(lineDir, zeroVector) || isTripleSame(goingDir, zeroVector) || t == numeric_limits<double>::infinity()) {
						for (int k = 0; k < 3; k++)
							point[k] = numeric_limits<double>::infinity();
						return;
					}

					for (int k = 0; k < 3; k++)
						point[k] = goingDir[k] * t + point[k];

					for (int m = 0; m < j; m++) {
						double distance = halfSpaceSet[m].ABCD[3];
						for (int k = 0; k < 3; k++)
							distance += halfSpaceSet[m].ABCD[k] * point[k];
						if (distance > 0) {
							for (int k = 0; k < 3; k++)
								point[k] = numeric_limits<double>::infinity();
							return;
						}
					}

				}
			}
		}
	}
}

void KerGen::findInitialPoint_4(double* point) {

	Vertex v = hostMeshptr->getVertex(0);
	for (int k = 0; k < 3; k++)
		point[k] = v.coords[k];

	point[0] = -19.4649;
	point[1] = 142.087;
	point[2] = 44.1074;

	for (int i = 0; i < hostMeshptr->getNumOfTris(); i++) {

		// construct the halfspace
		Triangle tri = hostMeshptr->getTriangle(i);
		HalfSpace halfSpace(hostMeshptr->getVertex(tri.corners[0]).coords, tri.normal, false);
		halfSpaceSet.push_back(halfSpace);

		// calculate the distance of current point to i'th plane
		double distance = halfSpaceSet[i].ABCD[3];
		for (int k = 0; k < 3; k++)
			distance += halfSpaceSet[i].ABCD[k] * point[k];
		// if distance > 0, then it means point does not satisfy the i'th halfspace
		if (distance > 3 * EPSILON) {
			// project the point onto i'th plane, now it is the new point
			for (int k = 0; k < 3; k++)
				point[k] = -halfSpaceSet[i].ABCD[k] * distance + point[k];

			bool found = false;
			double avg[3] = { 0, 0, 0 };
			// we made the previous point satisfy all the halfpaces upto i'th, check them for the new point
			for (int j = 0; j < i; j++) {
				// calculate the distance of the new point to j'th plane
				double distance = halfSpaceSet[j].ABCD[3];
				for (int k = 0; k < 3; k++)
					distance += halfSpaceSet[j].ABCD[k] * point[k];
				// if distance > 0, then it means that new point does not satisfy the j'th halfspace
				// but the previous point was satisfying j'th half-space
				// then we understand that j'th halfspace passes through between the previous point and the new point (current point)
				//if (distance > 3 * EPSILON) {
				double lineDir[3], goingDir[3];
				crossProduct(halfSpaceSet[j].ABCD, halfSpaceSet[i].ABCD, lineDir);
				crossProduct(lineDir, halfSpaceSet[i].ABCD, goingDir);
				Line goingLine(point, goingDir);
				double t = findLinePlaneIntersection(goingLine, halfSpaceSet[j]);

				double zeroVector[3] = { 0, 0, 0 };
				if (isTripleSame(lineDir, zeroVector) || isTripleSame(goingDir, zeroVector) || t == numeric_limits<double>::infinity()) {
					for (int k = 0; k < 3; k++)
						point[k] = numeric_limits<double>::infinity();
					return;
				}

				double newpoint[3];
				for (int k = 0; k < 3; k++) {
					newpoint[k] = goingDir[k] * t + point[k];
					avg[k] += newpoint[k];
				}

				bool not_yet = false;
				for (int m = 0; m < i; m++) {
					double distance = halfSpaceSet[m].ABCD[3];
					for (int k = 0; k < 3; k++)
						distance += halfSpaceSet[m].ABCD[k] * newpoint[k];
					if (distance > 3 * EPSILON) {
						not_yet = true;
						break;
					}
				}
				if (!not_yet) {
					for (int k = 0; k < 3; k++)
						point[k] = newpoint[k];
					found = true;
					break;
				}
			}
			//}
			if (!found) {
				for (int k = 0; k < 3; k++)
					avg[k] /= (i - 1);
				double t = findClosestValueSatisfiedByPoint(avg, halfSpaceSet);
				if (t > 3 * EPSILON) {
					for (int k = 0; k < 3; k++)
						point[k] = numeric_limits<double>::infinity();
					return;
				}
				else {
					for (int k = 0; k < 3; k++)
						point[k] = avg[k];
				}
			}
		}
	}

}

void KerGen::findInitialPoint_5(double* point) {

	Vertex v = hostMeshptr->getVertex(0);
	for (int k = 0; k < 3; k++)
		point[k] = v.coords[k];

	for (int i = 0; i < hostMeshptr->getNumOfTris(); i++) {

		// construct the halfspace
		Triangle tri = hostMeshptr->getTriangle(i);
		HalfSpace halfSpace(hostMeshptr->getVertex(tri.corners[0]).coords, tri.normal, false);
		halfSpaceSet.push_back(halfSpace);

		// calculate the distance of current point to i'th plane
		double distance = halfSpaceSet[i].ABCD[3];
		for (int k = 0; k < 3; k++)
			distance += halfSpaceSet[i].ABCD[k] * point[k];
		// if distance > 0, then it means point does not satisfy the i'th halfspace
		if (distance > EPSILON) {
			// project the point onto i'th plane, now it is the new point
			double newpoint[3];
			for (int k = 0; k < 3; k++)
				newpoint[k] = -halfSpaceSet[i].ABCD[k] * distance + point[k];

			// we made the previous point satisfy all the halfpaces upto i'th, check them for the new point
			vector<double> triedPoints;
			for (int j = 0; j < i; j++) {
				// calculate the distance of the new point to j'th plane
				double distance = halfSpaceSet[j].ABCD[3] + dotProduct(halfSpaceSet[j].ABCD, newpoint);
				// if distance > 0, then it means that new point does not satisfy the j'th halfspace
				// but the previous point was satisfying j'th half-space
				// then we understand that j'th halfspace passes through between the previous point and the new point (current point)
				if (distance > 1e-1 * EPSILON) {
					double lineDir[3], goingDir[3];
					crossProduct(halfSpaceSet[j].ABCD, halfSpaceSet[i].ABCD, lineDir);
					crossProduct(lineDir, halfSpaceSet[i].ABCD, goingDir);
					Line goingLine(newpoint, goingDir);
					double t = findLinePlaneIntersection(goingLine, halfSpaceSet[j]);

					bool is_loop = false;
					for (int t = 0; t < triedPoints.size(); t += 3) {
						if (abs(triedPoints[t] - newpoint[0]) < 1e-2 * EPSILON &&
							abs(triedPoints[t + 1] - newpoint[1]) < 1e-2 * EPSILON &&
							abs(triedPoints[t + 2] - newpoint[2]) < 1e-2 * EPSILON) {
							is_loop = true;
							break;
						}
					}

					double zeroVector[3] = { 0, 0, 0 };
					if (isTripleSame(lineDir, zeroVector) || isTripleSame(goingDir, zeroVector) || t == numeric_limits<double>::infinity() || is_loop) {
						for (int k = 0; k < 3; k++)
							point[k] = numeric_limits<double>::infinity();
						return;
					}

					for (int k = 0; k < 3; k++) {
						triedPoints.push_back(newpoint[k]);
						newpoint[k] = goingDir[k] * t + newpoint[k];
					}
					j = -1;
				}
			}

			triedPoints.clear();
			for (int k = 0; k < 3; k++)
				point[k] = newpoint[k];
		}
	}
}

void KerGen::filterRepetitions(double* distances, vector<double>& scalarsVector) {

	double _EPSILON = EPSILON * 3;
	vector<HalfSpace> temp_halfSpaceSet;

	for (int i = 0; i < halfSpaceSet.size(); i++) {
		bool isTheSame = false;
		for (int j = 0; j < temp_halfSpaceSet.size(); j++) {
			if (abs(distances[i] - scalarsVector[j]) < _EPSILON) {
				if (isTripleSame(halfSpaceSet[i].ABCD, temp_halfSpaceSet[j].ABCD)) {
					isTheSame = true;
					break;
				}
			}
		}
		if (!isTheSame) {
			temp_halfSpaceSet.push_back(halfSpaceSet[i]);
			scalarsVector.push_back(distances[i]);
		}
	}

	halfSpaceSet.clear();
	delete[] distances;
	distances = NULL;
	for (int j = 0; j < temp_halfSpaceSet.size(); j++)
		halfSpaceSet.push_back(temp_halfSpaceSet[j]);
	temp_halfSpaceSet.clear();

}

void KerGen::filterRepetitions2(double* distances, vector<double>& scalarsVector) {

	for (int i = 0; i < halfSpaceSet.size(); i++) {
		scalarsVector.push_back(distances[i]);
	}
	delete[] distances;
}

