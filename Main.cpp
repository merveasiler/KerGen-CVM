// @author Merve Asiler

#include "KernelComputation.h"
#include "SceneManager.h"
#include "CGALUtils.h"
#include "STLConverter.h"
#include "CommonUtils.h"

int main(int argc, char* argv[])
{

	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_DEBUG);
	{
		string command_type = "kergen_ideal";	// argv[1];
		//string shape_path = "D:/VS_Workspace/3D_Databases/DB-Basic/nonsymmetric.off"; // argv[2];
		string shape_path = "D:/VS_Workspace/3D_Databases/DB-Star-shaped-meshes/data/star.off";
		//string shape_path = "D:/VS_Workspace/3D_Databases/DB-ItalianStarShapes/Complex_Models/bot_eye.off";
		//string shape_path = "D:/VS_Workspace/3D_Databases/DB-ItalianStarShapes/Thingi/data/129975.off";
		//string shape_path = "D:/VS_Workspace/3D_Databases/DB-ItalianStarShapes/Refinements/vase/vase6.off";
		//string shape_path = "D:/VS_Workspace/3D_Databases/DB-StarCandidates/teddy.off";
		//string shape_path = "D:/VS_Workspace/3D_Databases/DB-Kids/0001.isometry.1.off";
		//string shape_path = "D:/VS_Workspace/3D_Databases/DB-Horse/1.obj";

		//string shape_path = "D:/VS_Workspace/3D_Databases/DB-Star-shaped-meshes/KernelResults-PolyhedronKernel/acorn_batch_polyhedronkernel_kernel.off";

		string source_path = "D:/VS_Workspace/3D_Databases/DB-Basic/nonsymmetric.off";
		string target_path = "D:/VS_Workspace/3D_Databases/DB-Basic/nonsymmetric_rotated.off";
		//string source_path = "D:/VS_Workspace/3D_Databases/DB-Star-shaped-meshes/19340_Star_v1.obj";
		//string target_path = "D:/VS_Workspace/3D_Databases/DB-Star-shaped-meshes/19340_Star_v1_rotated.off";
		//string source_path = "D:/VS_Workspace/3D_Databases/DB-ItalianStarShapes/Complex_Models/star.off";
		//string target_path = "D:/VS_Workspace/3D_Databases/DB-ItalianStarShapes/Complex_Models/star_rotated.off";

		//string shape_folder = "D:/VS_Workspace/3D_Databases/DB-Star-shaped-meshes/data";
		//string shape_folder = "D:/VS_Workspace/3D_Databases/DB-ItalianStarShapes/Thingi/data";
		string shape_folder = "D:/VS_Workspace/3D_Databases/DB-ItalianStarShapes/Refinements/data";
		//string shape_folder = "D:/VS_Workspace/3D_Databases/DB-Thingi/data";
		//string shape_folder = "D:/VS_Workspace/3D_Databases/DB-Princeton/data";
		//string shape_folder = "C:/Users/Merve/Desktop/problematics";

		//string output_folder = "D:/VS_Workspace/3D_Databases/DB-Star-shaped-meshes/KernelResults-KerGen-Ideal-2";
		//string output_folder = "D:/VS_Workspace/3D_Databases/DB-ItalianStarShapes/Thingi/KernelResults-KerGen";
		string output_folder = "D:/VS_Workspace/3D_Databases/DB-ItalianStarShapes/Refinements/KernelResults-KerGen-Ideal";
		//string output_folder = "D:/VS_Workspace/3D_Databases/DB-Thingi/KernelResults-KerGen";
		//string output_folder = "D:/VS_Workspace/3D_Databases/DB-Princeton/KernelResults-KerGen";
		//string output_folder = "C:/Users/Merve/Desktop/KernelResults";


		//string inpFileName = "D:/VS_Workspace/3D_Databases/DB-ItalianStarShapes/Thingi/KernelResults-KerTrack/KernelResults_batch_kertrack.txt";
		//string outFileName = "D:/VS_Workspace/3D_Databases/DB-ItalianStarShapes/Thingi/KernelResults-KerTrack/Statistics_kertrack";
		//string inpFileName = "D:/VS_Workspace/3D_Databases/DB-ItalianStarShapes/Refinements/KernelResults-CGAL/KernelResults_batch_kernel_by_cgal.txt";
		//string outFileName = "D:/VS_Workspace/3D_Databases/DB-ItalianStarShapes/Refinements/KernelResults-CGAL/Statistics_kernel_by_cgal";
		string inpFileName = "D:/VS_Workspace/3D_Databases/DB-Princeton/KernelResults-KerTrack/KernelResults_batch_kertrack.txt";
		string outFileName = "D:/VS_Workspace/3D_Databases/DB-Princeton/KernelResults-KerTrack/Statistics_kertrack";


		// DRAW:
		// Example:	draw C:/Users/Merve/3D_DATABASES/DB-Kids/method_alexa/11to15at0.5.off
		//			draw C:/Users/Merve/3D_DATABASES/DB-Camel/camel-02.obj
		if (command_type == "draw")
			drawMeshToScene(shape_path);

		// DRAW:
		// Example:	draw C:/Users/Merve/3D_DATABASES/DB_Kids/method_alexa/11to15at0.5.off
		//			draw C:/Users/Merve/3D_DATABASES/DB_Camel/camel-02.obj
		else if (command_type == "rotate")
			drawRotatedMeshToScene(shape_path);

		// COMPARE KERNEL RESULTS for CGAL & MDF-Ker-Plus & SDF-Ker-Plus
		else if (command_type == "experiment")
			doExperimentForPaper(shape_path);

		// COMPUTE KERNEL BY KERGEN_IDEAL
		// Example: kergen_ideal C:/Users/Merve/3D_DATABASES/DB-Star-shaped-meshes/Banana.obj
		else if (command_type == "kergen_ideal")
			ComputeKernel(shape_path, command_type);

		// COMPUTE KERNEL BY CGAL
		// Example: kernelByCGAL C:/Users/Merve/3D_DATABASES/DB-Star-shaped-meshes/Banana.obj
		else if (command_type == "kernel_by_cgal")
			ComputeKernel(shape_path, command_type);

		// COMPUTE KERNEL BY KERTRACK FOR ALL MESH FILES IN THE GIVEN FOLDER
		// Example: batch_kergen_ideal C:/Users/Merve/3D_DATABASES/DB-Star-shaped-meshes
		else if (command_type == "batch_kergen_ideal")
			ComputeBatchKernel(shape_folder, output_folder, command_type);

		// COMPUTE KERNEL BY CGAL FOR ALL MESH FILES IN THE GIVEN FOLDER
		// Example: batch_kernel_cgal C:/Users/Merve/3D_DATABASES/DB-Star-shaped-meshes
		else if (command_type == "batch_kernel_cgal")
			ComputeBatchKernel(shape_folder, output_folder, command_type);

		// FIND A KERNEL POINT MAXIMIZING A STATED COST FUNCTION by THIRD PARTY LIBRARY : SDLP
		// Example: sdlp C:/Users/Merve/3D_DATABASES/DB-Star-shaped-meshes/Banana.obj
		else if (command_type == "findkernelpoint_SDLP")
			FindKernelPoint_SDLP(shape_path);

		// APPLY SPHERICAL PARAMETRIZE by SENDING RAYS FROM A KERNEL POINT INSIDE THE MESH
		// Example: sp C:/Users/Merve/3D_DATABASES/DB-Star-shaped-meshes/Rock_6.obj
		else if (command_type == "sp")
			SphericalParametrize(shape_path);

		// APPLY SHAPE MORPH by INTERPOLATING VERTEX ANGLES AND DIRECTIONS FROM A KERNEL POINT
		// Example: sm C:/Users/Merve/3D_DATABASES/DB-Star-shaped-meshes/Rock_6.obj C:/Users/Merve/3D_DATABASES/DB-Star-shaped-meshes/liver.obj
		else if (command_type == "sm")
			ShapeMorphByKernel(source_path, target_path);

		// APPLY LINEAR INTERPOLATION on vertex positions
		// Example: lerp C:/Users/Merve/3D_DATABASES/DB-Star-shaped-meshes/Rock_6.obj C:/Users/Merve/3D_DATABASES/DB-Star-shaped-meshes/liver.obj
		else if (command_type == "lerp")
			ShapeMorphByLerp(source_path, target_path);

		// EXTRACT CONVEX HULL OF MESH
		// Example:	convexhull C:/Users/Merve/3D_DATABASES/DB_FaustRegistrations/tr_reg_000.off
		else if (command_type == "convexhull")
			computeConvexHull(shape_path);

		// CONVERT STL FILES TO OFF FILES
		// Example:	convertSTL C:/Users/Merve/3D_DATABASES/DB_FaustRegistrations/tr_reg_000.off
		else if (command_type == "convertSTL")
			convertFromSTLToOFF(shape_folder, output_folder);

		// FETCH NUMBER OF FACES AND TIME OBTAINED IN KERNEL COMPUTATION
		else if (command_type == "fetchStatistics")
			fetchStatistics(inpFileName, outFileName, true);

		// UNDEFINED:
		else
			cout << "Undefined Operation!" << endl;
	}
	_CrtDumpMemoryLeaks();

	return 0;

}
