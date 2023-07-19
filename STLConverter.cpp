// @author Merve Asiler

#include "STLReader.h"
#include "STLConverter.h"
#include "Mesh.h"

void convertFromSTLToOFF(string inputFolder, string outputFolder) {

    // Read folder, fecth mesh names
    int index = 0;
    path p(inputFolder);
    for (auto i = directory_iterator(p); i != directory_iterator(); i++)
    {
        if (!is_directory(i->path())) //we eliminate directories 
        {
            string meshName = i->path().filename().string();
            
            if (meshName.substr(meshName.length() - 3) != "stl")
                continue;

            try {
                stl_reader::StlMesh <float, unsigned int> stlmesh(inputFolder + "/" + meshName);
                
                std::ofstream outFile;
                outFile.open(outputFolder + "/" + meshName.substr(0, meshName.length() - 3) + "off");
                outFile << "OFF" << endl;

                // Write the number of vertices and triangles (dummy is the number of edges)
                outFile << stlmesh.num_vrts() << " " << stlmesh.num_tris() << " 0" << endl;

                // Write the vertices
                for (int j = 0; j < stlmesh.num_vrts(); j++)
                    outFile << stlmesh.vrt_coords(j)[0] << " " << stlmesh.vrt_coords(j)[1] << " " << stlmesh.vrt_coords(j)[2] << endl;

                // Write the triangles
                for (int j = 0; j < stlmesh.num_tris(); j++)
                    outFile << "3 " << stlmesh.tri_corner_ind(j, 0) << " " << stlmesh.tri_corner_ind(j, 1) << " " << stlmesh.tri_corner_ind(j, 2) << endl;

                outFile.close();
                
                cout << index++ << endl;
            }
            catch (std::exception& e) {
                e.what();
            }
            
        }
    }

}