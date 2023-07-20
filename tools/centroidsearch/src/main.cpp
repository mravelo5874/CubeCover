#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "polyscope/polyscope.h"
#include "polyscope/curve_network.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/render/engine.h"

#include "ReadFrameField.h"
#include "FrameField.h"
#include "TetMeshConnectivity.h"
#include "readMeshFixed.h"
#include "Search.h"

#include <random>
#include <iostream>
#include <fstream>

int main(int argc, char *argv[])
{
    // check arguments
    if (argc != 4)
    {
        std::cerr << "Usage: centroidsearch [int start idx] [.mesh file] [.fra file]" << std::endl;
        return -1;
    }

    int start_node = std::stoi(argv[1]);
    std::string path_mesh = argv[2];
    std::string path_fra = argv[3];

    Eigen::MatrixXd V;
    Eigen::MatrixXi T;
    Eigen::MatrixXi F;

    CubeCover::TetMeshConnectivity mesh;
    if (!CubeCover::readMESH(path_mesh, V, T, F))
        return -1;
    mesh = CubeCover::TetMeshConnectivity(T);

    Eigen::MatrixXd frames;
    Eigen::MatrixXi assignments;
    
    if (!CubeCover::readFrameField(path_fra, "", T, frames, assignments, true))
        return -1;
    
    CubeCover::FrameField* field;
    field = CubeCover::fromFramesAndAssignments(mesh, frames, assignments, true);
    if (!field)
        return -1;

    std::cout.flush();
    field->computeLocalAssignments();
    field->combAssignments();

    // perform search
    std::vector<Eigen::Vector2i> tree_traversal;
    centroidBFS(V, mesh, *field, start_node % mesh.nTets(), tree_traversal);

    // create curve network from tree_traversal
    Eigen::MatrixXd nodes;
    Eigen::MatrixXi edges;
    buildCurveNetwork(V, mesh, *field, nodes, edges, tree_traversal);

    std::cout << "nodes: " << nodes.size() << std::endl;
    std::cout << "edges: " << edges.size() << std::endl;

    // make a mesh out of all of the boundary faces
    int nbdry = 0;
    int nfaces = mesh.nFaces();
    for (int i = 0; i < nfaces; i++)
    {
        if (mesh.isBoundaryFace(i))
            nbdry++;
    }
    Eigen::MatrixXi bdryF(nbdry, 3);
    int curidx = 0;
    for (int i = 0; i < nfaces; i++)
    {
        if (mesh.isBoundaryFace(i))
        {
            for (int j = 0; j < 3; j++)
            {
                bdryF(curidx, j) = mesh.faceVertex(i, j);

            }
            // fix triangle orientations
            int tet = mesh.faceTet(i, 0);
            if (tet == -1)
            {
                std::swap(bdryF(curidx, 0), bdryF(curidx, 1));
            }
            curidx++;
        }
    }

    // polyscope stuff!
    {
        polyscope::init();
           
        // create curve network
        auto *curve_network = polyscope::registerCurveNetwork("search_network", nodes, edges);

        // color edges based on search depth
        std::vector<std::array<double, 3>> edge_colors(edges.size() / 2);
        polyscope::render::ValueColorMap colormap = polyscope::render::engine->getColorMap("viridis");

        for (size_t i = 0; i < edges.size() / 2; i++) {
            double val = (double)i / (double)(edges.size() / 2);
            glm::vec3 color = colormap.getValue(val);
            edge_colors[i] = {{color.r, color.g, color.b}};
        }
        polyscope::getCurveNetwork("search_network")->addEdgeColorQuantity("search_color", edge_colors);

        // surface mesh
        auto *psMesh = polyscope::registerSurfaceMesh("Boundary Mesh", V, bdryF);
        psMesh->setTransparency(0.1);
        psMesh->setSurfaceColor({ 0.5, 0.5, 0.5 });

        // visualize!
        polyscope::show();
    }
}


