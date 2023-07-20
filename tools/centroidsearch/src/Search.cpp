#include <Eigen/Geometry>

#include "Search.h"
#include "TetMeshConnectivity.h"
#include "FrameField.h"

#include <iostream>
#include <deque>
#include <string>
#include <map>

/* 
This function takes in a source id, and does breadth first search on the cetroids in the mesh.
*/
void centroidBFS(
    const Eigen::MatrixXd& V,
    const CubeCover::TetMeshConnectivity& mesh,
    const CubeCover::FrameField& field,
    int startTetId,
    std::vector<Eigen::Vector2i>& tree_traversal
)
{
    int num_tets = mesh.nTets();
    std::vector<int> tet_added_to_queue(num_tets, 0);
    std::deque<int> tet_queue;
    std::vector<Eigen::Vector2i> edge_queue; // intexed as tetid, tetedge 
    tree_traversal.clear();

    int cur_tet_id = startTetId;

    std::cout << "num_tets: " << num_tets << std::endl;
    std::cout << "starting tet: " << cur_tet_id << std::endl;

    // init search
    int init_edge_id = mesh.tetEdge(cur_tet_id, 0);
    Eigen::Vector2i init_edge = Eigen::Vector2i(mesh.edgeVertex(init_edge_id, 0), mesh.edgeVertex(init_edge_id, 1));

    tet_added_to_queue[cur_tet_id] += 1;
    tet_queue.push_back(cur_tet_id);

    int iter = 0;
    while (!tet_queue.empty()) 
    {
        cur_tet_id = tet_queue.front();
        tet_queue.pop_front();
        
        // iterate through all 4 potential tet neighboors
        for (int tetVert = 0; tetVert < 4; tetVert++)
        {
            int opp_tet_id = mesh.tetOppositeVertex(cur_tet_id, tetVert);
            int cur_face_id = mesh.tetFace(cur_tet_id, tetVert); // this might be wrong
            
            // make sure neighboor tet is vaild
            if (opp_tet_id > -1 && field.faceAssignment(cur_face_id).isIdentity())
            {
                // if tet is not in queue, add it
                if (tet_added_to_queue[opp_tet_id] == 0)
                {
                    tet_added_to_queue[opp_tet_id] = 1;
                    tet_queue.push_back(opp_tet_id);
                       
                    // add graph edge
                    Eigen::Vector2i graph_edge = Eigen::Vector2i(cur_tet_id, opp_tet_id);
                    tree_traversal.push_back(graph_edge);

                    // std::cout << "graph-edge: " << cur_tet_id << " -> " << opp_tet_id << std::endl;
                }
            }
        }
        iter++;
    }
}

/*
This function takes in a source id, and does depth first search on the cetroids in the mesh.
*/
void centroidDFS(
    const Eigen::MatrixXd& V,
    const CubeCover::TetMeshConnectivity& mesh,
    const CubeCover::FrameField& field,
    int startTetId,
    std::vector<Eigen::Vector2i>& tree_traversal
)
{
    int num_tets = mesh.nTets();
    std::vector<int> tet_added_to_queue(num_tets, 0);
    std::deque<int> tet_queue;
    std::vector<Eigen::Vector2i> edge_queue; // intexed as tetid, tetedge 
    tree_traversal.clear();

    int cur_tet_id = startTetId;

    std::cout << "num_tets: " << num_tets << std::endl;
    std::cout << "starting tet: " << cur_tet_id << std::endl;

    // init search
    int init_edge_id = mesh.tetEdge(cur_tet_id, 0);
    Eigen::Vector2i init_edge = Eigen::Vector2i(mesh.edgeVertex(init_edge_id, 0), mesh.edgeVertex(init_edge_id, 1));

    tet_added_to_queue[cur_tet_id] += 1;
    tet_queue.push_back(cur_tet_id);

    int iter = 0;
    while (!tet_queue.empty())
    {
        cur_tet_id = tet_queue.front();
        tet_queue.pop_front();

        // iterate through all 4 potential tet neighboors
        for (int tetVert = 0; tetVert < 4; tetVert++)
        {
            int opp_tet_id = mesh.tetOppositeVertex(cur_tet_id, tetVert);
            int cur_face_id = mesh.tetFace(cur_tet_id, tetVert); // this might be wrong

            // make sure neighboor tet is vaild
            if (opp_tet_id > -1 && field.faceAssignment(cur_face_id).isIdentity())
            {
                // if tet is not in queue, add it
                if (tet_added_to_queue[opp_tet_id] == 0)
                {
                    tet_added_to_queue[opp_tet_id] = 1;
                    tet_queue.push_back(opp_tet_id);

                    // add graph edge
                    Eigen::Vector2i graph_edge = Eigen::Vector2i(cur_tet_id, opp_tet_id);
                    tree_traversal.push_back(graph_edge);

                    // std::cout << "graph-edge: " << cur_tet_id << " -> " << opp_tet_id << std::endl;
                }
            }
        }
        iter++;
    }
}

/*
    This function takes in the result of a centroid search and creates a curve network for polyscope
*/
void buildCurveNetwork(
    const Eigen::MatrixXd& V,
    const CubeCover::TetMeshConnectivity& mesh,
    const CubeCover::FrameField& field,
    Eigen::MatrixXd& nodes, 
    Eigen::MatrixXi& edges,
    std::vector<Eigen::Vector2i>& tree_traversal
)
{
    // add all tet centroids
    int num_tets = mesh.nTets();
    nodes.resize(num_tets, 3);
    for (int i = 0; i < num_tets; i++)
    {
        Eigen::Vector3d cent(0, 0, 0);

        // iterate through all 4 tet vertex
        for (int tetVert = 0; tetVert < 4; tetVert++)
        {
            int vert_idx = mesh.tetVertex(i, tetVert);
            Eigen::Vector3d vert_pos = V.row(vert_idx).transpose();
            cent += vert_pos;
        }
        cent /= 4;
        nodes.row(i) = cent;
    }

    // transfer tree_traversal edges to matrix
    int num_edges = tree_traversal.size();
    edges.resize(tree_traversal.size(), 2);
    for (int i = 0; i < num_edges; i++)
    {
        edges(i, 0) = tree_traversal.at(i)[0];
        edges(i, 1) = tree_traversal.at(i)[1];
    }
}