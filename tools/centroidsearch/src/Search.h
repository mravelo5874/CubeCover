#ifndef SEARCH_H
#define SEARCH_H

#include <Eigen/Core>
#include <vector>

namespace CubeCover
{
    class TetMeshConnectivity;
    class FrameField;
};

void centroidBFS(
    const Eigen::MatrixXd& V,
    const CubeCover::TetMeshConnectivity& mesh,
    const CubeCover::FrameField& field,
    int startTetId,
    std::vector<Eigen::Vector2i>& tree_traversal
);

void centroidDFS(
    const Eigen::MatrixXd& V,
    const CubeCover::TetMeshConnectivity& mesh,
    const CubeCover::FrameField& field,
    int startTetId,
    std::vector<Eigen::Vector2i>& tree_traversal
);

void buildCurveNetwork(
    const Eigen::MatrixXd& V,
    const CubeCover::TetMeshConnectivity& mesh,
    const CubeCover::FrameField& field,
    Eigen::MatrixXd& nodes,
    Eigen::MatrixXi& edges,
    std::vector<Eigen::Vector2i>& tree_traversal
);

#endif