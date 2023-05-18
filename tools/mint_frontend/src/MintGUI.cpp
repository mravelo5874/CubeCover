#include <igl/file_dialog_open.h>
#include <igl/readOBJ.h>
#include <igl/bfs_orient.h>

#include "MintGUI.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/volume_mesh.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"

// #include <imgui/misc/cpp/imgui_stdlib.h>
// #include <polyscope/deps/imgui/imgui/imgui.h>
// #include <imgui/imgui.h>  
#include "TetMeshConnectivity.h"
#include "readMeshFixed.h"
#include "ReadMoments.h"

#include "FrameFieldVis.h"
#include "FrameField.h"
// #include "TetMeshConnectivity.h"
#include "ReadFrameField.h"
#include "SingularCurveNetwork.h"

#include <misc/cpp/imgui_stdlib.h>
#include <algorithm>
#include <cctype>
#include <string>


#include <nlohmann/json.hpp>

// for convenience
using json = nlohmann::json;


// This is c++14 experimental feature.  In c++17 >= this is part of STD.
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;


namespace MintFrontend
{


static void HelpMarker(const char* desc)
{
    ImGui::TextDisabled("(?)");
    if (ImGui::IsItemHovered())
    {
        ImGui::BeginTooltip();
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
        ImGui::TextUnformatted(desc);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}


	MintGUI::MintGUI()
	{
        path_mesh = new char[512];
        path_fra = new char[512];
        path_constraints = new char[512];
        path_outdir = new char[512];

        mesh = CubeCover::TetMeshConnectivity();
		moment_view_mode = Moments_To_Show::fourth;
        frame_field_view_mode = Frames_To_Show::frames;

        showBoundary = false;
        showInteriorTets = true;
        useSameColorRangeForAllMoments = false;


        const char *args[] = {"x^4", "x^3 y", "x^3 z", "x^2 y^2", "x^2 y z", 
                              "x^2 z^2", "x y^3", "x y^2 z", "x y z^2", "x z^3", 
                              "y^4", "y^3 z", "y^2 z^2", "y z^3", "z^4", 
							  "1", 
							  "x^2", "y^2", "z^2", 
			                  "xy", "xz", "yz"};
        std::vector<std::string> tmp(args, std::end(args));
        moment_labels = tmp;



        // exploded_spacing = 120.; // square

        exploded_spacing_inp = 1.;
        exploded_spacing_intern = 1.;

        color_range_min = 0.;
        color_range_max = 1.;

        sel_idx_mom = -1;
        sel_idx_fra = -1;

        cur_solver = Mint_Linear_Solver::exact;
        mint_mode = Mint_Integrability_Mode::free;
        shell_mode = Mint_Frame_Projection::offshell;
        metric_mode = Mint_Moment_Metric::four;

        transparancy_interior = .7;
        transparancy_boundary = .9;



        // Load config file here. 

        
	}


// TODO this is kinda buggy on linux and could work better...

char * fileSelectSubroutine()
{
    std::string picked_file = igl::file_dialog_open();
    return strdup(picked_file.c_str());
    // return picked_file.data();
}






//! Using only text manipulation, splits a full path into component file parts
FileParts MintGUI::fileparts(const std::string &fullpath)
{
    using namespace std;

    size_t idxSlash = fullpath.rfind("/");
    if (idxSlash == string::npos) {
        idxSlash = fullpath.rfind("\\");
    }
    size_t idxDot = fullpath.rfind(".");

    FileParts fp;
    if (idxSlash != string::npos && idxDot != string::npos) {
        fp.path = fullpath.substr(0, idxSlash + 1);
        fp.name = fullpath.substr(idxSlash + 1, idxDot - idxSlash - 1);
        fp.ext  = fullpath.substr(idxDot);
    } else if (idxSlash == string::npos && idxDot == string::npos) {
        fp.name = fullpath;
    } else if (/* only */ idxSlash == string::npos) {
        fp.name = fullpath.substr(0, idxDot);
        fp.ext  = fullpath.substr(idxDot);
    } else { // only idxDot == string::npos
        fp.path = fullpath.substr(0, idxSlash + 1);
        fp.name = fullpath.substr(idxSlash + 1);
    }
    return fp;
}


void MintGUI::show_base_mesh()
{
    std::cout << "show_base_mesh" << std::endl;
    clear_polyscope_state();
    polyscope::options::automaticallyComputeSceneExtents = true;

    auto tet_mesh = polyscope::registerTetMesh("tet_mesh", V, T)->setEdgeWidth(0.5)->setTransparency(.7);
    auto surf_mesh = polyscope::registerSurfaceMesh("surf_mesh", V, bdryF)->setEdgeWidth(1)->setTransparency(.7);

    rescale_structure(tet_mesh);
    rescale_structure(surf_mesh);

    // tet_mesh->resetTransform();
    // surf_mesh->resetTransform();
    
        // std::cout << "V: " << V.size() << std::endl;
        // std::cout << "T: " << T << std::endl;

    polyscope::view::resetCameraToHomeView();

	std::cout << "end show_base_mesh" << std::endl;
    //                 polyscope::state::boundingBox = 
    // std::tuple<glm::vec3, glm::vec3>{ {-1., -1., -1.}, {1., 1., 1.} };
}


void MintGUI::show_exploded_moments(MintFrontend::Moments_To_Show moment_view_mode)
{
	std::cout << "show_exploded_moments" << std::endl;

    if (M_curr.cols() != 22)
    {
        show_moments_all();
    }
    else
    {
	if (moment_view_mode == MintFrontend::fourth)
	{
		show_moments_4th();
	}
	if (moment_view_mode == MintFrontend::second)
	{
		show_moments_2nd();
	}
	if (moment_view_mode == MintFrontend::both)
	{
		show_moments_4th();
		show_moments_2nd();
	}
    }



	std::cout << "end show_exploded_moments" << std::endl;

}


void MintGUI::show_frame_field(MintFrontend::Frames_To_Show frame_field_view_mode)
{
    std::cout << "show_exploded_frame_field" << std::endl;


    clear_polyscope_state();
    set_base_mesh();
    set_frame_field();

	if (frame_field_view_mode == MintFrontend::Frames_To_Show::frames)
	{
		show_gl3_frame_field();
	}
	if (frame_field_view_mode == MintFrontend::Frames_To_Show::split_frames)
	{
		// show_moments_2nd();
	}
	if (frame_field_view_mode == MintFrontend::Frames_To_Show::split_moments)
	{
		// show_moments_4th();
		// show_moments_2nd();
	}
    if (frame_field_view_mode == MintFrontend::Frames_To_Show::split_difference)
	{
		// show_moments_4th();
		// show_moments_2nd();
	}
    



	std::cout << "end show_exploded_frame_field" << std::endl;
}


void MintGUI::show_gl3_frame_field()
{
   std::cout << "show_constraint_vals" << std::endl;
    

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

        auto *tetc = polyscope::registerPointCloud("Centroids", centroids);
        glm::vec3 dotcolor(0.1, 0.1, 0.1);
        tetc->setPointColor(dotcolor);
        tetc->setPointRadius(0.001);
        int vpf = framefieldvecs.size();
        for (int i = 0; i < vpf; i++)
        {
            std::stringstream ss;
            ss << "Frame Vector " << i;
            auto *vf = tetc->addVectorQuantity(ss.str(), framefieldvecs[i]);
            double mag = framefieldvecs[i].row(0).norm() +  framefieldvecs[i].row(1).norm() +  framefieldvecs[i].row(2).norm();
            std::cout << mag << std::endl;

            vf->setVectorColor({ dist(rng),dist(rng),dist(rng) });
            vf->setVectorRadius(0.001);
            vf->setEnabled(true);
        }

        rescale_structure(tetc);

        

        auto *green = polyscope::registerCurveNetwork("singularity(+1/4)", Pgreen, Egreen);
        green->setColor({ 0.0,1.0,0.0 });

        auto *blue = polyscope::registerCurveNetwork("singularity(-1/4)", Pblue, Eblue);
        blue->setColor({ 0.0,0.0,1.0 });

        auto *black = polyscope::registerCurveNetwork("singularity(other)", Pblack, Eblack);
        black->setColor({ 0.0,0.0,0.0 });

}





void MintGUI::set_frame_field()
{

    Eigen::MatrixXd frames;
    Eigen::MatrixXi assignments;
    if (!CubeCover::readFrameField(path_fra, "", T, frames, assignments, true))
        return ;

    // CubeCover::TetMeshConnectivity mesh(T);
    
    field = CubeCover::fromFramesAndAssignments(mesh, frames, assignments, true);
    if (!field)
        return ;

    // if (recomputeperms)
    // {
    std::cout << "No face assignments provided, recomputing: ";
    std::cout.flush();
    field->computeLocalAssignments();
    std::cout << "found " << field->nSingularEdges() << " singular edges" << std::endl;
    // }
    field->combAssignments();


    extractSingularCurveNetwork(V, mesh, *field, Pgreen, Egreen, Pblue, Eblue, Pblack, Eblack);

    buildFrameVectors(V, mesh, *field, 1.0, centroids, framefieldvecs);

/*{


    // visualize the seams

    std::vector<int> seamfaces;

    int ninverted = 0;
    int nnontrivial = 0;

    for (int i = 0; i < nfaces; i++)
    {
        if (!field->faceAssignment(i).isIdentity())
        {
            seamfaces.push_back(i);
            nnontrivial++;
        }

        if (field->faceAssignment(i).orientation() == -1)
            ninverted++;
    }

    std::cout << "Non-identity face assignments: " << nnontrivial << std::endl;
    if (ninverted > 0)
    {
        std::cout << "Warning: " << ninverted << " face assignments are orientation-reversing" << std::endl;
    }

    int nseamtris = seamfaces.size();

    Eigen::MatrixXd seamV(3 * nseamtris, 3);
    Eigen::MatrixXi seamF(nseamtris, 3);
    for (int i = 0; i < nseamtris; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            seamF(i, j) = 3 * i + j;
            seamV.row(3 * i + j) = V.row(mesh.faceVertex(seamfaces[i], j));
        }
    }


*/

}



void MintGUI::rescale_structure(polyscope::Structure* m)
{
    m->resetTransform();
    // m->centerBoundingBox();
    std::tuple<glm::vec3, glm::vec3> bbox = m->boundingBox();
    glm::vec3 bbox_diff = std::get<1>(bbox) - std::get<0>(bbox);
    double scale_fac = std::max(bbox_diff.x, bbox_diff.y);

    // std::cout << "scale fac " << scale_fac << std::endl;
    exploded_spacing_intern = scale_fac;
    double currScale = scale_fac;
    float s = static_cast<float>(1.0 / scale_fac);
    glm::mat4x4 newTrans = glm::scale(glm::mat4x4(1.0), glm::vec3{s, s, s});
    m->setTransform(newTrans * m->getTransform());
    m->centerBoundingBox();

}


///// HACKY :( 
void MintGUI::show_moments_all()
{

    std::cout << "show_constraint_vals" << std::endl;
    clear_polyscope_state();
    // polyscope::options::automaticallyComputeSceneExtents = true;
    // polyscope::state::lengthScale = polyscope::state::lengthScale * 1/5.;

    polyscope::options::automaticallyComputeSceneExtents = false;

    polyscope::state::boundingBox = std::tuple<glm::vec3, glm::vec3>{glm::vec3{-10,-5, -1.}, glm::vec3{6,2, 1.}};
    polyscope::state::lengthScale = 1.2;


    std::pair<double,double> colrng = std::pair<double,double>(color_range_min,color_range_max);

    int tot_cols = M_curr.cols();

    for( int i = 0; i < 3; i++)
    {
        for( int j = 0; j < 5; j++)
        {
            int cur_id = i*5 + j;
            if (cur_id >= tot_cols)
                break;

            std::cout << glm::to_string(std::get<0>(polyscope::state::boundingBox) ) << " bounding box " << glm::to_string(std::get<1>(polyscope::state::boundingBox) ) << " bounding box " << polyscope::state::lengthScale << " length scale " << std::endl;

            float exploded_spacing = exploded_spacing_intern * exploded_spacing_inp;
            glm::vec3 shift(-j * exploded_spacing, -i * exploded_spacing, 0);

            if (showInteriorTets)
            {
                // auto cur_mesh = polyscope::registerTetMesh(moment_labels.at(cur_id) + "_tm_" + std::to_string(1000000 + cur_id), V, T);
                auto cur_mesh = polyscope::registerTetMesh(std::to_string(10000 + cur_id) + "____" +moment_labels.at(cur_id), V, T);


                cur_mesh->setEdgeWidth(0.5)->setTransparency(transparancy_interior);
                rescale_structure(cur_mesh);
 

                cur_mesh->translate(shift);
                Eigen::VectorXd tmp_mvals = M_curr.block(0,cur_id,mesh.nTets(),cur_id+1);
                cur_mesh->addCellScalarQuantity("cur-moment", tmp_mvals)->setEnabled(true);
                if (useSameColorRangeForAllMoments)
                {
                    cur_mesh->addCellScalarQuantity("cur-moment", tmp_mvals)->setMapRange(colrng)->setEnabled(true);
                }
            }

            if (showBoundary)
            {
                // auto surf_mesh = polyscope::registerSurfaceMesh(moment_labels.at(cur_id) + "_sm_" + std::to_string(1000000 + cur_id), V, bdryF);
                auto surf_mesh = polyscope::registerSurfaceMesh(std::to_string(20000 + cur_id) + "____" + moment_labels.at(cur_id), V, bdryF);



                surf_mesh->setEdgeWidth(1)->setTransparency(transparancy_boundary);
                rescale_structure(surf_mesh);
                surf_mesh->translate(shift);

                // std::cout << M_curr.rows()-mesh.nTets() << " diff " << M_curr.rows()-mesh.nTets() - bdryF.rows() << std::endl;
                Eigen::VectorXd tmp_mvals = M_curr.block(mesh.nTets(),cur_id,bdryF.rows(),cur_id+1);
                surf_mesh->addFaceScalarQuantity("cur-moment", tmp_mvals)->setEnabled(true);
                if (useSameColorRangeForAllMoments)
                {
                    surf_mesh->addFaceScalarQuantity("cur-moment", tmp_mvals)->setMapRange(colrng)->setEnabled(true);
                }
            }


            // std::cout << "tet_mesh_" + std::to_string(1000 + cur_id) + " " <<  glm::to_string(cur_mesh->getTransform()) << std::endl;
        }
    }
    // polyscope::registerTetMesh("tet_mesh", V, T)->setEdgeWidth(0.5)->setTransparency(.7);
    // polyscope::registerSurfaceMesh("surf_mesh", V, bdryF)->setEdgeWidth(1)->setTransparency(.7);
    polyscope::options::automaticallyComputeSceneExtents = false;
    polyscope::state::lengthScale = polyscope::state::lengthScale * 3.;

    polyscope::view::resetCameraToHomeView();



}


void MintGUI::show_moments_4th()
{

    std::cout << "show_constraint_vals" << std::endl;
    clear_polyscope_state();
    polyscope::options::automaticallyComputeSceneExtents = false;
    // polyscope::state::lengthScale = polyscope::state::lengthScale * 1/5.;

    polyscope::state::boundingBox = std::tuple<glm::vec3, glm::vec3>{glm::vec3{-10,-5, -1.}, glm::vec3{6,2, 1.}};
    polyscope::state::lengthScale = 1.4;


    std::pair<double,double> colrng = std::pair<double,double>(color_range_min,color_range_max);

    for( int i = 0; i < 3; i++)
    {
        for( int j = 0; j < 5; j++)
        {
            int cur_id = i*5 + j;

            float exploded_spacing = exploded_spacing_intern * exploded_spacing_inp;
            glm::vec3 shift(-j * exploded_spacing, -i * exploded_spacing, 0);

            if (showInteriorTets)
            {
                // auto cur_mesh = polyscope::registerTetMesh(moment_labels.at(cur_id) + "_tm_" + std::to_string(1000000 + cur_id), V, T);
                auto cur_mesh = polyscope::registerTetMesh(std::to_string(10000 + cur_id) + "____" +moment_labels.at(cur_id), V, T);


                cur_mesh->setEdgeWidth(0.5)->setTransparency(.7);
                rescale_structure(cur_mesh);
                cur_mesh->translate(shift);
                Eigen::VectorXd tmp_mvals = M_curr.block(0,cur_id,mesh.nTets(),cur_id+1);
                cur_mesh->addCellScalarQuantity("cur-moment", tmp_mvals)->setEnabled(true);
                if (useSameColorRangeForAllMoments)
                {
                    cur_mesh->addCellScalarQuantity("cur-moment", tmp_mvals)->setMapRange(colrng)->setEnabled(true);
                }
            }

            if (showBoundary)
            {
                // auto surf_mesh = polyscope::registerSurfaceMesh(moment_labels.at(cur_id) + "_sm_" + std::to_string(1000000 + cur_id), V, bdryF);
                auto surf_mesh = polyscope::registerSurfaceMesh(std::to_string(20000 + cur_id) + "____" + moment_labels.at(cur_id), V, bdryF);



                surf_mesh->setEdgeWidth(1)->setTransparency(.9);
                rescale_structure(surf_mesh);
                surf_mesh->translate(shift);

                // std::cout << M_curr.rows()-mesh.nTets() << " diff " << M_curr.rows()-mesh.nTets() - bdryF.rows() << std::endl;
                Eigen::VectorXd tmp_mvals = M_curr.block(mesh.nTets(),cur_id,bdryF.rows(),cur_id+1);
                surf_mesh->addFaceScalarQuantity("cur-moment", tmp_mvals)->setEnabled(true);
                if (useSameColorRangeForAllMoments)
                {
                    surf_mesh->addFaceScalarQuantity("cur-moment", tmp_mvals)->setMapRange(colrng)->setEnabled(true);
                }
            }


            // std::cout << "tet_mesh_" + std::to_string(1000 + cur_id) + " " <<  glm::to_string(cur_mesh->getTransform()) << std::endl;
        }
    }
    // polyscope::registerTetMesh("tet_mesh", V, T)->setEdgeWidth(0.5)->setTransparency(.7);
    // polyscope::registerSurfaceMesh("surf_mesh", V, bdryF)->setEdgeWidth(1)->setTransparency(.7);
    polyscope::options::automaticallyComputeSceneExtents = false;
    polyscope::state::lengthScale = polyscope::state::lengthScale * 3.;

    polyscope::view::resetCameraToHomeView();
    // polyscope::state::boundingBox = 
    //     std::tuple<glm::vec3, glm::vec3>{ {-2.5, -1.5, -1.}, {2.5, 1.5, 1.} };
}


void MintGUI::show_moments_2nd()
{

	std::cout << "show_constraint_vals" << std::endl;
	clear_polyscope_state();
	// polyscope::options::automaticallyComputeSceneExtents = true;
	// // polyscope::state::lengthScale = polyscope::state::lengthScale * 1/5.;

    polyscope::options::automaticallyComputeSceneExtents = false;

    float exploded_spacing = exploded_spacing_intern * exploded_spacing_inp;
    polyscope::state::boundingBox = std::tuple<glm::vec3, glm::vec3>{glm::vec3{1,-2, -1.}*exploded_spacing, glm::vec3{3,1, 1.}*exploded_spacing};
    polyscope::state::lengthScale = 1.2;

    std::pair<double,double> colrng = std::pair<double,double>(color_range_min,color_range_max);

	for (int i = 0; i < 6; i++)
	{

			int cur_id = i + 16;

			double ang = i / 6. * 2. * 3.141562;
			glm::vec3 shift( (2 + std::cos(ang) ) * exploded_spacing, std::sin(ang) * exploded_spacing, 0);

			if (showInteriorTets)
			{
				// auto cur_mesh = polyscope::registerTetMesh(moment_labels.at(cur_id) + "_tm_" + std::to_string(1000000 + cur_id), V, T);
				auto cur_mesh = polyscope::registerTetMesh(std::to_string(10000 + cur_id) + "____" + moment_labels.at(cur_id), V, T);


				cur_mesh->setEdgeWidth(0.5)->setTransparency(.7);
				rescale_structure(cur_mesh);
				cur_mesh->translate(shift);
				Eigen::VectorXd tmp_mvals = M_curr.block(0, cur_id, mesh.nTets(), cur_id + 1);
                cur_mesh->addCellScalarQuantity("cur-moment", tmp_mvals)->setEnabled(true);
                if (useSameColorRangeForAllMoments)
                {
                    cur_mesh->addCellScalarQuantity("cur-moment", tmp_mvals)->setMapRange(colrng)->setEnabled(true);
                }
				
			}

			if (showBoundary)
			{
				// auto surf_mesh = polyscope::registerSurfaceMesh(moment_labels.at(cur_id) + "_sm_" + std::to_string(1000000 + cur_id), V, bdryF);
				auto surf_mesh = polyscope::registerSurfaceMesh(std::to_string(20000 + cur_id) + "____" + moment_labels.at(cur_id), V, bdryF);



				surf_mesh->setEdgeWidth(1)->setTransparency(.9);
				rescale_structure(surf_mesh);
				surf_mesh->translate(shift);

				// std::cout << M_curr.rows()-mesh.nTets() << " diff " << M_curr.rows()-mesh.nTets() - bdryF.rows() << std::endl;
				Eigen::VectorXd tmp_mvals = M_curr.block(mesh.nTets(), cur_id, bdryF.rows(), cur_id + 1);
				surf_mesh->addFaceScalarQuantity("cur-moment", tmp_mvals)->setEnabled(true);
                if (useSameColorRangeForAllMoments)
                {
                    surf_mesh->addFaceScalarQuantity("cur-moment", tmp_mvals)->setMapRange(colrng)->setEnabled(true);
                }
			}


			// std::cout << "tet_mesh_" + std::to_string(1000 + cur_id) + " " <<  glm::to_string(cur_mesh->getTransform()) << std::endl;
		
	}
	// polyscope::registerTetMesh("tet_mesh", V, T)->setEdgeWidth(0.5)->setTransparency(.7);
	// polyscope::registerSurfaceMesh("surf_mesh", V, bdryF)->setEdgeWidth(1)->setTransparency(.7);
	polyscope::options::automaticallyComputeSceneExtents = false;
	polyscope::state::lengthScale = polyscope::state::lengthScale * 3.;

	polyscope::view::resetCameraToHomeView();
	// polyscope::state::boundingBox = 
	//     std::tuple<glm::vec3, glm::vec3>{ {-2.5, -1.5, -1.}, {2.5, 1.5, 1.} };
}









void MintGUI::clear_polyscope_state()
{

    std::cout << "clear_polyscope_state" << std::endl;
    polyscope::removeStructure("tet_mesh", false);
    polyscope::removeStructure("surf_mesh", false);

	for (int cur_id = 0; cur_id < 22; cur_id++)
	{
		polyscope::removeStructure(std::to_string(10000 + cur_id) + "____" + moment_labels.at(cur_id), false);
		polyscope::removeStructure(std::to_string(20000 + cur_id) + "____" + moment_labels.at(cur_id), false);
	}

    // single frame field
    polyscope::removeStructure("singularity(+1/4)", false);
    polyscope::removeStructure("singularity(-1/4)", false);
    polyscope::removeStructure("singularity(other)", false);
    polyscope::removeStructure("Centroids", false);

	std::cout << "end_clear_polyscope_state" << std::endl;


}









void MintGUI::set_base_mesh()
{
    if (!CubeCover::readMESH(path_mesh, V, T, F))
    {
        polyscope::warning("Unable to load selected mesh");
        mesh = CubeCover::TetMeshConnectivity();
    }
    else{
        mesh = CubeCover::TetMeshConnectivity(T);
        // make boundary mesh out of volume mesh
            // make a mesh out of all of the boundary faces

        auto parts = fileparts(path_mesh);
        
        Eigen::MatrixXd bdryV;
        Eigen::MatrixXi tmp;

        igl::readOBJ(parts.path + parts.name + ".obj", bdryV, bdryF);

        igl::bfs_orient(bdryF, bdryF, tmp);

        std::cout << "num flipped face ids " << tmp.size() << std::endl;
        
        // int nbdry = 0;
        // int nfaces = mesh.nFaces();
        // for (int i = 0; i < nfaces; i++)
        // {
        //     if (mesh.isBoundaryFace(i))
        //         nbdry++;
        // }

        // bdryF.resize(nbdry, 3);
        // // std::cout << "nbdry" << nbdry << std::endl<< std::endl<< std::endl;
        // // std::cout << "bdryF " << bdryF.size() << std::endl;
        // int curidx = 0;
        // for (int i = 0; i < nfaces; i++)
        // {
        //     if (mesh.isBoundaryFace(i))
        //     {
        //         for (int j = 0; j < 3; j++)
        //         {
        //             bdryF(curidx, j) = mesh.faceVertex(i, j);
                    
        //         }
        //         // fix triangle orientations
        //         int tet = mesh.faceTet(i, 0);
        //         if (tet == -1)
        //         {
        //             std::swap(bdryF(curidx, 0), bdryF(curidx, 1));
        //         }
        //         curidx++;
        //     }
        // }

        
        show_base_mesh();
    }


}

void MintGUI::load_state_from_output_dir()
{
    folder_contents.clear();
    folder_contents_fra.clear();
    file_names.clear();
    file_names_fra.clear();
    file_names_mesh.clear();
    adj_folder_names.clear();

    const int size = strlen(path_outdir);
    char tmp_path[512];

    int n = sprintf (tmp_path, "%.*s", size-1, path_outdir);

    // strcpy(tmp_path, path_outdir);

    // if (path_outdir[size - 1] == '/')

    // path_outdir[size - 1] = '\0';
    // memset(tmp_path, '\0', sizeof(tmp_path)-1);
    std::cout << "Load_State_From_Output_Dir: " << tmp_path << "size: " << size << std::endl;
    FileParts fp_tmp = fileparts(tmp_path);

    std::cout << "orig_path: " << path_outdir << std::endl;
    std::cout << "path: " << fp_tmp.path << std::endl << "name: " << fp_tmp.name << std::endl; 


    for (const auto & entry : fs::directory_iterator(fp_tmp.path))
    {
        std::cout << entry.path() << std::endl;
        // std::string tmp = entry.path();
        // FileParts fp = fileparts(tmp);
        // if (fp.ext == ".mom"){
        //     folder_contents.push_back( entry.path() );
        //     file_names.push_back( fp.name );
        // }

        // if (fp.ext == ".mesh")
        // {
        //     const char* tmp_path_mesh = tmp.c_str();
        //     std::cout << tmp_path_mesh << std::endl;
        //     path_mesh = new char[512];
        //     strncpy(path_mesh, tmp_path_mesh, 512);
        // }
        //     // folder_contents.push_back( entry.path() );
    }


    for (const auto  &entry : fs::directory_iterator(path_outdir))
    {
        // std::cout << entry.path() << std::endl;
		std::string tmp =  entry.path().u8string();
        FileParts fp = fileparts(tmp);
        if (fp.ext == ".mom"){
            folder_contents.push_back( entry.path().u8string());
            file_names.push_back( fp.name );
        }

        if (fp.ext == ".fra"){
            folder_contents_fra.push_back( entry.path().u8string());
            file_names_fra.push_back( fp.name );
        }

        if (fp.ext == ".mesh")
        {
            const char* tmp_path_mesh = tmp.c_str();
            std::cout << tmp_path_mesh << std::endl;
            path_mesh = new char[512];
            strncpy(path_mesh, tmp_path_mesh, 512);
            file_names_mesh.push_back( fp.name );
        }
            // folder_contents.push_back( entry.path() );
    }

    sel_idx_mom = -1;
    sel_idx_fra = -1;

    std::sort(folder_contents.begin(), folder_contents.end());
    std::sort(folder_contents_fra.begin(), folder_contents_fra.end());
    std::sort(file_names.begin(), file_names.end());
    std::sort(file_names_fra.begin(), file_names_fra.end());

    // for (int i = 0; i < folder_contents.size(); i++)
    // {
    //     std::cout << folder_contents.at(i) << std::endl;

    // }

    set_base_mesh();
    show_base_mesh();

}



void MintGUI::gui_file_explorer_callback()
{
     ImGui::SetNextItemOpen(true, ImGuiCond_Once);

    if (ImGui::TreeNode("Contents of chosen dir"))
    {
             
        if (ImGui::TreeNode("Show Frames and Curls"))
        {
             
            for (int i = 0; i < file_names_fra.size(); i++)
            {
                if (ImGui::Selectable(file_names_fra.at(i).c_str(), sel_idx_fra == i))
                {
                    sel_idx_fra = i;
                    path_fra = new char[512];
                    strncpy(path_fra, folder_contents_fra.at(i).c_str(), 512);
                    // CubeCover::readMoments(path_constraints, M_curr, true);
                    // std::cout << path_constraints << std::endl;
                    show_frame_field(frame_field_view_mode);
					// show_exploded_moments(moment_view_mode);

                }

            }
            ImGui::TreePop();
        }


        ImGui::SetNextItemOpen(true, ImGuiCond_Once);

        if (ImGui::TreeNode("TODO: Select Mesh from directory"))
        {

            ImGui::TreePop();
        }
        
        ImGui::SetNextItemOpen(true, ImGuiCond_Once);

        if (ImGui::TreeNode("Select Moments to Visualize"))
        {
             

            for (int i = 0; i < file_names.size(); i++)
            {
                if (ImGui::Selectable(file_names.at(i).c_str(), sel_idx_mom == i))
                {
                    sel_idx_mom = i;
                    path_constraints = new char[512];
                    strncpy(path_constraints, folder_contents.at(i).c_str(), 512);
                    CubeCover::readMoments(path_constraints, M_curr, true);
                    // std::cout << path_constraints << std::endl;
					show_exploded_moments(moment_view_mode);

                }

            }




            ImGui::PushItemWidth(300);
            ImGui::InputTextWithHint("3", "path to exploded moments", path_constraints, 512);
            ImGui::PopItemWidth();
            ImGui::SameLine();
            HelpMarker("Choose moments to visualize in a different directory");
            ImGui::SameLine();

            if (ImGui::Button("Pick .mom")) {
                char* tmp_path_constraints = fileSelectSubroutine();
                path_constraints = new char[512];
                strncpy(path_constraints, tmp_path_constraints, 512);
				show_exploded_moments(moment_view_mode);

                sel_idx_mom = -1;

            }

            ImGui::TreePop();
        }


        ImGui::TreePop();

    }



    ///////////////////////////////////////////////////////////////////////////
}


void MintGUI::gui_run_mint_callback()
{
    if (ImGui::Button("Run Mint")) {
    // executes when button is pressed
    // directory_path = fileSelectSubroutine();
        polyscope::warning("NOT IMPLEMENTED");
    }

    ImGui::SameLine();
    HelpMarker("TODO: start mint call in seperate thread");
    ImGui::SameLine();




    if (ImGui::TreeNode("Exploded Solver options"))
    {
        ImGui::Text("Integrability Mode");
        if (ImGui::RadioButton("Free", mint_mode == Mint_Integrability_Mode::free))  
        { 
            mint_mode = Mint_Integrability_Mode::free;
        } ImGui::SameLine();
        if (ImGui::RadioButton("Mint (the main attraction!)", mint_mode == Mint_Integrability_Mode::mint))  
        { 
            mint_mode = Mint_Integrability_Mode::mint;
        } 


        ImGui::Text("Shell (local projection) Switch");
        if (ImGui::RadioButton("Off-Shell", shell_mode == Mint_Frame_Projection::offshell))  
        { 
            shell_mode = Mint_Frame_Projection::offshell;
        } ImGui::SameLine();
        if (ImGui::RadioButton("On-Shell", shell_mode == Mint_Frame_Projection::onshell))  
        { 
            shell_mode = Mint_Frame_Projection::onshell;
        } 


        ImGui::Text("Delta Metric");
        if (ImGui::RadioButton("Fourth Moments", metric_mode == Mint_Moment_Metric::four))  
        { 
            metric_mode = Mint_Moment_Metric::four;
        } ImGui::SameLine();
        if (ImGui::RadioButton("Second Moments", metric_mode == Mint_Moment_Metric::sec))  
        { 
            metric_mode = Mint_Moment_Metric::sec;
        } 
        if (ImGui::RadioButton("Scale Invariant: 4 + 2 \\oplus 2", metric_mode == Mint_Moment_Metric::four_plus_two_tensor_two))  
        { 
            metric_mode = Mint_Moment_Metric::four_plus_two_tensor_two;
        } ImGui::SameLine();
        if (ImGui::RadioButton("Scale Dependent: 4 + 2", metric_mode == Mint_Moment_Metric::four_plus_two))  
        { 
            metric_mode = Mint_Moment_Metric::four_plus_two;
        } 


        ImGui::Text(" \\ solver good, gmres for big models ");
    if (ImGui::RadioButton("Exact", cur_solver == Mint_Linear_Solver::exact))  { cur_solver = Mint_Linear_Solver::exact; } ImGui::SameLine();
    if (ImGui::RadioButton("GMRes", cur_solver == Mint_Linear_Solver::gmres))  { cur_solver = Mint_Linear_Solver::gmres; } 
    

	if (ImGui::Button("(re)Load Boundary")) {

		polyscope::warning("The chosen .bound file does not match the loaded mesh.");
		// executes when button is pressed
		// directory_path = fileSelectSubroutine();
		show_exploded_moments(moment_view_mode);
	}

	ImGui::SameLine();

	ImGui::SetNextItemOpen(true, ImGuiCond_Once);


	if (ImGui::TreeNode("Exploded Boundary Select options"))
	{

		// ImGui::Text("Start penzil.app to annotate model");

		// ImGui::Text("Load From File");

		if (ImGui::Button("Start penzil.app to annotate model")) {

	//		show_exploded_moments(moment_view_mode);

		}
		if (ImGui::Button("Load From File")) {

	//		show_exploded_moments(moment_view_mode);
		}

        if (ImGui::Button("Project Bound to closest GL(3) field")) {
    // executes when button is pressed
    // directory_path = fileSelectSubroutine();
    
        }




		ImGui::TreePop();

	}




	if (ImGui::Button("Compute and Load normal boundary")) {

		show_exploded_moments(moment_view_mode);

	}


        ImGui::TreePop();

    }

}

void MintGUI::gui_file_select_BROKEN()
{
    //    int rightWindowsWidth = 1000;
    	// ImGui::PushItemWidth(1000); // Make ui elements 100 pixels wide,
							   // instead of full width. Must have 
							   // matching PopItemWidth() below.

//	ImGui::InputInt("num points", &nPts);             // set a int variable
//	ImGui::InputFloat("param value", &anotherParam);  // set a float variable

//	if (ImGui::Button("run subroutine")) {
		// executes when button is pressed
//		mySubroutine();
//	}
    // std::string blah;
    // std::string blah2;
    // ImGui::InputText("Mesh Path", blah);
    // ImGui::InputText("Boundary Constraints Path", blah2);



    // std::cout << rightWindowsWidth << std::endl;


    ///////////////////////////////////////////////////////////////////////////

    ImGui::PushItemWidth(300);
    ImGui::InputTextWithHint("1", "abs mesh path, or use picker button.", path_mesh, 512);   
    ImGui::PopItemWidth();
    ImGui::SameLine();
    HelpMarker("Mesh Path");
    ImGui::SameLine();


    


    if (ImGui::Button("Pick .mesh")) {
    // executes when button is pressed
        char* tmp_path_mesh = fileSelectSubroutine();
        path_mesh = new char[512];
        strncpy(path_mesh, tmp_path_mesh, 512);
        auto cur_path_parts = fileparts(path_mesh);
        std::cout << cur_path_parts.ext << std::endl;

        std::string data = cur_path_parts.ext;
        std::transform(data.begin(), data.end(), data.begin(),
            [](unsigned char c){ return std::tolower(c); });


        std::cout << data << std::endl;
        if (data == ".mesh" )
        {
            set_base_mesh();
        }
        else 
        {
            polyscope::warning("Please pick a .mesh file to load.  Support for .obj coming eventually...");
        }
    }

    if (ImGui::Button("(re)load view mesh")) {
        auto cur_path_parts = fileparts(path_mesh);
        std::cout << cur_path_parts.ext << std::endl;

        std::string data = cur_path_parts.ext;
        std::transform(data.begin(), data.end(), data.begin(),
            [](unsigned char c){ return std::tolower(c); });
        if (data == ".mesh" )
        {
            set_base_mesh();
        }
        else 
        {
            polyscope::warning("Please pick a .mesh file to load.  Support for .obj coming eventually...");
        }
    }


    ImGui::SameLine();

    if (ImGui::Button("Run Full Pipeline")) {
        // auto cur_path_parts = fileparts(path_mesh);
        // std::cout << cur_path_parts.ext << std::endl;

        // std::string data = cur_path_parts.ext;
        // std::transform(data.begin(), data.end(), data.begin(),
        //     [](unsigned char c){ return std::tolower(c); });
        // if (data == ".mesh" )
        // {
        //     set_base_mesh();
        // }
        // else 
        // {
            polyscope::warning("TODO: Not implemented");
        // }
    }



    ///////////////////////////////////////////////////////////////////////////

        ImGui::PushItemWidth(300);
    ImGui::InputTextWithHint("2", "enter rel or abs path, or use picker.", path_outdir, 512);
    ImGui::PopItemWidth();


    ImGui::SameLine();
    HelpMarker("Specify directory to load a previous run, or to choose output for next run");
    ImGui::SameLine();


    if (ImGui::Button("Pick mint output dir")) {
        char* tmp_path_outdir = fileSelectSubroutine();
        strncpy(path_outdir, tmp_path_outdir, 512);
    }

 
    
    if (ImGui::Button("(re)load or create directory")) {

      polyscope::warning("This directory did not exist, creating it");
    // executes when button is pressed
    // directory_path = fileSelectSubroutine();


        //     char* tmp_path_mesh = fileSelectSubroutine();
        // strncpy(path_mesh, tmp_path_mesh, 512);
        // auto cur_path_parts = fileparts(path_mesh);
        // std::cout << cur_path_parts.ext << std::endl;

        // std::string data = cur_path_parts.ext;
        // std::transform(data.begin(), data.end(), data.begin(),
        //     [](unsigned char c){ return std::tolower(c); });


        // std::cout << data << std::endl;
        // if (data == ".mesh" )
        // {
        //     set_base_mesh();
        // }
        // else 
        // {
        //     polyscope::warning("Please pick a .mesh file to load.  Support for .obj coming eventually...");
        // }




    }
}

void MintGUI::gui_folder_explorer_callback()
{

}


void MintGUI::gui_main_control_panel_callback()
{
        ImGui::SetNextItemOpen(true, ImGuiCond_Once);

    if (ImGui::TreeNode("Exploded GUI options"))
    {
        if (ImGui::Checkbox("Use consistent heatmap", &useSameColorRangeForAllMoments)) { show_exploded_moments(moment_view_mode); }ImGui::SameLine();
        ImGui::PushItemWidth(100);
        ImGui::SliderFloat("m", &color_range_min, 0.0f, 1.0f, "min = %.3f");
        ImGui::PopItemWidth();ImGui::SameLine();
        ImGui::PushItemWidth(100);
        ImGui::SliderFloat("M", &color_range_max, 0.0f, 1.0f, "max = %.3f");
        ImGui::PopItemWidth();
 
 

        ImGui::PushItemWidth(100);
        
        ImGui::SliderFloat("Exploded Spacing", &exploded_spacing_inp, 0.0f, 150.0f, "ratio = %.3f");
        ImGui::PopItemWidth();

        ImGui::Text("Which moments to show");
        if (ImGui::RadioButton("4th", moment_view_mode == Moments_To_Show::fourth))  
        { 
            moment_view_mode = Moments_To_Show::fourth;
			show_exploded_moments(moment_view_mode);


            // polyscope::state::lengthScale = 5.;
        } ImGui::SameLine();
        if (ImGui::RadioButton("2nd", moment_view_mode == Moments_To_Show::second))  
		{ 
			moment_view_mode = Moments_To_Show::second; 
			show_exploded_moments(moment_view_mode);
		} ImGui::SameLine();
        if (ImGui::RadioButton("both", moment_view_mode == Moments_To_Show::both))  
		{ 
			moment_view_mode = Moments_To_Show::both; 
			show_exploded_moments(moment_view_mode);
		} 


        ImGui::Text("Which frames to show");
        if (ImGui::RadioButton("GL(3)", frame_field_view_mode == Frames_To_Show::frames))  
        { 
            frame_field_view_mode = Frames_To_Show::frames;
			show_frame_field(frame_field_view_mode);


            // polyscope::state::lengthScale = 5.;
        } ImGui::SameLine();
        if (ImGui::RadioButton("split frames", frame_field_view_mode == Frames_To_Show::split_frames))  
        { 
            frame_field_view_mode = Frames_To_Show::split_frames;
			show_frame_field(frame_field_view_mode);


            // polyscope::state::lengthScale = 5.;
        } ImGui::SameLine();
        if (ImGui::RadioButton("split_moments", frame_field_view_mode == Frames_To_Show::split_moments))  
        { 
            frame_field_view_mode = Frames_To_Show::split_moments;
			show_frame_field(frame_field_view_mode);


            // polyscope::state::lengthScale = 5.;
        } ImGui::SameLine();
        if (ImGui::RadioButton("split_difference", frame_field_view_mode == Frames_To_Show::split_difference))  
        { 
            frame_field_view_mode = Frames_To_Show::split_difference;
			show_frame_field(frame_field_view_mode);


            // polyscope::state::lengthScale = 5.;
        } 
        
        // if (ImGui::RadioButton("2nd", moment_view_mode == Moments_To_Show::second))  
		// { 
		// 	moment_view_mode = Moments_To_Show::second; 
		// 	show_exploded_moments(moment_view_mode);
		// } ImGui::SameLine();
        // if (ImGui::RadioButton("both", moment_view_mode == Moments_To_Show::both))  
		// { 
		// 	moment_view_mode = Moments_To_Show::both; 
		// 	show_exploded_moments(moment_view_mode);
		// } 




		if (ImGui::Checkbox("show boundary", &showBoundary)) { show_exploded_moments(moment_view_mode); }ImGui::SameLine();
        if (ImGui::Checkbox("show interior tets", &showInteriorTets)) { show_exploded_moments(moment_view_mode); }


        ImGui::TreePop();

    }
}



void MintGUI::gui_callback()
{
    ImGui::PushID("dummy");
    ImGui::Begin("d", nullptr, ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoSavedSettings);
        
  
    
    
    ImGui::End();
    ImGui::PopID();


    ImGui::PushID("folder_browser");
    
    ImGui::SetNextWindowPos(ImVec2(400, 20),ImGuiCond_Once);
    ImGui::SetNextWindowSize(ImVec2(1000, 0));

    ImGui::Begin("Folder Browser", nullptr, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings);

    // ImGui::SetNextWindowPos(ImVec2((400 ), 50), ImGuiCond_Once);
    // ImGui::SetNextWindowSize(ImVec2(1000, 0.));
    // ImGuiWindowFlags_NoMove

    //     // Main body of the Demo window starts here.
    // if (!ImGui::Begin("File Browser", nullptr, ImGuiWindowFlags_AlwaysAutoResize))
    // {
    //     // Early out if the window is collapsed, as an optimization.
    //     ImGui::End();
    //     return;
    // }
    gui_file_select_BROKEN();
  gui_folder_explorer_callback();



    ImGui::End();
    ImGui::PopID();





    ImGui::PushID("file_browser");
    
    ImGui::SetNextWindowPos(ImVec2(polyscope::view::windowWidth -(600 ), 600),ImGuiCond_Once);
    ImGui::SetNextWindowSize(ImVec2(1000, polyscope::view::windowHeight/2));

    ImGui::Begin("File Browser", nullptr, ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoSavedSettings);

    // ImGui::SetNextWindowPos(ImVec2((400 ), 50), ImGuiCond_Once);
    // ImGui::SetNextWindowSize(ImVec2(1000, 0.));
// ImGuiWindowFlags_NoMove

    //     // Main body of the Demo window starts here.
    // if (!ImGui::Begin("File Browser", nullptr, ImGuiWindowFlags_AlwaysAutoResize))
    // {
    //     // Early out if the window is collapsed, as an optimization.
    //     ImGui::End();
    //     return;
    // }
  gui_file_explorer_callback();



    ImGui::End();
    ImGui::PopID();




    ImGui::PushID("mint_user_callback");
    
    ImGui::SetNextWindowPos(ImVec2(polyscope::view::windowWidth - (500 ), 200), ImGuiCond_Once);
    ImGui::SetNextWindowSize(ImVec2(1000, 0.));

        // Main body of the Demo window starts here.
    // if (!ImGui::Begin("Mint Control Panel", nullptr, ImGuiWindowFlags_AlwaysAutoResize))
    // {
    //     // Early out if the window is collapsed, as an optimization.
    //     ImGui::End();
    //     return;
    // }

    ImGui::Begin("Mint Control Panel", nullptr, ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoSavedSettings);




  
   



    
  
gui_main_control_panel_callback();


gui_run_mint_callback();



	// ImGui::PopItemWidth();

    auto winsizehack = ImGui::GetWindowSize();


    ImGui::End();
    ImGui::PopID();

    ImGui::ShowDemoWindow(); 

    ImGui::SetWindowSize( winsizehack ); 
}

}









/*


    static ImGuiTextFilter filter;
        ImGui::Text("Filter usage:\n"
                    "  \"\"         display all lines\n"
                    "  \"xxx\"      display lines containing \"xxx\"\n"
                    "  \"xxx,yyy\"  display lines containing \"xxx\" or \"yyy\"\n"
                    "  \"-xxx\"     hide lines containing \"xxx\"");
        filter.Draw();


        
        // memset(lines, 0, sizeof(lines));
        
        // for (int i = 0; i < IM_ARRAYSIZE(lines); i++)
        //     if (filter.PassFilter(lines[i]))
        //         ImGui::BulletText("%s", lines[i]);




  if (ImGui::Button("run subroutine")) {
    // executes when button is pressed
    // directory_path = fileSelectSubroutine();
  }

	ImGui::SameLine();
	if (ImGui::Button("hi")) {
		polyscope::warning("hi");
	}
    */

