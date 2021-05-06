#include "precomp.h"
#include "prickliness.h"

#include "DrawPrickliness.h"

namespace Prickliness
{
    //Get a vector to a point inside a convex DCEL plane
    Vector_3 PricklinessArrangement::GetVectorInPlane(const Arrangement_2::Face_const_handle face) const
    {
        std::vector<Point_2> edgePoints;

        if (face->has_outer_ccb())
        {
            Arrangement_2::Ccb_halfedge_const_circulator currentHalfEdge = face->outer_ccb();

            //Select a point on 2 of the (non-fictitious) edges to calculate a point inside the face
            do
            {
                //Fictitious: https://doc.cgal.org/latest/Arrangement_on_surface_2/index.html#title32
                if (!currentHalfEdge->is_fictitious())
                {
                    if (currentHalfEdge->curve().is_ray())
                    {
                        //Edge is an infinite ray originating from a vertex on an unbounded face, push an arbitrary point on the ray
                        auto rayStart = currentHalfEdge->curve().ray().start();
                        auto rayVector = currentHalfEdge->curve().ray().direction().to_vector();

                        edgePoints.push_back(rayStart + (rayVector * Kernel::FT(0.5)));
                    }
                    else if (currentHalfEdge->curve().is_segment())
                    {
                        //Edge is a line segment bounded by two vertices, push the point in the middle of the line segment
                        auto heSegment = currentHalfEdge->curve().segment();
                        auto segmentVector = heSegment.end() - heSegment.start();

                        edgePoints.push_back(heSegment.start() + (segmentVector * Kernel::FT(0.5)));
                    }
                }

            } while (++currentHalfEdge != face->outer_ccb() && edgePoints.size() < 2);
        }

        //Get the point in between the two edge points and convert it to 3D
        //the arrangement plane is at z = 10
        Point_2 startFacePoint = edgePoints.at(0) + ((edgePoints.at(1) - edgePoints.at(0)) * Kernel::FT(0.5));
        Vector_3 faceVector(startFacePoint.x(), startFacePoint.y(), 10.0);

        return faceVector;
    }

    bool PricklinessArrangement::IsBoundaryVert(const Delaunay& DT, Delaunay::Face_circulator adjacentFaces) const
    {
        Delaunay::Face_circulator first = adjacentFaces;

        do
        {
            if (DT.is_infinite(adjacentFaces))
            {
                return true;
            }

        } while (++adjacentFaces != first);

        return false;
    }

    PricklinessArrangement::PricklinessArrangement(const Delaunay& DT)
    {
        Plane_3 arrangementPlane(Point_3(0.0, 0.0, 10.0), Vector_3(0.0, 0.0, 1.0));

        std::vector<Arr_X_monotone_curve_2> arrangementLines;
        std::vector<std::vector<Plane_3>> perpPlanes;
        std::vector<int> vertexDegrees;
        int currentVertIndex = 0;

        std::cout << "\t Constructing perpendicular planes" << std::endl;


        size_t DTVertCount = DT.number_of_vertices();
        //size_t percent = DTVertCount / 100;

        //Loop through all the finite vertices of the DT, the underlying structures adds a point at infinity to connect all the boundary verts together

        for (Del_finite_vertices_iterator vertIterator = DT.finite_vertices_begin(); vertIterator != DT.finite_vertices_end(); vertIterator++)
        {
            if (IsBoundaryVert(DT, vertIterator->incident_faces()))
            {
                continue;
            }


            //Add a perp plane vector for this vertex
            perpPlanes.push_back(std::vector<Plane_3>());

            //Store current vertex index for the result plane traversal
            PerpPlaneIndices vnpTriple;
            vnpTriple.vertexIndex = currentVertIndex;

            //Get next vertex and its neighbours iterator
            Point_3 currentVertex = vertIterator->point();
            Del_vertex_circulator neighbourCirculator = vertIterator->incident_vertices(), done(neighbourCirculator);


            int neighbourIndex = 0;
            if (neighbourCirculator != 0)
            {
                do
                {
                    //Skip the vertex at infinity that the triangulation DS uses to connect the boundary verts
                    if (!DT.is_infinite(neighbourCirculator))
                    {
                        //Get neighbour vertex and subtract current current vertex to move it towards the origin
                        Point_3 neighbour = neighbourCirculator->point();
                        Vector_3 neighbourToOrigin = neighbour - currentVertex;

                        //Construct a line going through the neighbour and current vector at the origin
                        //Flip the direction so it goes towards the current vector
                        //(it's important to make sure it's oriented this way because we need to check if we are on the local maxima side of the plane with a dot on the normal)
                        Line_3 lineThroughOrigin = Line_3(CGAL::ORIGIN, neighbourToOrigin).opposite();

                        //Store the index of the neighbour so we can use it to get it's perp plane
                        vnpTriple.planeIndex = perpPlanes.at(currentVertIndex).size();

                        //Construct a plane perpendiculat to the neighbour line at the origin
                        Plane_3 perpPlane = lineThroughOrigin.perpendicular_plane(CGAL::ORIGIN);
                        perpPlanes.at(currentVertIndex).push_back(perpPlane);

                        //Intersect the perp plane with the arrangement plane, store the resulting 2D line so we can mass add them in the arrangement
                        auto intersection = CGAL::intersection(arrangementPlane, perpPlane);

                        if (intersection)
                        {
                            if (Line_3* arrangementLine = boost::get<Line_3>(&*intersection))
                            {
                                Point_3 p3 = arrangementLine->point();
                                Vector_3 v3 = arrangementLine->to_vector();

                                Arr_data_line_2 arrangementLine2(Point_2(p3.x(), p3.y()), Vector_2(v3.x(), v3.y()));

                                arrangementLines.push_back(Arr_X_monotone_curve_2(arrangementLine2, vnpTriple));

                                //CGAL::insert(resultArrangement, Arr_X_monotone_curve_2(arrangementLine2, vnpTriple));
                            }
                            else
                            {
                                //This will never happen because the perp planes go through the origin and the arrangement plane lies above it.
                                std::cout << "Result of intersection between perpendicular plane and arrangement plane is a plane!" << std::endl;
                            }
                        }
                        else
                        {
                            //Parallel planes degenerate case: Ignore? Should never happen unless two vertices are right on top of each other
                            std::cout << "Empty intersection between perpendicular plane and arrangement plane!" << std::endl;
                        }

                        neighbourIndex++;
                    }
                } while (++neighbourCirculator != done);

                vertexDegrees.push_back(neighbourIndex); //Save degree this way because vert->degree() counts fictitious neighbours..
            }

            currentVertIndex++;

            //if (currentVertIndex % percent == 0)
            //{
            //    std::cout << currentVertIndex / percent << std::endl;
            //}

        }

        std::cout << "\t Constructing 2d arrangement (DCEL)" << std::endl;
        //Construct an arrangement (DCEL) that keeps track of the original lines the half-edges are created from
        //Arrangement_2 resultArrangement;
        CGAL::insert(resultArrangement, arrangementLines.begin(), arrangementLines.end());
        arrangementLines.clear();
        //FIFO queues for breadth first traversal of the arrangement
        std::queue<Arrangement_2::Face_handle> faceQueue;
        std::queue<std::vector<int>> abovePlaneCountsQueue;

        std::cout << "\t Determining local maxima for DCEL faces" << std::endl;
        if (resultArrangement.number_of_faces() > 1)
        {
            //Calculate the local maxima for the (random) starting face
            auto startingFace = resultArrangement.faces_begin();

            //Get a vector to a point in the face so we can check if we are above the perp planes
            Vector_3 startFaceVector = GetVectorInPlane(startingFace);

            //Construct the data structure keeping track of planes and the local maxima
            std::vector<int> abovePlanesCount;
            abovePlanesCount.reserve(perpPlanes.size());

            for (std::vector<Plane_3>& planeVector : perpPlanes)
            {
                int aboveCount = 0;
                for (Plane_3& plane : planeVector)
                {
                    //Count as above if perpendicular plane normal dot face > 0
                    if (CGAL::scalar_product(plane.orthogonal_vector(), startFaceVector) > Kernel::FT(0))
                    {
                        aboveCount++;
                    }
                }
                abovePlanesCount.push_back(aboveCount);
            }

            //Count for how many terrain vertices we are in the local maxima cone, this is the # local maxima for this arrangement face
            //We are in the cone if we are above all of its planes
            int localMaximaCount = 0;
            for (size_t i = 0; i < abovePlanesCount.size(); i++)
            {
                if (abovePlanesCount.at(i) == vertexDegrees.at(i))
                {
                    localMaximaCount++;
                }
            }

            startingFace->data().localMaximaCount = localMaximaCount;

            //Put the starting face and the counters for the planes on the queues for the breadth first iteration
            faceQueue.push(startingFace);
            abovePlaneCountsQueue.push(abovePlanesCount);

            //Breadth first search over the terrain:
            //Everytime we cross an edge check its corresponding perp plane(s) with a vector to the point in the next face:
            //if we are above a plane increment its corresponding above counter, decrement otherwise.
            //if the above counter is equal to the degree of the vertex after an increment we are within its local maxima cone and increment # local maxima
            //if the above counter is equal to the degree of the vertex before an decrement we have left the local maxima cone and decrement # local maxima
            while (!faceQueue.empty())
            {
                auto currentFace = faceQueue.front();
                faceQueue.pop();

                std::vector<int> currentAbovePlaneCounts = abovePlaneCountsQueue.front();
                abovePlaneCountsQueue.pop();


                //Iterate through the incident edges
                Arrangement_2::Ccb_halfedge_circulator currentHalfEdge = currentFace->outer_ccb();
                do
                {
                    //Fictitious: https://doc.cgal.org/latest/Arrangement_on_surface_2/index.html#title32
                    if (!currentHalfEdge->is_fictitious())
                    {
                        //Skip neighbour face if already visited
                        if (currentHalfEdge->twin()->face()->data().localMaximaCount != -1)
                        {
                            continue;
                        }

                        std::vector<int> neighbourAbovePlanesCounts = currentAbovePlaneCounts;
                        int localMaximaCount = currentFace->data().localMaximaCount;


                        Vector_3 faceVector = GetVectorInPlane(currentHalfEdge->twin()->face());

                        //This is a for loop because the DCEL stores the data in a list in case of the degenerate case where multiple lines overlap
                        for (PerpPlaneIndices ppIndices : currentHalfEdge->curve().data())
                        {
                            Plane_3 perpPlane = perpPlanes.at(ppIndices.vertexIndex).at(ppIndices.planeIndex);

                            //Which side of the perp plane did we go to?
                            if (CGAL::scalar_product(perpPlane.orthogonal_vector(), faceVector) > Kernel::FT(0))
                            {
                                //Above, increment above counter and check if inside local maxima cone
                                neighbourAbovePlanesCounts.at(ppIndices.vertexIndex)++;

                                if (neighbourAbovePlanesCounts.at(ppIndices.vertexIndex) == vertexDegrees.at(ppIndices.vertexIndex))
                                {
                                    localMaximaCount++;
                                }
                            }
                            else
                            {
                                //Below, check if we were inside the local maxima cone and decrement above counter
                                if (neighbourAbovePlanesCounts.at(ppIndices.vertexIndex) == vertexDegrees.at(ppIndices.vertexIndex))
                                {
                                    localMaximaCount--;
                                }
                                neighbourAbovePlanesCounts.at(ppIndices.vertexIndex)--;
                            }
                        }

                        currentHalfEdge->twin()->face()->data().localMaximaCount = localMaximaCount;

                        //Local maxima count for this face current most?
                        if (prickliness < localMaximaCount)
                        {
                            prickliness = localMaximaCount;
                            pricklyDirection = faceVector;
                        }

                        faceQueue.push(currentHalfEdge->twin()->face());
                        abovePlaneCountsQueue.push(neighbourAbovePlanesCounts);
                    }
                } while (++currentHalfEdge != currentFace->outer_ccb());
            }
        }

        //CGAL::draw(resultArrangement);
    }

    void PricklinessArrangement::Draw()
    {
        Arrangement_2::Vertex_const_iterator vit = resultArrangement.vertices_begin();

        Kernel::FT minX = vit->point().x();
        Kernel::FT minY = vit->point().y();
        Kernel::FT maxX = vit->point().x();
        Kernel::FT maxY = vit->point().y();

        for (; vit != resultArrangement.vertices_end(); ++vit) {

            if (vit->point().x() < minX) { minX = vit->point().x(); }
            if (vit->point().y() < minY) { minY = vit->point().y(); }
            if (vit->point().x() > maxX) { maxX = vit->point().x(); }
            if (vit->point().y() > maxY) { maxY = vit->point().y(); }
        }

        Point_2 minPoint(minX, minY);
        Point_2 maxPoint(maxX, maxY);


        minPoint -= Vector_2(10, 10);
        maxPoint += Vector_2(10, 10);
        std::cout << "Min:" << minPoint << " Max:" << maxPoint << std::endl;

        Point_2 topLeft(minPoint.x(), maxPoint.y());
        Point_2 bottomRight(maxPoint.x(), minPoint.y());

        std::vector<Arr_traits_2::Curve_2> boundarybox;
        boundarybox.push_back(Arr_traits_2::Curve_2(topLeft, minPoint));
        boundarybox.push_back(Arr_traits_2::Curve_2(minPoint, bottomRight));
        boundarybox.push_back(Arr_traits_2::Curve_2(bottomRight, maxPoint));
        boundarybox.push_back(Arr_traits_2::Curve_2(maxPoint, topLeft));

        CGAL::insert(resultArrangement, boundarybox.begin(), boundarybox.end());

        std::cout << "Drawing Prickliness (" << prickliness << "): " << resultArrangement.number_of_faces() << " Faces, " << resultArrangement.number_of_edges() << " Edges, " << resultArrangement.number_of_vertices() << " Vertices." << std::endl;

        DrawPrickliness(resultArrangement, prickliness, pricklyDirection);
    }

    void PricklinessArrangement::WriteToIpe(const std::filesystem::path& filePath, bool centerOnOrigin)
    {
        //Merge faces
        Arrangement_2::Edge_iterator edgeIt;
        std::vector<Arrangement_2::Halfedge_handle> redundantEdges;
        std::vector<int> pricklyValues;
        for (edgeIt = resultArrangement.edges_begin(); edgeIt != resultArrangement.edges_end(); ++edgeIt)
        {
            Arrangement_2::Halfedge_handle hEdge = edgeIt;

            if (hEdge->face()->data().localMaximaCount == hEdge->twin()->face()->data().localMaximaCount)
            {
                int prickly = hEdge->face()->data().localMaximaCount;
                redundantEdges.push_back(hEdge);
                pricklyValues.push_back(hEdge->face()->data().localMaximaCount);
            }
        }

        int currIndex = 0;
        for (Arrangement_2::Halfedge_handle hEdge : redundantEdges)
        {
            Arrangement_2::Face_handle fh = resultArrangement.remove_edge(hEdge, true, true);
            fh->data().localMaximaCount = pricklyValues.at(currIndex);
            ++currIndex;
        }

        //Remove redundant vertices
        std::vector<Arrangement_2::Vertex_handle> vertexHandles;
        for (Arrangement_2::Vertex_handle vertexHandle : resultArrangement.vertex_handles())
        {
            vertexHandles.push_back(vertexHandle);
        }
        for (Arrangement_2::Vertex_handle vertexHandle : vertexHandles)
        {
            CGAL::remove_vertex(resultArrangement, vertexHandle);
        }
        vertexHandles.clear();

        //Add bounding box
#if 0  //Boundary box variable size
        Arrangement_2::Vertex_const_iterator vit = resultArrangement.vertices_begin();

        Kernel::FT minX = vit->point().x();
        Kernel::FT minY = vit->point().y();
        Kernel::FT maxX = vit->point().x();
        Kernel::FT maxY = vit->point().y();

        for (; vit != resultArrangement.vertices_end(); ++vit) {

            if (vit->point().x() < minX) { minX = vit->point().x(); }
            if (vit->point().y() < minY) { minY = vit->point().y(); }
            if (vit->point().x() > maxX) { maxX = vit->point().x(); }
            if (vit->point().y() > maxY) { maxY = vit->point().y(); }
        }

        Point_2 minPoint(minX, minY);
        Point_2 maxPoint(maxX, maxY);


        minPoint -= Vector_2(1, 1);
        maxPoint += Vector_2(1, 1);
#else //Static size

        Point_2 minPoint;
        Point_2 maxPoint;

        if (centerOnOrigin)
        {
            //square centered on prickly direction
            minPoint = CGAL::ORIGIN;
            maxPoint = CGAL::ORIGIN;
        }
        else
        {
            //square centered on prickly direction
            minPoint = Point_2(pricklyDirection.x(), pricklyDirection.y());
            maxPoint = Point_2(pricklyDirection.x(), pricklyDirection.y());
        }

        minPoint -= Vector_2(4, 4);
        maxPoint += Vector_2(4, 4);

#endif // boundary box

        std::cout << "Min:" << minPoint << " Max:" << maxPoint << std::endl;

        Point_2 topLeft(minPoint.x(), maxPoint.y());
        Point_2 bottomRight(maxPoint.x(), minPoint.y());

        std::vector<Arr_traits_2::Curve_2> boundarybox;
        boundarybox.push_back(Arr_traits_2::Curve_2(topLeft, minPoint));
        boundarybox.push_back(Arr_traits_2::Curve_2(minPoint, bottomRight));
        boundarybox.push_back(Arr_traits_2::Curve_2(bottomRight, maxPoint));
        boundarybox.push_back(Arr_traits_2::Curve_2(maxPoint, topLeft));

        CGAL::insert(resultArrangement, boundarybox.begin(), boundarybox.end());

        const std::string ipeFileStart = "<?xml version=\"1.0\"?>\n"
            "<!DOCTYPE ipe SYSTEM \"ipe.dtd\">\n"
            "<ipe version = \"70206\" creator=\"Ipe 7.2.7\">\n"
            "<page>\n"
            "<layer name = \"alpha\"/>\n"
            "<view layers = \"alpha\" active=\"alpha\"/>\n";

        const std::string ipeFileEnd = "</page>\n</ipe>";

        Arrangement_2::Face_const_iterator fit;

        //Find lowest and highest values in current view (This double loop stuff is kind of ugly but it's only for some basic i/o so whatever)
        int currentViewMaximaMax = 0;
        int currentViewMaximaMin = prickliness;

        for (fit = resultArrangement.faces_begin(); fit != resultArrangement.faces_end(); ++fit)
        {
            //Only export finite faces
            if (!fit->is_unbounded() && !fit->is_fictitious() && fit->has_outer_ccb())
            {
                if (fit->data().localMaximaCount != -1)
                {
                    auto edgeCirculator = fit->outer_ccb();
                    auto startEdgeCirc = fit->outer_ccb();

                    bool inside = true;

                    do
                    {
                        if (edgeCirculator->source()->point().x() < minPoint.x() ||
                            edgeCirculator->source()->point().x() > maxPoint.x() ||
                            edgeCirculator->source()->point().y() < minPoint.y() ||
                            edgeCirculator->source()->point().y() > maxPoint.y())
                        {
                            inside = false;
                            break;
                        }

                    } while (++edgeCirculator != startEdgeCirc);

                    //Skip faces outside of bounding box
                    if (!inside)
                    {
                        continue;
                    }

                    if (currentViewMaximaMax < fit->data().localMaximaCount)
                    {
                        currentViewMaximaMax = fit->data().localMaximaCount;
                    }

                    if (currentViewMaximaMin > fit->data().localMaximaCount)
                    {
                        currentViewMaximaMin = fit->data().localMaximaCount;
                    }
                }
            }
        }

        std::cout << "Min prickly in view: " << currentViewMaximaMin << std::endl;
        std::cout << "Max prickly in view: " << currentViewMaximaMax << std::endl;

        std::ostringstream fileName;
        fileName << "prickliness_";

        if (centerOnOrigin)
        {
            fileName << "O";
        }
        else
        {
            fileName << "P";
        }

        fileName << "_" << currentViewMaximaMin << "_to_" << currentViewMaximaMax << ".ipe";

        std::ofstream ipeFileStream(filePath / fileName.str());

        if (!ipeFileStream.is_open())
        {
            std::cout << "\t Ipe file path invalid or not openable.." << std::endl;
            return;
        }

        ipeFileStream << ipeFileStart;


        //Export the arrangement faces.
        for (fit = resultArrangement.faces_begin(); fit != resultArrangement.faces_end(); ++fit)
        {
            //Only export finite faces
            if (!fit->is_unbounded() && !fit->is_fictitious() && fit->has_outer_ccb())
            {

                auto edgeCirculator = fit->outer_ccb();
                auto startEdgeCirc = fit->outer_ccb();

                bool inside = true;

                do
                {
                    if (edgeCirculator->source()->point().x() < minPoint.x() ||
                        edgeCirculator->source()->point().x() > maxPoint.x() ||
                        edgeCirculator->source()->point().y() < minPoint.y() ||
                        edgeCirculator->source()->point().y() > maxPoint.y())
                    {
                        inside = false;
                        break;
                    }

                } while (++edgeCirculator != startEdgeCirc);

                //Skip faces outside of bounding box
                if (!inside)
                {
                    continue;
                }




                //Determine face color
                if (fit->data().localMaximaCount != -1)
                {
                    //double shade = (((double)fit->data().localMaximaCount / (double)prickliness));
                    //double shade = (((double)(fit->data().localMaximaCount - currentViewMaximaMin) / (double)(currentViewMaximaMax - currentViewMaximaMin)));


                    ////Green -> Red (https://stackoverflow.com/a/7947812)
                    //double r = std::fmin(2.0 * shade, 1.0);
                    //double g = std::fmin(2.0 * (1.0 - shade), 1.0);
                    //double b = 0.0;

                    ////Yellow -> Red
                    //double r = 1.0;
                    //double g = 1.0 - shade;
                    //double b = 0;

                    //PricklyColor color = GetColorForPercentage((double)(fit->data().localMaximaCount - currentViewMaximaMin) / (double)(currentViewMaximaMax - currentViewMaximaMin));
                    PricklyColor color = GetColorForPercentage((double)(fit->data().localMaximaCount) / (double)(currentViewMaximaMax));


                    //Polygon opening tag, set both stroke and fill to pickly shade so we don't get the wierd lines in between the faces
                    ipeFileStream << "<path stroke = \"" << color.r << " " << color.g << " " << color.b << "\" fill=\"" << color.r << " " << color.g << " " << color.b << "\">\n";
                }
                else
                {
                    //White face for prickliless (I just made that word)faces
                    ipeFileStream << "<path stroke = \"1 1 1\" fill=\"1 1 1\">\n";
                }


                edgeCirculator = fit->outer_ccb();
                startEdgeCirc = fit->outer_ccb();

                //Circle along the edges and add the vertices to the ipe polygon
                bool first = true;
                do
                {
                    ipeFileStream << (edgeCirculator->source()->point().x() * Kernel::FT(10)) << " " << (edgeCirculator->source()->point().y() * Kernel::FT(10));

                    if (first)
                    {
                        ipeFileStream << " m\n";
                        first = false;
                    }
                    else
                    {
                        ipeFileStream << " l\n";
                    }
                } while (++edgeCirculator != startEdgeCirc);

                ipeFileStream << "h\n";

                //Polygon closing tag
                ipeFileStream << "</path>\n";

            }
        }


        ipeFileStream << ipeFileEnd;

    }

    //Thank you https://stackoverflow.com/a/7128796
    PricklinessArrangement::PricklyColor PricklinessArrangement::GetColorForPercentage(double percent) const
    {
        static std::vector<std::pair<double, PricklyColor>> percentColors
        {
            { 0.0,   { 1.0,	        1.0,	     1.0 }},
            { 0.125, { 0.929411765,	0.97254902,	 0.694117647 }},
            { 0.25,	 { 0.780392157,	0.91372549,	 0.705882353}},
            { 0.375, { 0.498039216,	0.803921569, 0.733333333}},
            { 0.5,	 { 0.254901961,	0.71372549,	 0.768627451}},
            { 0.625, { 0.11372549,	0.568627451, 0.752941176}},
            { 0.75,	 { 0.133333333,	0.368627451, 0.658823529}},
            { 0.875, { 0.145098039,	0.203921569, 0.580392157}},
            { 1.0,	 { 0.031372549,	0.11372549,	 0.345098039}}
        };

        int i = 1;
        for (; i < percentColors.size() - 1; i++)
        {
            if (percent < percentColors.at(i).first)
            {
                break;
            }
        }

        std::pair<double, PricklyColor> lower = percentColors.at(i - 1);
        std::pair<double, PricklyColor> upper = percentColors.at(i);
        double range = upper.first - lower.first;
        double rangePercent = (percent - lower.first) / range;
        double percentLower = 1 - rangePercent;
        double percentUpper = rangePercent;

        PricklyColor color;
        color.r = lower.second.r * percentLower + upper.second.r * percentUpper;
        color.g = lower.second.g * percentLower + upper.second.g * percentUpper;
        color.b = lower.second.b * percentLower + upper.second.b * percentUpper;

        return color;
    }
}