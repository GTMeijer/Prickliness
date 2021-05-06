#pragma once

namespace Prickliness
{
    //2D arrangements
    struct PerpPlaneIndices
    {
        bool operator==(const PerpPlaneIndices& rhs) const
        {
            return (vertexIndex == rhs.vertexIndex && planeIndex == rhs.planeIndex);
        }

        int vertexIndex = -1;
        int planeIndex = -1;
    };

    //Struct wrapping the localMaximaCount per face so we can instantiate it with -1
    struct FaceCounter
    {
        int localMaximaCount = -1;
    };

    typedef CGAL::Arr_linear_traits_2<Kernel> Arr_traits_2;
    //typedef CGAL::Arr_segment_traits_2<Kernel> Arr_segment_traits_2;
    //Use the consolidated curve data traits to add a data field to the lines, if they overlap these traits will automatically append them in a list
    //This is used to keep track of the accompanying perp plane and terrain vertex
    typedef CGAL::Arr_consolidated_curve_data_traits_2<Arr_traits_2, PerpPlaneIndices> Arr_segment_data_traits_2;
    //Construct a DCEL that adds auxiliary data to every face (local maxima count)
    typedef CGAL::Arr_face_extended_dcel<Arr_segment_data_traits_2, FaceCounter> DCEL;
    typedef CGAL::Arrangement_2<Arr_segment_data_traits_2, DCEL> Arrangement_2;

    typedef Arr_segment_data_traits_2::X_monotone_curve_2 Arr_X_monotone_curve_2;
    typedef Arr_segment_data_traits_2::Line_2 Arr_data_line_2;
    typedef Arr_segment_data_traits_2::Curve_2 Arr_data_curve_2;
    typedef Arrangement_2::Face_handle DCELFaceHandle;

    class PricklinessArrangement
    {
    public:

        PricklinessArrangement(const Delaunay& DT);

        void Draw();
        void WriteToIpe(const std::filesystem::path& filePath, bool centerOnOrigin);

        int prickliness = 0;
        Vector_3 pricklyDirection;

        void WriteToFile(const std::filesystem::path& filePath) const;
        void ReadFromFile(const std::filesystem::path& filePath) const;

    private:

        //Get a vector to a point inside a convex DCEL plane
        Vector_3 GetVectorInPlane(const Arrangement_2::Face_const_handle face) const;

        bool IsBoundaryVert(const Delaunay& DT, Delaunay::Face_circulator adjacentFaces) const;

        Arrangement_2 resultArrangement;

        struct PricklyColor
        {
            double r;
            double g;
            double b;
        };

        PricklyColor GetColorForPercentage(double percent) const;
    };
}