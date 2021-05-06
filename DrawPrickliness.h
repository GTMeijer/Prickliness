#ifndef CGAL_DRAW_PRICKLINESS_H
#define CGAL_DRAW_PRICKLINESS_H

#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef CGAL_USE_BASIC_VIEWER

//#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Random.h>

namespace Prickliness
{

    // Viewer class for Arr
    class PricklyArrangementViewerQt : public CGAL::Basic_viewer_qt
    {
        typedef CGAL::Basic_viewer_qt Base;

    public:
        /// Construct the viewer.
        /// @param a_a the arrangement to view
        /// @param title the title of the window
        /// @param anofaces if true, do not draw faces (faces are not computed; this can be
        ///        usefull for very big object where this time could be long)
        PricklyArrangementViewerQt(QWidget* parent,
            const Arrangement_2& a_arr,
            int maxPrickliness,
            Vector_3 pricklyDirection,
            const char* title = "Prickliness Arrangement Viewer",
            bool anofaces = false) :
            // First draw: vertices; edges, faces; multi-color; no inverse normal
            Base(parent, title, true, true, true, false, false),
            arr(a_arr),
            m_nofaces(anofaces),
            maxPrickliness(maxPrickliness),
            pricklyDirection(pricklyDirection)
        {
            compute_elements();
        }

    protected:

        void compute_face(Arrangement_2::Face_const_handle fh)
        {
            if (fh->is_unbounded())
            {
                return;
            }

            int shade = 255;
            if (fh->data().localMaximaCount != -1)
            {
                shade = 255 - (int)(255.0 * ((double)fh->data().localMaximaCount / (double)maxPrickliness));
            }

            CGAL::Color c(shade, shade, shade);

            face_begin(c);
            
            auto edgeCirculator = fh->outer_ccb();
            auto startEdgeCirc = fh->outer_ccb();

            //Circle along the edges and add the vertices to the draw face
            do
            {
                add_point_in_face(edgeCirculator->source()->point());
            } while (++edgeCirculator != startEdgeCirc);


            face_end();
        }
        void compute_elements()
        {
            clear();

            // Draw the arrangement faces.
            Arrangement_2::Face_const_iterator fit;
            for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit)
            {
                //Only draw finite faces
                if (!fit->is_unbounded() && !fit->is_fictitious() && fit->has_outer_ccb())
                {
                    if (fit->data().localMaximaCount != -1)
                    {
                        compute_face(fit);
                    }
                }
            }

            // Draw the arrangement edges.
            CGAL::Color edgeColor(0, 0, 0);
            Arrangement_2::Edge_const_iterator eit;
            for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit)
            {

                if (!eit->is_fictitious() && !eit->curve().is_ray())
                {
                    if (eit->face()->data().localMaximaCount != -1 || eit->twin()->face()->data().localMaximaCount != -1)
                    {
                        add_segment(eit->source()->point(), eit->target()->point(), edgeColor);
                    }
                }
            }

            //// Draw the arrangement vertices.
            //Arrangement_2::Vertex_const_iterator vit;
            //for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
            //{
            //    add_point(vit->point());
            //}
        }

        virtual void keyPressEvent(QKeyEvent* e)
        {
            Base::keyPressEvent(e);
        }

    protected:
        const Arrangement_2& arr;
        bool m_nofaces;
        int maxPrickliness;
        Vector_3 pricklyDirection;
    };


    inline void DrawPrickliness(const Arrangement_2& a_arr,
        int maxPrickliness,
        Vector_3 pricklyDirection,
        const char* title = "Prickliness Arrangement Viewer",
        bool nofill = false)
    {
#if defined(CGAL_TEST_SUITE)
        bool cgal_test_suite = true;
#else
        bool cgal_test_suite = qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

        if (!cgal_test_suite)
        {
            int argc = 1;
            const char* argv[2] = { "arr_viewer","\0" };
            QApplication app(argc, const_cast<char**>(argv));
            PricklyArrangementViewerQt mainwindow(app.activeWindow(), a_arr, maxPrickliness, pricklyDirection, title);
            mainwindow.show();
            app.exec();
        }
    }

} // End namespace CGAL

#endif // CGAL_DRAW_PRICKLINESS_H

#endif // CGAL_DRAW_LCC_H