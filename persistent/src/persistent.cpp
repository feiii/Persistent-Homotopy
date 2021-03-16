#include "HandleTunnelLoop.h"

namespace DartLib
{

/*!
 *	constructor of the CHandleTunnelLoop
 */
CHandleTunnelLoop::CHandleTunnelLoop() : m_pMesh(NULL) {}

CHandleTunnelLoop::CHandleTunnelLoop(M* pMesh) { set_mesh(pMesh); }

void CHandleTunnelLoop::set_mesh(M* pMesh)
{
    m_pMesh = pMesh;

    for (M::VertexIterator viter(m_pMesh); !viter.end(); viter++)
    {
        M::CVertex* pV = *viter;
        pV->pair() = NULL;
    }

    for (M::EdgeIterator eiter(m_pMesh); !eiter.end(); eiter++)
    {
        M::CEdge* pE = *eiter;
        pE->pair() = NULL;
    }

    for (M::FaceIterator fiter(m_pMesh); !fiter.end(); fiter++)
    {
        M::CFace* pF = *fiter;
        pF->pair() = NULL;
    }

    _extract_boundary_surface();
    _extract_interior_volume();
};

/*! extract the boundry surface, define the filtration
* calculate genus
*/
void CHandleTunnelLoop::_extract_boundary_surface()
{
    for (M::FaceIterator fiter(m_pMesh); !fiter.end(); fiter++)
    {
        M::CFace* pF = *fiter;
        if (m_pMesh->boundary(pF))
        {
            m_boundary_faces.insert(pF);
        }
    }

    for (auto pF : m_boundary_faces)
    {
        for (M::FaceEdgeIterator feiter(pF); !feiter.end(); ++feiter)
        {
            M::CEdge* pE = *feiter;
            m_boundary_edges.insert(pE);
        }
    }

    for (auto pE : m_boundary_edges)
    {
        M::CVertex* pV = m_pMesh->edge_vertex(pE, 0);
        M::CVertex* pW = m_pMesh->edge_vertex(pE, 1);

        m_boundary_vertices.insert(pV);
        m_boundary_vertices.insert(pW);
    }

    int euler_number = m_boundary_vertices.size() + m_boundary_faces.size() - m_boundary_edges.size();
    m_genus = (2 - euler_number) / 2;
    std::cout << "Genus of the Boundary Mesh is " << m_genus << std::endl;
    std::cout << "Number of boundary vertices is  " << m_boundary_vertices.size() << std::endl;
    std::cout << "Number of boundary edges is " << m_boundary_edges.size() << std::endl;
    std::cout << "Number of boundary faces is " << m_boundary_faces.size() << std::endl;
}

/*!
 *	extract interior volume, define the filtration
 *  id boundary vertices, edges and faces first. Then to interior vertices, edges and faces.
 *  group vertices, edges, faces together (no matter they are boundary or interior) and id them sequentially.
 */
void CHandleTunnelLoop::_extract_interior_volume()
{
    int vid = 1;
    for (auto pV : m_boundary_vertices)
    {
        pV->idx() = vid++;
    }

    for (M::VertexIterator viter(m_pMesh); !viter.end(); viter++)
    {
        M::CVertex* pV = *viter;
        if (pV->idx() > 0)
            continue;

        pV->idx() = vid++;
        m_inner_vertices.insert(pV);
    }

    int eid = 1;
    for (auto pE : m_boundary_edges)
    {
        pE->idx() = eid++;
    }

    for (M::EdgeIterator eiter(m_pMesh); !eiter.end(); eiter++)
    {
        M::CEdge* pE = *eiter;
        if (pE->idx() > 0)
            continue;

        pE->idx() = eid++;
        m_inner_edges.insert(pE);
    }

    int fid = 1;
    for (auto pF : m_boundary_faces)
    {
        pF->idx() = fid++;
    }

    for (M::FaceIterator fiter(m_pMesh); !fiter.end(); fiter++)
    {
        M::CFace* pF = *fiter;
        if (pF->idx() > 0)
            continue;

        pF->idx() = fid++;
        m_inner_faces.insert(pF);
    }
}
/*!
*   pair the interior simplices 
*/
void CHandleTunnelLoop::interior_volume_pair()
{
    _pair(m_inner_vertices);
    std::cout << "finished inner vertices pair " << std::endl;
    _pair(m_inner_edges);
    std::cout << "finished inner edges pair " << std::endl;
    _pair(m_inner_faces);
    std::cout << "finished inner faces pair " << std::endl;

    std::cout << "After pairing the Interior Volume :" << std::endl;

    std::vector<M::CEdge*> handle_loops;
    for (auto eiter = m_generators.begin(); eiter != m_generators.end(); eiter++)
    {
        M::CEdge* pE = *eiter;
        if (pE->generator() && pE->pair() == NULL)
        {
            std::cout << "Generator Edge " << pE->idx() << std::endl;
        }
        else
        {
            handle_loops.push_back(pE);
            std::cout << "Killer Face " << pE->pair()->idx() << std::endl;
        }
    }
    
    for (size_t i = 0; i < handle_loops.size(); i++)
    {
        M::CFace* pF = handle_loops[i]->pair();
        _mark_loop(pF);
    }
}

CMyTMesh::CDart* CHandleTunnelLoop::surface_reflection(M::CEdge* e){
    M::CDart* pD = e->dart();
    //M::CFace* pf = pD->cell(2); //original face
    M::CDart* pD0 = pD;
    M::CFace* res = NULL;

    do{
        pD0 = pD0->beta(2);
        if (pD0->beta(3) != NULL){
            pD0 = pD0->beta(3);
        }
        res = pD0->cell(2);
    } while(!m_pMesh->boundary(res));
    if (pD == pD0 || res->visited()){
        return NULL;
    }
    else{
        return pD0;
    }
}

/*
Tunnel and Handle loop algo will give the generators of homology/homotopy groups.
Here, we use Van Kampen theorem to calculate the relationship of generators. From CW complex point
of view, if a surface is obtained by attaching D^2 along some 1-cells (via j:D2->S), 
if some of the 1-cells are generators, then the image of these 1-cells in j(S^1) are relations.

With an existing triangulation, we attach D^2 along the 1-cells to obtain target surface S. For any
edge e not in the 1-cell, the D^2 must attach in the way of ee^-1=1, tivial. So we iterate all the 
faces, and mark the visited faces. If an edge is shared by two faces, then the opposing faces should 
have opposite direction to cancel the shared edges. And if an edge has no unvisited surfaces sharing it,
we sharp the edge.
*/

void CHandleTunnelLoop::_van_kampen(M::CFace* face){ 
    std::cout << "use Van Kampen theorem to draw generator relations" << std::endl;

    std::vector<M::CFace*> F;
    std::vector<M::CEdge*> E;
    std::vector<M::CFace*> visited;
    std::vector<M::CEdge*> sharped;

    F.push_back(face);
    visited.push_back(face);
    M::CEdge* e = NULL;
    for (M::FaceEdgeIterator feiter(face); !feiter.end(); ++feiter){
        M::CEdge* pE = *feiter;
        E.push_back(pE);
    }

    while (!F.empty()){ 
        //new faces and edges obtained in this iteration
        std::vector<M::CFace*> temp_F;
        std::vector<M::CEdge*> temp_E;

        for (auto eiter = E.begin(); eiter != E.end(); eiter++){
            M::CEdge* e = *eiter;
            for (M::EdgeFaceIterator efiter(e); !efiter.end(); efiter++ ){
                M::CFace* f = *efiter;
                if (std::find(F.begin(), F.end(), f)!= F.end() || !m_pMesh->boundary(f)) continue;
                else if (std::find(visited.begin(), visited.end(), f) == visited.end()){
                    //face reflected along e that is not visited or in the previous iteration.
                    temp_F.push_back(f);

                    for (M::FaceEdgeIterator feiter1(f); !feiter1.end(); ++feiter1){
                        M::CEdge* e1 = *feiter1;
                        if (std::find(E.begin(), E.end(), e1) == E.end()){
                            //edges of the reflected face that are not in previous iteratoin or double counted.
                            //count vertices valences, for those doesn't belong to the circle, remove.
                            M::CVertex* v0 = m_pMesh->edge_vertex(e1, 0);
                            M::CVertex* v1 = m_pMesh->edge_vertex(e1, 1);      
                            v0->valence() += 1;
                            v1->valence() += 1;   

                            temp_E.push_back(e1);      
                        }
                    }
                }
                else{
                    //e->sharp() = true;
                    sharped.push_back(e);
                }
            }
        }

        F = temp_F;
        E = temp_E;

        //reset all vertex valence to 0.
        for (auto eiter1 = E.begin(); eiter1 != E.end(); eiter1++){
            M::CVertex* v0 = m_pMesh->edge_vertex(e1, 0);
            M::CVertex* v1 = m_pMesh->edge_vertex(e1, 1);     

            v0->valence() = 0;
            v1->valence() = 0;          
        }

        std::cout << "length of F = " << F.size() << std::endl;
        std::cout << "length of E = " << E.size() << std::endl;
        std::cout << "length of sharped = " << sharped.size() << std::endl;

        visited.insert(visited.end(), temp_F.begin(), temp_F.end());

        if (sharped.size() > 0){
            for (auto iter = E.begin(); iter != E.end(); iter++){
                M::CEdge* te = *iter;
                te->sharp() = true;
            }
            break;
        }
    }
}

/*!
*   pair simplices on the boundary surface
*/
void CHandleTunnelLoop::boundary_surface_pair()
{
    _pair(m_boundary_vertices);
    std::cout << "finished boundary vertices pair " << std::endl;
    _pair(m_boundary_edges);
    std::cout << "finished boundary edges pair " << std::endl;
    _pair(m_boundary_faces);
    std::cout << "finished boundary faces pair " << std::endl;

    std::cout << "After Pairing the boundary surface: " << std::endl;

    for (auto eiter = m_boundary_edges.begin(); eiter != m_boundary_edges.end(); eiter++)
    {
        M::CEdge* pE = *eiter;
        if (pE->generator() && pE->pair() == NULL)
        {
            std::cout << "Generator Edge " << pE->idx() << std::endl;
            m_generators.insert(pE);
        }
    }

    for (auto fiter = m_boundary_faces.begin(); fiter != m_boundary_faces.end(); fiter++){
        M::CFace* pF = *fiter;
        if (pF->generator()){
            std::cout << "generator face after boundary paring = " << pF->idx() << std::endl;
            _van_kampen(pF);
        }
    }
}

/*!
 * pair vertices
 */
void CHandleTunnelLoop::_pair(std::set<M::CVertex*>& vertices)
{
    for (auto viter = vertices.begin(); viter != vertices.end(); viter++)
    {
        M::CVertex* pV = *viter;
        //label the generators
        pV->generator() = true;
    }
};

/*!
 *	pair edges;
 */
void CHandleTunnelLoop::_pair(std::set<M::CEdge*>& edges)
{

    for (auto eiter = edges.begin(); eiter != edges.end(); eiter++)
    {
        M::CEdge* pE = *eiter;
        //std::cout << ".";
        Cycle<M::CVertex, Compare<M::CVertex>> vcycle;

        M::CVertex* v  = NULL;
        M::CVertex* v1 = m_pMesh->edge_vertex(pE, 0);
        M::CVertex* v2 = m_pMesh->edge_vertex(pE, 1);    

        vcycle.add(v1);
        vcycle.add(v2);

        v  = vcycle.head();

        while (v != NULL && v->pair() != NULL && !vcycle.empty()){
            M::CEdge* e    = v->pair();
            M::CVertex* v1 = m_pMesh->edge_vertex(e, 0);
            M::CVertex* v2 = m_pMesh->edge_vertex(e, 1);

            vcycle.add(v1);
            vcycle.add(v2);

            v = vcycle.head();
        }

        if (!vcycle.empty()){
            v->pair() = pE;
        }
        else{
            pE->generator() = true;
        }            
    }
};

/*!
 *	pair faces
 */
void CHandleTunnelLoop::_pair(std::set<M::CFace*>& faces)
{
    for (auto fiter = faces.begin(); fiter != faces.end(); fiter++)
    {
        M::CFace* pF = *fiter;
        //std::cout << "-";
        Cycle<M::CEdge, Compare<M::CEdge>> ecycle;

        M::CEdge* e = NULL;
        for (M::FaceEdgeIterator feiter(pF); !feiter.end(); ++feiter){
            M::CEdge* pE = *feiter;
            ecycle.add(pE);
        }
        e = ecycle.head(); //need to be youngest
        while (e != NULL && !e->generator()){ //need to be positive
            ecycle.add(e);
            e = ecycle.head();
        }        

        //loop if youngest generator is paired and ecycle not empty
        while (e != NULL && e->pair() != NULL && !ecycle.empty()){ 
            M::CFace* f = e->pair();
            //change orientation of the face
            for (M::FaceEdgeIterator feiter(f); !feiter.end(); ++feiter){
                M::CEdge* e1 = *feiter;
                ecycle.add(e1);
            }

            e = ecycle.head(); //need to be youngest
            while (e != NULL && !e->generator()){ //need to be positive
                ecycle.add(e);
                e = ecycle.head();
            }
        }

        if (!ecycle.empty()){
            e->pair() = pF;       
        }
        else{
            pF->generator() = true;
        }
    }
};

//pair boundary faces with halfedges.
// void CHandleTunnelLoop::_ordered_pair(std::set<M::CFace*>& faces)
// {
//     std::cout << "number of boundary faces = " << faces.size() << std::endl;
//     for (auto fiter = faces.begin(); fiter != faces.end(); fiter++)
//     {
//         M::CFace* pF = *fiter;
//         std::cout << "face index = " << pF->idx() << std::endl;

//         ordered_Cycle<M::CEdge, Compare<M::CEdge>> ecycle;

//         M::CEdge* e = NULL;
//         for (M::FaceEdgeIterator feiter(pF); !feiter.end(); ++feiter){
//             M::CEdge* pE = *feiter;
//             ecycle.add(pE);
//         }
//         ecycle._clean_cycle();
//         e = ecycle.head(); //need to be youngest
//         while (e != NULL && !e->generator()){ //need to be positive
//             ecycle.remove(e);
//             e = ecycle.head();
//         }        

//         //loop if youngest generator is paired and ecycle not empty
//         while (e != NULL && e->pair() != NULL && !ecycle.empty()){ 
//             std::cout << "edge index = " << e->idx() << std::endl;

//             M::CFace* f = e->pair();
//             //order the boundary: change direction of e, then go iterate the boundary in that direction.
//             M::CDart* D1 = e->dart();
//             M::CVertex* v1 = D1->cell(0);
//             std::cout << "e0 = " << v1->idx() << std::endl;
//             D1 = D1->beta(2);
//             M::CEdge* e1 = D1->cell(1);
//             M::CVertex* v2 = D1->cell(0);
//             std::cout << "e0 = " << v2->idx() << std::endl;
//             //std::cout << "e1 edge index = " << e->idx() << std::endl;
//             ecycle.insert( e, e1 );
//             M::CDart* D2 = e1->dart();
//             D2 = D2->beta(1);
//             M::CEdge* e2 = D2->cell(1);
//             M::CVertex* v3 = D2->cell(0);
//             std::cout << "e0 = " << v3->idx() << std::endl;

//             while (e2 != e){
//                 //std::cout << "e2 edge index = " << e->idx() << std::endl;
//                 ecycle.insert(e, e2);
//                 M::CDart* D2 = e1->dart();
//                 D2 = D2->beta(1);
//                 M::CEdge* e2 = D2->cell(1);
//             } 

//             ecycle._clean_cycle();
//             e = ecycle.head(); //need to be youngest
//             while (e != NULL && !e->generator()){ //need to be positive
//                 ecycle.remove(e);
//                 e = ecycle.head();
//             }
//         }

//         if (!ecycle.empty()){
//             e->pair() = pF;       
//         }
//         else{
//             pF->generator() = true;
//         }
//     }
// };

/*!
 *	mark the handle and tunnel loops as sharp edges
 *  my logic: start from the original face, extract all its edges. For each positive edge,
 *  find the face it pairs with and repeat. Every edge only pairs with one face, and every face
 *  only kills one edge.
 */
void CHandleTunnelLoop::_mark_loop(M::CFace* face)
{
    //std::queue<M::CFace*> section;
    // std::set<M::CEdge*> visited;
    // if (m_inner_faces.find(face) != m_inner_faces.end()){
    //     visited_faces.insert(face);
    // }
        Cycle<M::CEdge, Compare<M::CEdge>> ecycle;

        M::CEdge* e = NULL;
        for (M::FaceEdgeIterator feiter(face); !feiter.end(); ++feiter){
            M::CEdge* pE = *feiter;
            ecycle.add(pE);
        }
        e = ecycle.head(); //need to be youngerst
        while (e != NULL && !e->generator()){ //need to be positive
            ecycle.add(e);
            e = ecycle.head();
        }        

        while (e != NULL && e->pair() != NULL && !ecycle.empty()){
            std::cout << "head id = " << e->idx() << std::endl;

            M::CFace* f = e->pair();
            for (M::FaceEdgeIterator feiter(f); !feiter.end(); ++feiter){
                M::CEdge* e1 = *feiter;
                ecycle.add(e1);
            }

            e = ecycle.head(); //need to be youngest
            while (e != NULL && (!e->generator() || m_pMesh->boundary(e) )){ //need to be positive
                if(m_pMesh->boundary(e)) e->sharp() = true;
                ecycle.add(e);
                e = ecycle.head();
            }
        }

    // while (!ecycle.empty()){
    //     M::CEdge* e = ecycle.head();
    //     if (m_pMesh->boundary(e)) e->sharp() = true;
    //     ecycle.add(e);
    // }
};


void CHandleTunnelLoop::write_m(const std::string& output)
{
    std::fstream _os(output, std::fstream::out);

    if (_os.fail())
    {
        fprintf(stderr, "Error is opening file %s\n", output);
        return;
    }

    M::CBoundary boundary(m_pMesh);
    const auto& surface = boundary.boundary_surface();
    
    std::set<M::CVertex*> vSet;
    for (auto pF : surface)
    {
        for (M::FaceVertexIterator fviter(pF); !fviter.end(); ++fviter)
        {
            M::CVertex* pV = *fviter;
            vSet.insert(pV);
        }
    }

    for (auto pV : vSet)
    {
        CPoint& p = pV->point();
        _os << "Vertex " << pV->idx();
        _os << " " << p[0] << " " << p[1] << " " << p[2] << "\n";
    }

    for (auto pF : surface)
    {
        _os << "Face " << pF->idx();
        for (M::FaceVertexIterator fviter(pF); !fviter.end(); ++fviter)
        {
            M::CVertex* pV = *fviter;
            _os << " " << pV->idx();
        }
        _os << "\n";
    }

    for (auto pE : m_boundary_edges)
    {
        M::CVertex* pv1 = m_pMesh->edge_vertex(pE, 0);
        M::CVertex* pv2 = m_pMesh->edge_vertex(pE, 1);

        if (pE->sharp())
        {
            _os << "Edge " << pv1->idx() << " " << pv2->idx() << " ";
            _os << "{sharp}" << std::endl;
        }
    }

    _os.close();
}

void CHandleTunnelLoop::exact_boundary(S& surface) 
{
    M::CBoundary boundary(m_pMesh);
    const auto& boundary_surface = boundary.boundary_surface();

    std::set<M::CVertex*> vSet;
    for (auto pF : boundary_surface)
    {
        for (M::FaceVertexIterator fviter(pF); !fviter.end(); ++fviter)
        {
            M::CVertex* pV = *fviter;
            vSet.insert(pV);
        }
    }

    std::vector<CPoint> pts;
    int vid = 1;
    for (auto pV : vSet)
    {
        pV->idx() = vid++;
        pts.push_back(pV->point());
    }

    std::vector<int> indices;
    for (auto pF : boundary_surface)
    {
        for (M::FaceVertexIterator fviter(pF); !fviter.end(); ++fviter)
        {
            M::CVertex* pV = *fviter;
            indices.push_back( pV->idx() - 1 );
        }
    }

    surface.load(pts, indices);
}

void CHandleTunnelLoop::prune()
{
    bool isContinue = false;
    do
    {
        _prune();
        isContinue = _shrink_triangles();
    } while (isContinue);
}

bool CHandleTunnelLoop::_shrink_triangles()
{
    int count = 0;

    for (auto pF : m_boundary_faces)
    {
        int nSharp = 0;
        std::vector<M::CEdge*> edges;
        for (M::FaceEdgeIterator feiter(pF); !feiter.end(); ++feiter)
        {
            M::CEdge* pE = *feiter;
            edges.push_back(pE);
            if (pE->sharp())
                nSharp++;
        }

        /*
               .------.          .      .
                \    /   ----->   \    / 
                 \  /   (unsharp)  \  /
                  \/                \/
         */
        if (nSharp == 3)
        {
            edges[0]->sharp() = false;
            count++;
        }

        /*
               .      .          .------.
                \    /   ----->   
                 \  /    (switch  
                  \/   sharp edges) 
         */
        if (nSharp == 2)
        {
            std::vector<M::CVertex*> verts;
            M::CVertex* pCommonVertex = NULL;
            for (int i = 0; i < 3; ++i)
            {
                if (edges[i]->sharp())
                {
                    for (int j = 0; j < 2; ++j)
                    {
                        M::CVertex* pV = m_pMesh->edge_vertex(edges[i], j);
                        if (std::find(verts.begin(), verts.end(), pV) == verts.end())
                            verts.push_back(pV);
                        else
                        {
                            pCommonVertex = pV;
                            break;
                        }
                    }
                }
                if (pCommonVertex)
                    break;
            }

            if (pCommonVertex && pCommonVertex->valence() == 2)
            {
                for (int i = 0; i < 3; ++i)
                    edges[i]->sharp() = !edges[i]->sharp();
                count++;
            }
        }
    }

    return count > 0;
}

void CHandleTunnelLoop::_prune() //remove branches, try to form circles
{
    std::set<M::CVertex*> vSet;

    for (auto pE : m_boundary_edges)
    {
        if (pE->sharp())
        {
            for (int i = 0; i < 2; ++i)
            {
                M::CVertex* pV = m_pMesh->edge_vertex(pE, i);
                pV->valence() = 0;
                vSet.insert(pV);
            }
        }
    }

    for (auto pE : m_boundary_edges)
    {
        if (pE->sharp())
        {
            for (int i = 0; i < 2; ++i)
            {
                M::CVertex* pV = m_pMesh->edge_vertex(pE, i);
                pV->valence() += 1;
            }
        }
    }

    std::queue<M::CVertex*> vQueue;
    for (auto pV : vSet)
    {
        if (pV->valence() == 1)
            vQueue.push(pV);
    }

    while (!vQueue.empty())
    {
        M::CVertex* pV = vQueue.front();
        vQueue.pop();

        for (M::VertexEdgeIterator veiter(m_pMesh, pV); !veiter.end(); ++veiter)
        {
            M::CEdge* pE = *veiter;
            if (m_pMesh->boundary(pE) && pE->sharp())
            {
                M::CVertex* pW = pV == m_pMesh->edge_vertex(pE, 0) ? m_pMesh->edge_vertex(pE, 1) : m_pMesh->edge_vertex(pE, 0);
                pW->valence() -= 1;

                pE->sharp() = false;

                if (pW->valence() == 1)
                    vQueue.push(pW);

                break;
            }
        }
    }
}

}