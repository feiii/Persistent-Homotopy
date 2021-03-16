#ifndef _DARTLIB_HANDLE_TUNNEL_LOOP_H_
#define _DARTLIB_HANDLE_TUNNEL_LOOP_H_

#include <stdio.h>
#include <list>
#include <set>
#include <queue>
#include <fstream>
#include <iostream>
#include <limits>

#include "MyTMesh.h"
#include "MyMesh.h"

namespace DartLib
{

class CHandleTunnelLoop
{
  public:
    using M = CMyTMesh;
    using S = CMyMesh;

    CHandleTunnelLoop();

    CHandleTunnelLoop(M* pMesh);

    void set_mesh(M* pMesh);

    void boundary_surface_pair();

    void interior_volume_pair();

    std::set<M::CEdge*>& boundary_edges() { return m_boundary_edges; };

    void write_m(const std::string& output);
    
    void exact_boundary(S& surface);

    void prune();
  protected:
    void _extract_boundary_surface();

    void _extract_interior_volume();

    M::CDart* surface_reflection(M::CEdge* edge);

    void _pair(std::set<M::CVertex*>& vertices);
    void _pair(std::set<M::CEdge*>& edges);
    void _pair(std::set<M::CFace*>& faces);
    void _ordered_pair(std::set<M::CFace*>& faces);

    //add the following to compute Van Kempen relationship
    void _van_kampen(M::CFace* faces);

    void _mark_loop(M::CFace* face);
    void _mark_loop(M::CEdge* edge);

    bool _shrink_triangles();
    void _prune();
    
  protected:
    M* m_pMesh;

    int m_genus;

    /* boundary surface */
    std::set<M::CVertex*> m_boundary_vertices;
    std::set<M::CEdge*>   m_boundary_edges;
    std::set<M::CFace*>   m_boundary_faces;

    /* interior volume */
    std::set<M::CVertex*> m_inner_vertices;
    std::set<M::CEdge*>   m_inner_edges;
    std::set<M::CFace*>   m_inner_faces;

    std::set<M::CEdge*>   m_generators;
};

template <typename T>
class Compare
{
  public:
    bool operator()(T*& t1, T*& t2)
    {
        if (t1->idx() < t2->idx())
            return true;
        return false;
    }
};

template <class T, class C>
class Cycle
{
  public:
    T* head()
    {
        T* phead = NULL;
        int max_id = -1;
        typename std::vector<T*>::iterator iter;

        for (typename std::vector<T*>::iterator viter = m_cycle.begin(); viter != m_cycle.end(); viter++)
        {
            T* pv = *viter;
            // if (pv->idx() > max_id && pv->generator())
            if (pv->idx() > max_id)
            {
                phead = pv;
                max_id = pv->idx();
            }
        }

        // if (phead != NULL) m_cycle.remove(phead);
        return phead;
    }

    void add(T* pw)
    {
        auto viter = find(m_cycle.begin(), m_cycle.end(), pw);
        if (viter != m_cycle.end())
            m_cycle.erase(viter);
        else
            m_cycle.push_back(pw);
    };

    bool empty() { return m_cycle.empty(); };

    void print()
    {
        for (typename std::vector<T*>::iterator viter = m_cycle.begin(); viter != m_cycle.end(); viter++)
        {
            T* pv = *viter;
            std::cout << pv->idx() << " ";
        }
        std::cout << std::endl;
    }

    std::vector<T*> _output(){ 
        return m_cycle;
    }    

  protected:
    std::vector<T*> m_cycle;
};

//New container for persistence homotopy
template <class T, class C>
class ordered_Cycle
{
  public:
    using M = CMyTMesh;
    
    T* head()
    {
        T* phead = NULL;
        int max_id = -std::numeric_limits<int>::max();
        typename std::vector<T*>::iterator iter;

        for (typename std::vector<T*>::iterator viter = m_cycle.begin(); viter != m_cycle.end(); viter++)
        {
            T* pv = *viter;
            // if (pv->idx() > max_id && pv->generator())
            if (pv->idx() > max_id)
            {
                phead = pv;
                max_id = pv->idx();
            }
        }
        // if (phead != NULL) m_cycle.remove(phead);
        return phead;
    }

    void add(T* pw){
        m_cycle.push_back(pw);
    }

    void remove(T* pw){
        auto viter = find(m_cycle.begin(), m_cycle.end(), pw);
        if (viter != m_cycle.end()){
            m_cycle.erase(viter);
        }
        else{
            std::cout << "pw doesn't exist, can't remove" << std::endl;
        }
    }

    void insert(T* pv, std::vector<T*> boundary ) //insert boudary edges after pv
    {
        auto viter = find(m_cycle.begin(), m_cycle.end(), pv); //viter is the index, pw the value
        
        if (viter != m_cycle.end()){
            ++viter;
            m_cycle.insert( viter, boundary.begin(), boundary.end() );  
            return ;
        }
        else{
            std::cout << "pv doesn't exist, can't insert" << std::endl;
            //m_cylce.push_back(pw);
        }
    };

    bool empty() { return m_cycle.empty(); };

    int length() { return m_cycle.size(); };

    void print()
    {
        for (typename std::vector<T*>::iterator viter = m_cycle.begin(); viter != m_cycle.end(); viter++)
        {
            T* pv = *viter;
            std::cout << pv->idx() << " ";
        }
        std::cout << std::endl;
    }

    void _clean_cycle(){
        for (int i = 1; i != m_cycle.size();){
            M::CEdge* curr = m_cycle[i];
            M::CEdge* prev = m_cycle[i-1];
            //std::cout << "this index = " << curr->idx() << std::endl;
            if (curr != NULL && prev != NULL && curr == prev){
                m_cycle.erase(m_cycle.begin()+i-1, m_cycle.begin()+i+1);
                i -= 1;;
            }
            else{
                i += 1;
            }
        }            
    }

    std::vector<T*> _output(){ 
        return m_cycle;
    }

  protected:
    std::vector<T*> m_cycle;
};

};     // namespace DartLib
#endif //!_DARTLIB_HANDLE_TUNNEL_LOOP_H_