#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

// Integer evidence ONLY. No FLOAT moments, observation factors or state indices.
// Edges mean N(a)-N(b); each independent edge carries WL then IF soft noise.
// These variances are the research output regularisation model, not a claim
// of statistically independent integer decisions or calibrated fix probability.
struct ZhangUserHoldNode
{
    std::string name, identity, domain;
    bool operator==(const ZhangUserHoldNode& other) const
    { return name==other.name && identity==other.identity && domain==other.domain; }
};

struct ZhangUserHoldEdge
{
    ZhangUserHoldNode a,b;
    std::int64_t wl=0,l1=0;
};

inline std::int64_t zhangUserHoldChecked(__int128 value)
{
    if (value<std::numeric_limits<std::int64_t>::min() || value>std::numeric_limits<std::int64_t>::max())
        throw std::overflow_error("USER_HOLD_INTEGER_OVERFLOW");
    return static_cast<std::int64_t>(value);
}

inline std::int64_t zhangUserHoldInteger(double value)
{
    if (!std::isfinite(value) || value!=std::floor(value) || std::abs(value)>9007199254740991.0)
        throw std::runtime_error("USER_HOLD_NONEXACT_INPUT_INTEGER");
    return static_cast<std::int64_t>(value);
}

struct ZhangUserHoldProjection
{
    std::vector<ZhangUserHoldEdge> edges;
    Eigen::MatrixXd transform; // edge-path coefficients, WL/IF interleaved
    Eigen::MatrixXd noise;
};

class ZhangUserConditionalHold
{
public:
    std::vector<ZhangUserHoldEdge> edges;
    Eigen::MatrixXd noise=Eigen::MatrixXd(0,0);

    // Return an exact path a-b in the stored forest; no rational-rowspace or
    // approximate numerical rank test can authorise an integer here.
    bool path(const ZhangUserHoldNode& a,const ZhangUserHoldNode& b,
        std::vector<int>& weights,std::int64_t& wl,std::int64_t& l1) const
    {
        weights.assign(edges.size(),0); wl=0; l1=0;
        if (a.domain!=b.domain || a.domain.empty()) return false;
        std::vector<ZhangUserHoldNode> visited{a};
        std::vector<std::vector<int>> paths{weights};
        for (std::size_t next=0;next<visited.size();++next)
        {
            if (visited[next]==b)
            {
                weights=paths[next];
                __int128 w=0,n=0;
                for (std::size_t i=0;i<edges.size();++i)
                { w+=(__int128)weights[i]*edges[i].wl; n+=(__int128)weights[i]*edges[i].l1; }
                wl=zhangUserHoldChecked(w); l1=zhangUserHoldChecked(n); return true;
            }
            for (std::size_t i=0;i<edges.size();++i)
            {
                ZhangUserHoldNode node; int sign=0;
                if (edges[i].a==visited[next]) { node=edges[i].b; sign=1; }
                else if (edges[i].b==visited[next]) { node=edges[i].a; sign=-1; }
                if (!sign || std::find(visited.begin(),visited.end(),node)!=visited.end()) continue;
                visited.push_back(node); auto p=paths[next]; p[i]=sign; paths.push_back(std::move(p));
            }
        }
        return false;
    }

    // Caller puts its preferred current reference first. Missing old reference
    // nodes are algebraically eliminated; remaining pairs keep the FULL path
    // noise covariance, including off-diagonal terms from a shared old edge.
    ZhangUserHoldProjection project(const std::vector<ZhangUserHoldNode>& active) const
    {
        ZhangUserHoldProjection out;
        std::vector<std::vector<int>> weights;
        std::vector<bool> covered(active.size(),false);
        for (std::size_t root=0;root<active.size();++root)
        {
            if (covered[root]) continue;
            covered[root]=true;
            for (std::size_t i=root+1;i<active.size();++i)
            {
                if (covered[i]) continue;
                std::vector<int> p; std::int64_t wl,l1;
                if (!path(active[i],active[root],p,wl,l1)) continue;
                covered[i]=true; out.edges.push_back({active[i],active[root],wl,l1}); weights.push_back(std::move(p));
            }
        }
        out.transform=Eigen::MatrixXd::Zero(2*out.edges.size(),2*edges.size());
        for (std::size_t i=0;i<weights.size();++i)
            for (std::size_t j=0;j<edges.size();++j)
            {
                out.transform(2*i,2*j)=weights[i][j];
                out.transform(2*i+1,2*j+1)=weights[i][j];
            }
        out.noise=out.transform*noise*out.transform.transpose();
        return out;
    }

    std::vector<ZhangUserHoldNode> nodes() const
    {
        std::vector<ZhangUserHoldNode> result;
        for (const auto& e:edges) for (const auto& node:{e.a,e.b})
            if (std::find(result.begin(),result.end(),node)==result.end()) result.push_back(node);
        return result;
    }

    // Retire a replaced physical arc/state/product chart, but first preserve
    // every consequence amongst the unaffected old nodes. Absence alone is
    // NOT a slip. A returning node must match its complete certified identity.
    int retireReplaced(const std::vector<ZhangUserHoldNode>& observed)
    {
        const auto old=nodes(); std::vector<ZhangUserHoldNode> retain; int retired=0;
        for (const auto& node:old)
        {
            const auto replacement=std::find_if(observed.begin(),observed.end(),
                [&](const auto& n){return n.name==node.name;});
            if (replacement!=observed.end() && !(node==*replacement)) ++retired;
            else retain.push_back(node);
        }
        if (retired)
        { auto p=project(retain); edges=std::move(p.edges); noise=std::move(p.noise); }
        return retired;
    }

    // All-or-nothing candidate union. Redundant identical constraints add no
    // rows/noise. Conflicts cannot replace previously accepted integer evidence.
    bool merge(const std::vector<ZhangUserHoldEdge>& candidate,double wlVariance,double ifVariance,int& added)
    {
        auto proposed=*this; added=0;
        for (const auto& edge:candidate)
        {
            if (edge.a.domain.empty() || edge.a.domain!=edge.b.domain || edge.a==edge.b ||
                edge.a.identity.empty() || edge.b.identity.empty() ||
                !(wlVariance>0) || !(ifVariance>0)) return false;
            std::vector<int> p; std::int64_t wl,l1;
            if (proposed.path(edge.a,edge.b,p,wl,l1))
            { if (wl!=edge.wl || l1!=edge.l1) { added=0; return false; } continue; }
            const int size=proposed.noise.rows();
            proposed.noise.conservativeResize(size+2,size+2);
            proposed.noise.bottomRows(2).setZero(); proposed.noise.rightCols(2).setZero();
            proposed.noise(size,size)=wlVariance; proposed.noise(size+1,size+1)=ifVariance;
            proposed.edges.push_back(edge); ++added;
        }
        *this=std::move(proposed); return true;
    }
};
