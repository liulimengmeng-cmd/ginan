#pragma once

#include <string>
#include <yaml-cpp/yaml.h>

// Map misses must not allocate undefined descendants in the source YAML's
// shared memory holder. Keep every existing non-map's legacy lookup semantics.
inline YAML::Node zhangYamlReadOnlyChildPreservingFallback(
    YAML::Node parent, const std::string& key)
{
    if (!parent.IsDefined())
        return YAML::Node(YAML::NodeType::Undefined);
    if (parent.IsMap())
    {
        const YAML::Node& readonly = parent;
        const YAML::Node child = readonly[key];
        // reset(Zombie) throws before IsDefined can be tested. Normalize first
        // to preserve tryGetFromYaml's as<T>/IsSequence/IsScalar fallback.
        if (!child.IsDefined())
            return YAML::Node(YAML::NodeType::Undefined);
        return child;
    }
    // Preserve Sequence-to-Map/Null mutations, original Mark and Scalar
    // BadSubscript. Do not Clone or approximate sequence string indexing.
    return parent[key];
}
