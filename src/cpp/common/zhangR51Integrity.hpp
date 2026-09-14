#pragma once
#include "common/zhangIntegerAudit.hpp"
#include <cstdlib>
#include <cstring>
#include "common/eigenIncluder.hpp"

inline bool zhangR51Enabled() {
    static const bool enabled=[] {const char* p=std::getenv("ZHANG_R51_ENABLE");return p && std::string(p)=="1";}();
    return enabled;
}
struct ZhangR51ConflictWitness {
    bool valid=false,feasible=false;
    std::string kind="INVALID_INPUT";
    ZhangExactInteger modulus=0,rhs=0,remainder=0;
    ZhangExactVector combination;
};
inline ZhangR51ConflictWitness zhangR51ExactConflictWitness(
    ZhangExactMatrix matrix,
    ZhangExactVector values
)
{
    ZhangR51ConflictWitness out;
    const std::size_t rowCount = matrix.size();
    if (values.size()!=rowCount) return out;
    const std::size_t columnCount = matrix.empty()
        ? 0
        : matrix.front().size();
    for (const auto& row : matrix)
    {
        if (row.size() != columnCount)
        {
            return {};
        }
    }
    if (rowCount == 0) {out.valid=true;out.feasible=true;return out;}
    ZhangExactMatrix leftTransform = zhangExactIdentityMatrix(rowCount);
    auto swapRows = [&](std::size_t left, std::size_t right)
    {
        std::swap(matrix[left], matrix[right]);
        std::swap(values[left], values[right]);
        std::swap(leftTransform[left], leftTransform[right]);
    };
    auto addRowMultiple = [&](std::size_t destination,
                              std::size_t source,
                              const ZhangExactInteger& multiplier)
    {
        for (std::size_t column = 0; column < columnCount; column++)
        {
            matrix[destination][column] += multiplier * matrix[source][column];
        }
        values[destination] += multiplier * values[source];
        for (std::size_t c=0;c<rowCount;++c)
            leftTransform[destination][c] += multiplier * leftTransform[source][c];
    };
    auto swapColumns = [&](std::size_t left, std::size_t right)
    {
        for (auto& row : matrix)
        {
            std::swap(row[left], row[right]);
        }
    };
    auto addColumnMultiple = [&](std::size_t destination,
                                 std::size_t source,
                                 const ZhangExactInteger& multiplier)
    {
        for (auto& row : matrix)
        {
            row[destination] += multiplier * row[source];
        }
    };

    std::size_t rank = 0;
    while (rank < rowCount && rank < columnCount)
    {
        std::size_t selectedRow = rowCount;
        std::size_t selectedColumn = columnCount;
        for (std::size_t row = rank; row < rowCount; row++)
        {
            for (std::size_t column = rank; column < columnCount; column++)
            {
                if (matrix[row][column] != 0 &&
                    (selectedRow == rowCount ||
                     zhangExactAbs(matrix[row][column]) <
                         zhangExactAbs(matrix[selectedRow][selectedColumn])))
                {
                    selectedRow = row;
                    selectedColumn = column;
                }
            }
        }
        if (selectedRow == rowCount)
        {
            break;
        }
        swapRows(rank, selectedRow);
        swapColumns(rank, selectedColumn);

        while (true)
        {
            bool restart = false;
            for (std::size_t row = rank + 1; row < rowCount; row++)
            {
                if (matrix[row][rank] == 0)
                {
                    continue;
                }
                ZhangExactInteger quotient =
                    matrix[row][rank] / matrix[rank][rank];
                addRowMultiple(row, rank, -quotient);
                if (matrix[row][rank] != 0)
                {
                    swapRows(row, rank);
                }
                restart = true;
                break;
            }
            if (restart)
            {
                continue;
            }
            for (std::size_t column = rank + 1;
                 column < columnCount;
                 column++)
            {
                if (matrix[rank][column] == 0)
                {
                    continue;
                }
                ZhangExactInteger quotient =
                    matrix[rank][column] / matrix[rank][rank];
                addColumnMultiple(column, rank, -quotient);
                if (matrix[rank][column] != 0)
                {
                    swapColumns(column, rank);
                }
                restart = true;
                break;
            }
            if (!restart)
            {
                break;
            }
        }
        rank++;
    }

    out.valid=true;
    for (std::size_t i=0;i<rowCount;++i) {
        const ZhangExactInteger d=i<rank?zhangExactAbs(matrix[i][i]):ZhangExactInteger(0);
        if ((d==0 && values[i]!=0) || (d!=0 && values[i]%d!=0)) {
            out.feasible=false;out.modulus=d;out.rhs=values[i];
            out.remainder=d==0?values[i]:ZhangExactInteger(values[i]%d);
            out.combination=leftTransform[i];
            out.kind=d==0?"RATIONAL_AFFINE_CONTRADICTION":"INTEGER_DIVISIBILITY_CONTRADICTION";
            return out;
        }
    }
    out.feasible=true;out.kind="FEASIBLE";return out;
}

// Binary audit fingerprint avoids formatting a multi-GiB full covariance.
// The receipt is consumed in the same writer invocation; numerical candidate
// root equality remains enforced by the upstream delivery moment contract.
inline std::string zhangR51NumericRoot(const VectorXd& x,const MatrixXd& p) {
    std::uint64_t h=1469598103934665603ULL;
    auto add=[&](const char* bytes,std::size_t n){for(std::size_t i=0;i<n;++i){h^=static_cast<unsigned char>(bytes[i]);h*=1099511628211ULL;}};
    add(reinterpret_cast<const char*>(x.data()),x.size()*sizeof(double));
    add(reinterpret_cast<const char*>(p.data()),p.size()*sizeof(double));
    return std::to_string(x.size())+"x"+std::to_string(p.rows())+"x"+std::to_string(p.cols())+":"+std::to_string(h);
}
