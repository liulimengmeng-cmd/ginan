// Dense golden bodies copied from the frozen pre-P4 source, not wrappers
// calling the optimized HNF. Independent enumeration below checks congruence.
#include "common/zhangIntegerProductGainFrontier.hpp"
#include "common/zhangR51PhysicalImage.hpp"
#include <iostream>
#include <random>
#include <stdexcept>
#include <numeric>
namespace {
void need(bool condition,const char* reason) {if(!condition)throw std::runtime_error(reason);}
inline ZhangExactRowHnf denseGoldenHnf(
    ZhangExactMatrix rows,
    ZhangExactVector values = {},
    bool trackRowTransform = false
)
{
    ZhangExactRowHnf result;
    if (rows.empty())
    {
        return result;
    }
    const std::size_t columns = rows.front().size();
    for (const auto& row : rows)
    {
        zhangExactWorkCheckpoint("zhangIntegerAudit_LOOP");
        if (row.size() != columns)
        {
            result.consistent = false;
            return result;
        }
    }
    if (values.empty())
    {
        values.resize(rows.size());
    }
    if (values.size() != rows.size())
    {
        result.consistent = false;
        return result;
    }
    ZhangExactMatrix rowTransform;
    if (trackRowTransform)
    {
        rowTransform = zhangExactIdentityMatrix(rows.size());
    }

    auto addRowMultiple = [&](std::size_t destination,
                              std::size_t source,
                              const ZhangExactInteger& multiplier)
    {
        for (std::size_t column = 0; column < columns; column++)
        {
            zhangExactWorkCheckpoint("zhangIntegerAudit_LOOP");
            rows[destination][column] += multiplier * rows[source][column];
        }
        values[destination] += multiplier * values[source];
        if (trackRowTransform)
        {
            for (std::size_t column = 0;
                 column < rowTransform.size(); column++)
            {
                zhangExactWorkCheckpoint("zhangIntegerAudit_LOOP");
                rowTransform[destination][column] +=
                    multiplier * rowTransform[source][column];
            }
        }
    };
    auto negateRow = [&](std::size_t row)
    {
        for (auto& value : rows[row])
        {
            zhangExactWorkCheckpoint("zhangIntegerAudit_LOOP");
            value = -value;
        }
        values[row] = -values[row];
        if (trackRowTransform)
        {
            for (auto& value : rowTransform[row])
            {
                zhangExactWorkCheckpoint("zhangIntegerAudit_LOOP");
                value = -value;
            }
        }
    };

    std::size_t pivotRow = 0;
    for (std::size_t column = 0;
         column < columns && pivotRow < rows.size();
         column++)
    {
        zhangExactWorkCheckpoint("zhangIntegerAudit_LOOP");
        std::size_t selected = rows.size();
        for (std::size_t row = pivotRow; row < rows.size(); row++)
        {
            zhangExactWorkCheckpoint("zhangIntegerAudit_LOOP");
            if (rows[row][column] != 0 &&
                (selected == rows.size() ||
                 zhangExactAbs(rows[row][column]) <
                     zhangExactAbs(rows[selected][column])))
            {
                selected = row;
                // A nonzero integer cannot have absolute value below one.
                // The old strict-< scan also chose this FIRST unit pivot.
                if (zhangExactAbs(rows[row][column]) == 1) break;
            }
        }
        if (selected == rows.size())
        {
            continue;
        }
        std::swap(rows[pivotRow], rows[selected]);
        std::swap(values[pivotRow], values[selected]);
        if (trackRowTransform)
        {
            std::swap(rowTransform[pivotRow], rowTransform[selected]);
        }

        bool reduced = false;
        do
        {
            reduced = true;
            for (std::size_t row = pivotRow + 1; row < rows.size(); row++)
            {
                zhangExactWorkCheckpoint("zhangIntegerAudit_LOOP");
                if (rows[row][column] == 0)
                {
                    continue;
                }
                ZhangExactInteger quotient =
                    rows[row][column] / rows[pivotRow][column];
                addRowMultiple(row, pivotRow, -quotient);
                if (rows[row][column] != 0)
                {
                    std::swap(rows[row], rows[pivotRow]);
                    std::swap(values[row], values[pivotRow]);
                    if (trackRowTransform)
                    {
                        std::swap(
                            rowTransform[row], rowTransform[pivotRow]);
                    }
                }
                reduced = false;
                break;
            }
        }
        while (!reduced);

        if (rows[pivotRow][column] < 0)
        {
            negateRow(pivotRow);
        }
        for (std::size_t row = 0; row < pivotRow; row++)
        {
            zhangExactWorkCheckpoint("zhangIntegerAudit_LOOP");
            if (rows[row][column] == 0)
            {
                continue;
            }
            ZhangExactInteger quotient = zhangExactFloorDivide(
                rows[row][column],
                rows[pivotRow][column]
            );
            addRowMultiple(row, pivotRow, -quotient);
        }
        pivotRow++;
    }

    for (std::size_t row = 0; row < rows.size(); row++)
    {
        zhangExactWorkCheckpoint("zhangIntegerAudit_LOOP");
        bool nonzero = std::any_of(
            rows[row].begin(),
            rows[row].end(),
            [](const auto& value) { return value != 0; }
        );
        if (!nonzero)
        {
            result.consistent &= values[row] == 0;
            continue;
        }
        result.basis.push_back(std::move(rows[row]));
        result.values.push_back(std::move(values[row]));
        if (trackRowTransform)
        {
            result.rowTransform.push_back(std::move(rowTransform[row]));
        }
    }
    return result;
}

inline ZhangExactAffineIntegerQuotient denseGoldenAffine(
	const ZhangExactMatrix& deterministicRows,
	const ZhangExactVector& deterministicValues,
	int ambientDimension,
    ZhangExactQuotientWork work=ZhangExactQuotientWork::FULL_SEARCH_COORDINATES)
{
    ZhangPhaseTimer affineTimer("EXACT_AFFINE_TOTAL");
    ZhangR48ExactTimer feasibilityTimer("R48_EXACT_FEASIBILITY",deterministicRows,ambientDimension);
	ZhangExactAffineIntegerQuotient result;
	result.ambientDimension = ambientDimension;
	if (ambientDimension <= 0 ||
		deterministicRows.size() != deterministicValues.size())
	{
		result.failureReason = "AFFINE_QUOTIENT_DIMENSION_MISMATCH";
		return result;
	}
	for (const auto& row : deterministicRows)
	{
	    zhangExactWorkCheckpoint("zhangProductRelationSolver_LOOP");
		if (row.size() != static_cast<std::size_t>(ambientDimension))
		{
			result.failureReason = "AFFINE_QUOTIENT_INVALID_ROW";
			return result;
		}
	}
	const auto hnf = [&]() {
        ZhangPhaseTimer timer("EXACT_AFFINE_ROW_HNF");
        return denseGoldenHnf(deterministicRows, deterministicValues);
    }();
	if (!hnf.consistent)
	{
		result.failureReason = "AFFINE_QUOTIENT_INTEGER_CONFLICT";
		return result;
	}
	result.deterministicBasis = hnf.basis;
	result.deterministicValues = hnf.values;
	result.deterministicRank = static_cast<int>(hnf.basis.size());
	if (result.deterministicRank > ambientDimension)
	{
		result.failureReason = "AFFINE_QUOTIENT_OVER_RANK";
		return result;
	}
	if (work == ZhangExactQuotientWork::FEASIBILITY_ONLY)
	{
		if (result.deterministicRank == 0)
		{
			result.valid = true;
			result.failureReason = "NONE";
			return result;
		}
		ZhangExactMatrix columns(ambientDimension,
			ZhangExactVector(result.deterministicRank));
		for (int column = 0; column < ambientDimension; ++column)
		for (int row = 0; row < result.deterministicRank; ++row) {
			zhangExactWorkCheckpoint("AFFINE_FEASIBILITY_TRANSPOSE");
			columns[column][row] = result.deterministicBasis[row][column];
		}
		const auto feasibility = zhangIntegerRowLatticeContains(
			columns, result.deterministicValues, false);
		result.smithInvariants = feasibility.smithInvariants;
		result.valid = feasibility.contained;
		result.failureReason = result.valid ? "NONE" :
			"AFFINE_QUOTIENT_NO_INTEGER_ORIGIN";
		return result;
	}

	if (result.deterministicRank == 0)
	{
		result.particularSolution = ZhangExactVector(ambientDimension);
		if (work != ZhangExactQuotientWork::PARTICULAR_ONLY)
			result.kernelBasis = zhangExactIdentityMatrix(ambientDimension);
        if(work==ZhangExactQuotientWork::FULL_SEARCH_COORDINATES)
		    result.quotientProjector = zhangExactIdentityMatrix(ambientDimension);
		result.quotientRank = ambientDimension;
		result.valid = true;
		result.failureReason = "NONE";
		return result;
	}

	// Solve D z0=d exactly by treating the columns of D as generators of Z^r.
	ZhangExactMatrix columnGenerators(ambientDimension,
		ZhangExactVector(result.deterministicRank));
	for (int column = 0; column < ambientDimension; column++)
	for (int row = 0; row < result.deterministicRank; row++) {
		zhangExactWorkCheckpoint("AFFINE_PARTICULAR_TRANSPOSE");
		columnGenerators[column][row] = result.deterministicBasis[row][column];
	}
	const auto particular = [&]() {
        ZhangPhaseTimer timer("EXACT_AFFINE_SMITH_ORIGIN");
        // One Smith decomposition proves divisibility, rank and invariants.
        // H and H^T have the same nonzero Smith factors: decomposing H^T
        // separately for a zero membership query proves nothing additional.
        // Replay the SAME column-operation tape for this single RHS instead
        // of storing an ambientDimension^2 right-transform matrix.
        return zhangIntegerRowLatticeContains(
            columnGenerators, result.deterministicValues, true, true);
    }();
    result.smithInvariants = particular.smithInvariants;
    if (particular.rank != result.deterministicRank)
    {
        result.failureReason = "AFFINE_QUOTIENT_RANK_MISMATCH";
        return result;
    }
	if (!particular.contained ||
		particular.combination.size() !=
			static_cast<std::size_t>(ambientDimension))
	{
		result.failureReason = "AFFINE_QUOTIENT_NO_INTEGER_ORIGIN";
		return result;
	}
	result.particularSolution = particular.combination;
	if (work == ZhangExactQuotientWork::PARTICULAR_ONLY)
	{
		result.quotientRank = ambientDimension-result.deterministicRank;
		result.valid = true;
		result.failureReason = "NONE";
		return result;
	}

    feasibilityTimer.finish();
    {
        ZhangPhaseTimer kernelTimer("EXACT_AFFINE_KERNEL");
        ZhangR48ExactTimer timer("R48_EXACT_KERNEL",result.deterministicBasis,ambientDimension);
	    result.kernelBasis = zhangExactIntegerKernel(result.deterministicBasis, ambientDimension);
    }
	result.quotientRank = static_cast<int>(result.kernelBasis.size());
	if (result.deterministicRank + result.quotientRank != ambientDimension)
	{
		result.failureReason = "AFFINE_QUOTIENT_KERNEL_RANK_MISMATCH";
		return result;
	}

    if(work==ZhangExactQuotientWork::PARTICULAR_AND_KERNEL) {
        result.valid=true;result.failureReason="NONE";return result;
    }
	// Find an exact integer left inverse L K^T=I.  The saturated kernel is
	// primitive, so every quotient unit vector must lie in the column lattice
	// of K.
	ZhangExactMatrix kernelColumns(ambientDimension,
		ZhangExactVector(result.quotientRank));
	for (int column = 0; column < ambientDimension; column++)
	for (int row = 0; row < result.quotientRank; row++)
		kernelColumns[column][row] = result.kernelBasis[row][column];
    const auto inverseRows=[&](){
        ZhangR48ExactTimer timer("R48_EXACT_LEFT_INVERSE",kernelColumns);
        return zhangIntegerRowLatticeContainsBatch(kernelColumns,zhangExactIdentityMatrix(result.quotientRank));
    }();
	result.quotientProjector.reserve(result.quotientRank);
	for (int quotient = 0; quotient < result.quotientRank; quotient++)
	{
	    zhangExactWorkCheckpoint("zhangProductRelationSolver_LOOP");
		const auto& inverseRow=inverseRows[quotient];
		if (!inverseRow.contained ||
			inverseRow.combination.size() !=
				static_cast<std::size_t>(ambientDimension))
		{
			result.failureReason = "AFFINE_QUOTIENT_LEFT_INVERSE_FAILED";
			return result;
		}
		result.quotientProjector.push_back(inverseRow.combination);
	}
	const auto identityAudit = zhangExactMultiply(
		result.quotientProjector,
		[&]()
		{
			ZhangExactMatrix transpose(ambientDimension,
				ZhangExactVector(result.quotientRank));
			for (int row = 0; row < result.quotientRank; row++)
			for (int column = 0; column < ambientDimension; column++)
				transpose[column][row] = result.kernelBasis[row][column];
			return transpose;
		}());
	if (identityAudit != zhangExactIdentityMatrix(result.quotientRank))
	{
		result.failureReason = "AFFINE_QUOTIENT_LEFT_INVERSE_AUDIT_FAILED";
		return result;
	}
	result.valid = true;
	result.failureReason = "NONE";
	return result;
}

inline ZhangPrimitiveImageCoordinates denseGoldenImage(const ZhangExactMatrix& rows)
{
	ZhangPrimitiveImageCoordinates result;
	if (rows.empty() || rows.front().empty()) return result;
	const int m = rows.size(), n = rows.front().size();
	ZhangExactMatrix columns(n, ZhangExactVector(m));
	for (int i=0; i<m; ++i)
	{
	    zhangExactWorkCheckpoint("zhangProductRelationSolver_LOOP");
		if (rows[i].size() != n) return result;
		for (int j=0; j<n; ++j) columns[j][i]=rows[i][j];
	}
	const auto image = denseGoldenHnf(columns);
	if (!image.consistent) return result;
	const int r = image.basis.size();
	result.primitiveRows.assign(r, ZhangExactVector(n));
	result.imageGenerators.assign(m, ZhangExactVector(r));
	for (int i=0;i<m;++i) for(int k=0;k<r;++k)
		result.imageGenerators[i][k]=image.basis[k][i];
    const auto memberships=zhangIntegerRowLatticeContainsBatch(image.basis,columns);
	for (int j=0;j<n;++j)
	{
	    zhangExactWorkCheckpoint("zhangProductRelationSolver_LOOP");
		const auto& membership = memberships[j];
		if (!membership.contained || membership.combination.size() != r) return result;
		for (int k=0;k<r;++k) result.primitiveRows[k][j]=membership.combination[k];
	}
	if (r && (!zhangExactPrimitiveRowLattice(result.primitiveRows, n) ||
		zhangExactMultiply(result.imageGenerators, result.primitiveRows) != rows)) return result;
	result.valid=true;
	return result;
}

inline ZhangR47ProductSearchFrame denseGoldenFrame(
    const ZhangExactMatrix& targets,const ZhangExactMatrix& conditioners,
    const ZhangExactVector& values,int dimension,
    ZhangProductFrameWork work=ZhangProductFrameWork::INTEGER_PROJECTOR)
{
    ZhangPhaseTimer frameTimer("PRODUCT_FRAME_TOTAL");
    zhangCaptureProductFrame(targets,conditioners,values,dimension,static_cast<int>(work));
    ZhangR48ExactTimer productTimer("R48_PRODUCT_IMAGE",targets,dimension);
    ZhangR47ProductSearchFrame out;
    if(targets.empty() || conditioners.size()!=values.size() ||
       !zhangExactRectangularMatrix(targets,dimension) ||
       !zhangExactRectangularMatrix(conditioners,dimension)) return out;
    for(int c=0;c<dimension;++c)
    {
        zhangExactWorkCheckpoint("zhangR47ProductDomain_LOOP");
        bool used=false;
        for(const auto& row:targets) used |= row[c]!=0;
        for(const auto& row:conditioners) used |= row[c]!=0;
        if(used) out.columns.push_back(c);
    }
    if(out.columns.empty()) {
        // No physical column remains.  Zero rows are feasible only when every
        // right-hand side is zero; do not let this fast path certify 0 = b.
        for(const auto& value:values) if(value!=0) {
            out.reason="AFFINE_QUOTIENT_INTEGER_CONFLICT";return out;
        }
        out.affine=std::make_shared<const ZhangExactAffineIntegerQuotient>();
        out.valid=true;out.reason="ZERO_PRODUCT_IMAGE";return out;
    }
    auto compact=[&](const ZhangExactMatrix& input)
    {
        ZhangExactMatrix result(input.size(),ZhangExactVector(out.columns.size()));
        for(std::size_t r=0;r<input.size();++r) for(std::size_t c=0;c<out.columns.size();++c) {
            zhangExactWorkCheckpoint("PRODUCT_DOMAIN_COMPACT_COPY");
            result[r][c]=input[r][out.columns[c]];
        }
        return result;
    };
    const auto t=[&](){ZhangPhaseTimer timer("PRODUCT_FRAME_TARGET_COMPACT");return compact(targets);}();
    const auto h=[&](){ZhangPhaseTimer timer("PRODUCT_FRAME_CONDITIONER_COMPACT");return compact(conditioners);}();
    if(work==ZhangProductFrameWork::INTEGER_PROJECTOR)
    {
        ZhangPhaseTimer timer("PRODUCT_FRAME_TARGET_HNF");
        out.targetRank=denseGoldenHnf(t).basis.size();
    }
    const auto quotientWork=work==ZhangProductFrameWork::IMAGE_GENERATORS_ONLY
        ? ZhangExactQuotientWork::PARTICULAR_ONLY
        : ZhangExactQuotientWork::PARTICULAR_AND_KERNEL;
    out.affine=std::make_shared<const ZhangExactAffineIntegerQuotient>(
        denseGoldenAffine(h,values,out.columns.size(),quotientWork));
    if(!out.affine->valid) {out.reason=out.affine->failureReason;return out;}
    out.conditionerRank=out.affine->deterministicRank;
    if(!out.affine->quotientRank) {out.valid=true;out.reason="PRODUCT_FULLY_DETERMINED";return out;}
    ZhangExactMatrix k,m;
    if(work==ZhangProductFrameWork::IMAGE_GENERATORS_ONLY)
    {
        ZhangR48ExactTimer timer("R48_PRODUCT_IMAGE_PROJECTED_KERNEL",
            out.affine->deterministicBasis,out.columns.size());
        m=zhangExactProjectedKernelColumns(
            out.affine->deterministicBasis,t,out.columns.size());
        if(m.empty() || m.front().size()!=
            static_cast<std::size_t>(out.affine->quotientRank))
        {out.reason="PRODUCT_PROJECTED_KERNEL_FAILED";return out;}
    }
    else
    {
        ZhangPhaseTimer timer("PRODUCT_FRAME_KERNEL_PROJECT");
        k.assign(out.columns.size(),ZhangExactVector(out.affine->quotientRank));
        for(int r=0;r<out.affine->quotientRank;++r)
        for(int c=0;c<out.columns.size();++c)
            k[c][r]=out.affine->kernelBasis[r][c];
        m=zhangExactMultiply(t,k);
    }
    const auto image=[&](){ZhangPhaseTimer timer("PRODUCT_FRAME_IMAGE_HNF");return denseGoldenImage(m);}();
    if(!image.valid) {out.reason="PRODUCT_IMAGE_COORDINATES_FAILED";return out;}
    out.imageGenerators=image.imageGenerators;
    out.searchRank=image.primitiveRows.size();
    // R51 constructs and audits its own rational left inverse of this exact
    // image. This path keeps T*V during column operations and never stores
    // the full physical kernel or the unused network integer lift.
    if(work==ZhangProductFrameWork::IMAGE_GENERATORS_ONLY) {
        out.valid=true;
        out.reason=out.searchRank?"EXACT_PRODUCT_IMAGE_QUOTIENT":"PRODUCT_FULLY_DETERMINED";
        return out;
    }
    // Solve only the product-image coordinate RHS, never construct the
    // entire network quotient left inverse. K is saturated, so these rows
    // have exact integer lifts. Audit each lift before exposing it to ILS.
    const auto lifts=[&]() {
        ZhangPhaseTimer liftTimer("PRODUCT_FRAME_INTEGER_LIFT");
        ZhangR48ExactTimer timer("R48_PRODUCT_IMAGE_LIFT",k);
        return zhangIntegerRowLatticeContainsBatch(k,image.primitiveRows);
    }();
    ZhangExactMatrix projected;
    for(std::size_t i=0;i<lifts.size();++i) {
        zhangExactWorkCheckpoint("zhangR47ProductDomain_LOOP");
        if(!lifts[i].contained || lifts[i].combination.size()!=out.columns.size()) {
            out.reason="PRODUCT_IMAGE_INTEGER_LIFT_FAILED";return out;
        }
        projected.push_back(lifts[i].combination);
    }
    if(!projected.empty() && zhangExactMultiply(projected,k)!=image.primitiveRows) {
        out.reason="PRODUCT_IMAGE_INTEGER_LIFT_AUDIT_FAILED";return out;
    }
    for(const auto& row:projected)
    {
        zhangExactWorkCheckpoint("zhangR47ProductDomain_LOOP");
        ZhangExactVector full(dimension); ZhangExactInteger offset=0;
        for(std::size_t c=0;c<out.columns.size();++c)
        {
            zhangExactWorkCheckpoint("zhangR47ProductDomain_LOOP");full[out.columns[c]]=row[c];offset-=row[c]*out.affine->particularSolution[c];}
        out.projector.push_back(std::move(full));out.offsets.push_back(offset);
    }
    out.valid=true;out.reason=out.searchRank?"EXACT_PRODUCT_IMAGE_QUOTIENT":"PRODUCT_FULLY_DETERMINED";
    return out;
}

void compareHnf(const ZhangExactMatrix& rows,const ZhangExactVector& values,bool witnesses)
{
    const auto golden=denseGoldenHnf(rows,values,witnesses);
    ZhangExactRowHnfWorkMetrics metrics;
    const auto actual=zhangExactRowHermiteNormalForm(rows,values,witnesses,&metrics);
    need(actual.consistent==golden.consistent,"HNF consistency changed");
    need(actual.basis==golden.basis,"HNF basis/order differs from dense golden");
    need(actual.values==golden.values,"HNF RHS differs from dense golden");
    need(actual.rowTransform==golden.rowTransform,"integer row witness/order changed");
    need(metrics.inputRows==rows.size(),"per-call input row metric");
    need(metrics.coefficientUpdates<=metrics.denseCoefficientUpdates,"sparse update bound");
    if(witnesses && !rows.empty() && !actual.rowTransform.empty())
    {
        need(zhangExactMultiply(actual.rowTransform,rows)==actual.basis,"independent coefficient pullback");
        ZhangExactVector rhs=values.empty()?ZhangExactVector(rows.size()):values;
        for(std::size_t i=0;i<actual.rowTransform.size();++i)
        {
            ZhangExactInteger v=0;
            for(std::size_t j=0;j<rhs.size();++j)v+=actual.rowTransform[i][j]*rhs[j];
            need(v==actual.values[i],"independent RHS pullback");
        }
    }
}
void compareFrame(const ZhangExactMatrix& targets,const ZhangExactMatrix& h,
    const ZhangExactVector& values,int dimension,ZhangProductFrameWork work)
{
    const auto golden=denseGoldenFrame(targets,h,values,dimension,work);
    const auto actual=zhangR47CompileProductSearchFrame(targets,h,values,dimension,work);
    need(actual.valid==golden.valid && actual.reason==golden.reason,"frame outcome changed");
    need(actual.columns==golden.columns && actual.targetRank==golden.targetRank &&
         actual.conditionerRank==golden.conditionerRank && actual.searchRank==golden.searchRank,
         "frame catalogue/rank order changed");
    need(actual.imageGenerators==golden.imageGenerators && actual.projector==golden.projector &&
         actual.offsets==golden.offsets,"frame image/projector/offset coordinates changed");
    if(!actual.valid || !actual.affine || !golden.affine)return;
    const auto& a=*actual.affine;const auto& b=*golden.affine;
    need(a.valid==b.valid && a.particularSolution==b.particularSolution &&
         a.deterministicBasis==b.deterministicBasis && a.deterministicValues==b.deterministicValues &&
         a.kernelBasis==b.kernelBasis && a.quotientProjector==b.quotientProjector &&
         a.smithInvariants==b.smithInvariants,"complete affine golden changed");
    const ZhangExactMatrix numeric=zhangExactIdentityMatrix(targets.size());
    const auto ag=zhangR51PhysicalImageFromFrame(targets,numeric,actual);
    const auto bg=zhangR51PhysicalImageFromFrame(targets,numeric,golden);
    need(ag.valid==bg.valid && ag.reason==bg.reason && ag.particularTarget==bg.particularTarget &&
         ag.generators==bg.generators && ag.projector==bg.projector && ag.offsets==bg.offsets,
         "physical image target/generators/projector/offset changed");
}
}
int main(int argc,char** argv)
{
    std::size_t hnfCases=0,imageCases=0;
    std::mt19937 rng(0x50443431);
    compareHnf({}, {}, true);++hnfCases;
    compareHnf({{0,0},{0,0}}, {0,1}, true);++hnfCases;
    compareHnf({{1,0},{0}}, {0,0}, true);++hnfCases;
    compareHnf({{1,0},{0,1}}, {0}, true);++hnfCases;
    for(int trial=0;trial<500;++trial)
    {
        const int dimension=1+rng()%11,count=1+rng()%13;
        ZhangExactMatrix h(count,ZhangExactVector(dimension));
        ZhangExactVector x(dimension),values(count);
        for(auto& a:x)a=int(rng()%11)-5;
        for(int r=0;r<count;++r)for(int c=0;c<dimension;++c)
        {
            if(rng()%4==0)h[r][c]=int(rng()%17)-8;
            if(trial%17==0 && r==0 && c==0)h[r][c]=ZhangExactInteger(1)<<103;
            values[r]+=h[r][c]*x[c];
        }
        if(count>1 && trial%3==0){h.back()=h.front();values.back()=values.front();}
        if(count>1 && trial%7==0){h.back()=h.front();values.back()=values.front()+1;}
        for(bool witnesses:{false,true}){compareHnf(h,values,witnesses);++hnfCases;}
        std::reverse(h.begin(),h.end());std::reverse(values.begin(),values.end());
        compareHnf(h,values,true);++hnfCases;
    }
    // Force source-pivot replacement and both forward/back elimination.
    compareHnf({{21,0,0,7},{15,0,2,0},{6,3,0,0},{0,0,1,2}}, {7,3,4,5},true);++hnfCases;
    ZhangExactMatrix sparse(64,ZhangExactVector(2048));ZhangExactVector rhs(64);
    for(int r=0;r<64;++r){sparse[r][r]=1;sparse[r][100+r]=2;sparse[r][2047]=-1;rhs[r]=r-32;}
    ZhangExactRowHnfWorkMetrics sparseMetrics;
    const auto sparseActual=zhangExactRowHermiteNormalForm(sparse,rhs,true,&sparseMetrics);
    const auto sparseGolden=denseGoldenHnf(sparse,rhs,true);
    need(sparseActual.basis==sparseGolden.basis && sparseActual.values==sparseGolden.values &&
         sparseActual.rowTransform==sparseGolden.rowTransform,"sparse large-width golden");
    // Add shared first-column entries to require actual row additions.
    for(int r=1;r<64;++r)sparse[r][0]=r%3+1;
    zhangExactRowHermiteNormalForm(sparse,rhs,false,&sparseMetrics);
    need(sparseMetrics.rowAdditions>0 && sparseMetrics.sourceIndexBuilds>0,"sparse fixture did not execute additions");
    need(sparseMetrics.coefficientUpdates*16<sparseMetrics.denseCoefficientUpdates,"sparse arithmetic did not avoid dense updates");
    compareHnf(sparse,rhs,true);++hnfCases;
    for(int trial=0;trial<120;++trial)
    {
        const int d=1+rng()%5,rows=rng()%5,targetsCount=1+rng()%3;
        ZhangExactMatrix h(rows,ZhangExactVector(d)),t(targetsCount,ZhangExactVector(d));
        ZhangExactVector x(d),v(rows);
        for(auto& a:x)a=int(rng()%7)-3;
        for(int r=0;r<rows;++r)for(int c=0;c<d;++c){h[r][c]=int(rng()%7)-3;v[r]+=h[r][c]*x[c];}
        for(auto& row:t)for(auto& a:row)a=int(rng()%7)-3;
        if(trial%2==0 && d>1){for(auto& row:t)row.back()=0;for(auto& row:h)row.back()=0;}
        if(rows>1 && trial%5==0){h.back()=h.front();v.back()=v.front();}
        if(rows>1 && trial%11==0){h.back()=h.front();v.back()=v.front()+1;}
        for(auto work:{ZhangProductFrameWork::INTEGER_PROJECTOR,ZhangProductFrameWork::IMAGE_GENERATORS_ONLY})
        {compareFrame(t,h,v,d,work);++imageCases;}
    }
    // Independent complete-set oracle: 2a-b=1 => b odd, and target 2c even.
    const auto parity=zhangR51PhysicalImage({{0,1,0},{0,0,2}},{{2,-1,0}},{1},{{1,0},{0,1}},3);
    need(parity.valid && parity.generators==ZhangExactMatrix{{2,0},{0,2}},"finite-index target lattice was saturated");
    need(parity.particularTarget[0]%2!=0 && parity.particularTarget[1]%2==0,"parity affine origin");
    for(int u=-9;u<=9;++u)for(int v=-8;v<=8;++v)
    {
        const bool oracle=(u%2!=0 && v%2==0);
        const bool image=((ZhangExactInteger(u)-parity.particularTarget[0])%2==0 &&
                          (ZhangExactInteger(v)-parity.particularTarget[1])%2==0);
        need(oracle==image,"enumerated complete integer image differs");
        if(oracle){const int a=(u+1)/2,c=v/2;need(2*a-u==1 && 2*c==v,"independent integer reconstruction");}
    }
    need(!zhangR51PhysicalImage({{1}},{{2}},{1},{{1}},1).valid,"2N=1 was accepted");
    const auto feasible=zhangExactAffineIntegerQuotient({{2}},{4},1,ZhangExactQuotientWork::PARTICULAR_ONLY);
    need(feasible.valid && feasible.particularSolution==ZhangExactVector{2},"integer feasibility vs rational entailment");
    need(!zhangIntegerRowLatticeContains({{2}},{1},false).contained,"rational consequence became integer-row witness");
    const auto unknown=zhangR51PhysicalImage({{0,1}},{{1,0}},{4},{{1}},2);
    need(unknown.valid && unknown.generators==ZhangExactMatrix{{1}} && unknown.projector.size()==1,
         "unknown physical variable was replaced by zero");
    // Both the full-width target view and an appended unknown physical column
    // exercise production verified-domain borrowing with unchanged coordinates.
    const auto domain=zhangR47CompileProductSearchFrame({{1,0}},{{1,-2}},{1},2,
        ZhangProductFrameWork::IMAGE_GENERATORS_ONLY);
    need(domain.valid,"verified fixture domain");
    for(const auto& targets:{ZhangExactMatrix{{1,0},{0,1}},ZhangExactMatrix{{1,0,0},{0,0,1}}})
    {
        const int d=targets.front().size();
        auto h=ZhangExactMatrix{{1,-2}};h.front().resize(d);
        const auto old=denseGoldenFrame(targets,h,{1},d,ZhangProductFrameWork::IMAGE_GENERATORS_ONLY);
        const auto view=zhangR47CompileTargetOnVerifiedDomain(targets,domain.affine,d);
        need(view.valid && view.columns==old.columns && view.imageGenerators==old.imageGenerators,
             "verified target identity view changed exact image");
        const auto numeric=zhangExactIdentityMatrix(targets.size());
        const auto a=zhangR51PhysicalImageFromFrame(targets,numeric,view);
        const auto b=zhangR51PhysicalImageFromFrame(targets,numeric,old);
        need(a.valid && b.valid && a.particularTarget==b.particularTarget && a.generators==b.generators &&
             a.projector==b.projector && a.offsets==b.offsets,"verified target PAR coordinates changed");
        ++imageCases;
    }
    std::ostringstream log;
    {ZhangPhaseTimer timer(log,"P4_TEST");
        zhangExactRowHermiteNormalForm({{1,0},{0,1}},{0,0});
        zhangR47CompileProductSearchFrame({{1,0},{0,1}},{{1,-2}},{1},2,
            ZhangProductFrameWork::IMAGE_GENERATORS_ONLY);
        zhangR47CompileTargetOnVerifiedDomain({{1,0},{0,1}},domain.affine,2);
    }
    need(log.str().find("ZHANG_EXACT_HNF_INPUT_RESULT input_rows=2 columns=2 input_nonzeros=2")!=std::string::npos,
         "per-input HNF metrics missing");
    need(log.str().find("output_nonzeros=2")!=std::string::npos,"HNF residual metrics missing");
    need(log.str().find("identity_borrowed=1 compact_copied_cells=0")!=std::string::npos,
         "full-width exact matrix borrowing not exercised");
    need(log.str().find("call_sequence=")!=std::string::npos &&
         log.str().find("phase=P4_TEST")!=std::string::npos,"HNF call scope correlation missing");
    ZhangExactRowHnfWorkMetrics bigMetrics;
    zhangExactRowHermiteNormalForm({{-(ZhangExactInteger(1)<<103)-5,0}}, {0},false,&bigMetrics);
    need(bigMetrics.statisticsAvailable && bigMetrics.inputMaxBits==104 && bigMetrics.outputMaxBits==104,
         "no-copy signed-magnitude bit length incorrect");
    // Finite deterministic budgets retain the old dense checkpoint sequence;
    // interrupted work cannot return a partial exact result as success.
    std::size_t cancelled=0;
    for(std::uint64_t limit=1;limit<250;limit+=7)
    {
        ZhangExactWorkBudget oldBudget,newBudget;
        oldBudget.maximumCheckpoints=newBudget.maximumCheckpoints=limit;
        ZhangExactRowHnf oldResult,newResult;
        const bool before=zhangRunWithinExactWorkBudget(oldBudget,[&](){return denseGoldenHnf(sparse,rhs,true);},oldResult);
        const bool after=zhangRunWithinExactWorkBudget(newBudget,[&](){return zhangExactRowHermiteNormalForm(sparse,rhs,true);},newResult);
        need(before==after && oldBudget.checkpoints==newBudget.checkpoints,
             "finite budget checkpoint sequence changed");
        if(!after)
        {
            ++cancelled;
            need(newBudget.exhausted && newResult.basis.empty() && newResult.rowTransform.empty(),
                 "finite budget exposed partial HNF");
        }
    }
    need(cancelled>0,"cancellation fixture did not cancel");
    ZhangExactWorkBudget expired;expired.deadline=ZhangExactWorkBudget::Clock::now()-std::chrono::seconds(1);
    ZhangExactRowHnf expiredResult;
    need(!zhangRunWithinExactWorkBudget(expired,[&](){return zhangExactRowHermiteNormalForm(sparse,rhs,true);},expiredResult) &&
         expired.exhausted && expiredResult.basis.empty(),"expired wall budget exposed partial HNF");
    std::cout<<"P4_EXACT_IMAGE_COPY_SPARSE_CONTRACT PASS hnf_cases="<<hnfCases
        <<" image_coordinate_cases="<<imageCases<<" sparse_coefficient_updates="<<sparseMetrics.coefficientUpdates
        <<" old_dense_updates="<<sparseMetrics.denseCoefficientUpdates<<" cancellations="<<cancelled<<"\n";
    if(argc>1 && std::string(argv[1])=="--benchmark")
    {
        // Kernel-only paired microtiming; never infer whole-network speedup
        // or an epoch RSS peak from these synthetic, in-process fixtures.
        auto benchmark=[&](const char* name,const ZhangExactMatrix& h,const ZhangExactVector& v)
        {
            std::vector<double> oldSeconds,newSeconds,profiledSeconds;
            for(int repetition=0;repetition<8;++repetition)
            {
                ZhangExactRowHnf golden,fresh,profiled;
                auto timed=[&](auto&& operation,std::vector<double>& times)
                {
                    const auto start=std::chrono::steady_clock::now();auto result=operation();
                    times.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());
                    return result;
                };
                auto oldWork=[&](){return denseGoldenHnf(h,v,true);};
                auto newWork=[&](){return zhangExactRowHermiteNormalForm(h,v,true);};
                if(repetition%2==0){golden=timed(oldWork,oldSeconds);fresh=timed(newWork,newSeconds);}
                else {fresh=timed(newWork,newSeconds);golden=timed(oldWork,oldSeconds);}
                ZhangExactRowHnfWorkMetrics metrics;
                profiled=timed([&](){return zhangExactRowHermiteNormalForm(h,v,true,&metrics);},profiledSeconds);
                need(golden.consistent==fresh.consistent && golden.basis==fresh.basis &&
                     golden.values==fresh.values && golden.rowTransform==fresh.rowTransform,
                     "benchmark old/new HNF mismatch");
                need(profiled.consistent==fresh.consistent && profiled.basis==fresh.basis &&
                     profiled.values==fresh.values && profiled.rowTransform==fresh.rowTransform,
                     "profiling changed HNF");
            }
            auto emit=[&](const char* mode,std::vector<double> times)
            {
                const double mean=std::accumulate(times.begin(),times.end(),0.0)/times.size();
                std::sort(times.begin(),times.end());
                std::cout<<"P4_HNF_KERNEL_MICROBENCH fixture="<<name<<" mode="<<mode
                    <<" input_rows="<<h.size()<<" columns="<<h.front().size()<<" repetitions="<<times.size()
                    <<" mean_seconds="<<std::setprecision(17)<<mean<<" min_seconds="<<times.front()
                    <<" median_seconds="<<(times[3]+times[4])/2<<" max_seconds="<<times.back()
                    <<" complete_basis_rhs_rowtransform_equal=1 network_speedup_proven=0\n";
            };
            emit("frozen_dense",oldSeconds);emit("sparse_no_metrics",newSeconds);emit("sparse_metrics",profiledSeconds);
        };
        benchmark("sparse_64x2048",sparse,rhs);
        ZhangExactMatrix dense(24,ZhangExactVector(32));ZhangExactVector denseRhs(24);
        for(int r=0;r<24;++r)for(int c=0;c<32;++c){dense[r][c]=int((r*7+c*11+r*c)%13)-6;denseRhs[r]+=dense[r][c]*(c%5-2);}
        benchmark("dense_24x32",dense,denseRhs);
        ZhangExactMatrix big(24,ZhangExactVector(256));ZhangExactVector bigRhs(24);
        for(int r=0;r<24;++r)
        {
            big[r][r]=(ZhangExactInteger(1)<<103)+r+1;big[r][128+r]=r+1;big[r][255]=-7;
            if(r)big[r][0]=r+3;
            bigRhs[r]=big[r][r]*(r%3-1)+big[r][128+r];
        }
        benchmark("sparse_103bit_24x256",big,bigRhs);
    }
}
