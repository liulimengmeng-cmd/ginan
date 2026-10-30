#include "common/zhangUserResidualQuotient.hpp"
#include "common/zhangR47ProductDomain.hpp"

#include <cstdlib>
#include <iostream>
#include <set>
#include <sstream>
#include <vector>

static unsigned assertions = 0;
#define CHECK(condition) do { ++assertions; if (!(condition)) { \
    std::cerr << "FAIL line " << __LINE__ << ": " #condition "\n"; return EXIT_FAILURE; \
} } while (false)

// Independent reference oracle: literal small-integer enumeration and ordinary
// dot products. It does not invoke HNF, lattice membership, affine quotient,
// numerical conditioning, or any other production core function.
static bool satisfies(const ZhangExactMatrix& rows, const ZhangExactVector& values,
                      const ZhangExactVector& point)
{
    for (std::size_t r = 0; r < rows.size(); ++r)
    {
        ZhangExactInteger sum = 0;
        for (std::size_t c = 0; c < point.size(); ++c) sum += rows[r][c] * point[c];
        if (sum != values[r]) return false;
    }
    return true;
}

static bool independentAffineEnumeration(const ZhangUserResidualQuotient& frame,
                                        const ZhangExactMatrix& h,
                                        const ZhangExactVector& rhs)
{
    const int ambient = frame.exact.ambientDimension;
    const int free = frame.exact.quotientRank;
    const int width = 7;
    long combinations = 1;
    for (int j = 0; j < ambient; ++j) combinations *= width;
    for (long code = 0; code < combinations; ++code)
    {
        long residue = code;
        ZhangExactVector n(ambient);
        for (int j = 0; j < ambient; ++j) { n[j] = residue % width - 3; residue /= width; }
        if (!satisfies(h, rhs, n)) continue;
        ZhangExactVector q(free);
        for (int j = 0; j < free; ++j)
        for (int k = 0; k < ambient; ++k)
            q[j] += frame.exact.quotientProjector[j][k] *
                (n[k] - frame.exact.particularSolution[k]);
        auto recovered = frame.exact.particularSolution;
        for (int j = 0; j < free; ++j)
        for (int k = 0; k < ambient; ++k)
            recovered[k] += frame.exact.kernelBasis[j][k] * q[j];
        if (recovered != n) return false;
    }
    combinations = 1;
    for (int j = 0; j < free; ++j) combinations *= width;
    for (long code = 0; code < combinations; ++code)
    {
        long residue = code;
        auto n = frame.exact.particularSolution;
        for (int j = 0; j < free; ++j)
        {
            const int value = residue % width - 3;
            residue /= width;
            for (int k = 0; k < ambient; ++k) n[k] += frame.exact.kernelBasis[j][k] * value;
        }
        if (!satisfies(h, rhs, n)) return false;
    }
    return true;
}

#ifdef EIGEN_USE_BLAS
extern "C" int __wrap_main()
#else
int main()
#endif
{
    VectorXd mean(2); mean << 0.25, 0.8;
    MatrixXd covariance(2, 2); covariance << 4, 1, 1, 9;
    const auto originalMean = mean;
    const auto originalCovariance = covariance;
    auto empty = zhangBuildUserResidualQuotient(mean, covariance, {}, {});
    CHECK(empty.valid);
    CHECK(empty.exact.deterministicRank == 0 && empty.exact.quotientRank == 2);
    CHECK((empty.quotientMean - mean).norm() == 0);
    CHECK((empty.quotientCovariance - covariance).norm() == 0);
    CHECK(independentAffineEnumeration(empty, {}, {}));

    const ZhangExactMatrix parityRows{{2, 1}};
    const ZhangExactVector parityRhs{1};
    auto parity = zhangBuildUserResidualQuotient(mean, covariance, parityRows, parityRhs);
    CHECK(parity.valid);
    CHECK(parity.exact.deterministicRank == 1 && parity.exact.quotientRank == 1);
    CHECK(independentAffineEnumeration(parity, parityRows, parityRhs));
    CHECK((mean - originalMean).norm() == 0);
    CHECK((covariance - originalCovariance).norm() == 0);
    // Hand-calculated Gaussian equality update: H=[2,1], S=29,
    // P H'=[9,11]', residual=-0.3. This reference does not call a conditioner.
    VectorXd referenceMean(2); referenceMean << 0.25 - 2.7 / 29, 0.8 - 3.3 / 29;
    MatrixXd referenceCovariance(2, 2);
    referenceCovariance << 4 - 81.0 / 29, 1 - 99.0 / 29,
                           1 - 99.0 / 29, 9 - 121.0 / 29;
    CHECK((parity.conditionalMean - referenceMean).norm() < 1e-12);
    CHECK((parity.conditionalCovariance - referenceCovariance).norm() < 1e-12);
    double referenceQMean = 0, referenceQCovariance = 0;
    for (int c = 0; c < 2; ++c)
    {
        const double l = parity.exact.quotientProjector[0][c].convert_to<double>();
        referenceQMean += l * (referenceMean(c) -
            parity.exact.particularSolution[c].convert_to<double>());
        for (int d = 0; d < 2; ++d)
            referenceQCovariance += l * referenceCovariance(c, d) *
                parity.exact.quotientProjector[0][d].convert_to<double>();
    }
    CHECK(std::abs(parity.quotientMean(0) - referenceQMean) < 1e-12);
    CHECK(std::abs(parity.quotientCovariance(0, 0) - referenceQCovariance) < 1e-12);
    CHECK(parity.minimumQuotientEigenvalue > 0);
    // The integer image has y odd, not arbitrary y in Z. Enumerating the
    // production origin/kernel must never introduce the inadmissible even y.
    for (int q = -7; q <= 7; ++q)
    {
        const auto y = parity.exact.particularSolution[1] + parity.exact.kernelBasis[0][1] * q;
        CHECK((y - 1) % 2 == 0);
    }
    // Independently audit the actual target-image interface used by the
    // shared WL/L1 search. The target y is 1+2Z, never a saturated Z image.
    const auto image = zhangR47CompileProductSearchFrame({{0, 1}}, parityRows, parityRhs, 2);
    CHECK(image.valid && image.searchRank == 1);
    CHECK(image.imageGenerators.size() == 1 && image.imageGenerators[0].size() == 1);
    CHECK(image.imageGenerators[0][0] == 2 || image.imageGenerators[0][0] == -2);
    const auto targetOrigin = image.affine.particularSolution[1];
    for (int x = -5; x <= 5; ++x)
    {
        const ZhangExactVector n{x, 1 - 2*x};
        ZhangExactInteger searchCoordinate = image.offsets[0];
        for (int c=0;c<2;++c) searchCoordinate += image.projector[0][c]*n[c];
        CHECK(targetOrigin + image.imageGenerators[0][0]*searchCoordinate == n[1]);
    }
    // Accepted integer search coordinate u=2 means projector*N=2-offset.
    // Reversing the offset sign would describe a different affine set.
    const auto acceptedImageUnion = zhangAuditUserIntegerUnion(
        {parityRows[0], image.projector[0]}, {1, 2-image.offsets[0]}, 2);
    CHECK(acceptedImageUnion.valid && acceptedImageUnion.rank == 2);
    for (int x=-8;x<=8;++x)
    {
        const ZhangExactVector n{x,1-2*x};
        const bool admitted = satisfies(acceptedImageUnion.rows,acceptedImageUnion.values,n);
        CHECK(admitted == (n[1] == targetOrigin + 2*image.imageGenerators[0][0]));
    }

    // Temporary candidate mask: N_masked remains a live latent integer.
    // Projecting from the FULL H preserves odd-only live candidates. Dropping
    // N_masked equality rows before projection would wrongly admit even ones.
    const auto maskedImage = zhangProjectUserCandidateIntegerImage(parity.exact, {{0,1}});
    CHECK(maskedImage.valid && maskedImage.candidateDimension == 1);
    CHECK(maskedImage.fixedRank == 0 && maskedImage.freeImageRank == 1);
    CHECK(maskedImage.fixedRows.empty() && maskedImage.fixedValues.empty());
    CHECK(maskedImage.imageGenerators[0][0] == 2 || maskedImage.imageGenerators[0][0] == -2);
    for (int live=-7;live<=7;++live)
    {
        bool exists = false;
        for(int masked=-8;masked<=8;++masked)
            if (2*masked+live==1) exists = true;
        const bool imageContains = (ZhangExactInteger(live)-maskedImage.particularTarget[0])%
            maskedImage.imageGenerators[0][0] == 0;
        CHECK(imageContains == exists);
        const auto admission = zhangAuditUserIntegerUnion({{2,1},{0,1}},{1,live},2);
        CHECK(admission.valid == exists);
    }
    // Free-domain particular origin 0 is NOT evidence that an unknown live
    // integer equals zero. The candidate image is all Z and its fixed rank 0.
    const auto unknownImage = zhangProjectUserCandidateIntegerImage(empty.exact,{{0,1}});
    CHECK(unknownImage.valid && unknownImage.fixedRank==0 && unknownImage.freeImageRank==1);
    CHECK(unknownImage.fixedRows.empty());
    CHECK(unknownImage.imageGenerators[0][0]==1 || unknownImage.imageGenerators[0][0]==-1);
    CHECK(zhangProjectUserCandidateIntegerImage(parity.exact,{}).valid);
    CHECK(!zhangProjectUserCandidateIntegerImage(parity.exact,{{0,1},{0,2}}).valid);
    CHECK(!zhangProjectUserCandidateIntegerImage(parity.exact,{{1}}).valid);
    const auto insufficientDomain = zhangExactAffineIntegerQuotient(
        {{2,1}},{1},2,ZhangExactQuotientWork::FEASIBILITY_ONLY);
    CHECK(insufficientDomain.valid);
    CHECK(!zhangProjectUserCandidateIntegerImage(insufficientDomain,{{0,1}}).valid);
    // Pure unimodular rechart u=masked+live, v=live. Then H=[2,-1]
    // instead of [2,1]. Both candidate images must describe precisely odd Z.
    MatrixXd rechart(2,2); rechart << 1,1,0,1;
    MatrixXd inverseRechart(2,2); inverseRechart << 1,-1,0,1;
    const VectorXd rechartedMean = rechart*mean;
    const MatrixXd rechartedCovariance = rechart*covariance*rechart.transpose();
    const auto rechartedDomain = zhangBuildUserResidualQuotient(
        rechartedMean,rechartedCovariance,{{2,-1}},{1});
    CHECK(rechartedDomain.valid);
    CHECK((inverseRechart*rechartedDomain.conditionalMean-parity.conditionalMean).norm()<1e-12);
    CHECK((inverseRechart*rechartedDomain.conditionalCovariance*inverseRechart.transpose()-
        parity.conditionalCovariance).norm()<1e-12);
    const auto rechartedImage = zhangProjectUserCandidateIntegerImage(rechartedDomain.exact,{{0,1}});
    CHECK(rechartedImage.valid && rechartedImage.fixedRank==maskedImage.fixedRank);
    CHECK(rechartedImage.freeImageRank==maskedImage.freeImageRank);
    CHECK(rechartedImage.imageGenerators==maskedImage.imageGenerators);
    CHECK((rechartedImage.particularTarget[0]-maskedImage.particularTarget[0])%2==0);
    for(int live=-7;live<=7;++live)
    {
        bool oldExists=false,newExists=false;
        for(int coordinate=-8;coordinate<=8;++coordinate)
        {
            oldExists=oldExists || 2*coordinate+live==1;
            newExists=newExists || 2*coordinate-live==1;
        }
        CHECK(oldExists==newExists);
    }

    auto noInteger = zhangBuildUserResidualQuotient(mean, covariance, {{2, 0}}, {1});
    CHECK(!noInteger.valid);
    CHECK(noInteger.reason == "AFFINE_QUOTIENT_NO_INTEGER_ORIGIN");
    CHECK(!zhangBuildUserResidualQuotient(mean, covariance, {{0, 0}}, {1}).valid);
    CHECK(!zhangBuildUserResidualQuotient(mean, covariance, {{1, 0}, {1, 0}}, {1, 2}).valid);
    CHECK(!zhangBuildUserResidualQuotient(mean, covariance, {{1}}, {1}).valid);
    CHECK(!zhangBuildUserResidualQuotient(mean, covariance, {{1, 0}}, {}).valid);

    auto nonprimitive = zhangBuildUserResidualQuotient(mean, covariance, {{2, 0}}, {4});
    CHECK(nonprimitive.valid);
    CHECK(std::abs(nonprimitive.conditionalMean(0) - 2) < 1e-12);
    CHECK(independentAffineEnumeration(nonprimitive, {{2, 0}}, {4}));
    // User L1 prerequisite is exact logical constancy of the candidate WL,
    // after full-domain integer feasibility, not integer-row membership. A
    // held nonprimitive 2*WL=6 really fixes WL=3 and has no free target image.
    const auto nonprimitiveWl = zhangBuildUserResidualQuotient(
        mean,covariance,{{2,0}},{6});
    CHECK(nonprimitiveWl.valid);
    const auto constantWlImage = zhangProjectUserCandidateIntegerImage(
        nonprimitiveWl.exact,{{1,0}});
    CHECK(constantWlImage.valid && constantWlImage.fixedRank==1);
    CHECK(constantWlImage.freeImageRank==0 && constantWlImage.particularTarget[0]==3);
    CHECK(constantWlImage.imageGenerators.size()==1 && constantWlImage.imageGenerators[0].empty());
    CHECK(satisfies(constantWlImage.fixedRows,constantWlImage.fixedValues,{3}));
    CHECK(!satisfies(constantWlImage.fixedRows,constantWlImage.fixedValues,{2}));
    // 2*WL=1 is rationally solvable but has no integer solution. It must not
    // authorize L1 merely because a rational elimination reports a constant.
    CHECK(!noInteger.exact.valid);
    CHECK(!zhangProjectUserCandidateIntegerImage(noInteger.exact,{{1,0}}).valid);
    auto duplicate = zhangBuildUserResidualQuotient(mean, covariance,
        {{-4, -2}, {2, 1}, {6, 3}}, {-2, 1, 3});
    CHECK(duplicate.valid);
    CHECK(duplicate.exact.deterministicRank == 1);
    CHECK((duplicate.conditionalMean - parity.conditionalMean).norm() < 1e-12);
    CHECK((duplicate.conditionalCovariance - parity.conditionalCovariance).norm() < 1e-12);
    CHECK(independentAffineEnumeration(duplicate, parityRows, parityRhs));

    VectorXd mean3(3); mean3 << .2, -.1, .6;
    MatrixXd p3 = MatrixXd::Identity(3, 3);
    const ZhangExactMatrix mixedRows{{2, 1, 1}, {0, 2, -2}};
    const ZhangExactVector mixedValues{2, 0};
    auto mixed = zhangBuildUserResidualQuotient(mean3, p3, mixedRows, mixedValues);
    CHECK(mixed.valid);
    CHECK(mixed.exact.quotientRank == 1);
    CHECK(independentAffineEnumeration(mixed, mixedRows, mixedValues));
    auto mixedReordered = zhangBuildUserResidualQuotient(mean3, p3,
        {{0, -2, 2}, {2, 3, -1}, {2, 1, 1}}, {0, 2, 2});
    CHECK(mixedReordered.valid);
    CHECK((mixed.conditionalMean - mixedReordered.conditionalMean).norm() < 1e-12);
    CHECK((mixed.conditionalCovariance - mixedReordered.conditionalCovariance).norm() < 1e-12);
    CHECK(independentAffineEnumeration(mixedReordered, mixedRows, mixedValues));
    // Candidate equality rank is independent from retained full-domain rank.
    // 2*m+a=1, a+b=3 fixes one candidate direction, but their free image still
    // has an odd a / even b congruence and must not be saturated.
    const auto maskedMixedDomain = zhangBuildUserResidualQuotient(mean3,p3,
        {{2,1,0},{0,1,1}},{1,3});
    CHECK(maskedMixedDomain.valid);
    const auto candidateMixed = zhangProjectUserCandidateIntegerImage(
        maskedMixedDomain.exact,{{0,1,0},{0,0,1}});
    CHECK(candidateMixed.valid && candidateMixed.candidateDimension==2);
    CHECK(candidateMixed.fixedRank==1 && candidateMixed.freeImageRank==1);
    CHECK(candidateMixed.fixedRows.size()==1 && candidateMixed.ambientFixedRows.size()==1);
    CHECK(satisfies(candidateMixed.fixedRows,candidateMixed.fixedValues,{1,2}));
    CHECK(!satisfies(candidateMixed.fixedRows,candidateMixed.fixedValues,{1,3}));
    CHECK(satisfies(candidateMixed.ambientFixedRows,candidateMixed.fixedValues,{0,1,2}));
    for(int m=-4;m<=4;++m)
    {
        const ZhangExactVector point{m,1-2*m,2+2*m};
        CHECK(satisfies(candidateMixed.ambientFixedRows,candidateMixed.fixedValues,point));
        CHECK((point[1]-candidateMixed.particularTarget[0])%
            candidateMixed.imageGenerators[0][0]==0);
    }

    // A nonprimitive quotient equation retains its divisibility check.
    CHECK(!zhangMapUserQuotientDecision(parity, {{2}}, {3}).valid);
    auto mapped = zhangMapUserQuotientDecision(parity, {{2}}, {2});
    CHECK(mapped.valid && mapped.independentAddedRank == 1);
    CHECK(mapped.united.rank == 2 && mapped.united.quotientDimension == 0);
    const auto integerAtOne = [&]() {
        auto n = parity.exact.particularSolution;
        for (int c = 0; c < 2; ++c) n[c] += parity.exact.kernelBasis[0][c];
        return n;
    }();
    CHECK(satisfies(mapped.ambientRows, mapped.ambientValues, integerAtOne));
    CHECK(satisfies(mapped.united.rows, mapped.united.values, integerAtOne));
    auto fixed = zhangBuildUserResidualQuotient(mean, covariance,
        mapped.united.rows, mapped.united.values);
    CHECK(fixed.valid && fixed.exact.quotientRank == 0);
    CHECK(fixed.quotientMean.size() == 0 && fixed.quotientCovariance.rows() == 0);
    CHECK(fixed.conditionalCovariance.norm() < 1e-12);
    const auto whollyFixedCandidate = zhangProjectUserCandidateIntegerImage(fixed.exact,{{1,0},{0,1}});
    CHECK(whollyFixedCandidate.valid && whollyFixedCandidate.fixedRank==2);
    CHECK(whollyFixedCandidate.freeImageRank==0 && whollyFixedCandidate.fixedRows.size()==2);
    CHECK(satisfies(whollyFixedCandidate.fixedRows,whollyFixedCandidate.fixedValues,integerAtOne));
    auto redundantDecision = zhangMapUserQuotientDecision(fixed, {{}}, {0});
    CHECK(redundantDecision.valid && redundantDecision.independentAddedRank == 0);
    CHECK(!zhangMapUserQuotientDecision(fixed, {{}}, {1}).valid);

    // Residual PD, already-conditioned held direction, and absence of jitter.
    VectorXd deterministicMean(2); deterministicMean << 3, .4;
    MatrixXd singularP = MatrixXd::Zero(2, 2); singularP(1, 1) = 2;
    auto previouslyConditioned = zhangBuildUserResidualQuotient(
        deterministicMean, singularP, {{1, 0}}, {3});
    CHECK(previouslyConditioned.valid && previouslyConditioned.exact.quotientRank == 1);
    CHECK(previouslyConditioned.conditionalCovariance(0, 0) == 0);
    CHECK(std::abs(previouslyConditioned.quotientCovariance(0, 0) - 2) < 1e-12);
    CHECK(!zhangBuildUserResidualQuotient(deterministicMean, singularP, {{1, 0}}, {4}).valid);
    auto uncertifiedNull = zhangBuildUserResidualQuotient(deterministicMean, singularP, {}, {});
    CHECK(!uncertifiedNull.valid);
    CHECK(uncertifiedNull.reason == "RESIDUAL_QUOTIENT_NOT_POSITIVE_DEFINITE");
    MatrixXd negativeP = MatrixXd::Identity(2, 2); negativeP(0, 0) = -1;
    CHECK(!zhangBuildUserResidualQuotient(mean, negativeP, {}, {}).valid);
    auto allHeld = zhangBuildUserResidualQuotient(mean, covariance, {{1, 1}, {1, -1}}, {1, 3});
    CHECK(allHeld.valid && allHeld.exact.quotientRank == 0);
    CHECK(std::abs(allHeld.conditionalMean(0) - 2) < 1e-12);
    CHECK(std::abs(allHeld.conditionalMean(1) + 1) < 1e-12);
    CHECK(independentAffineEnumeration(allHeld, {{1, 1}, {1, -1}}, {1, 3}));
    CHECK(!zhangAuditUserIntegerUnion({{2, 0}}, {1}, 2).valid);
    CHECK(zhangAuditUserIntegerUnion({}, {}, 0).valid);
    CHECK(!zhangAuditUserIntegerUnion({{}}, {1}, 0).valid);
    CHECK(zhangBuildUserResidualQuotient(VectorXd(0), MatrixXd(0, 0), {}, {}).valid);
    VectorXd largeOrigin(1); largeOrigin << 0;
    CHECK(!zhangBuildUserResidualQuotient(largeOrigin, MatrixXd::Identity(1, 1),
        {{1}}, {ZhangExactInteger(1) << 54}).valid);

    // Reference loss does not require throwing away all rows touching the
    // reference: (ref-a)=7, (ref-b)=3 retain b-a=4 via their difference.
    const ZhangExactMatrix starRows{{1, -1, 0}, {1, 0, -1}};
    const ZhangExactVector starValues{7, 3};
    auto restricted = zhangRestrictUserIntegerLedger(starRows, starValues, {1, 2});
    CHECK(restricted.valid && restricted.rows.size() == 1);
    CHECK(restricted.discardedRank == 1 && restricted.knownFullProjection);
    CHECK(satisfies(restricted.rows, restricted.values, {-2, 2}));
    CHECK(!satisfies(restricted.rows, restricted.values, {-2, 3}));
    CHECK(restricted.sourceCombinations.size() == restricted.rows.size());
    for (std::size_t r=0;r<restricted.rows.size();++r)
    {
        ZhangExactVector reconstructed(3);
        ZhangExactInteger rhs=0;
        for (std::size_t source=0;source<starRows.size();++source)
        {
            const auto coefficient=restricted.sourceCombinations[r][source];
            rhs+=coefficient*starValues[source];
            for(int c=0;c<3;++c) reconstructed[c]+=coefficient*starRows[source][c];
        }
        CHECK(reconstructed[0]==0);
        CHECK(reconstructed[1]==restricted.rows[r][0] && reconstructed[2]==restricted.rows[r][1]);
        CHECK(rhs==restricted.values[r]);
    }
    for (int a = -4; a <= 4; ++a)
    for (int b = -4; b <= 4; ++b)
    {
        bool existential = false;
        for (int reference = -12; reference <= 12; ++reference)
            if (satisfies(starRows, starValues, {reference, a, b})) existential = true;
        CHECK(satisfies(restricted.rows, restricted.values, {a, b}) == existential);
    }
    // A slip on a retires that physical coordinate but leaves ref-b intact.
    auto slippedA = zhangRestrictUserIntegerLedger(starRows, starValues, {0, 2});
    CHECK(slippedA.valid && slippedA.rows.size() == 1 && slippedA.knownFullProjection);
    CHECK(satisfies(slippedA.rows, slippedA.values, {4, 1}));
    CHECK(!satisfies(slippedA.rows, slippedA.values, {4, 2}));
    // Changing column order without retiring anything is a complete rechart.
    auto reordered = zhangRestrictUserIntegerLedger(starRows, starValues, {2, 0, 1});
    CHECK(reordered.valid && reordered.discardedRank == 0 && reordered.knownFullProjection);
    CHECK(satisfies(reordered.rows, reordered.values, {2, 5, -2}));
    CHECK(!satisfies(reordered.rows, reordered.values, {2, 5, -1}));
    // Explicit limitation: eliminating old from 2*old+live=1 leaves parity.
    // Equality-only restriction returns no equality and MUST NOT claim full
    // projection equivalence; caller may retain old as a latent ledger column.
    auto parityRetirement = zhangRestrictUserIntegerLedger({{2, 1}}, {1}, {1});
    CHECK(parityRetirement.valid && parityRetirement.rows.empty());
    CHECK(parityRetirement.equalityConsequencesOnly && !parityRetirement.knownFullProjection);
    CHECK(parityRetirement.discardedRank == 1);
    CHECK(satisfies(parityRetirement.rows, parityRetirement.values, {0}));
    bool oldExistsForEven = false;
    for (int old = -5; old <= 5; ++old)
        if (satisfies({{2, 1}}, {1}, {old, 0})) oldExistsForEven = true;
    CHECK(!oldExistsForEven);
    CHECK(!zhangRestrictUserIntegerLedger(starRows, starValues, {1, 1}).valid);
    CHECK(!zhangRestrictUserIntegerLedger(starRows, starValues, {-1}).valid);
    CHECK(!zhangRestrictUserIntegerLedger({{2}}, {1}, {0}).valid);
    auto noColumns = zhangRestrictUserIntegerLedger(starRows, starValues, {});
    CHECK(noColumns.valid && noColumns.rows.empty() && noColumns.discardedRank == 2);
    auto noRows = zhangRestrictUserIntegerLedger({}, {}, {1, 0}, 2);
    CHECK(noRows.valid && noRows.knownFullProjection);
    std::cout << "USER_RESIDUAL_QUOTIENT_CONTRACT assertions=" << assertions
              << " status=PASS independent_oracle=ENUMERATION_AND_HAND_GAUSSIAN\n";
    return EXIT_SUCCESS;
}
