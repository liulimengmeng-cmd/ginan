#include "ambres/GNSSambres.hpp"
#include <algorithm>
#include <math.h>
#include <iterator>
#include <limits>
#include <tuple>

#define LOG_PI 1.14472988584940017
#define SQRT2 1.41421356237309510
#define AMB_RANG 10

bool   AR_VERBO      = false;
double FIXED_AMB_VAR = 1e-8;

ReceiverAmbiguityTransform buildReceiverAmbiguityIntegerTransform(
    const GinAR_mtx&        ambiguityResolution,
    const map<E_Sys, bool>& receiverAmbiguityPivot
)
{
    using GroupKey = std::tuple<string, E_Sys, int>;

    ReceiverAmbiguityTransform result;
    const int columnCount = ambiguityResolution.aflt.size();

    map<GroupKey, vector<int>> groups;
    for (const auto& [localIndex, key] : ambiguityResolution.ambmap)
    {
        if (localIndex < 0 || localIndex >= columnCount)
        {
            continue;
        }
        groups[{key.str, key.Sat.sys, key.num}].push_back(localIndex);
    }

    int rowCount = 0;
    bool anySingleDifferenceEnabled = false;
    for (const auto& [group, members] : groups)
    {
        const E_Sys system = std::get<1>(group);
        const auto enabled = receiverAmbiguityPivot.find(system);
        const bool singleDifference =
            enabled != receiverAmbiguityPivot.end() && enabled->second;
        anySingleDifferenceEnabled = anySingleDifferenceEnabled || singleDifference;
        rowCount += singleDifference
            ? std::max(0, static_cast<int>(members.size()) - 1)
            : static_cast<int>(members.size());
    }

    if (!anySingleDifferenceEnabled)
    {
        result.matrix = MatrixXd::Identity(columnCount, columnCount);
    }
    else
    {
        result.matrix = MatrixXd::Zero(rowCount, columnCount);
    }

    int row = 0;
    for (const auto& [group, members] : groups)
    {
        const auto& [receiver, system, observation] = group;
        const auto enabled = receiverAmbiguityPivot.find(system);
        const bool singleDifference =
            enabled != receiverAmbiguityPivot.end() && enabled->second;

        ReceiverAmbiguityDatum datum;
        datum.receiver = receiver;
        datum.system = system;
        datum.observation = observation;
        datum.memberCount = members.size();
        datum.singleDifferenced = singleDifference;

        if (!singleDifference)
        {
            result.identityGroupCount++;
            result.groups.push_back(datum);
            if (anySingleDifferenceEnabled)
            {
                for (int member : members)
                {
                    result.matrix(row++, member) = 1;
                }
            }
            continue;
        }

        if (members.size() < 2)
        {
            result.droppedSingletonGroupCount++;
            if (!members.empty())
            {
                datum.pivot = ambiguityResolution.ambmap.at(members.front()).Sat;
            }
            result.groups.push_back(datum);
            continue;
        }

        const int pivot = *std::min_element(
            members.begin(),
            members.end(),
            [&](int left, int right)
            {
                const double leftVariance = ambiguityResolution.Paflt(left, left);
                const double rightVariance = ambiguityResolution.Paflt(right, right);
                const bool leftFinite = std::isfinite(leftVariance);
                const bool rightFinite = std::isfinite(rightVariance);
                if (leftFinite != rightFinite)
                {
                    return leftFinite;
                }
                if (leftFinite && leftVariance != rightVariance)
                {
                    return leftVariance < rightVariance;
                }
                return ambiguityResolution.ambmap.at(left).Sat <
                       ambiguityResolution.ambmap.at(right).Sat;
            }
        );

        datum.pivot = ambiguityResolution.ambmap.at(pivot).Sat;
        result.singleDifferencedGroupCount++;
        result.groups.push_back(datum);

        for (int member : members)
        {
            if (member == pivot)
            {
                continue;
            }
            result.matrix(row, member) = +1;
            result.matrix(row, pivot) = -1;
            row++;
        }
    }

    return result;
}

bool mapIntegerAmbiguityConstraintsToOriginalState(
    GinAR_mtx&                        integerAmbiguityResolution,
    const ReceiverAmbiguityTransform& integerTransform,
    const map<int, KFKey>&            originalAmbiguityMap
)
{
    const int resolvedCombinationCount = integerAmbiguityResolution.Ztrs.rows();
    const int integerCoordinateCount    = integerTransform.matrix.rows();
    const int originalAmbiguityCount    = integerTransform.matrix.cols();

    if (integerAmbiguityResolution.Ztrs.cols() != integerCoordinateCount ||
        integerAmbiguityResolution.zfix.size() != resolvedCombinationCount ||
        originalAmbiguityMap.size() != originalAmbiguityCount)
    {
        return false;
    }

    for (int column = 0; column < originalAmbiguityCount; column++)
    {
        if (originalAmbiguityMap.find(column) == originalAmbiguityMap.end())
        {
            return false;
        }
    }

    integerAmbiguityResolution.Ztrs =
        (integerAmbiguityResolution.Ztrs * integerTransform.matrix).eval();
    integerAmbiguityResolution.ambmap = originalAmbiguityMap;

    return true;
}

ConditionalIntegerComplement buildConditionalIntegerComplement(
    const GinAR_mtx& inputAmbiguities,
    const GinAR_mtx& acceptedPartialResolution
)
{
    ConditionalIntegerComplement result;
    const int inputCount = inputAmbiguities.aflt.size();
    const int fixedCount = acceptedPartialResolution.zfix.size();
    const int complementCount = inputCount - fixedCount;

    if (inputCount <= 0 || inputAmbiguities.Paflt.rows() != inputCount ||
        inputAmbiguities.Paflt.cols() != inputCount)
    {
        result.diagnosticStatus = "INVALID_INPUT_DIMENSIONS";
        return result;
    }
    if (fixedCount <= 0 || complementCount <= 0 ||
        acceptedPartialResolution.Ztrs.rows() != fixedCount ||
        acceptedPartialResolution.Ztrs.cols() != inputCount ||
        acceptedPartialResolution.fullDecorrelatedTransform.rows() != inputCount ||
        acceptedPartialResolution.fullDecorrelatedTransform.cols() != inputCount)
    {
        result.diagnosticStatus = "INVALID_PARTIAL_RESOLUTION_DIMENSIONS";
        return result;
    }

    const MatrixXd& fullTransform =
        acceptedPartialResolution.fullDecorrelatedTransform;
    const MatrixXd expectedFixedTransform = fullTransform.bottomRows(fixedCount);
    if (!acceptedPartialResolution.Ztrs.isApprox(expectedFixedTransform, 1e-12))
    {
        // Common-set LAMBDA may select a non-contiguous row subset.  It needs
        // a separate integer-complement construction and is not silently
        // treated as the standard partial-LAMBDA tail.
        result.diagnosticStatus = "FIXED_ROWS_NOT_DECORRELATED_TAIL";
        return result;
    }

    const VectorXd transformedFloat = fullTransform * inputAmbiguities.aflt;
    const MatrixXd transformedCovariance =
        fullTransform * inputAmbiguities.Paflt * fullTransform.transpose();
    const MatrixXd complementFixedCovariance = transformedCovariance.block(
        0,
        complementCount,
        complementCount,
        fixedCount
    );
    const MatrixXd fixedCovariance = transformedCovariance.bottomRightCorner(
        fixedCount,
        fixedCount
    );
    const LDLT<MatrixXd> fixedSolver(fixedCovariance);
    if (fixedSolver.info() != Eigen::Success || !fixedSolver.isPositive())
    {
        result.diagnosticStatus = "FIXED_COVARIANCE_NOT_POSITIVE_DEFINITE";
        return result;
    }

    const VectorXd fixedInnovation =
        acceptedPartialResolution.zfix - transformedFloat.tail(fixedCount);
    result.ambiguityResolution.aflt =
        transformedFloat.head(complementCount) +
        complementFixedCovariance * fixedSolver.solve(fixedInnovation);
    result.ambiguityResolution.Paflt =
        transformedCovariance.topLeftCorner(complementCount, complementCount) -
        complementFixedCovariance *
            fixedSolver.solve(complementFixedCovariance.transpose());
    result.ambiguityResolution.Paflt =
        0.5 * (result.ambiguityResolution.Paflt +
               result.ambiguityResolution.Paflt.transpose()).eval();
    result.transformToInputCoordinates = fullTransform.topRows(complementCount);

    if (!result.ambiguityResolution.aflt.allFinite() ||
        !result.ambiguityResolution.Paflt.allFinite())
    {
        result.ambiguityResolution = {};
        result.transformToInputCoordinates.resize(0, 0);
        result.diagnosticStatus = "NONFINITE_CONDITIONAL_SOLUTION";
        return result;
    }

    result.diagnosticStatus = "CONDITIONAL_INTEGER_COMPLEMENT_READY";
    return result;
}

/** Probability of error (assuming normal distribution) */
double round_perr(
    double dx,  ///< Distance between value and mean
    double var  ///< Variance
)
{
    if (var < 1e-20)
        return 0;

    double p0   = 0;
    double fact = -0.25 / var;

    for (int i = 1; i < AMB_RANG; i++)
    {
        p0 += exp((i + 2 * dx) * i * fact);
        p0 += exp((i - 2 * dx) * i * fact);
    }

    return p0 / (p0 + 1);
}

/** Simple integer Rounding */
int simple_round(
    Trace&     trace,  ///< Debug trace
    GinAR_mtx& mtrx,   ///< Reference to structure containing float values and covariance
    GinAR_opt  opt     ///< Object containing processing options
)
{
    MatrixXd P    = mtrx.Paflt;
    VectorXd ret  = mtrx.aflt;
    int      namb = ret.size(), nfix = 0;

    mtrx.Ztrs.resize(0, 0);
    mtrx.zfix.resize(0);

    if (namb <= 0)
    {
        mtrx.diagnosticStatus = "NO_ELIGIBLE_AMBIGUITIES";
        return 0;
    }

    double ratthr = 1 / (opt.ratthr + 1);
    double sucthr = 1 - pow(opt.sucthr, 1.0 / namb);
    tracepdeex(4, trace, "\n#ARES_RND Using integer rounding ... %.4e  %.4f", sucthr, ratthr);

    vector<int> zind;
    vector<int> xind;
    xind.reserve(namb);

    for (int i = 0; i < namb; i++)
    {
        xind.push_back(i);
        double dv   = ret(i) - ROUND(ret(i));
        double perr = round_perr(dv, P(i, i));

        if (fabs(dv) < ratthr && perr < sucthr)
        {
            ret(i) = ret(i) - dv;
            nfix++;
            zind.push_back(i);
        }
    }

    MatrixXd Z = MatrixXd::Identity(namb, namb);

    if (nfix == 0)
    {
        mtrx.diagnosticStatus = "ROUND_VALIDATION_FAILED";
        return 0;
    }

    mtrx.Ztrs = Z(zind, xind);
    mtrx.zfix = ret(zind);
    mtrx.selectedDecorrelatedAmbiguityCount = nfix;
    mtrx.diagnosticStatus = "RESOLVED_BY_ROUNDING";

    return nfix;
}

/** Iterative Rounding */
int interat_round(
    Trace&     trace,  ///< Debug trace
    GinAR_mtx& mtrx,   ///< Reference to structure containing float values and covariance
    GinAR_opt& opt     ///< Object containing processing options
)
{
    MatrixXd P    = mtrx.Paflt;
    VectorXd x    = mtrx.aflt;
    int      namb = x.size();

    mtrx.Ztrs.resize(0, 0);
    mtrx.zfix.resize(0);

    if (namb <= 0)
    {
        mtrx.diagnosticStatus = "NO_ELIGIBLE_AMBIGUITIES";
        return 0;
    }

    double sucthr = 1 - pow(opt.sucthr, 1.0 / namb);
    double ratthr = 1 / (opt.ratthr + 1);

    vector<int> zind;
    vector<int> xind;
    xind.reserve(namb);

    for (int i = 0; i < namb; i++)
        xind.push_back(i);

    MatrixXd I = MatrixXd::Identity(namb, namb);
    MatrixXd Ztrs;
    VectorXd xfix;
    VectorXd dvvct = VectorXd::Zero(namb);

    int nfix = 0;
    int nnew = 0;

    for (int iter = 0; iter < opt.nitr; iter++)
    {
        zind.clear();
        nnew = 0;

        for (int i = 0; i < namb; i++)
        {
            double dv   = x(i) - ROUND(x(i));
            double perr = round_perr(dv, P(i, i));

            if ((fabs(dv) < ratthr) && (perr < sucthr))
            {
                nnew++;
                zind.push_back(i);
                dvvct(i) = dv;
            }
        }

        if (nnew <= nfix)
            break;

        VectorXd dx = dvvct(zind);

        Ztrs = I(zind, xind);
        xfix = x(zind) - dx;

        MatrixXd Psel = P(zind, zind) + (FIXED_AMB_VAR * I(zind, zind));
        MatrixXd K    = P(xind, zind) * Psel.inverse();
        MatrixXd S    = K * P(zind, xind);
        x             = x - K * dx;
        P             = P - S;
        nfix          = nnew;
    }

    if (nfix == 0)
    {
        mtrx.diagnosticStatus = "ITERATIVE_ROUND_VALIDATION_FAILED";
        return 0;
    }

    mtrx.Ztrs = Ztrs;
    mtrx.zfix = xfix;
    mtrx.selectedDecorrelatedAmbiguityCount = nfix;
    mtrx.diagnosticStatus = "RESOLVED_BY_ITERATIVE_ROUNDING";

    return nfix;
}

bool LTDL_factorization(
    GinAR_mtx& mtrx  ///< Reference to structure containing float values and covariance
)
{
    int      n = mtrx.aflt.size();
    MatrixXd P = mtrx.Paflt;

    MatrixXd L = MatrixXd::Zero(n, n);
    VectorXd D = VectorXd::Zero(n);

    for (int i = n - 1; i >= 0; i--)
    {
        if (P(i, i) <= 0)
            return false;

        D(i)                    = P(i, i);
        double a                = sqrt(P(i, i));
        L.block(i, 0, 1, i + 1) = P.block(i, 0, 1, i + 1) / a;
        for (int j = 0; j < i; j++)
            P.block(j, 0, 1, j + 1) -= L(i, j) * L.block(i, 0, 1, j + 1);
        L.block(i, 0, 1, i + 1) /= L(i, i);
    }

    mtrx.Ltrs = L;
    mtrx.Dtrs = D;
    return true;
}

/** Lambda decorrelation (trough Z transform) */
int Ztrans_reduction(
    Trace&     trace,  ///< Debug trace
    GinAR_mtx& mtrx    ///< Reference to structure containing float values and covariance
)
{
    int n = mtrx.aflt.size();
    if (n < 1)
        return -1;

    if (LTDL_factorization(mtrx) == false)
    {
        tracepdeex(
            1,
            trace,
            "WARNING: LD decomposition error, ambiguity matrix may not positive definite\n"
        );
        return -1;
    }

    VectorXd x = mtrx.aflt;
    VectorXd D = mtrx.Dtrs;
    MatrixXd L = mtrx.Ltrs;
    MatrixXd Z = MatrixXd::Identity(n, n);

    if (AR_VERBO)
    {
        trace << std::setprecision(8);
        trace << "\n"
              << "x =" << "\n"
              << x.transpose() << "\n";
        trace << "\n"
              << "Px=" << "\n"
              << mtrx.Paflt << "\n";
        trace << "\n"
              << "Lx=" << "\n"
              << L << "\n";
        trace << "\n"
              << "Dx=" << "\n"
              << D.transpose() << "\n";
    }

    int k = n - 2;
    int j = n - 2;
    while (j >= 0)
    {
        if (j <= k)
            for (int i = j + 1; i < n; i++)
            {
                double mu = ROUND(L(i, j));
                if (mu != 0)
                {
                    L.col(j) -= mu * L.col(i);
                    Z.col(j) -= mu * Z.col(i);
                }
            }

        double del = D(j) + L(j + 1, j) * L(j + 1, j) * D(j + 1);
        if ((del + 1E-6) < D(j + 1))
        {
            double eta = D(j) / del;
            double lam = D(j + 1) * L(j + 1, j) / del;

            D(j)     = eta * D(j + 1);
            D(j + 1) = del;

            MatrixXd a0             = L.block(j, 0, 1, j);
            MatrixXd a1             = L.block(j + 1, 0, 1, j);
            L.block(j, 0, 1, j)     = a1 - L(j + 1, j) * a0;
            L.block(j + 1, 0, 1, j) = lam * a1 + eta * a0;
            L(j + 1, j)             = lam;

            VectorXd Ltmp                       = L.block(j + 2, j, n - j - 2, 1);
            L.block(j + 2, j, n - j - 2, 1)     = L.block(j + 2, j + 1, n - j - 2, 1);
            L.block(j + 2, j + 1, n - j - 2, 1) = Ltmp;

            VectorXd Ztmp = Z.col(j);
            Z.col(j)      = Z.col(j + 1);
            Z.col(j + 1)  = Ztmp;

            k = j;
            j = n - 2;
        }
        else
            j--;
    }

    mtrx.Ztrs = Z.transpose();
    mtrx.zflt = mtrx.Ztrs * x;
    mtrx.Ltrs = L;
    mtrx.Dtrs = D;

    if (AR_VERBO)
    {
        trace << std::setprecision(8);
        trace << "\n"
              << "z =" << "\n"
              << mtrx.zflt.transpose() << "\n";
        trace << "\n"
              << "Zt=" << "\n"
              << mtrx.Ztrs << "\n";
        trace << "\n"
              << "Lz=" << "\n"
              << L << "\n";
        trace << "\n"
              << "Dz=" << "\n"
              << D.transpose() << "\n";
    }

    return n;
}

/** Integer bootstrapping */
// int integer_bootst(
// 	Trace& trace,		///< Debug trace
// 	GinAR_mtx& mtrx,	///< Reference to structure containing float values and covariance
// 	GinAR_opt opt)		///< Object containing processing options
// {
// 	int info = Ztrans_reduction(trace, mtrx);

// 	if (info < 0)
// 		return 0;

// 	GinAR_mtx mtrx2;
// 	mtrx2.aflt = mtrx.zflt;

// 	MatrixXd Z   = mtrx.Ztrs;
// 	mtrx2.Paflt	 = Z*mtrx.Paflt*Z.transpose();

// 	int nfix = interat_round (trace, mtrx2, opt);

// 	mtrx.Ztrs = mtrx2.Ztrs * Z;
// 	mtrx.zfix = mtrx2.zfix;

// 	return nfix;
// }

int integer_bootst(
    Trace&     trace,  ///< Debug trace
    GinAR_mtx& mtrx,   ///< Reference to structure containing float values and covariance
    GinAR_opt& opt     ///< Object containing processing options
)
{
    LDLT<MatrixXd> ldlt_;
    ldlt_.compute(mtrx.Paflt);

    if (ldlt_.isPositive() == false)
    {
        mtrx.diagnosticStatus = "COVARIANCE_NOT_POSITIVE_DEFINITE";
        tracepdeex(
            1,
            trace,
            "WARNING: LD decomposition error, ambiguity matrix may not positive definite\n"
        );
        return 0;
    }

    MatrixXd L_ = ldlt_.matrixL();
    auto     tr = ldlt_.transpositionsP();

    int      siz = mtrx.aflt.size();
    MatrixXd I0  = MatrixXd::Identity(siz, siz);
    MatrixXd Zt  = tr * I0;
    VectorXd z_  = tr * mtrx.aflt;

    for (int j = siz - 2; j >= 0; j--)
    {
        for (int i = j + 1; i < siz; i++)
        {
            double mu = ROUND(L_(i, j));

            if (mu != 0)
            {
                L_.row(i) -= mu * L_.row(j);
                Zt.row(i) -= mu * Zt.row(j);
                z_(i) -= mu * z_(j);
            }
        }
    }

    MatrixXd Pz = Zt * mtrx.Paflt * Zt.transpose();

    GinAR_mtx mtrx2;
    mtrx2.aflt  = z_;
    mtrx2.Paflt = Pz;

    if (AR_VERBO)
    {
        trace << "\n"
              << "x_=" << "\n"
              << mtrx.aflt.transpose() << "\n";
        trace << "\n"
              << "Px=" << "\n"
              << mtrx.Paflt << "\n";
        trace << "\n"
              << "Zt=" << "\n"
              << Zt << "\n";
        trace << "\n"
              << "z_=" << "\n"
              << z_.transpose() << "\n";
        trace << "\n"
              << "Pz=" << "\n"
              << Pz << "\n";
    }

    int nfix = interat_round(trace, mtrx2, opt);

    mtrx.Ztrs = mtrx2.Ztrs * Zt;
    mtrx.zfix = mtrx2.zfix;
    mtrx.diagnosticStatus = mtrx2.diagnosticStatus;
    mtrx.selectedDecorrelatedAmbiguityCount =
        mtrx2.selectedDecorrelatedAmbiguityCount;

    return nfix;
}

/** Lambda algorithm and its variations (ILQ, Common set, BIE) */
int lambda_search(
    Trace&     trace,  ///< Debug trace
    GinAR_mtx& mtrx,   ///< Reference to structure containing float values and covariance
    GinAR_opt  opt     ///< Object containing processing options
)
{
    int info = Ztrans_reduction(trace, mtrx);

    if (info < 0)
    {
        mtrx.diagnosticStatus = "DECORRELATION_FAILED";
        tracepdeex(2, trace, "\n Matrix decorrelation failed ... ");
        return 0;
    }

    int nmax = mtrx.Dtrs.size();
    mtrx.fullDecorrelatedTransform = mtrx.Ztrs;
    int k    = nmax - 1;
    int kmax = k;

    double succ = erf(sqrt(1 / (8 * mtrx.Dtrs(k--))));
    mtrx.bootstrappedSuccessRate = succ;

    if (succ < opt.sucthr)
    {
        mtrx.selectedDecorrelatedAmbiguityCount = 1;
        mtrx.diagnosticStatus = "SUCCESS_RATE_BELOW_THRESHOLD";
        return 0;
    }

    int zsiz = 1;

    while (k >= 0)
    {
        const double nextSucc = succ * erf(sqrt(1 / (8 * mtrx.Dtrs(k--))));
        if (nextSucc < opt.sucthr)
            break;
        succ = nextSucc;
        zsiz++;
    }
    mtrx.bootstrappedSuccessRate = succ;
    mtrx.selectedDecorrelatedAmbiguityCount = zsiz;
    if (zsiz < opt.minimumDecorrelatedAmbiguityCount)
    {
        mtrx.diagnosticStatus = "INSUFFICIENT_DECORRELATED_AMBIGUITIES";
        return 0;
    }

    int kmin = kmax - zsiz + 1;

    // Keep equal-distance candidates distinct.  A tie is meaningful for the
    // ratio test (ratio == 1) and must not be overwritten by map::operator[].
    multimap<double, VectorXd> zfixList;

    // Plain LAMBDA needs only the best candidate and the ratio test needs
    // exactly the best two.  Common-set and BIE modes use the configured pool,
    // with two as their minimum useful size.
    const int candidateLimit = opt.mode == E_ARmode::LAMBDA
        ? 1
        : opt.mode == E_ARmode::LAMBDA_ALT ? 2 : std::max(2, opt.nset);
    const bool retainCutoffTies =
        opt.mode == E_ARmode::LAMBDA_AL2 || opt.mode == E_ARmode::LAMBDA_BIE;
    const auto distanceTolerance = [](double distance)
    {
        return 64 * std::numeric_limits<double>::epsilon() *
               std::max(1.0, std::abs(distance));
    };

    MatrixXd L    = mtrx.Ltrs;
    VectorXd D    = mtrx.Dtrs;
    VectorXd zflt = mtrx.zflt;

    VectorXd dist = VectorXd::Zero(nmax);
    VectorXd zadj = VectorXd::Zero(nmax);
    VectorXd zfix = VectorXd::Zero(nmax);
    VectorXd zdif = VectorXd::Zero(nmax);
    VectorXd step = VectorXd::Zero(nmax);

    k              = kmax;
    zadj(k)        = zflt(k);
    zfix(k)        = ROUND(zadj(k));
    zdif(k)        = zadj(k) - zfix(k);
    step(k)        = zdif(k) < 0 ? -1 : 1;
    bool   search  = true;
    double maxdist = 1e99;
    int    ncand   = 0;

    while (search)
    {
        double newdist = dist(k) + zdif(k) * zdif(k) / D(k);

        if (newdist <= maxdist + distanceTolerance(maxdist))
        {
            if (k != kmin)
            {
                k--;
                dist(k) = newdist;

                zadj(k) = zflt(k);
                for (int j = k + 1; j < nmax; j++)
                    zadj(k) -= zdif(j) * L(j, k);

                zfix(k) = ROUND(zadj(k));
                zdif(k) = zadj(k) - zfix(k);
                step(k) = zdif(k) < 0 ? -1 : 1;
            }
            else
            {
                VectorXd zcut     = zfix.tail(zsiz);
                zfixList.emplace(newdist, zcut);
                ncand             = zfixList.size();
                double maxd       = newdist * opt.ratthr;

                if (ncand > 1 && maxd < maxdist)
                    maxdist = maxd;

                if (ncand >= candidateLimit)
                {
                    auto cutoff = zfixList.begin();
                    std::advance(cutoff, candidateLimit - 1);
                    const double cutoffDistance = cutoff->first;

                    if (retainCutoffTies)
                    {
                        const double cutoffTolerance = distanceTolerance(cutoffDistance);
                        for (auto it = std::next(cutoff); it != zfixList.end();)
                        {
                            if (it->first <= cutoffDistance + cutoffTolerance)
                            {
                                ++it;
                            }
                            else
                            {
                                it = zfixList.erase(it);
                            }
                        }
                    }
                    else
                    {
                        zfixList.erase(std::next(cutoff), zfixList.end());
                    }

                    maxd = cutoffDistance;

                    if (maxd < maxdist)
                        maxdist = maxd;

                    ncand = zfixList.size();
                }

                zfix(kmin) += step(kmin);
                zdif(kmin) = zadj(kmin) - zfix(kmin);
                step(kmin) = -step(kmin) + (step(kmin) < 0 ? 1 : -1);
            }
        }
        else
        {
            if (k == kmax)
                break;
            else
            {
                k++;
                zfix(k) += step(k);
                zdif(k) = zadj(k) - zfix(k);
                step(k) = -step(k) + (step(k) < 0 ? 1 : -1);
            }
        }
    }

    if (zfixList.size() < 1)
    {
        mtrx.diagnosticStatus = "NO_INTEGER_CANDIDATES";
        return 0;
    }

    double   mindist = zfixList.begin()->first;
    VectorXd zfix0   = zfixList.begin()->second;
    mtrx.zfix        = zfix0;
    MatrixXd Z       = mtrx.Ztrs.bottomRows(zsiz);
    mtrx.Ztrs        = Z;
    mtrx.integerCandidateCount = zfixList.size();
    mtrx.bestSquaredNorm = mindist;

    if (zfixList.size() > 1)
    {
        auto secondCandidate = std::next(zfixList.begin());
        mtrx.secondSquaredNorm = secondCandidate->first;
        if (mindist > 0)
        {
            mtrx.solutionRatio = mtrx.secondSquaredNorm / mindist;
        }
        else if (mtrx.secondSquaredNorm == 0)
        {
            mtrx.solutionRatio = 1;
        }
        else
        {
            mtrx.solutionRatio = std::numeric_limits<double>::max();
        }
    }

    switch (opt.mode)
    {
        case E_ARmode::LAMBDA:
            mtrx.diagnosticStatus = "RESOLVED_WITH_SUCCESS_RATE_ONLY";
            return zfix0.size();

        case E_ARmode::LAMBDA_ALT:
        {
            if (mtrx.secondSquaredNorm < 0)
            {
                mtrx.diagnosticStatus = "INSUFFICIENT_RATIO_CANDIDATES";
                return 0;
            }
            else if (mtrx.solutionRatio < opt.ratthr)
            {
                mtrx.diagnosticStatus = "RATIO_BELOW_THRESHOLD";
                return 0;
            }
            else
            {
                mtrx.diagnosticStatus = "RESOLVED_RATIO_ACCEPTED";
                return zfix0.size();
            }
        }

        case E_ARmode::LAMBDA_AL2:
        {
            for (auto& [dis, fixvec] : zfixList)
            {
                if ((dis / mindist) > opt.ratthr)
                    break;

                for (int l = 0; l < zfix0.size(); l++)
                {
                    if (zfix0(l) == -99999.5)
                        continue;

                    if (zfix0(l) != fixvec(l))
                        zfix0(l) = -99999.5;
                }
            }

            vector<int> zind;
            for (int k = 0; k < zfix0.size(); k++)
                if (zfix0(k) != -99999.5)
                    zind.push_back(k);
            tracepdeex(2, trace, "... %d ambiguties in common\n", zind.size());

            vector<int> xind;
            for (int k = 0; k < nmax; k++)
                xind.push_back(k);

            mtrx.zfix = zfix0(zind);
            mtrx.Ztrs = Z(zind, xind);

            mtrx.diagnosticStatus = zind.empty()
                ? "NO_COMMON_INTEGER_COMBINATIONS"
                : "RESOLVED_COMMON_SET";

            return zind.size();
        }

        case E_ARmode::LAMBDA_BIE:
        {
            double acum = 0;

            for (auto& [dis, fixvec] : zfixList)
            {
                double fct = exp(-0.5 * (dis - mindist));
                acum += fct;
            }

            VectorXd zbie = VectorXd::Zero(zsiz);

            for (auto& [dis, fixvec] : zfixList)
            {
                double fct = exp(-0.5 * (dis - mindist)) / acum;
                if (AR_VERBO)
                    trace << "\n"
                          << "BIE Candidate found:" << fixvec.transpose() << ";   dist= " << dis
                          << ";   fact= " << fct;
                zbie += fct * fixvec;
            }

            mtrx.zfix = zbie;
            mtrx.diagnosticStatus = "BIE_WEIGHTED_ESTIMATE";

            return zbie.size();
        }
    }

    return 0;
}

/** Ambiguity resolution function for Ginan */
int GNSS_AR(
    Trace&     trace,  ///< Debug trace
    GinAR_mtx& mtrx,   ///< Reference to structure containing float values and covariance
    GinAR_opt  opt     ///< Object containing processing options
)
{
    mtrx.diagnosticStatus = "NOT_RUN";
    mtrx.selectedDecorrelatedAmbiguityCount = 0;
    mtrx.integerCandidateCount = 0;
    mtrx.bootstrappedSuccessRate = -1;
    mtrx.bestSquaredNorm = -1;
    mtrx.secondSquaredNorm = -1;
    mtrx.solutionRatio = -1;
    mtrx.fullDecorrelatedTransform.resize(0, 0);

    switch (opt.mode)
    {
        case E_ARmode::OFF:
            mtrx.diagnosticStatus = "MODE_OFF";
            return 0;
        case E_ARmode::ROUND:
            return simple_round(trace, mtrx, opt);
        case E_ARmode::ITER_RND:
            return interat_round(trace, mtrx, opt);
        case E_ARmode::BOOTST:
            return integer_bootst(trace, mtrx, opt);
        case E_ARmode::LAMBDA:
            return lambda_search(trace, mtrx, opt);
        case E_ARmode::LAMBDA_ALT:
            return lambda_search(trace, mtrx, opt);
        case E_ARmode::LAMBDA_AL2:
            return lambda_search(trace, mtrx, opt);
        case E_ARmode::LAMBDA_BIE:
            return lambda_search(trace, mtrx, opt);
            // default:							tracepdeex(1, trace, "\n AR mode not supported \n");
    }

    return 0;
}
