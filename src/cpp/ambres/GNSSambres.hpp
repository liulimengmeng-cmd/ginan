#pragma once

#include "common/algebra.hpp"
#include "common/common.hpp"
#include "common/constants.hpp"
#include "common/eigenIncluder.hpp"
#include "common/observations.hpp"
#include "common/receiver.hpp"
#include "common/satSys.hpp"
#include "common/trace.hpp"

extern double FIXED_AMB_VAR;
extern bool   AR_VERBO;

struct GinAR_mtx
{
    map<int, KFKey> ambmap;
    VectorXd        aflt;
    MatrixXd        Paflt;

    MatrixXd Ztrs;
    MatrixXd Ltrs;
    VectorXd Dtrs;

    // Full unimodular transform produced by LAMBDA before partial ambiguity
    // selection.  Ztrs may later contain only the accepted bottom rows.
    MatrixXd fullDecorrelatedTransform;

    VectorXd zflt;
    VectorXd zfix;

    VectorXd afix;
    MatrixXd Pafix;

    string diagnosticStatus = "NOT_RUN";
    int    selectedDecorrelatedAmbiguityCount = 0;
    int    integerCandidateCount = 0;
    double bootstrappedSuccessRate = -1;
    double bestSquaredNorm = -1;
    double secondSquaredNorm = -1;
    double solutionRatio = -1;
};

struct ConditionalIntegerComplement
{
    GinAR_mtx ambiguityResolution;
    MatrixXd transformToInputCoordinates;
    string diagnosticStatus = "NOT_RUN";
};

struct GinAR_opt
{
    string           recv;
    map<E_Sys, bool> sys_solve;

    bool     endu = false;
    E_ARmode mode = E_ARmode::OFF;  /* AR mode */

    int nset = 0;                   /* candidate set size for lambda */
    int nitr = 3;                   /* number of iterations for iter_rnd */
    int minimumDecorrelatedAmbiguityCount = 3;

    double MIN_Elev_prc = D2R * 10; /* min elevation for processing */
    double MIN_Elev_AR  = D2R * 15; /* min elevation for AR */
    double MIN_Elev_piv = D2R * 20; /* min elevation for pivot */

    double sucthr = 0.9999;         /* success rate threshold */
    double ratthr = 3;              /* ratio test threshold */

    bool   clear_old_amb = false;
    int    Max_Hold_epc  = 0;   /* max hold (epoch) */
    double Max_Hold_tim  = 600; /* max hold (seconds) */
};

struct ReceiverAmbiguityDatum
{
    string receiver;
    E_Sys  system = E_Sys::NONE;
    int    observation = 0;
    SatSys pivot;
    int    memberCount = 0;
    bool   singleDifferenced = false;
};

struct ReceiverAmbiguityTransform
{
    MatrixXd matrix;
    vector<ReceiverAmbiguityDatum> groups;
    int singleDifferencedGroupCount = 0;
    int identityGroupCount = 0;
    int droppedSingletonGroupCount = 0;
};

ReceiverAmbiguityTransform buildReceiverAmbiguityIntegerTransform(
    const GinAR_mtx&        ambiguityResolution,
    const map<E_Sys, bool>& receiverAmbiguityPivot
);

bool mapIntegerAmbiguityConstraintsToOriginalState(
    GinAR_mtx&                        integerAmbiguityResolution,
    const ReceiverAmbiguityTransform& integerTransform,
    const map<int, KFKey>&            originalAmbiguityMap
);

ConditionalIntegerComplement buildConditionalIntegerComplement(
    const GinAR_mtx& inputAmbiguities,
    const GinAR_mtx& acceptedPartialResolution
);

int GNSS_AR(Trace& trace, GinAR_mtx& mtrx, GinAR_opt opt);
