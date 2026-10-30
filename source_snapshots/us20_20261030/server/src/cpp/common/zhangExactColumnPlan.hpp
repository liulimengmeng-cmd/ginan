#pragma once

#include "common/zhangIntegerAudit.hpp"

#include <stdexcept>
#include <atomic>

/** A bounded, resumable version of the column operations used by
 * zhangExactProjectedKernelColumns. This object is computational state, NOT
 * an integer certificate: the owner must already have verified its affine
 * domain. Row operations are private; only the unimodular column operations
 * are retained, never the full ambient square transform or physical kernel.
 *
 * Cancellation is allowed while scanning or constructing a scratch row /
 * column. A mutation and its journal entry commit together with no budget
 * checkpoint in between. Thus a cancelled object is always restartable. */
struct ZhangExactColumnPlanLimitExceeded {};

struct ZhangExactColumnOperation
{
    enum class Kind { SWAP, ADD } kind = Kind::SWAP;
    std::size_t destination = 0, source = 0;
    ZhangExactInteger factor = 0;
};

struct ZhangExactColumnPlanLimits
{
    std::size_t maximumMatrixCells = 16000000;
    std::size_t maximumMatrixBytes = 512 * 1024 * 1024;
    std::size_t maximumOperations = 1000000;
    std::size_t maximumJournalBytes = 128 * 1024 * 1024;
    unsigned maximumCoefficientBits = 16384;
};

class ZhangExactColumnPlan
{
    struct PoolReservation
    {
        inline static std::atomic<std::size_t> total{0};
        std::size_t bytes=0;
        void resize(std::size_t desired)
        {
            if(desired<=bytes) {total.fetch_sub(bytes-desired);bytes=desired;return;}
            const auto increase=desired-bytes;
            auto observed=total.load();
            do
            {
                if(observed>768*1024*1024 || increase>768*1024*1024-observed)
                    throw ZhangExactColumnPlanLimitExceeded{};
            } while(!total.compare_exchange_weak(observed,observed+increase));
            bytes=desired;
        }
        ~PoolReservation() {total.fetch_sub(bytes);}
        PoolReservation()=default;
        PoolReservation(const PoolReservation&)=delete;
        PoolReservation& operator=(const PoolReservation&)=delete;
    } reservation;
    enum class Phase { PIVOT, ELIMINATE_ROWS, ELIMINATE_COLUMNS, COMPLETE };
    ZhangExactMatrix matrix;
    Phase phase = Phase::PIVOT;
    std::size_t scanRow = 0, scanColumn = 0;
    std::size_t selectedRow = 0, selectedColumn = 0;
    bool pivotReady = false;
    ZhangExactColumnPlanLimits limits;
    std::size_t matrixBytes=0;

    static ZhangExactInteger magnitude(const ZhangExactInteger& x)
    { return x < 0 ? -x : x; }

    void checkCoefficient(const ZhangExactInteger& x) const
    {
        if (x != 0 && boost::multiprecision::msb(magnitude(x)) >= limits.maximumCoefficientBits)
            throw ZhangExactColumnPlanLimitExceeded{};
    }

    static std::size_t integerBytes(const ZhangExactInteger& x)
    {return sizeof(x)+(x==0?0:(boost::multiprecision::msb(magnitude(x))+8)/8);}

    void checkMatrixBytes(std::size_t candidate) const
    {if(candidate>limits.maximumMatrixBytes) throw ZhangExactColumnPlanLimitExceeded{};}

    void appendOperation(ZhangExactColumnOperation operation)
    {
        checkCoefficient(operation.factor);
        const std::size_t bits = operation.factor == 0 ? 0
            : boost::multiprecision::msb(magnitude(operation.factor)) + 1;
        const std::size_t bytes = sizeof(operation) + (bits + 7) / 8;
        if (operations.size() >= limits.maximumOperations ||
            journalBytes > limits.maximumJournalBytes ||
            bytes > limits.maximumJournalBytes - journalBytes)
            throw ZhangExactColumnPlanLimitExceeded{};
        operations.push_back(std::move(operation));
        journalBytes += bytes;
    }

    void addRow(std::size_t dst,std::size_t src,const ZhangExactInteger& factor)
    {
        ZhangExactVector scratch(columns);
        std::size_t candidateBytes=matrixBytes;
        for (std::size_t c=0;c<columns;++c)
        {
            zhangExactWorkCheckpoint("COLUMN_PLAN_ROW_SCRATCH");
            scratch[c]=matrix[dst][c]+factor*matrix[src][c];
            checkCoefficient(scratch[c]);
            candidateBytes-=integerBytes(matrix[dst][c]);
            candidateBytes+=integerBytes(scratch[c]);
        }
        checkMatrixBytes(candidateBytes);
        zhangExactWorkCheckpoint("COLUMN_PLAN_ROW_COMMIT");
        reservation.resize(candidateBytes);
        matrix[dst].swap(scratch);
        matrixBytes=candidateBytes;
    }

    void addColumn(std::size_t dst,std::size_t src,const ZhangExactInteger& factor)
    {
        ZhangExactVector scratch(rows);
        std::size_t candidateBytes=matrixBytes;
        for (std::size_t r=0;r<rows;++r)
        {
            zhangExactWorkCheckpoint("COLUMN_PLAN_COLUMN_SCRATCH");
            scratch[r]=matrix[r][dst]+factor*matrix[r][src];
            checkCoefficient(scratch[r]);
            candidateBytes-=integerBytes(matrix[r][dst]);
            candidateBytes+=integerBytes(scratch[r]);
        }
        checkMatrixBytes(candidateBytes);
        zhangExactWorkCheckpoint("COLUMN_PLAN_COLUMN_COMMIT");
        reservation.resize(candidateBytes);
        appendOperation({ZhangExactColumnOperation::Kind::ADD,dst,src,factor});
        for (std::size_t r=0;r<rows;++r) matrix[r][dst].swap(scratch[r]);
        matrixBytes=candidateBytes;
    }

public:
    const std::size_t rows,columns;
    std::size_t rank=0,journalBytes=0,advanceCalls=0,sourceNonzeros=0,sourceMaxBits=0;
    std::vector<ZhangExactColumnOperation> operations;

    explicit ZhangExactColumnPlan(const ZhangExactMatrix& source,
        std::size_t emptyColumnCount=0,ZhangExactColumnPlanLimits bounds={})
        : limits(bounds),rows(source.size()),
          columns(source.empty()?emptyColumnCount:source.front().size())
    {
        if (columns && rows > limits.maximumMatrixCells/columns)
            throw ZhangExactColumnPlanLimitExceeded{};
        for (const auto& row:source)
        {
            if (row.size()!=columns) throw std::invalid_argument("COLUMN_PLAN_SHAPE");
            for (const auto& x:row)
            {
                checkCoefficient(x);matrixBytes+=integerBytes(x);checkMatrixBytes(matrixBytes);
                if(x!=0)
                {
                    ++sourceNonzeros;
                    sourceMaxBits=std::max(sourceMaxBits,
                        std::size_t(boost::multiprecision::msb(magnitude(x))+1));
                }
            }
        }
        // The single bounded copy is not interruptible, exactly like an
        // Eigen/library allocation. No partially copied object is published.
        reservation.resize(matrixBytes);
        matrix=source;
        selectedRow=rows;selectedColumn=columns;
        if (!rows || !columns) phase=Phase::COMPLETE;
    }

    bool complete() const {return phase==Phase::COMPLETE;}
    std::size_t residentMatrixBytes() const {return matrixBytes;}
    static std::size_t totalResidentMatrixBytes() {return PoolReservation::total.load();}

    void advance()
    {
        ++advanceCalls;
        while (!complete())
        {
            if (phase==Phase::PIVOT)
            {
                if (rank>=rows || rank>=columns)
                {phase=Phase::COMPLETE;break;}
                // Indices persist across timeout, including large zero scans.
                while (!pivotReady && scanRow<rows)
                {
                    while (scanColumn<columns)
                    {
                        zhangExactWorkCheckpoint("COLUMN_PLAN_PIVOT_SCAN");
                        const auto& x=matrix[scanRow][scanColumn];
                        if (x!=0 && (selectedRow==rows ||
                            magnitude(x)<magnitude(matrix[selectedRow][selectedColumn])))
                        {
                            selectedRow=scanRow;selectedColumn=scanColumn;
                            if (magnitude(x)==1) pivotReady=true;
                        }
                        ++scanColumn;
                        if (pivotReady) break;
                    }
                    if (!pivotReady) {++scanRow;scanColumn=rank;}
                }
                if (selectedRow==rows) {phase=Phase::COMPLETE;break;}
                // No cancellation between the initial row swap and column
                // swap: otherwise retrying the row swap would undo it.
                zhangExactWorkCheckpoint("COLUMN_PLAN_PIVOT_COMMIT");
                appendOperation({ZhangExactColumnOperation::Kind::SWAP,rank,selectedColumn,0});
                std::swap(matrix[rank],matrix[selectedRow]);
                for (auto& row:matrix) std::swap(row[rank],row[selectedColumn]);
                phase=Phase::ELIMINATE_ROWS;scanRow=rank+1;
            }
            if (phase==Phase::ELIMINATE_ROWS)
            {
                while (scanRow<rows)
                {
                    zhangExactWorkCheckpoint("COLUMN_PLAN_ROW_SCAN");
                    if (matrix[scanRow][rank]==0) {++scanRow;continue;}
                    const ZhangExactInteger factor=-matrix[scanRow][rank]/matrix[rank][rank];
                    addRow(scanRow,rank,factor);
                    if (matrix[scanRow][rank]!=0) std::swap(matrix[scanRow],matrix[rank]);
                    scanRow=rank+1;
                }
                phase=Phase::ELIMINATE_COLUMNS;scanColumn=rank+1;
            }
            if (phase==Phase::ELIMINATE_COLUMNS)
            {
                bool restarted=false;
                while (scanColumn<columns)
                {
                    zhangExactWorkCheckpoint("COLUMN_PLAN_COLUMN_SCAN");
                    if (matrix[rank][scanColumn]==0) {++scanColumn;continue;}
                    const ZhangExactInteger factor=-matrix[rank][scanColumn]/matrix[rank][rank];
                    addColumn(scanColumn,rank,factor);
                    // Complete the Euclidean step without a checkpoint after
                    // its committed add, so cancellation cannot repeat it.
                    if (matrix[rank][scanColumn]!=0)
                    {
                        // A resource-limit exception discards the whole plan;
                        // only budget cancellation is a resumable condition.
                        appendOperation({ZhangExactColumnOperation::Kind::SWAP,scanColumn,rank,0});
                        for (auto& row:matrix) std::swap(row[scanColumn],row[rank]);
                    }
                    phase=Phase::ELIMINATE_ROWS;scanRow=rank+1;
                    restarted=true;break;
                }
                if (restarted) continue;
                ++rank;phase=Phase::PIVOT;scanRow=rank;scanColumn=rank;
                selectedRow=rows;selectedColumn=columns;pivotReady=false;
            }
        }
        if (complete()) {ZhangExactMatrix{}.swap(matrix);matrixBytes=0;reservation.resize(0);}
    }
};

/** Target-local replay state. All extra columns are genuinely free integer
 * variables, appended after the source-domain columns in the original order. */
class ZhangExactColumnProjection
{
    ZhangExactMatrix transformed;
    ZhangExactMatrix partialImage;
    std::size_t nextOperation=0;
    std::size_t outputRow=0,outputColumn=0;
    std::size_t transformedBytes=0;
    std::size_t outputBytes=0;
    bool finished=false,outputInitialized=false;
    static std::size_t integerBytes(const ZhangExactInteger& x)
    {
        if(x==0) return sizeof(x);
        const ZhangExactInteger magnitude=x<0?-x:x;
        const auto bits=boost::multiprecision::msb(magnitude)+1;
        if(bits>16384) throw ZhangExactColumnPlanLimitExceeded{};
        return sizeof(x)+(bits+7)/8;
    }
    static void checkBytes(std::size_t bytes)
    {if(bytes>64*1024*1024) throw ZhangExactColumnPlanLimitExceeded{};}
public:
    ZhangExactMatrix image;
    explicit ZhangExactColumnProjection(ZhangExactMatrix targets)
        :transformed(std::move(targets))
    {
        for(const auto& row:transformed) for(const auto& x:row)
        {transformedBytes+=integerBytes(x);checkBytes(transformedBytes);}
    }
    std::size_t completedOperations() const {return nextOperation;}
    std::size_t residentBytes() const {return transformedBytes+outputBytes;}
    bool complete() const {return finished;}
    void advance(const ZhangExactColumnPlan& plan)
    {
        if (finished) return;
        if (!plan.complete()) throw std::logic_error("COLUMN_PLAN_INCOMPLETE");
        const std::size_t width=transformed.empty()?plan.columns:transformed.front().size();
        for (const auto& row:transformed)
            if (row.size()!=width || width<plan.columns)
                throw std::invalid_argument("COLUMN_PROJECTION_SHAPE");
        while (nextOperation<plan.operations.size())
        {
            const auto& operation=plan.operations[nextOperation];
            zhangExactWorkCheckpoint("COLUMN_PLAN_TARGET_REPLAY");
            if (operation.kind==ZhangExactColumnOperation::Kind::SWAP)
            {
                for (auto& row:transformed)
                    std::swap(row[operation.destination],row[operation.source]);
            }
            else
            {
                ZhangExactVector scratch(transformed.size());
                std::size_t candidateBytes=transformedBytes;
                for (std::size_t r=0;r<transformed.size();++r)
                {
                    zhangExactWorkCheckpoint("COLUMN_PLAN_TARGET_SCRATCH");
                    scratch[r]=transformed[r][operation.destination]+
                        operation.factor*transformed[r][operation.source];
                    candidateBytes-=integerBytes(transformed[r][operation.destination]);
                    candidateBytes+=integerBytes(scratch[r]);
                }
                checkBytes(candidateBytes);
                zhangExactWorkCheckpoint("COLUMN_PLAN_TARGET_COMMIT");
                for (std::size_t r=0;r<transformed.size();++r)
                    transformed[r][operation.destination].swap(scratch[r]);
                transformedBytes=candidateBytes;
            }
            ++nextOperation;
        }
        if (!outputInitialized)
        {
            partialImage.assign(transformed.size(),ZhangExactVector(width-plan.rank));
            outputBytes=transformed.size()*(width-plan.rank)*sizeof(ZhangExactInteger);
            outputColumn=plan.rank;outputInitialized=true;
        }
        while (outputRow<transformed.size())
        {
            while (outputColumn<width)
            {
                zhangExactWorkCheckpoint("COLUMN_PLAN_IMAGE_COPY");
                partialImage[outputRow][outputColumn-plan.rank]=transformed[outputRow][outputColumn];
                outputBytes+=integerBytes(transformed[outputRow][outputColumn])-sizeof(ZhangExactInteger);
                ++outputColumn;
            }
            ++outputRow;outputColumn=plan.rank;
        }
        zhangExactWorkCheckpoint("COLUMN_PLAN_IMAGE_COMMIT");
        image.swap(partialImage);ZhangExactMatrix{}.swap(transformed);transformedBytes=0;finished=true;
    }
};
