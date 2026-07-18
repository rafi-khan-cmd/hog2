#include "RubiksCubeCorners.h"
#include "TemplateAStar.h"
#include <algorithm>
#include <chrono>
#include <cfloat>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <vector>




class FixedCorner : public RubiksCorner {
public:
    bool GoalTest(const RubiksCornerState &node,
                  const RubiksCornerState &goal) const override {
        for (int loc = 0; loc < 8; loc++) {
            if (node.GetCubeInLoc(loc) != goal.GetCubeInLoc(loc))
                return false;
        }
        for (int cube = 0; cube < 8; cube++) {
            if (node.GetCubeOrientation(cube) != goal.GetCubeOrientation(cube))
                return false;
        }
        return true;
    }
};

static RubiksCornerState GetCornerDual(const RubiksCornerState &s) {
    RubiksCornerState dual;
    for (int loc = 0; loc < 8; loc++) {
        int cube = s.GetCubeInLoc(loc);
        int orient = (int)s.GetCubeOrientation(cube);
        dual.SetCubeInLoc(cube, loc);
        dual.SetCubeOrientation(loc, (3 - orient) % 3);
    }
    return dual;
}

class AnyGoalCornerHeuristic : public Heuristic<RubiksCornerState> {
public:
    explicit AnyGoalCornerHeuristic(const RubikCornerPDB &pdb) : pdb(pdb) {}

    double HCost(const RubiksCornerState &node,
                 const RubiksCornerState &target) const override {
        static const RubiksCornerState goal;
        RubiksCornerState relative(node);

        int targetLocOfCube[8];
        int targetRotationOfCube[8];
        for (int loc = 0; loc < 8; loc++) {
            int cube = target.GetCubeInLoc(loc);
            targetLocOfCube[cube] = loc;
            targetRotationOfCube[cube] = target.GetCubeOrientation(cube);
        }

        for (int loc = 0; loc < 8; loc++) {
            int cube = node.GetCubeInLoc(loc);
            int relativeCube = targetLocOfCube[cube];
            int nodeOrient = node.GetCubeOrientation(cube);
            int targetOrient = targetRotationOfCube[cube];
            relative.SetCubeInLoc(loc, relativeCube);
            relative.SetCubeOrientation(relativeCube, (3 - targetOrient + nodeOrient) % 3);
        }

        return pdb.HCost(relative, goal);
    }

private:
    const RubikCornerPDB &pdb;
};


void RunPracticeExperiment()
{
    FixedCorner Corner;
    RubiksCornerState Solved;
    Solved.Reset();
    RubiksCornerState SolvedClean;
    SolvedClean.Reset();
    constexpr uint64_t kPDB8Size = 88179840;


    std::vector<int> Corners6 = {0, 1, 2, 3, 4, 5};
    RubikCornerPDB pdb6(&Corner, Solved, Corners6);
    printf("Building 6-corner PDB...\n");
    pdb6.BuildPDB(Solved, std::thread::hardware_concurrency());
    printf("Done.\n");

    AnyGoalCornerHeuristic Heuristic(pdb6);

    std::vector<uint64_t> BucketStarts = {0};
    constexpr uint64_t BucketSizeToRun = 100000;

    for (size_t k = 0; k < BucketStarts.size(); k++){
        uint64_t TestStartingIndex = BucketStarts[k];
        uint64_t TestEndingIndex = TestStartingIndex + BucketSizeToRun;
        std::vector<uint8_t> pdb8(kPDB8Size, 255);
        std::vector<uint8_t> UniqueNodeExpansionTracker(kPDB8Size, 0);

        TemplateAStar<
        RubiksCornerState,
        RubiksCornersAction,
        RubiksCorner> AStar;

        uint64_t Contradictions = 0;
        uint64_t TargetFill = 0;
        uint64_t NodesExpanded = 0;
        uint64_t NodesSeen = 0;
        uint64_t ClosedListFillsCount = 0;
        /*
        uint64_t OpenListItemsMyCount = 0;
        uint64_t NumOpenListItems = 0;
        uint64_t NumClosedListItems = 0;
        */
        uint64_t ClosedListItemMyCount = 0;
        uint64_t UnsafeClosedSkips = 0;
        uint64_t UniqueNodesExpandedAcrossSearches = 0;
        uint64_t TheNumberOfSearches = 0;
        uint64_t SafeClosedInBucketItems = 0;
        uint64_t SafeClosedOutOfBucketItems = 0;
        uint64_t UnsafeClosedInBucketItems = 0;
        uint64_t UnsafeOutOfBucketItems = 0;
        uint64_t DuplicateStatesEncountered = 0;
        uint64_t DualFillCount = 0;
        uint64_t DualHashOutsideBucket = 0;
        uint64_t DualHashInsideBucket = 0;
        AStar.SetHeuristic(&Heuristic);
        AStar.SetReopenNodes(true);
        uint64_t SolvedHash = RubikCornerPDB::GetStateHash(SolvedClean);
        pdb8[SolvedHash] = 0;
        uint8_t GoalFill = 0;
        if (SolvedHash >= TestStartingIndex && SolvedHash < TestEndingIndex) {
            GoalFill = 1;
        }
        auto StartTime = std::chrono::high_resolution_clock::now();
        uint64_t CurrentFilledEntries = TargetFill + ClosedListFillsCount + GoalFill + DualFillCount;


        for (uint64_t i = TestStartingIndex; i < TestEndingIndex && CurrentFilledEntries < (TestEndingIndex - TestStartingIndex); i++){
            if (pdb8[i] != 255){
                continue;
            }
            RubiksCornerState target;
            RubikCornerPDB::GetStateFromHash(target, i);
            std::vector<RubiksCornerState> Path;
            AStar.GetPath(&Corner, SolvedClean, target, Path);
            TheNumberOfSearches = TheNumberOfSearches + 1;
            NodesExpanded = NodesExpanded + AStar.GetNodesExpanded();
            NodesSeen = NodesSeen + AStar.GetNodesTouched();
            uint64_t TotalItems = AStar.GetNumItems();

            uint8_t Distance;

            if (Path.empty()) {
                if (i == RubikCornerPDB::GetStateHash(SolvedClean)) {
                    Distance = 0;
                } else {
                    continue;
                }
            } else {
                Distance = Path.size() - 1;
            }

            if (pdb8[i] == 255){
                pdb8[i] = Distance;
                TargetFill = TargetFill + 1;
            }
            else if (pdb8[i] != 255 && pdb8[i] != Distance){
                Contradictions = Contradictions + 1;
            }
            /*uint64_t NumOpenItems = AStar.GetNumOpenItems();*/
            for (uint64_t j = 0; j < TotalItems; j++){
                auto Item = AStar.GetItem(j);
                if (Item.where == kClosedList){
                    uint64_t TheHash = RubikCornerPDB::GetStateHash(Item.data);
                    ClosedListItemMyCount = ClosedListItemMyCount + 1;
                    if (UniqueNodeExpansionTracker[TheHash] == 0){
                        UniqueNodeExpansionTracker[TheHash] = 1;
                        UniqueNodesExpandedAcrossSearches = UniqueNodesExpandedAcrossSearches + 1;
                    }
                    if (TheHash >= TestStartingIndex && TheHash < TestEndingIndex){
                            if (pdb8[TheHash] == 255){
                                pdb8[TheHash] = Item.g;
                                ClosedListFillsCount = ClosedListFillsCount + 1;
                                SafeClosedInBucketItems = SafeClosedInBucketItems + 1;
                        }
                            else{
                                DuplicateStatesEncountered = DuplicateStatesEncountered + 1;
                                if (pdb8[TheHash] != Item.g){
                                    Contradictions = Contradictions + 1;
                                }
                            }
                            
                    }
                    else{
                        SafeClosedOutOfBucketItems = SafeClosedOutOfBucketItems + 1;

                        }
                    RubiksCornerState DualState = GetCornerDual(Item.data);
                    uint64_t DualHash = RubikCornerPDB::GetStateHash(DualState);
                    if (DualHash >= TestStartingIndex && DualHash < TestEndingIndex){
                        DualHashInsideBucket = DualHashInsideBucket + 1;
                        if (pdb8[DualHash] == 255){
                            pdb8[DualHash] = Item.g;
                            DualFillCount = DualFillCount + 1;
                        }
                        else if (pdb8[DualHash] != 255){
                            if (pdb8[DualHash] != Item.g){
                                Contradictions = Contradictions + 1;
                            }
                        }
                    }
                    else {
                        DualHashOutsideBucket = DualHashOutsideBucket + 1;
                    }
                
                }
                }
                CurrentFilledEntries = TargetFill + ClosedListFillsCount + GoalFill + DualFillCount;
                if (TheNumberOfSearches % 500 == 0) {
                    auto CheckpointTime = std::chrono::high_resolution_clock::now();
                    double CheckpointElapsed = std::chrono::duration<double>(CheckpointTime - StartTime).count();
                    double CheckpointOverlap = UniqueNodesExpandedAcrossSearches == 0 ? 0 : static_cast<double>(NodesExpanded) / UniqueNodesExpandedAcrossSearches;
                    double CheckpointInBucketRate = ClosedListItemMyCount == 0 ? 0 : static_cast<double>(SafeClosedInBucketItems) / ClosedListItemMyCount;
                    double CheckpointOutBucketRate = ClosedListItemMyCount == 0 ? 0 : static_cast<double>(SafeClosedOutOfBucketItems) / ClosedListItemMyCount;
                    double CheckpointClosedFillRate = ClosedListItemMyCount == 0 ? 0 : static_cast<double>(ClosedListFillsCount) / ClosedListItemMyCount;
                    double CheckpointFillsPerSearch = TheNumberOfSearches == 0 ? 0 : static_cast<double>(CurrentFilledEntries) / TheNumberOfSearches;
                    printf("progress: searches=%" PRIu64 " filled=%" PRIu64 " bucket=%" PRIu64 " i=%" PRIu64 " target=%" PRIu64 " closed=%" PRIu64 " dual=%" PRIu64 " contradictions=%" PRIu64 " nodesExpanded=%" PRIu64 " nodesTouched=%" PRIu64 " closedItems=%" PRIu64 " uniqueExpanded=%" PRIu64 " overlap=%.6f duplicates=%" PRIu64 " safeIn=%" PRIu64 " safeOut=%" PRIu64 " dualIn=%" PRIu64 " dualOut=%" PRIu64 " unsafe=%" PRIu64 " inRate=%.6f outRate=%.6f closedFillRate=%.6f fillsPerSearch=%.6f elapsed=%.6f\n",
                           TheNumberOfSearches, CurrentFilledEntries, TestEndingIndex - TestStartingIndex, i,
                           TargetFill, ClosedListFillsCount, DualFillCount, Contradictions,
                           NodesExpanded, NodesSeen, ClosedListItemMyCount, UniqueNodesExpandedAcrossSearches,
                           CheckpointOverlap, DuplicateStatesEncountered, SafeClosedInBucketItems, SafeClosedOutOfBucketItems,
                           DualHashInsideBucket, DualHashOutsideBucket, UnsafeClosedSkips,
                           CheckpointInBucketRate, CheckpointOutBucketRate, CheckpointClosedFillRate,
                           CheckpointFillsPerSearch, CheckpointElapsed);
                }
            }
        
        auto EndTime = std::chrono::high_resolution_clock::now();
        double ElapsedSeconds = std::chrono::duration<double>(EndTime - StartTime).count();
        uint64_t BucketSize = TestEndingIndex - TestStartingIndex;
        uint64_t TotalFilledEntries = TargetFill + ClosedListFillsCount + GoalFill + DualFillCount;
        double OverlapCount = UniqueNodesExpandedAcrossSearches == 0 ? 0 : static_cast<double>(NodesExpanded) / UniqueNodesExpandedAcrossSearches;
        double InBucketClosedRate = ClosedListItemMyCount == 0 ? 0 : static_cast<double>(SafeClosedInBucketItems) / ClosedListItemMyCount;
        double OutOfBucketClosedRate = ClosedListItemMyCount == 0 ? 0 : static_cast<double>(SafeClosedOutOfBucketItems) / ClosedListItemMyCount;
        double ClosedFillRate = ClosedListItemMyCount == 0 ? 0 : static_cast<double>(ClosedListFillsCount) / ClosedListItemMyCount;
        double FillsPerSearch = TheNumberOfSearches == 0 ? 0 : static_cast<double>(TotalFilledEntries) / TheNumberOfSearches;

        printf("progress: searches=%" PRIu64 " filled=%" PRIu64 " bucket=%" PRIu64 " i=%" PRIu64 " target=%" PRIu64 " closed=%" PRIu64 " dual=%" PRIu64 " contradictions=%" PRIu64 " nodesExpanded=%" PRIu64 " nodesTouched=%" PRIu64 " closedItems=%" PRIu64 " uniqueExpanded=%" PRIu64 " overlap=%.6f duplicates=%" PRIu64 " safeIn=%" PRIu64 " safeOut=%" PRIu64 " dualIn=%" PRIu64 " dualOut=%" PRIu64 " unsafe=%" PRIu64 " inRate=%.6f outRate=%.6f closedFillRate=%.6f fillsPerSearch=%.6f elapsed=%.6f\n",
               TheNumberOfSearches, TotalFilledEntries, BucketSize, TestEndingIndex,
               TargetFill, ClosedListFillsCount, DualFillCount, Contradictions,
               NodesExpanded, NodesSeen, ClosedListItemMyCount, UniqueNodesExpandedAcrossSearches,
               OverlapCount, DuplicateStatesEncountered, SafeClosedInBucketItems, SafeClosedOutOfBucketItems,
               DualHashInsideBucket, DualHashOutsideBucket, UnsafeClosedSkips,
               InBucketClosedRate, OutOfBucketClosedRate, ClosedFillRate,
               FillsPerSearch, ElapsedSeconds);

        printf("Starting  %" PRIu64"\n", k);
        printf("Bucket Range: [%" PRIu64 ", %" PRIu64 ")\n", TestStartingIndex, TestEndingIndex);
        printf("Bucket Size k: %" PRIu64 "\n", BucketSize);
        printf("Ranking Type: POLYNOMIAL_RANK / lexicographic-style full 8-corner ranking\n");
        /*printf("Ranking Type: LINEAR_RANK / Myrvold-Ruskey-style full 8-corner ranking\n");*/
        printf("The Number of Searches Run: %" PRIu64 "\n", TheNumberOfSearches);
        printf("Contradictions: %" PRIu64 "\n", Contradictions);
        printf("Target Fills: %" PRIu64 "\n", TargetFill);
        printf("Total Filled Entries: %" PRIu64 "\n", TotalFilledEntries);
        printf("Nodes Expanded: %" PRIu64 "\n", NodesExpanded);
        printf("Nodes Touched: %" PRIu64 "\n", NodesSeen);
        printf("Closed List Fills: %" PRIu64 "\n", ClosedListFillsCount);
        printf("Dual List Fills: %" PRIu64 "\n", DualFillCount);
        printf("Dual Hash Inside Bucket: %" PRIu64 "\n", DualHashInsideBucket);
        printf("Dual Hash Outside Bucket: %" PRIu64 "\n", DualHashOutsideBucket);    
        printf("Closed List Items My Count: %" PRIu64 "\n", ClosedListItemMyCount);
        printf("Unsafe Closed Skips: %" PRIu64 "\n", UnsafeClosedSkips);
        printf("Unique Nodes Expanded Across Searches: %" PRIu64 "\n", UniqueNodesExpandedAcrossSearches);
        printf("Overlap Count: %.6f\n", OverlapCount);
        printf("Duplicate States Encountered: %" PRIu64 "\n", DuplicateStatesEncountered);
        printf("Safe Closed In-Bucket Items: %" PRIu64 "\n", SafeClosedInBucketItems);
        printf("Safe Closed Out-Of-Bucket Items: %" PRIu64 "\n", SafeClosedOutOfBucketItems);
        printf("Unsafe Closed In-Bucket Items: %" PRIu64 "\n", UnsafeClosedInBucketItems);
        printf("Unsafe Out Of Bucket Items: %" PRIu64 "\n", UnsafeOutOfBucketItems);
        printf("In-Bucket Closed Rate: %.6f\n", InBucketClosedRate);
        printf("Out-Of-Bucket Closed Rate: %.6f\n", OutOfBucketClosedRate);
        printf("Closed Fill Rate: %.6f\n", ClosedFillRate);
        printf("Fills Per Search: %.6f\n", FillsPerSearch);
        printf("Elapsed Time: %.6f seconds\n", ElapsedSeconds);

    }

}
