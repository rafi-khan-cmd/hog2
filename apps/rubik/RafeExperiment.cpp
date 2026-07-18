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

namespace {

constexpr uint64_t kPDB8Size = 88179840;
#if 0 // Reserved for the disabled staged PDB tests.
constexpr size_t kGeneratedTargetLimit = 250000;
constexpr int kGeneratedTargetDepth = 6;
constexpr uint64_t kArbitraryScanLimit = 2000000;
constexpr uint64_t kArbitrarySearchLimit = 1000;
constexpr uint64_t kArbitraryExpansionCap = 2000000;
constexpr uint64_t kContradictionPrintLimit = 10;
#endif
constexpr size_t kDiagnosticStateLimit = 2000;
constexpr size_t kDiagnosticTargetLimit = 200;
constexpr size_t kDiagnosticAdmissibilityPairs = 200;
constexpr uint64_t kDiagnosticAStarCap = 100000;
constexpr uint64_t kDiagnosticPrintLimit = 10;
constexpr bool kRunTable8OverlapDiagnostic = true;
constexpr uint64_t kTable8BucketStart = 40000000;
constexpr uint64_t kTable8BucketEnd = 40100000;
constexpr size_t kTable8Instances = 10;
constexpr size_t kTable8MaxGoals = 128;

#if 0 // Used only by the disabled staged PDB tests.
struct TestTarget {
    RubiksCornerState state;
    uint64_t hash;
    uint64_t expectedDepth;
};
#endif

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

class AnyGoalCornerHeuristic : public Heuristic<RubiksCornerState> {
public:
    explicit AnyGoalCornerHeuristic(const RubikCornerPDB &pdb) : pdb(pdb) {}

    double HCost(const RubiksCornerState &node,
                 const RubiksCornerState &target) const override {
        RubiksCornerState relative;
        relative.Reset();

        int nodeLocOfCube[8];
        for (int loc = 0; loc < 8; loc++) {
            nodeLocOfCube[node.GetCubeInLoc(loc)] = loc;
        }

        for (int loc = 0; loc < 8; loc++) {
            int targetCube = target.GetCubeInLoc(loc);
            relative.SetCubeInLoc(loc, nodeLocOfCube[targetCube]);
        }

        for (int loc = 0; loc < 8; loc++) {
            int targetCube = target.GetCubeInLoc(loc);
            int relativeCube = nodeLocOfCube[targetCube];
            int nodeOrient = (int)node.GetCubeOrientation(targetCube);
            int targetOrient = (int)target.GetCubeOrientation(targetCube);
            relative.SetCubeOrientation(relativeCube, (targetOrient - nodeOrient + 3) % 3);
        }

        return pdb.HCost(relative, relative);
    }

private:
    const RubikCornerPDB &pdb;
};

#if 0 // Full-cube parity helper; not a reachability filter in HOG2's corner-only model.
bool HasOddPermutationParity(const RubiksCornerState &state)
{
    int inversions = 0;
    for (int i = 0; i < 8; i++) {
        int a = state.GetCubeInLoc(i);
        for (int j = i + 1; j < 8; j++) {
            if (a > state.GetCubeInLoc(j))
                inversions++;
        }
    }
    return (inversions % 2) != 0;
}
#endif

bool RunAStar(TemplateAStar<RubiksCornerState, RubiksCornersAction, RubiksCorner> &astar,
              FixedCorner &corner,
              const RubiksCornerState &start,
              const RubiksCornerState &target,
              std::vector<RubiksCornerState> &path,
              uint64_t expansionCap)
{
    if (expansionCap == 0) {
        astar.GetPath(&corner, start, target, path);
        return !path.empty();
    }

    if (!astar.InitializeSearch(&corner, start, target, path))
        return true;

    while (!astar.DoSingleSearchStep(path)) {
        if (astar.GetNodesExpanded() >= expansionCap)
            return false;
    }
    return true;
}

#if 0 // Used only by the disabled normal PDB bucket fill and staged tests.
uint64_t ComputeDistanceToIndex(FixedCorner &corner,
                                AnyGoalCornerHeuristic &heuristic,
                                const RubiksCornerState &solved,
                                uint64_t index)
{
    RubiksCornerState target;
    RubikCornerPDB::GetStateFromHash(target, index);

    TemplateAStar<RubiksCornerState, RubiksCornersAction, RubiksCorner> astar;
    astar.SetHeuristic(&heuristic);
    astar.SetReopenNodes(true);

    std::vector<RubiksCornerState> path;
    astar.GetPath(&corner, solved, target, path);
    if (path.empty())
        return UINT64_MAX;
    return path.size() - 1;
}

uint64_t FillFromClosedList(TemplateAStar<RubiksCornerState, RubiksCornersAction, RubiksCorner> &astar,
                            std::vector<uint8_t> &pdb8,
                            uint64_t &closedCount,
                            uint64_t &contradictions,
                            uint64_t &unsafeSkipped)
{
    uint64_t newFills = 0;
    closedCount = 0;
    unsafeSkipped = 0;

    double minOpenG = DBL_MAX;
    for (unsigned int i = 0; i < astar.GetNumOpenItems(); i++)
        minOpenG = std::min(minOpenG, astar.GetOpenItem(i).g);

    for (int i = 0; i < astar.GetNumItems(); i++) {
        const auto &item = astar.GetItem(i);
        if (item.where != kClosedList)
            continue;

        closedCount++;
        if (item.g > minOpenG) {
            unsafeSkipped++;
            continue;
        }

        uint64_t h = RubikCornerPDB::GetStateHash(item.data);
        uint8_t g = (uint8_t)item.g;
        if (pdb8[h] == 0xFF) {
            pdb8[h] = g;
            newFills++;
        } else if (pdb8[h] != g) {
            contradictions++;
            if (contradictions <= kContradictionPrintLimit) {
                printf("CONTRADICTION: hash=%" PRIu64 " old=%u new=%u\n",
                       h, (unsigned)pdb8[h], (unsigned)g);
            }
        }
    }
    return newFills;
}

bool FillTargetFromPath(std::vector<uint8_t> &pdb8,
                        uint64_t targetHash,
                        uint64_t pathLength,
                        uint64_t &newFills,
                        uint64_t &contradictions)
{
    uint8_t g = (uint8_t)pathLength;
    if (pdb8[targetHash] == 0xFF) {
        pdb8[targetHash] = g;
        newFills++;
        return true;
    }
    if (pdb8[targetHash] != g) {
        contradictions++;
        if (contradictions <= kContradictionPrintLimit) {
            printf("CONTRADICTION: target hash=%" PRIu64 " old=%u new=%u\n",
                   targetHash, (unsigned)pdb8[targetHash], (unsigned)g);
        }
        return false;
    }
    return true;
}
#endif

void GenerateReachableSamples(FixedCorner &corner,
                              const RubiksCornerState &solved,
                              std::vector<RubiksCornerState> &states,
                              std::vector<uint64_t> &hashes,
                              size_t limit)
{
    states.clear();
    hashes.clear();
    states.push_back(solved);
    hashes.push_back(RubikCornerPDB::GetStateHash(solved));

    std::vector<RubiksCornerState> frontier;
    frontier.push_back(solved);
    while (!frontier.empty() && states.size() < limit) {
        std::vector<RubiksCornerState> nextFrontier;
        for (const auto &parent : frontier) {
            std::vector<RubiksCornerState> children;
            corner.GetSuccessors(parent, children);
            for (const auto &child : children) {
                uint64_t hash = RubikCornerPDB::GetStateHash(child);
                if (std::find(hashes.begin(), hashes.end(), hash) != hashes.end())
                    continue;
                states.push_back(child);
                hashes.push_back(hash);
                nextFrontier.push_back(child);
                if (states.size() >= limit)
                    break;
            }
            if (states.size() >= limit)
                break;
        }
        frontier.swap(nextFrontier);
    }
}

void RunHeuristicDiagnostics(FixedCorner &corner,
                             const RubiksCornerState &solved,
                             AnyGoalCornerHeuristic &heuristic)
{
    std::vector<RubiksCornerState> samples;
    std::vector<uint64_t> sampleHashes;
    GenerateReachableSamples(corner, solved, samples, sampleHashes, kDiagnosticStateLimit);

    size_t targetCount = std::min(samples.size(), kDiagnosticTargetLimit);

    uint64_t consistencyChecks = 0;
    uint64_t consistencyFailures = 0;
    for (size_t t = 0; t < targetCount; t++) {
        const RubiksCornerState &target = samples[t];
        for (size_t n = 0; n < samples.size(); n++) {
            double hNode = heuristic.HCost(samples[n], target);
            std::vector<RubiksCornerState> children;
            corner.GetSuccessors(samples[n], children);
            for (const auto &child : children) {
                double edgeCost = corner.GCost(samples[n], child);
                double hChild = heuristic.HCost(child, target);
                consistencyChecks++;
                if (hNode > edgeCost + hChild + 0.000001) {
                    consistencyFailures++;
                    if (consistencyFailures <= kDiagnosticPrintLimit) {
                        printf("CONSISTENCY FAILURE: target=%" PRIu64
                               " node=%" PRIu64 " child=%" PRIu64
                               " hNode=%.0f edge=%.0f hChild=%.0f\n",
                               sampleHashes[t], sampleHashes[n],
                               RubikCornerPDB::GetStateHash(child),
                               hNode, edgeCost, hChild);
                    }
                }
            }
        }
    }

    TemplateAStar<RubiksCornerState, RubiksCornersAction, RubiksCorner> astar;
    astar.SetHeuristic(&heuristic);
    astar.SetReopenNodes(true);

    uint64_t admissibilityChecks = 0;
    uint64_t admissibilityFailures = 0;
    uint64_t cappedAdmissibilitySearches = 0;
    uint64_t searchesWithReopens = 0;
    uint64_t reopenedItems = 0;
    uint64_t totalExpanded = 0;

    for (size_t i = 0; i < samples.size() && admissibilityChecks < kDiagnosticAdmissibilityPairs; i++) {
        for (size_t j = 0; j < targetCount && admissibilityChecks < kDiagnosticAdmissibilityPairs; j++) {
            std::vector<RubiksCornerState> path;
            bool completed = RunAStar(astar, corner, samples[i], samples[j], path, kDiagnosticAStarCap);
            totalExpanded += astar.GetNodesExpanded();
            if (!completed) {
                cappedAdmissibilitySearches++;
                continue;
            }

            uint64_t searchReopenedItems = 0;
            for (int itemIndex = 0; itemIndex < astar.GetNumItems(); itemIndex++) {
                if (astar.GetItem(itemIndex).reopened)
                    searchReopenedItems++;
            }
            if (searchReopenedItems != 0)
                searchesWithReopens++;
            reopenedItems += searchReopenedItems;

            uint64_t pathLength = path.empty() ? 0 : path.size() - 1;
            double hStart = heuristic.HCost(samples[i], samples[j]);
            admissibilityChecks++;
            if (hStart > pathLength + 0.000001) {
                admissibilityFailures++;
                if (admissibilityFailures <= kDiagnosticPrintLimit) {
                    printf("ADMISSIBILITY FAILURE: from=%" PRIu64
                           " to=%" PRIu64 " h=%.0f path=%" PRIu64 "\n",
                           sampleHashes[i], sampleHashes[j], hStart, pathLength);
                }
            }
        }
    }

    printf("Heuristic diagnostic complete.\n");
    printf("sampleStates=%zu targetStates=%zu\n", samples.size(), targetCount);
    printf("consistencyChecks=%" PRIu64 " consistencyFailures=%" PRIu64 "\n",
           consistencyChecks, consistencyFailures);
    printf("admissibilityChecks=%" PRIu64 " admissibilityFailures=%" PRIu64
           " cappedAdmissibilitySearches=%" PRIu64 "\n",
           admissibilityChecks, admissibilityFailures, cappedAdmissibilitySearches);
    printf("reopenDiagnostics: searchesWithReopens=%" PRIu64
           " reopenedItems=%" PRIu64 " totalExpanded=%" PRIu64 "\n",
           searchesWithReopens, reopenedItems, totalExpanded);
}

void RunTable8OverlapDiagnostic(FixedCorner &corner,
                                const RubiksCornerState &solved,
                                AnyGoalCornerHeuristic &heuristic)
{
    const size_t goalCounts[] = {2, 4, 8, 16, 32, 64, 128};
    const size_t requiredTargets = kTable8Instances * kTable8MaxGoals;
    std::vector<RubiksCornerState> targets;
    std::vector<uint64_t> baseHashes;
    std::srand(12345);

    while (targets.size() < requiredTargets) {
        uint64_t hash = kTable8BucketStart +
                        (std::rand() % (kTable8BucketEnd - kTable8BucketStart));
        if (std::find(baseHashes.begin(), baseHashes.end(), hash) != baseHashes.end())
            continue;

        RubiksCornerState base;
        RubikCornerPDB::GetStateFromHash(base, hash);
        std::vector<RubiksCornerState> localTargets;
        std::vector<uint64_t> localHashes;
        GenerateReachableSamples(corner, base, localTargets, localHashes, kTable8MaxGoals);
        if (localTargets.size() != kTable8MaxGoals)
            continue;

        targets.insert(targets.end(), localTargets.begin(), localTargets.end());
        baseHashes.push_back(hash);
    }

    TemplateAStar<RubiksCornerState, RubiksCornersAction, RubiksCorner> astar;
    astar.SetHeuristic(&heuristic);
    astar.SetReopenNodes(true);

    std::vector<uint8_t> nodesSeen(kPDB8Size, 0);
    std::vector<uint64_t> touchedNodes;

    printf("Table 8 close-goal overlap diagnostic: baseBucket=[%" PRIu64 ", %" PRIu64
           ") instances=%zu localTargetsPerInstance=%zu\n",
           kTable8BucketStart, kTable8BucketEnd,
           kTable8Instances, kTable8MaxGoals);

    for (size_t goalCount : goalCounts) {
        double ratioSum = 0.0;
        uint64_t totalExpansions = 0;
        uint64_t totalUniqueExpansions = 0;
        uint64_t completedInstances = 0;
        uint64_t emptyPaths = 0;

        for (size_t instance = 0; instance < kTable8Instances; instance++) {
            touchedNodes.clear();
            uint64_t instanceExpansions = 0;
            bool complete = true;
            size_t firstTarget = instance * kTable8MaxGoals;

            for (size_t goal = 0; goal < goalCount; goal++) {
                std::vector<RubiksCornerState> path;
                astar.GetPath(&corner, solved, targets[firstTarget + goal], path);
                if (path.empty()) {
                    emptyPaths++;
                    complete = false;
                    break;
                }

                instanceExpansions += astar.GetNodesExpanded();
                for (int item = 0; item < astar.GetNumItems(); item++) {
                    const auto &itemDetail = astar.GetItem(item);
                    if (itemDetail.where != kClosedList)
                        continue;

                    uint64_t hash = RubikCornerPDB::GetStateHash(itemDetail.data);
                    if (nodesSeen[hash] == 0) {
                        nodesSeen[hash] = 1;
                        touchedNodes.push_back(hash);
                    }
                }
            }

            uint64_t instanceUniqueExpansions = touchedNodes.size();
            for (uint64_t hash : touchedNodes)
                nodesSeen[hash] = 0;

            if (!complete || instanceUniqueExpansions == 0)
                continue;

            totalExpansions += instanceExpansions;
            totalUniqueExpansions += instanceUniqueExpansions;
            ratioSum += (double)instanceExpansions / (double)instanceUniqueExpansions;
            completedInstances++;
        }

        double averageRatio = completedInstances == 0 ? 0.0 :
            ratioSum / (double)completedInstances;
        double aggregateRatio = totalUniqueExpansions == 0 ? 0.0 :
            (double)totalExpansions / (double)totalUniqueExpansions;
        printf("table8 k=%zu instances=%" PRIu64 " emptyPaths=%" PRIu64
               " totalExpanded=%" PRIu64 " totalUniqueExpanded=%" PRIu64
               " averageRatio=%.6f aggregateRatio=%.6f\n",
               goalCount, completedInstances, emptyPaths,
               totalExpansions, totalUniqueExpansions,
               averageRatio, aggregateRatio);
    }
}

} // namespace

void RunExperiment()
{
    FixedCorner corner;
    RubiksCornerState solved;
    solved.Reset();
    RubiksCornerState solvedClean;
    solvedClean.Reset();

    std::vector<int> corners6 = {0, 1, 2, 3, 4, 5};
    RubikCornerPDB pdb6(&corner, solved, corners6);
    printf("Building 6-corner PDB...\n");
    pdb6.BuildPDB(solved, std::thread::hardware_concurrency());
    printf("Done.\n");

    AnyGoalCornerHeuristic heuristic(pdb6);
    RunHeuristicDiagnostics(corner, solvedClean, heuristic);
    if (kRunTable8OverlapDiagnostic) {
        RunTable8OverlapDiagnostic(corner, solvedClean, heuristic);
        return;
    }

#if 0 // Normal 8-corner PDB bucket fill; inactive while the overlap diagnostic is enabled.
    std::vector<uint8_t> pdb8(kPDB8Size, 0xFF);
    pdb8[RubikCornerPDB::GetStateHash(solvedClean)] = 0;

    TemplateAStar<RubiksCornerState, RubiksCornersAction, RubiksCorner> astar;
    astar.SetHeuristic(&heuristic);
    astar.SetReopenNodes(true);
    printf("Experiment setup complete; staged testing code is commented out.\n");

    std::vector<uint8_t> nodes_already_seen(kPDB8Size,0);
    uint64_t unique_nodes_expanded = 0;
    uint64_t bucket_start = 40000000;
    uint64_t bucket_end = 40100000;
    uint64_t contradictions = 0;
    uint64_t already_filled_skips = 0;
    uint64_t parity_skips = 0;
    uint64_t searches = 0;
    uint64_t empty_paths = 0;
    uint64_t new_fills = 0;
    uint64_t target_fills = 0;
    uint64_t closed_fills = 0;
    uint64_t matching_revisits = 0;
    uint64_t closed_items_seen = 0;
    uint64_t in_bucket_closed_items = 0;
    uint64_t unsafe_closed_skips = 0;
    uint64_t total_nodes_expanded = 0;
    uint64_t total_items_seen = 0;
    uint64_t min_fills_per_search = UINT64_MAX;
    uint64_t max_fills_per_search = 0;
    std::vector<uint64_t> written_by(bucket_end - bucket_start, UINT64_MAX);
    std::vector<uint64_t> written_search(bucket_end - bucket_start, 0);
    std::vector<uint8_t> written_kind(bucket_end - bucket_start, 0);

    printf("Starting bucket fill: [%" PRIu64 ", %" PRIu64 ")\n", bucket_start, bucket_end);
    auto bucket_start_time = std::chrono::steady_clock::now();

    for (int i = bucket_start; i < bucket_end; i++){
        if ((i - bucket_start) % 1000 == 0) {
            printf("bucket progress: i=%d searches=%" PRIu64
                   " newFills=%" PRIu64 " contradictions=%" PRIu64
                   " paritySkips=%" PRIu64 "\n",
                   i, searches, new_fills, contradictions, parity_skips);
        }
        if (pdb8[i] == 255) {
            RubiksCornerState target;
            RubikCornerPDB::GetStateFromHash(target, i);
            if (HasOddPermutationParity(target) == false) {
                std::vector<RubiksCornerState> path;
                astar.GetPath(&corner, solvedClean, target, path);
                if (path.empty()) {
                    empty_paths++;
                    continue;
                }
                uint64_t fills_before_search = new_fills;
                uint64_t target_distance = path.size() - 1;
                pdb8[i] = target_distance;
                written_by[i - bucket_start] = i;
                written_search[i - bucket_start] = searches + 1;
                written_kind[i - bucket_start] = 1;
                searches++;
                new_fills++;
                target_fills++;
                uint64_t nodes_expanded = astar.GetNodesExpanded();
                total_nodes_expanded += nodes_expanded;
                uint64_t total_states_seen = astar.GetNumItems();
                total_items_seen += total_states_seen;
                if (searches <= 20 || searches % 1000 == 0) {
                    printf("search %" PRIu64 ": target=%d distance=%" PRIu64
                           " expanded=%" PRIu64 " items=%" PRIu64 "\n",
                           searches, i, target_distance, nodes_expanded, total_states_seen);
                }
                double min_open_g = DBL_MAX;
                for (unsigned int open_index = 0; open_index < astar.GetNumOpenItems(); open_index++) {
                    min_open_g = std::min(min_open_g, astar.GetOpenItem(open_index).g);
                }
                for (int j = 0; j < total_states_seen; j++){
                    auto item_detail = astar.GetItem(j);
                    if (item_detail.where == kClosedList){
                        closed_items_seen++;
                        uint64_t index = RubikCornerPDB::GetStateHash(item_detail.data);
                        if (nodes_already_seen[index] == 0){
                            nodes_already_seen[index] = 1;
                            unique_nodes_expanded = unique_nodes_expanded + 1;
                        }
                        if (item_detail.g > min_open_g) {
                            unsafe_closed_skips++;
                            continue;
                        }
                        /*uint64_t index = RubikCornerPDB::GetStateHash(item_detail.data);*/
                        uint64_t g = item_detail.g;
                        if (index < bucket_start || index >= bucket_end) {
                            continue;
                        }

                        in_bucket_closed_items++;
                        if (pdb8[index] == 255){
                            pdb8[index] = g;
                            written_by[index - bucket_start] = i;
                            written_search[index - bucket_start] = searches;
                            written_kind[index - bucket_start] = 2;
                            new_fills++;
                            closed_fills++;
                        }
                        else if (pdb8[index] != g){
                            contradictions = contradictions + 1;
                            if (contradictions <= 20) {
                                uint64_t true_distance = ComputeDistanceToIndex(corner, heuristic, solvedClean, index);
                                const char *kind = written_kind[index - bucket_start] == 1 ? "target" : "closed";
                                printf("CONTRADICTION: index=%" PRIu64
                                       " old=%u new=%" PRIu64 " true=%" PRIu64
                                       " firstWrittenBy=%" PRIu64 " firstKind=%s firstSearch=%" PRIu64
                                       " currentTarget=%d currentSearch=%" PRIu64 "\n",
                                       index, (unsigned)pdb8[index], g, true_distance,
                                       written_by[index - bucket_start], kind,
                                       written_search[index - bucket_start], i, searches);
                            }
}
                        else {
                            matching_revisits++;
                        }
                    }
                }
                uint64_t fills_this_search = new_fills - fills_before_search;
                min_fills_per_search = std::min(min_fills_per_search, fills_this_search);
                max_fills_per_search = std::max(max_fills_per_search, fills_this_search);

            }
            else {
                parity_skips++;
            }
        }
        else {
            already_filled_skips++;
        }
    }

    auto bucket_end_time = std::chrono::steady_clock::now();
    double elapsed_seconds = std::chrono::duration<double>(bucket_end_time - bucket_start_time).count();
    uint64_t final_bucket_filled = 0;
    for (uint64_t i = bucket_start; i < bucket_end; i++) {
        if (pdb8[i] != 255)
            final_bucket_filled++;
    }
    double average_fills_per_search = searches == 0 ? 0.0 : (double)new_fills / (double)searches;
    if (searches == 0)
        min_fills_per_search = 0;

    double paper_expanded_ratio = 0.0;
    if (unique_nodes_expanded != 0) {
        paper_expanded_ratio =
            (double)total_nodes_expanded / (double)unique_nodes_expanded;
    }
    

    printf("Bucket fill complete.\n");
    printf("bucket=[%" PRIu64 ", %" PRIu64 ") searches=%" PRIu64
           " newFills=%" PRIu64 " contradictions=%" PRIu64 "\n",
           bucket_start, bucket_end, searches, new_fills, contradictions);
    printf("targetFills=%" PRIu64 " closedFills=%" PRIu64
           " finalBucketFilled=%" PRIu64 " bucketSize=%" PRIu64 "\n",
           target_fills, closed_fills, final_bucket_filled, bucket_end - bucket_start);
    printf("alreadyFilledSkips=%" PRIu64 " paritySkips=%" PRIu64
           " emptyPaths=%" PRIu64 "\n",
           already_filled_skips, parity_skips, empty_paths);
    printf("closedItemsSeen=%" PRIu64 " inBucketClosedItems=%" PRIu64
           " unsafeClosedSkips=%" PRIu64 " matchingRevisits=%" PRIu64 "\n",
           closed_items_seen, in_bucket_closed_items, unsafe_closed_skips, matching_revisits);
    printf("totalNodesExpanded=%" PRIu64 " totalItemsSeen=%" PRIu64 "\n",
           total_nodes_expanded, total_items_seen);
    printf("fills/search: min=%" PRIu64 " max=%" PRIu64 " average=%.3f\n",
           min_fills_per_search, max_fills_per_search, average_fills_per_search);
    printf("uniqueNodesExpanded=%" PRIu64 " paperExpandedRatio=%.6f\n",
           unique_nodes_expanded, paper_expanded_ratio);
    printf("elapsedSeconds=%.2f\n", elapsed_seconds);
#endif



#if 0
    std::vector<TestTarget> targets;
    std::vector<uint64_t> targetHashes;
    targetHashes.push_back(RubikCornerPDB::GetStateHash(solvedClean));

    auto addTarget = [&](const RubiksCornerState &state, uint64_t depth) {
        uint64_t hash = RubikCornerPDB::GetStateHash(state);
        if (std::find(targetHashes.begin(), targetHashes.end(), hash) == targetHashes.end()) {
            targets.push_back({state, hash, depth});
            targetHashes.push_back(hash);
        }
    };

    std::vector<RubiksCornerState> frontier;
    frontier.push_back(solvedClean);
    for (int depth = 1; depth <= kGeneratedTargetDepth; depth++) {
        std::vector<RubiksCornerState> nextFrontier;
        for (const auto &parent : frontier) {
            std::vector<RubiksCornerState> children;
            corner.GetSuccessors(parent, children);
            for (const auto &child : children) {
                uint64_t hash = RubikCornerPDB::GetStateHash(child);
                if (std::find(targetHashes.begin(), targetHashes.end(), hash) == targetHashes.end())
                    nextFrontier.push_back(child);
                addTarget(child, depth);
                if (targets.size() >= kGeneratedTargetLimit)
                    break;
            }
            if (targets.size() >= kGeneratedTargetLimit)
                break;
        }
        frontier.swap(nextFrontier);
        if (targets.size() >= kGeneratedTargetLimit)
            break;
    }

    uint64_t searches = 0;
    uint64_t totalExpansions = 0;
    uint64_t totalNewFills = 1;
    uint64_t contradictions = 0;
    uint64_t unsafeSkipped = 0;
    std::vector<uint64_t> fillsPerSearch;

    uint64_t heuristicOverestimates = 0;
    uint64_t distanceMismatches = 0;
    uint64_t generatedDepthMismatches = 0;

    printf("Starting generated-depth test with %zu targets up to depth %d.\n",
           targets.size(), kGeneratedTargetDepth);

    for (const auto &targetInfo : targets) {
        std::vector<RubiksCornerState> path;
        bool completed = RunAStar(astar, corner, solvedClean, targetInfo.state, path, 0);

        uint64_t closedCount = 0;
        uint64_t searchUnsafeSkipped = 0;
        uint64_t newFills = completed ? FillFromClosedList(astar, pdb8, closedCount, contradictions, searchUnsafeSkipped) : 0;
        uint64_t pathLength = path.empty() ? 0 : path.size() - 1;
        if (completed)
            FillTargetFromPath(pdb8, targetInfo.hash, pathLength, newFills, contradictions);
        uint64_t expansions = astar.GetNodesExpanded();
        double hStart = heuristic.HCost(solvedClean, targetInfo.state);
        totalExpansions += expansions;
        totalNewFills += newFills;
        unsafeSkipped += searchUnsafeSkipped;
        fillsPerSearch.push_back(newFills);
        searches++;

        printf("generated search %" PRIu64 ": targetHash=%" PRIu64
               " expectedDepth=%" PRIu64 " h=%.0f path=%" PRIu64 " targetPDB=%u"
               " expanded=%" PRIu64 " closed=%" PRIu64 " unsafeSkipped=%" PRIu64
               " newFills=%" PRIu64 "\n",
               searches, targetInfo.hash, targetInfo.expectedDepth, hStart, pathLength,
               (unsigned)pdb8[targetInfo.hash], expansions, closedCount,
               searchUnsafeSkipped, newFills);

        if (hStart > pathLength) {
            heuristicOverestimates++;
            printf("WARNING: heuristic overestimate at hash=%" PRIu64
                   " h=%.0f path=%" PRIu64 "\n",
                   targetInfo.hash, hStart, pathLength);
        }
        if (pathLength != targetInfo.expectedDepth) {
            generatedDepthMismatches++;
            printf("WARNING: generated-depth mismatch at hash=%" PRIu64
                   " expected=%" PRIu64 " path=%" PRIu64 "\n",
                   targetInfo.hash, targetInfo.expectedDepth, pathLength);
        }
        if (pdb8[targetInfo.hash] != pathLength) {
            distanceMismatches++;
            printf("WARNING: target distance mismatch at hash=%" PRIu64
                   " path=%" PRIu64 " pdb=%u\n",
                   targetInfo.hash, pathLength, (unsigned)pdb8[targetInfo.hash]);
        }
    }

    uint64_t arbitraryScanned = 0;
    uint64_t arbitraryOddSkipped = 0;
    uint64_t arbitraryFilledSkipped = 0;
    uint64_t arbitraryCompleted = 0;
    uint64_t arbitraryCapped = 0;

    printf("Starting arbitrary-hash probe: scan up to %" PRIu64
           " indices, complete up to %" PRIu64
           " searches, cap=%" PRIu64 " expansions/search.\n",
           kArbitraryScanLimit, kArbitrarySearchLimit, kArbitraryExpansionCap);

    for (uint64_t x = 0; x < kPDB8Size && arbitraryScanned < kArbitraryScanLimit &&
                         arbitraryCompleted < kArbitrarySearchLimit; x++) {
        arbitraryScanned++;

        if (pdb8[x] != 0xFF) {
            arbitraryFilledSkipped++;
            continue;
        }

        RubiksCornerState target;
        RubikCornerPDB::GetStateFromHash(target, x);
        if (HasOddPermutationParity(target)) {
            arbitraryOddSkipped++;
            continue;
        }

        std::vector<RubiksCornerState> path;
        bool completed = RunAStar(astar, corner, solvedClean, target, path, kArbitraryExpansionCap);
        uint64_t expansions = astar.GetNodesExpanded();
        double hStart = heuristic.HCost(solvedClean, target);

        if (!completed) {
            arbitraryCapped++;
            printf("arbitrary target=%" PRIu64 " h=%.0f capped after %" PRIu64
                   " expansions\n", x, hStart, expansions);
            continue;
        }

        uint64_t closedCount = 0;
        uint64_t searchUnsafeSkipped = 0;
        uint64_t newFills = FillFromClosedList(astar, pdb8, closedCount, contradictions, searchUnsafeSkipped);
        uint64_t pathLength = path.empty() ? 0 : path.size() - 1;
        FillTargetFromPath(pdb8, x, pathLength, newFills, contradictions);
        totalExpansions += expansions;
        totalNewFills += newFills;
        unsafeSkipped += searchUnsafeSkipped;
        fillsPerSearch.push_back(newFills);
        arbitraryCompleted++;
        searches++;

        printf("arbitrary search %" PRIu64 ": target=%" PRIu64
               " h=%.0f path=%" PRIu64 " targetPDB=%u expanded=%" PRIu64
               " closed=%" PRIu64 " unsafeSkipped=%" PRIu64
               " newFills=%" PRIu64 "\n",
               arbitraryCompleted, x, hStart, pathLength, (unsigned)pdb8[x],
               expansions, closedCount, searchUnsafeSkipped, newFills);

        if (hStart > pathLength) {
            heuristicOverestimates++;
            printf("WARNING: heuristic overestimate at arbitrary hash=%" PRIu64
                   " h=%.0f path=%" PRIu64 "\n", x, hStart, pathLength);
        }
        if (pdb8[x] != pathLength) {
            distanceMismatches++;
            printf("WARNING: arbitrary target distance mismatch at hash=%" PRIu64
                   " path=%" PRIu64 " pdb=%u\n", x, pathLength, (unsigned)pdb8[x]);
        }
    }

    std::sort(fillsPerSearch.begin(), fillsPerSearch.end());
    uint64_t minFills = fillsPerSearch.empty() ? 0 : fillsPerSearch.front();
    uint64_t maxFills = fillsPerSearch.empty() ? 0 : fillsPerSearch.back();

    printf("Staged test complete.\n");
    printf("searches=%" PRIu64 " generatedTargets=%zu arbitraryCompleted=%" PRIu64
           " arbitraryCapped=%" PRIu64 "\n",
           searches, targets.size(), arbitraryCompleted, arbitraryCapped);
    printf("arbitraryScanned=%" PRIu64 " arbitraryOddSkipped=%" PRIu64
           " arbitraryFilledSkipped=%" PRIu64 "\n",
           arbitraryScanned, arbitraryOddSkipped, arbitraryFilledSkipped);
    printf("totalExpansions=%" PRIu64 " totalNewFills=%" PRIu64
           " unsafeSkipped=%" PRIu64
           " contradictions=%" PRIu64 " heuristicOverestimates=%" PRIu64
           " generatedDepthMismatches=%" PRIu64
           " distanceMismatches=%" PRIu64 "\n",
           totalExpansions, totalNewFills, unsafeSkipped, contradictions,
           heuristicOverestimates, generatedDepthMismatches, distanceMismatches);
    printf("fills/search: min=%" PRIu64 " max=%" PRIu64 "\n",
           minFills, maxFills);
#endif
}
