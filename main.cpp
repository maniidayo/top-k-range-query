#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include <chrono>
#include <random>
#include <limits>   
#include <time.h>
#include <math.h>
#include <assert.h>
#include <malloc.h>
#include <vector>
#include <iomanip>
#include <sys/resource.h>

#include "getopt.h"
#include "def_global.h"
#include "containers/relation.h"
#include "indices/hint_m.h"
#include "indices/WIT.h"
#include "indices/IT.h"
#include "indices/SAIT.h"

using namespace std;

int BUCKET_SIZE;
int bucketnum;

int MEMCHECK = 0;

int isWeighted = 1;
ifstream file("/home/lee/top-k/samples/RENFE.csv");
ofstream outfile("R_RENFE.csv");

RunSettings settings;

HierarchicalIndex* idxR;

void printMemoryUsage() {
    std::ifstream status_file("/proc/self/status");
    std::string line;
    while (std::getline(status_file, line)) {
        if (line.find("VmSize:") != std::string::npos || line.find("VmRSS:") != std::string::npos) {
            outfile << line << endl;
        }
    }
}

int randKey(int floor, int ceiling) {
    int range = ceiling - floor;
    return floor + range * ((double) rand() / (double) (RAND_MAX + 1.0));
}

Interval randomQuery(int maxStart, int maxStop, double A, double B) {
    int length;
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> dist(A, B);
    // poisson_distribution<> dist(1000);
    do {
        length = static_cast<int>(dist(gen));
    } while (length < 0 || length > maxStop);
    int start = randKey(0, maxStop - length);
    int stop = start + length;
    int weight = 0;
    return Interval(0, start, stop, 0);
}

bool compareIntervals(const Interval& a, const Interval& b) {
    return a.weight > b.weight;
}

bool compareRecords(const Record& a, const Record& b) {
    return a.weight > b.weight;
}



/* Building Section */





vector<Interval> BuildBruteForce(vector<Interval> intervals) {
    typedef chrono::high_resolution_clock Clock;
    typedef chrono::milliseconds milliseconds;
  
    Clock::time_point t0, t1;
    t0 = Clock::now();    
    sort(intervals.begin(), intervals.end(), compareIntervals);
    t1 = Clock::now();

    milliseconds ms1 = chrono::duration_cast<milliseconds>(t1 - t0);
    outfile << ms1.count() << endl;

    return intervals;
}

IT BuildITree(vector<Interval> intervals) {
    typedef chrono::high_resolution_clock Clock;
    typedef chrono::milliseconds milliseconds;
  
    Clock::time_point t0, t1;

    t0 = Clock::now();
    IT tree(intervals);
    t1 = Clock::now();    

    milliseconds ms1 = chrono::duration_cast<milliseconds>(t1 - t0);
    outfile << ms1.count() << endl;

    return tree;
}

void BuildHINT(vector<Interval> intervals, int m) {

    typedef chrono::high_resolution_clock Clock;
    typedef chrono::milliseconds milliseconds;
  
    Clock::time_point t0, t1;
    
    Relation R;

    // Intervals -> Record (from HINT)

    t0 = Clock::now();

    // free(idxR);

    size_t sum = 0;

    for (auto interval : intervals)
    {        
        R.emplace_back(Record(interval.id, interval.start, interval.end, interval.weight));

        R.gstart = std::min(R.gstart, interval.start);
        R.gend   = std::max(R.gend  , interval.end);
        R.longestRecord = std::max(R.longestRecord, interval.end - interval.start + 1);
        sum += interval.end - interval.start;
    }

    R.avgRecordExtent = (float)sum/R.size();
    settings.maxBits = int(log2(R.gend - R.gstart)+1);    
    settings.numBits = determineOptimalNumBitsForHINT_M(R, 0.1);

    idxR = new HINT_M_SubsSort(R, settings.numBits, settings.maxBits);

    t1 = Clock::now();

    milliseconds ms1 = chrono::duration_cast<milliseconds>(t1 - t0);
    outfile << ms1.count() << endl;
}

WIT BuildWITree(vector<Interval> intervals) {
    typedef chrono::high_resolution_clock Clock;
    typedef chrono::milliseconds milliseconds;
  
    Clock::time_point t0, t1;

    t0 = Clock::now();
    WIT tree(intervals);
    t1 = Clock::now();    

    milliseconds ms1 = chrono::duration_cast<milliseconds>(t1 - t0);
    outfile << ms1.count() << endl;

    return tree;
}

vector<WIT> BuildWITreeBucket(vector<Interval> intervals) {    
    typedef chrono::high_resolution_clock Clock;
    typedef chrono::milliseconds milliseconds;
  
    Clock::time_point t0, t1;
    vector<WIT> intervalForest;

    t0 = Clock::now();

    int intervalNum = (int)intervals.size();
    sort(intervals.begin(), intervals.end(), compareIntervals);

    auto it = intervals.begin();
    bucketnum = intervalNum / BUCKET_SIZE;
    int remainder  = intervalNum % BUCKET_SIZE;
    vector<Interval> temp;

    for (int i = 0; i < bucketnum; ++i) {
        for (int j = 0; j < BUCKET_SIZE; ++j) {
            temp.push_back(std::move(*(it++)));
        }
        intervalForest.emplace_back(temp);
        temp.clear();
    }

    if (remainder != 0) {
        for (int j = 0; j < remainder; ++j) {
            temp.push_back(std::move(*(it++)));
        }
        intervalForest.emplace_back(temp);
        bucketnum++;
    }

    t1 = Clock::now();

    milliseconds ms1 = chrono::duration_cast<milliseconds>(t1 - t0);
    outfile << ms1.count() << endl;

    return intervalForest;
}

SAIT BuildSAIT(vector<Interval> intervals) {
    typedef chrono::high_resolution_clock Clock;
    typedef chrono::milliseconds milliseconds;
  
    Clock::time_point t0, t1;

    t0 = Clock::now();
    SAIT tree(intervals);
    t1 = Clock::now();    

    milliseconds ms1 = chrono::duration_cast<milliseconds>(t1 - t0);
    outfile << ms1.count() << endl;

    return tree;
}




/* Query Section */





vector<vector<Interval>> QueryBruteForce(vector<Interval>& intervals, const vector<Interval>& queries, int k) {   
    typedef chrono::high_resolution_clock Clock;
    typedef chrono::milliseconds milliseconds;
  
    Clock::time_point t0, t1;

    vector<vector<Interval>> bruteforceResults;

    t0 = Clock::now();
    for (const auto& q : queries) {
        vector<Interval> results;
        int count = 0;
        for (const auto& i : intervals) {
            if (i.end >= q.start && i.start <= q.end) {
                results.push_back(i);
                if (++count == k) break;
            }
        }
        bruteforceResults.push_back(results);
    }
    t1 = Clock::now();

    milliseconds ms1 = chrono::duration_cast<milliseconds>(t1 - t0);
    outfile << ms1.count() << ",";
    return bruteforceResults;
}

vector<vector<Interval>> QueryITree(IT& ITree, const vector<Interval>& queries, int k) {
    typedef chrono::high_resolution_clock Clock;
    typedef chrono::milliseconds milliseconds;
  
    Clock::time_point t0, t1;
    vector<vector<Interval>> ITResults;

    t0 = Clock::now();

    for (const auto& q : queries) {
        auto results = ITree.query(q.start, q.end, k);
        sort(results.begin(), results.end(), compareIntervals);  // 정렬 추가
        ITResults.push_back(results);
    }

    t1 = Clock::now();

    milliseconds ms1 = chrono::duration_cast<milliseconds>(t1 - t0);
    outfile << ms1.count() << ",";
    return ITResults;
}

vector<vector<Interval>> QueryITreeMW(IT& ITree, const vector<Interval>& queries, int k) {
    typedef chrono::high_resolution_clock Clock;
    typedef chrono::milliseconds milliseconds;
  
    Clock::time_point t0, t1;
    vector<vector<Interval>> ITMWResults;

    t0 = Clock::now();

    for (const auto& q : queries) {
        auto results = ITree.queryMW(q.start, q.end, k);
        sort(results.begin(), results.end(), compareIntervals);  // 정렬 추가
        ITMWResults.push_back(results);
    }

    t1 = Clock::now();

    milliseconds ms1 = chrono::duration_cast<milliseconds>(t1 - t0);
    outfile << ms1.count() << ",";
    return ITMWResults;
}

vector<vector<Interval>> QueryWITree(WIT& WITree, const vector<Interval>& queries, int k) {
    typedef chrono::high_resolution_clock Clock;
    typedef chrono::milliseconds milliseconds;
  
    Clock::time_point t0, t1;
    vector<vector<Interval>> WITResults;

    t0 = Clock::now();

    for (const auto& q : queries) {
        auto results = WITree.queryTopK2(q.start, q.end, k);
        sort(results.begin(), results.end(), compareIntervals);  // 정렬 추가
        WITResults.push_back(results);
    }

    t1 = Clock::now();

    milliseconds ms1 = chrono::duration_cast<milliseconds>(t1 - t0);
    
    outfile << ms1.count() << "," ;
    return WITResults;    
}

vector<vector<Interval>> QueryWITreeBucket(vector<WIT>& WIForest, const vector<Interval>& queries, int k) {
    typedef chrono::high_resolution_clock Clock;
    typedef chrono::milliseconds milliseconds;
  
    Clock::time_point t0, t1;
    vector<vector<Interval>> WITBResults;

    t0 = Clock::now();

    for (const auto& q : queries) {
        vector<Interval> combinedResults;
        int bucketIndex = 0, kk = k;
        do {
            auto results = WIForest[bucketIndex++].queryTopK2(q.start, q.end, kk);
            combinedResults.insert(combinedResults.end(), std::make_move_iterator(results.begin()), std::make_move_iterator(results.end()));
            kk -= results.size();
        } while (bucketIndex < bucketnum && kk > 0);
        WITBResults.push_back(combinedResults);
    }

    t1 = Clock::now();

    milliseconds ms1 = chrono::duration_cast<milliseconds>(t1 - t0);
    outfile << ms1.count() << "," ;
    return WITBResults;
}

vector<vector<Interval>> QuerySAIT(SAIT& SAITree, const vector<Interval>& queries, int k) {
    typedef chrono::high_resolution_clock Clock;
    typedef chrono::milliseconds milliseconds;
  
    Clock::time_point t0, t1;
    vector<vector<Interval>> SAITResults;

    t0 = Clock::now();

    for (const auto& q : queries) {
        auto results = SAITree.query(q.start, q.end, k);
        sort(results.begin(), results.end(), compareIntervals);  // 정렬 추가
        SAITResults.push_back(results);
    }

    t1 = Clock::now();

    milliseconds ms1 = chrono::duration_cast<milliseconds>(t1 - t0);
    outfile << ms1.count();
    return SAITResults;
}

vector<vector<Interval>> QueryHINT(const vector<Interval>& queries, int k) {
    typedef chrono::high_resolution_clock Clock;
    typedef chrono::milliseconds milliseconds;
  
    Clock::time_point t0, t1;

    vector<vector<Interval>> HINTResults;
    vector<vector<Record>> HINTResultsOriginal;

    vector<Record> queryresult;
    
    int count = 0;

    t0 = Clock::now();

    for (const auto& q : queries) {

        queryresult = idxR->executeBottomUp_gOverlaps(RangeQuery(++count, q.start, q.end), k);

        sort(queryresult.begin(), queryresult.end(), compareRecords);

        vector<Record> topKResult;

        if (queryresult.size() > k){
            topKResult.insert(topKResult.end(), queryresult.begin(), queryresult.begin() + k);
        } else {
            topKResult.insert(topKResult.end(), queryresult.begin(), queryresult.end());
        }

        HINTResultsOriginal.emplace_back(topKResult);
    }

    t1 = Clock::now();

    int iter = 0;

    for (auto results : HINTResultsOriginal) {
        vector<Interval> intervals;
        for (auto result : results) {
            intervals.emplace_back(Interval(result.id, result.start, result.end, result.weight));
        }
        HINTResults.emplace_back(intervals);
    }

    milliseconds ms1 = chrono::duration_cast<milliseconds>(t1 - t0);
    outfile << ms1.count() << "," ;

    return HINTResults;
}



/* Test Section */





void TEST(int k, vector<Interval>& intervals,
            vector<Interval>& queries, 
            vector<Interval>& sortedIntervals,
            IT& ITree,
            WIT& WITree,
            vector<WIT>& WIForest,
            SAIT& SAITree) {

    auto BFResults      = QueryBruteForce(sortedIntervals, queries, k);
    auto ITResults      = QueryITree(ITree, queries, k);
    auto ITMWResults    = QueryITreeMW(ITree, queries, k);
    auto HINTResults    = QueryHINT(queries, k);
    auto WITResults     = QueryWITree(WITree, queries, k);   
    auto WITBResults    = QueryWITreeBucket(WIForest, queries, k);   
    auto SAITResults    = QuerySAIT(SAITree, queries, k);   

    for (size_t i = 0; i < queries.size(); ++i) {
        auto& BFTopK = BFResults[i];
        auto& ITTopK = ITResults[i];
        auto& ITMWTopK = ITMWResults[i];
        auto& HINTTopK = HINTResults[i];
        auto& WITTopK = WITResults[i];
        auto& WITBTopK = WITBResults[i];
        auto& SAITTopK = SAITResults[i];

        int actualk = (int)BFResults[i].size();

        ITTopK.resize(actualk);
        ITMWTopK.resize(actualk);
        HINTTopK.resize(actualk);
        WITTopK.resize(actualk);
        WITBTopK.resize(actualk);
        SAITTopK.resize(actualk);

        sort(ITTopK.begin(), ITTopK.end(), compareIntervals);
        sort(ITMWTopK.begin(), ITMWTopK.end(), compareIntervals);
        sort(HINTTopK.begin(), HINTTopK.end(), compareIntervals);
        sort(WITTopK.begin(), WITTopK.end(), compareIntervals);
        sort(WITBTopK.begin(), WITBTopK.end(), compareIntervals);
        sort(SAITTopK.begin(), SAITTopK.end(), compareIntervals);
         
        for (int j = 0; j < actualk; ++j) {
            auto weight1 = BFTopK[j].weight;
            auto weight2 = ITTopK[j].weight;
            auto weight3 = ITMWTopK[j].weight;
            auto weight4 = HINTTopK[j].weight;
            auto weight5 = WITTopK[j].weight;
            auto weight6 = WITBTopK[j].weight;
            auto weight7 = SAITTopK[j].weight;

            if(weight1 != weight2 || 
               weight1 != weight3 || 
               weight1 != weight4 || 
               weight1 != weight5 || 
               weight1 != weight6 || 
               weight1 != weight7) {

                cout << queries[i] << endl;                
                cout << BFTopK << endl;
                cout << ITTopK << endl;
                cout << ITMWTopK << endl;
                cout << HINTTopK << endl;
                cout << WITTopK << endl;
                cout << WITBTopK << endl;
                cout << SAITTopK << endl;
            }

            assert(weight1 == weight2);
            assert(weight1 == weight3);
            assert(weight1 == weight4);
            assert(weight1 == weight5);
            assert(weight1 == weight6);
            assert(weight1 == weight7);
        }
    }
}

int main() {    
    settings.init();
    settings.method = "hint_m";

    outfile << "n,k,m,BruteForce,IT,ITMW,HINTm,WIT,WIT+Bucket,SAIT" << endl;
    outfile << "0,";

    vector<Interval> intervals;
    vector<Interval> queries;

    int count = 0;

    int min = 10000000000, max = 0, sum = 0;

    double m = 0.01; // query extent

    random_device           rd;
    mt19937                 gen(rd());
    poisson_distribution<>  distrib(100000);

    if(true == file.fail()) {
        cout << "no file found" << endl;
        return 1;
    }

    string line;
    string cell;

    cout << "csv start" << endl;
    
    if(isWeighted == 0) {
        while (getline(file, line)) {
            count++;

            stringstream lineStream(line);

            int start, end, weight;

            if (getline(lineStream, cell, ',')) 
                start = stoi(cell);

            if (getline(lineStream, cell, ',')) 
                end = stoi(cell);

            min = std::min(start, min);            
            max = std::max(end, max);

            weight = distrib(gen);

            intervals.push_back(Interval(count, start, end, weight));

            sum += (end - start);
        }
    }

    else {
        while (getline(file, line)) {
            count++;

            stringstream lineStream(line);

            int start, end, weight;

            if (getline(lineStream, cell, ',')) 
                start = stoi(cell);

            if (getline(lineStream, cell, ',')) 
                end = stoi(cell);

            if (getline(lineStream, cell, ',')) 
                weight = stoi(cell);

            min = std::min(start, min);            
            max = std::max(end, max);

            intervals.push_back(Interval(count, start, end, weight));

            sum += (end - start);
        }
    }

    // Construction


    int lengthAVG = (sum / count);

    BUCKET_SIZE = sqrt(count);

    printMemoryUsage();

    vector<Interval> sortedIntervals = BuildBruteForce(intervals);    
    cout << "BF build finished" << endl;

    printMemoryUsage();

    IT ITree = BuildITree(intervals);    
    cout << "IT build finished" << endl;

    printMemoryUsage();

    BuildHINT(intervals, m * 100);
    cout << "HINTm build finished" << endl;
    
    printMemoryUsage();    

    WIT WITree = BuildWITree(intervals);    
    cout << "WIT build finished" << endl;

    printMemoryUsage();    

    vector<WIT> WIForest = BuildWITreeBucket(intervals);    
    cout << "WITB build finished" << endl;

    printMemoryUsage();

    SAIT SAITree = BuildSAIT(intervals);    
    cout << "SAIT build finished" << endl;

    printMemoryUsage();    

    if(MEMCHECK == 1){
        return 0;
    }
      
    int queryNum = 10000;
    
    for (int i = 0; i < queryNum; ++i) {            
        queries.push_back(randomQuery(max, max, m * max, 0.1));
    }
    
    // Effect of k

    
    for (int k = 5; k < 50; k+=5) {
        outfile << intervals.size() << ",";
        outfile << k << ",";
        outfile << m << ",";
        TEST(k, intervals, queries, sortedIntervals, ITree, WITree, WIForest, SAITree);
        outfile << endl;
    }

    double temp1[5] = {0.001, 0.005, 0.01, 0.05, 0.1};
    int k = 5;

    // Effect of query extent

    
    for (int i = 0; i < 5; i++) {
        m = temp1[i];

        queries.clear();

        for (int i = 0; i < queryNum; ++i) {            
            queries.push_back(randomQuery(max, max, m * max, 0.1));
        }

        outfile << intervals.size() << ",";
        outfile << k << ",";
        outfile << m << ",";
        TEST(k, intervals, queries, sortedIntervals, ITree, WITree, WIForest, SAITree);
        outfile << endl;
    }

    double temp2[4] = {0.2, 0.4, 0.6, 0.8};
    m = 0.01;
    k = 5;


    // Effect of dataset size


    for (int i = 0; i < 4; i++) {

        int DataSize = temp2[i] * intervals.size();        

        vector<Interval> intervals_new;

        sample(intervals.begin(), intervals.end(), back_inserter(intervals_new), DataSize, mt19937{rd()});

        max = 0;
        sum = 0;

        for(auto interval : intervals_new) {
            int start = interval.start;
            int end = interval.end;
            max = (end > max) ? end : max;
            sum += (end - start);
        }

        queries.clear();

        for (int i = 0; i < queryNum; ++i) {            
            queries.push_back(randomQuery(max, max, m * max, 0.1));
        }

        printMemoryUsage();

        vector<Interval> sortedIntervals = BuildBruteForce(intervals_new);    
        cout << "BF build finished" << endl;

        printMemoryUsage();

        IT ITree = BuildITree(intervals_new);    
        cout << "IT build finished" << endl;

        printMemoryUsage();

        BuildHINT(intervals_new, m * 100);
        cout << "HINTm build finished" << endl;

        printMemoryUsage();    

        WIT WITree = BuildWITree(intervals_new);    
        cout << "WIT build finished" << endl;

        printMemoryUsage();    

        vector<WIT> WIForest = BuildWITreeBucket(intervals_new);    
        cout << "WITB build finished" << endl;

        printMemoryUsage();

        SAIT SAITree = BuildSAIT(intervals_new);    
        cout << "SAIT build finished" << endl;

        printMemoryUsage();    

        outfile << intervals_new.size() << ",";
        outfile << k << ",";
        outfile << m << ",";
        TEST(k, intervals, queries, sortedIntervals, ITree, WITree, WIForest, SAITree);
        outfile << endl;
    }
  
    /*
    for (BUCKET_SIZE = 5000; BUCKET_SIZE < 100000; BUCKET_SIZE += 5000) {
        outfile << BUCKET_SIZE << ",";
        TEST2(6, intervals, queries, BUCKET_SIZE);
        outfile << endl;
    }
    */

    return 0;
}