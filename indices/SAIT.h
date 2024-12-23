// SAIT.h

#ifndef SAIT_H
#define SAIT_H

#include <vector>
#include <queue>
#include <algorithm>
#include <list>
#include <iostream>
#include "Interval.h"

using namespace std;

struct SegmentTreeNode {
    int startIndex; 
    int endIndex;
    list<Interval> intervals;
    SegmentTreeNode* left;
    SegmentTreeNode* right; 

    SegmentTreeNode(int start, int end) : startIndex(start), endIndex(end), left(nullptr), right(nullptr) {}
};

class SegmentTree {
public:
    SegmentTreeNode* root;

    SegmentTree() : root(nullptr) {}

    SegmentTree(const vector<Interval>& intervals) {
        root = buildSegmentTree(intervals, 0, intervals.size() - 1);
    }

private:
    SegmentTreeNode* buildSegmentTree(const vector<Interval>& intervals, int start, int end) {
        if (start > end) return nullptr;

        SegmentTreeNode* node = new SegmentTreeNode(start, end);

        if (start == end) {
            node->intervals.push_back(intervals[start]);
        } else {
            int mid = start + (end - start) / 2;
            node->left = buildSegmentTree(intervals, start, mid);
            node->right = buildSegmentTree(intervals, mid + 1, end);

            node->intervals = mergeSortedIntervals(node->left->intervals, node->right->intervals);
        }

        return node;
    }

    list<Interval> mergeSortedIntervals(const list<Interval>& leftIntervals, const list<Interval>& rightIntervals) {
        list<Interval> mergedIntervals;
        auto itLeft = leftIntervals.begin();
        auto itRight = rightIntervals.begin();

        while (itLeft != leftIntervals.end() && itRight != rightIntervals.end()) {
            if (itLeft->weight > itRight->weight) {
                mergedIntervals.push_back(*itLeft);
                ++itLeft;
            } else {
                mergedIntervals.push_back(*itRight);
                ++itRight;
            }
        }

        while (itLeft != leftIntervals.end()) {
            mergedIntervals.push_back(*itLeft);
            ++itLeft;
        }
        while (itRight != rightIntervals.end()) {
            mergedIntervals.push_back(*itRight);
            ++itRight;
        }

        return mergedIntervals;
    }
};

struct SAITNode {
    int center;                                                     
    int maxWeight;

    vector<Interval>    intervalsInsideWeight;                      

    vector<Interval>    intervalsInsideStartList;
    vector<Interval>    intervalsInsideEndList;
    vector<Interval>    intervalsOutsideStartList;
    vector<Interval>    intervalsOutsideEndList;

    SegmentTree         intervalsInsideStartTree;         
    SegmentTree         intervalsInsideEndTree;           
    SegmentTree         intervalsOutsideStartTree;    
    SegmentTree         intervalsOutsideEndTree;       

    SAITNode* left;                                                 
    SAITNode* right;                                                

    SAITNode(int c) : center(c), left(nullptr), right(nullptr) {};
};

class SAIT {
public:
    SAITNode* root = nullptr;

    SAIT(const vector<Interval>& intervals) {
        auto StartList = intervals;
        sort(StartList.begin(), StartList.end(), [](const Interval& a, const Interval& b) {
            return a.start < b.start;
        });

        auto EndList = intervals;
        sort(EndList.begin(), EndList.end(), [](const Interval& a, const Interval& b) {
            return a.end > b.end;
        });

        auto WeightList = intervals;
        sort(WeightList.begin(), WeightList.end(), [](const Interval& a, const Interval& b) {
            return a.weight > b.weight;
        });
        
        root = buildSAIT(StartList, EndList, WeightList, 1);
    }

    SAIT(const SAIT&) = delete;
    SAIT& operator=(const SAIT&) = delete;

    SAIT(SAIT&& other) noexcept : root(other.root) {
        other.root = nullptr;
    }
    SAIT& operator=(SAIT&& other) noexcept {
        if (this != &other) {
            root = other.root;
            other.root = nullptr;
        }
        return *this;
    }
    vector<Interval> query(int start, int end, int k) {
        priority_queue<Interval, vector<Interval>, CompareWeight> minHeap;
        query(root, start, end, k, minHeap);

        vector<Interval> result;
        while (!minHeap.empty()) {
            result.push_back(minHeap.top());
            minHeap.pop();
        }
        reverse(result.begin(), result.end());
        return result;
    }

private:
    SAITNode* buildSAIT(const vector<Interval>& StartList, const vector<Interval>& EndList, const vector<Interval>& WeightList, int flag) {
        if (StartList.empty()) {
            return nullptr;
        }

        // Compute midpoints
        vector<int> midpoints(StartList.size());
        for (size_t i = 0; i < StartList.size(); ++i) {
            midpoints[i] = StartList[i].start + (StartList[i].end - StartList[i].start) / 2;
        }

        size_t mid = midpoints.size() / 2;
        nth_element(midpoints.begin(), midpoints.begin() + mid, midpoints.end());
        int median = midpoints[mid];
       
        // Create node
        SAITNode* node = new SAITNode(median);
        
        // OutsideStartList and OutsideStartTree
        node->intervalsOutsideStartList = StartList;
        if(flag == 0) node->intervalsOutsideStartTree = SegmentTree(node->intervalsOutsideStartList);

        // OutsideEndList and OutsideEndTree
        node->intervalsOutsideEndList = EndList;
        if(flag == 0) node->intervalsOutsideEndTree = SegmentTree(node->intervalsOutsideEndList);

        // Partition intervals
        vector<Interval> leftStartList;
        vector<Interval> leftEndList;
        vector<Interval> leftWeightList;

        vector<Interval> rightStartList;
        vector<Interval> rightEndList;
        vector<Interval> rightWeightList;

        for (const auto& interval : StartList) {
            if (interval.end < median) {
                leftStartList.push_back(interval);
            } else if (interval.start > median) {
                rightStartList.push_back(interval);
            } else {
                node->intervalsInsideStartList.push_back(interval);
            }
        }

        for (const auto& interval : EndList) {
            if (interval.end < median) {
                leftEndList.push_back(interval);
            } else if (interval.start > median) {
                rightEndList.push_back(interval);
            } else {
                node->intervalsInsideEndList.push_back(interval);
            }
        }

        for (const auto& interval : WeightList) {
            if (interval.end < median) {
                leftWeightList.push_back(interval);
            } else if (interval.start > median) {
                rightWeightList.push_back(interval);
            } else {
                node->intervalsInsideWeight.push_back(interval);
            }
        }

        node->intervalsInsideStartTree = SegmentTree(node->intervalsInsideStartList);
        node->intervalsInsideEndTree = SegmentTree(node->intervalsInsideEndList);

        // Recursively build left and right subtrees
        node->left = buildSAIT(leftStartList, leftEndList, leftWeightList, 0);
        node->right = buildSAIT(rightStartList, rightEndList, rightWeightList, 0);

        return node;
    }

    int findUpperBound(const vector<Interval>& intervals, int target) {
        int left = 0;
        int right = intervals.size() - 1;
        int result = -1;

        while (left <= right) {
            int mid = left + (right - left) / 2;

            if (intervals[mid].start <= target) {
                result = mid;  // This index satisfies the condition
                left = mid + 1;  // Search in the right half for a possible larger index
            } else {
                right = mid - 1;  // Search in the left half
            }
        }
        
        return result;
    }

    int findLowerBound(const vector<Interval>& intervals, int target) {
        int left = 0;
        int right = intervals.size() - 1;
        int result = -1;

        while (left <= right) {
            int mid = left + (right - left) / 2;

            if (intervals[mid].end >= target) {
                result = mid;  // This index satisfies the condition
                left = mid + 1;  // Search in the right half for a possible larger index
            } else {
                right = mid - 1;  // Search in the left half
            }
        }

        
        return result;
    }
    
    void collectRelevantNodes(SegmentTreeNode* node, int index1, int index2, vector<SegmentTreeNode*>& nodeList) {
        int temp;
        if (!node) {
            return;
        } else if (node->startIndex == node->endIndex || (node->startIndex == index1 && node->endIndex == index2)) {
            nodeList.push_back(node);
        } else if (node->left->endIndex < index2) {            
            nodeList.push_back(node->left);      
            collectRelevantNodes(node->right, node->right->startIndex, index2, nodeList);
        } else if (node->left->endIndex >= index2) {
            collectRelevantNodes(node->left, index1, index2, nodeList);
        }
        return;
    }

    void searchSegmentTree(SegmentTree& tree, int index, int k, priority_queue<Interval, vector<Interval>, CompareWeight>& minHeap) {
        vector<SegmentTreeNode*> nodeList;
        collectRelevantNodes(tree.root, 0, index, nodeList);
        
        /*int temp;
        cout << "FINAL LIST" << endl;
        cout << "from 0 to " << index << endl;
        for(auto node : nodeList) {
            cout << node->startIndex << " and " << node->endIndex << endl;
        }
        cin >> temp;*/

        struct HeapNode {
            Interval interval;
            SegmentTreeNode* node;
            list<Interval>::iterator it;
            bool operator<(const HeapNode& other) const {
                return interval.weight < other.interval.weight;
            }
        };

        priority_queue<HeapNode> nodeHeap;

        for (SegmentTreeNode* node : nodeList) {
            if (!node->intervals.empty()) {
                nodeHeap.push({*(node->intervals.begin()), node, node->intervals.begin()});
            }
        }

        while (!nodeHeap.empty()) {
            HeapNode top = nodeHeap.top();
            nodeHeap.pop();

            if ((int)minHeap.size() < k) {
                minHeap.push(top.interval);
            } else if (top.interval.weight > minHeap.top().weight) {
                minHeap.pop();
                minHeap.push(top.interval);
            } else {
                break;
            }

            ++(top.it);
            if (top.it != top.node->intervals.end()) {
                nodeHeap.push({*(top.it), top.node, top.it});
            }
        }
    }

    /* void searchSegmentTree(SegmentTree& tree, int index, int k, priority_queue<Interval, vector<Interval>, CompareWeight>& minHeap) {
        vector<SegmentTreeNode*> nodeList;
        collectRelevantNodes(tree.root, index, nodeList);

        struct HeapNode {
            Interval interval;
            SegmentTreeNode* node;
            list<Interval>::iterator it;
            bool operator<(const HeapNode& other) const {
                return interval.weight < other.interval.weight;
            }
        };

        priority_queue<HeapNode> NodeHeap; // Max-Heap

        for (SegmentTreeNode* node : nodeList) {
            if (!node->intervals.empty()) {
                NodeHeap.push({*(node->intervals.begin()), node, node->intervals.begin()});
            }
        }

        while (!NodeHeap.empty()) {
            HeapNode top = NodeHeap.top();
            NodeHeap.pop();

            if ((int)minHeap.size() < k) {
                minHeap.push(top.interval);
            } else if (top.interval.weight > minHeap.top().weight) {                
                minHeap.pop();
                minHeap.push(top.interval);
            } else {
                break;
            }

            ++(top.it);
            if (top.it != top.node->intervals.end()) {
                NodeHeap.push({*(top.it), top.node, top.it});
            }
        }
    } */

    void query(SAITNode* node, int query_start, int query_end, int k, priority_queue<Interval, vector<Interval>, CompareWeight>& minHeap) {
        int targetIndex;
        
        if (!node) return;

        if (query_end < node->center) {
            targetIndex = findUpperBound(node->intervalsInsideStartList, query_end);
            if (targetIndex != -1) {
                searchSegmentTree(node->intervalsInsideStartTree, targetIndex, k, minHeap);
            }
            if (node->left) {
                query(node->left, query_start, query_end, k, minHeap);
            }
        } else if (node->center < query_start) {
            targetIndex = findLowerBound(node->intervalsInsideEndList, query_start);
            if (targetIndex != -1) {
                searchSegmentTree(node->intervalsInsideEndTree, targetIndex, k, minHeap);
            }
            if (node->right) {
                query(node->right, query_start, query_end, k, minHeap);
            }
        } else {
            for (const auto& interval : node->intervalsInsideWeight) {
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                } else {
                    break;
                }
            }
            if (node->left) {
                int targetIndex = findLowerBound(node->left->intervalsOutsideEndList, query_start);
                if (targetIndex != -1) {
                    searchSegmentTree(node->left->intervalsOutsideEndTree, targetIndex, k, minHeap);
                }
            }
            if (node->right) {
                int targetIndex = findUpperBound(node->right->intervalsOutsideStartList, query_end);
                if (targetIndex != -1) {
                    searchSegmentTree(node->right->intervalsOutsideStartTree, targetIndex, k, minHeap);
                }
            }
        }
        return;
    }
};

#endif
