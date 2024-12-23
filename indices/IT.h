// IT.h

#ifndef IT_H
#define IT_H

#include <vector>
#include <queue>
#include <algorithm>
#include <list>
#include <iostream>
#include "Interval.h"

using namespace std;

struct ITNode {
    int center;
    int MW;

    vector<Interval>    intervalsStartList;
    vector<Interval>    intervalsEndList;  

    ITNode* left;                                                 
    ITNode* right;                                                

    ITNode(int c) : center(c), MW(0), left(nullptr), right(nullptr) {};
};

class IT {
public:
    ITNode* root = nullptr;

    IT(const vector<Interval>& intervals) {
        auto StartList = intervals;
        sort(StartList.begin(), StartList.end(), [](const Interval& a, const Interval& b) {
            return a.start < b.start;
        });

        auto EndList = intervals;
        sort(EndList.begin(), EndList.end(), [](const Interval& a, const Interval& b) {
            return a.end > b.end;
        });
        
        root = buildIT(StartList, EndList);
    }

    IT(const IT&) = delete;
    IT& operator=(const IT&) = delete;

    IT(IT&& other) noexcept : root(other.root) {
        other.root = nullptr;
    }
    IT& operator=(IT&& other) noexcept {
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

    vector<Interval> queryMW(int start, int end, int k) {
        priority_queue<Interval, vector<Interval>, CompareWeight> minHeap;
        queryMW(root, start, end, k, minHeap);

        vector<Interval> result;
        while (!minHeap.empty()) {
            result.push_back(minHeap.top());
            minHeap.pop();
        }
        reverse(result.begin(), result.end());
        return result;
    }

private:
    ITNode* buildIT(const vector<Interval>& StartList, const vector<Interval>& EndList) {
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
        ITNode* node = new ITNode(median);
        
        for (auto interval : StartList) {
            node->MW = node->MW > interval.weight ? node->MW : interval.weight;
        }

        // Partition intervals
        vector<Interval> leftStartList;
        vector<Interval> leftEndList;

        vector<Interval> rightStartList;
        vector<Interval> rightEndList;

        for (const auto& interval : StartList) {
            if (interval.end < median) {
                leftStartList.push_back(interval);
            } else if (interval.start > median) {
                rightStartList.push_back(interval);
            } else {
                node->intervalsStartList.push_back(interval);
            }
        }

        for (const auto& interval : EndList) {
            if (interval.end < median) {
                leftEndList.push_back(interval);
            } else if (interval.start > median) {
                rightEndList.push_back(interval);
            } else {
                node->intervalsEndList.push_back(interval);
            }
        }

        // Recursively build left and right subtrees
        node->left = buildIT(leftStartList, leftEndList);
        node->right = buildIT(rightStartList, rightEndList);

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

    void query(ITNode* node, int query_start, int query_end, int k, priority_queue<Interval, vector<Interval>, CompareWeight>& minHeap) {
        int targetIndex;
        
        if (!node) return;

        if (query_end < node->center) {

            targetIndex = findUpperBound(node->intervalsStartList, query_end);

            if (targetIndex != -1) {
                for (int i = 0; i <= targetIndex; i++) {
                    Interval interval = node->intervalsStartList[i];

                    if ((int)minHeap.size() < k) {
                        minHeap.push(interval);
                    } else if (interval.weight > minHeap.top().weight) {                
                        minHeap.pop();
                        minHeap.push(interval);
                    }
                }
            }                
            
            if (node->left) {
                query(node->left, query_start, query_end, k, minHeap);
            }

        } else if (node->center < query_start) {

            targetIndex = findLowerBound(node->intervalsEndList, query_start);

            if (targetIndex != -1) {
                for (int i = 0; i <= targetIndex; i++) {
                    Interval interval = node->intervalsEndList[i];

                    if ((int)minHeap.size() < k) {
                        minHeap.push(interval);
                    } else if (interval.weight > minHeap.top().weight) {                
                        minHeap.pop();
                        minHeap.push(interval);
                    }
                }
            }

            if (node->right) {
                query(node->right, query_start, query_end, k, minHeap);
            }

        } else {

            for (auto interval : node->intervalsStartList) {
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                }
            }

            if (node->left) {
                query(node->left, query_start, query_end, k, minHeap);
            }
            if (node->right) {
                query(node->right, query_start, query_end, k, minHeap);
            }
        }

        return;
    }


    void queryMW(ITNode* node, int query_start, int query_end, int k, priority_queue<Interval, vector<Interval>, CompareWeight>& minHeap) {
        int targetIndex;
        
        if (!node) return;

        if (query_end < node->center) {

            targetIndex = findUpperBound(node->intervalsStartList, query_end);

            if (targetIndex != -1) {
                for (int i = 0; i <= targetIndex; i++) {
                    Interval interval = node->intervalsStartList[i];

                    if ((int)minHeap.size() < k) {
                        minHeap.push(interval);
                    } else if (interval.weight > minHeap.top().weight) {                
                        minHeap.pop();
                        minHeap.push(interval);
                    }
                }
            }                
            
            if (node->left && ((int)minHeap.size() < k || node->left->MW > minHeap.top().weight)) {
                query(node->left, query_start, query_end, k, minHeap);
            }

        } else if (node->center < query_start) {

            targetIndex = findLowerBound(node->intervalsEndList, query_start);

            if (targetIndex != -1) {
                for (int i = 0; i <= targetIndex; i++) {
                    Interval interval = node->intervalsEndList[i];

                    if ((int)minHeap.size() < k) {
                        minHeap.push(interval);
                    } else if (interval.weight > minHeap.top().weight) {                
                        minHeap.pop();
                        minHeap.push(interval);
                    }
                }
            }

            if (node->right && ((int)minHeap.size() < k || node->right->MW > minHeap.top().weight)) {
                queryMW(node->right, query_start, query_end, k, minHeap);
            }

        } else {

            for (auto interval : node->intervalsStartList) {
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                }
            }

            if (node->left && ((int)minHeap.size() < k || node->left->MW > minHeap.top().weight)) {
                queryMW(node->left, query_start, query_end, k, minHeap);
            }
            if (node->right && ((int)minHeap.size() < k || node->right->MW > minHeap.top().weight)) {
                queryMW(node->right, query_start, query_end, k, minHeap);
            }
        }

        return;
    }


};

#endif
