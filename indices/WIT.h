#ifndef WIT_H
#define WIT_H

#include <vector>
#include <algorithm>
#include <queue>
#include "Interval.h"

using namespace std;

struct WITNode {
    Interval interval;
    int min;
    int max;
    int MW;
    WITNode* left;
    WITNode* right;
    WITNode(const Interval& i)
    {
        interval = i;
        min = i.start;
        max = i.end;
        MW = i.weight;
        left = nullptr;
        right = nullptr;
    }
};

class WIT {
public:
    WIT(const vector<Interval>& intervals) {
        for (const auto& interval : intervals) {
            roots.push_back(interval);
        }
        sort(roots.begin(), roots.end(), [](const Interval& a, const Interval& b) {
            return a.start < b.start;
        });
        build();
    }

    // 복사 생성자와 복사 대입 연산자 삭제
    WIT(const WIT&) = delete;
    WIT& operator=(const WIT&) = delete;

    // 이동 생성자와 이동 대입 연산자 정의 (필요한 경우)
    WIT(WIT&& other) noexcept : roots(std::move(other.roots)), root(other.root) {
        other.root = nullptr;
    }
    WIT& operator=(WIT&& other) noexcept {
        if (this != &other) {
            deleteTree(root);
            roots = std::move(other.roots);
            root = other.root;
            other.root = nullptr;
        }
        return *this;
    }
   
    ~WIT() {
        deleteTree(root);
    }

    vector<Interval> query(int start, int end) {
        vector<Interval> result;
        query(root, start, end, result);
        return result;
    }

    vector<Interval> queryTopK(int start, int end, int k) {
        priority_queue<Interval, vector<Interval>, CompareWeight> minHeap;
        queryTopK(root, start, end, k, minHeap);

        vector<Interval> result;
        while (!minHeap.empty()) {
            result.push_back(minHeap.top());
            minHeap.pop();
        }
        reverse(result.begin(), result.end());
        return result;
    }

    vector<Interval> queryTopK_H(int start, int end, int k) {
        priority_queue<Interval, vector<Interval>, CompareWeight> minHeap;        
        priority_queue<WITNode*, vector<WITNode*>, CompareMW> MWHeap;
        
        queryTopK_H(root, start, end, k, minHeap, MWHeap);

        while(!MWHeap.empty()){
            auto target = MWHeap.top();
            if(target->max >= start && target->min <= end){
                if((int)minHeap.size() < k || target->MW > minHeap.top().weight)
                    queryTopK_H(target, start, end, k, minHeap, MWHeap);
            }
            MWHeap.pop();
        }

        vector<Interval> result;
        while (!minHeap.empty()) {
            result.push_back(minHeap.top());
            minHeap.pop();
        }
        reverse(result.begin(), result.end());
        return result;
    }

private:
    vector<Interval> roots;
    WITNode* root = nullptr;    

    struct CompareMW {
        bool operator()(const WITNode* a, const WITNode* b) {
            return a->MW < b->MW;
        }
    };

    void build() {
        root = build(roots, 0, roots.size() - 1);
    }

    WITNode* build(const vector<Interval>& intervals, int start, int end) {
        if (start > end) {
            return nullptr;
        }
        int mid = start + (end - start) / 2;
        WITNode* node = new WITNode(intervals[mid]);
        node->left = build(intervals, start, mid - 1);
        node->right = build(intervals, mid + 1, end);
        if (node->left) {            
            node->min = min(node->min, node->left->min);
            node->max = max(node->max, node->left->max);
            if(node->MW > node->left->MW) {
                // node->MW = node->MW;
            } else if(node->MW == node->left->MW) {
                // node->MW = node->MW;=
            } else {
                node->MW = node->left->MW;      
            }
        }
        if (node->right) {           
            node->min = min(node->min, node->right->min);
            node->max = max(node->max, node->right->max);
            if(node->MW > node->right->MW) {
                // node->MW = node->MW;
            } else if(node->MW == node->right->MW) {
                // node->MW = node->MW;
            } else {
                node->MW = node->right->MW;
            }
        }
        return node;
    }
    
    void query(WITNode* node, int query_start, int query_end, vector<Interval>& result) {
        if (!node) {
            return;
        }
        if (node->interval.start <= query_end && node->interval.end >= query_start) {
            result.push_back(node->interval);
        }
        if (node->left && node->left->max >= query_start && node->left->min <= query_end) {
            query(node->left, query_start, query_end, result);
        }
        if (node->right && node->right->max >= query_start && node->right->min <= query_end) {
            query(node->right, query_start, query_end, result);
        }
    }

    void queryTopK(WITNode* node, int query_start, int query_end, int k, priority_queue<Interval, vector<Interval>, CompareWeight>& minHeap) {
        if (!node) {
            return;
        }
        if (node->interval.start <= query_end && node->interval.end >= query_start) {
            if ((int)minHeap.size() < k) {
                minHeap.push(node->interval);
            } else if (node->interval.weight > minHeap.top().weight) {
                minHeap.pop();
                minHeap.push(node->interval);
            }
        }
        if (node->left && node->left->max >= query_start && node->left->min <= query_end && ((int)minHeap.size() < k || node->left->MW > minHeap.top().weight)) {
            queryTopK(node->left, query_start, query_end, k, minHeap);
        }
        if (node->right && node->right->max >= query_start && node->right->min <= query_end && ((int)minHeap.size() < k || node->right->MW > minHeap.top().weight)) {
            queryTopK(node->right, query_start, query_end, k, minHeap);
        }
    }

    void queryTopK_H(WITNode* node, int query_start, int query_end, int k,
                        priority_queue<Interval, vector<Interval>, CompareWeight>& minHeap,                         
                        priority_queue<WITNode*, vector<WITNode*>, CompareMW>& MWHeap) {
        if (!node) {
            return;
        }
        if (node->interval.start <= query_end && node->interval.end >= query_start) {
            if ((int)minHeap.size() < k) {
                minHeap.push(node->interval);
            } else if (node->interval.weight > minHeap.top().weight) {
                minHeap.pop();
                minHeap.push(node->interval);
            }
        }
        if (node->left) {
            MWHeap.push(node->left);
        }
        if (node->right) {
            MWHeap.push(node->right);
        }
    }

    void deleteTree(WITNode* node) {
        if (!node) return;
        deleteTree(node->left);
        deleteTree(node->right);
        delete node;
    }
};

#endif