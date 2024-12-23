#ifndef INTERVAL_H
#define INTERVAL_H

#include <iostream>
#include <vector>

struct Interval {
    int id;
    int start;
    int end;
    int weight = 0;
    Interval(int id, int s, int e, int w) : id(id), start(s), end(e), weight(w) {}
    Interval() : id(0),  start(0), end(0), weight(0) {}

    bool operator==(const Interval& other) const {
        return start == other.start && end == other.end && weight == other.weight;
    }
};

std::ostream& operator<<(std::ostream& os, const Interval& interval) {
    os << "Interval(" << interval.start << ", " << interval.end << ", " << interval.weight << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<Interval>& intervals) {
    os << "[";
    for (size_t i = 0; i < intervals.size(); ++i) {
        os << intervals[i];
        if (i != intervals.size() - 1) {
            os << ", ";
        }
    }
    os << "]";
    return os;
}

struct CompareWeight {
    bool operator()(const Interval& a, const Interval& b) {
        return a.weight > b.weight;  // Lightweight intervals should have higher priority
    }
};

#endif