/******************************************************************************
 * Project:  hint
 * Purpose:  Indexing interval data
 * Author:   Panagiotis Bouros, pbour@github.io
 * Author:   Nikos Mamoulis
 ******************************************************************************
 * Copyright (c) 2020 - 2022
 *
 * All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 ******************************************************************************/

#ifndef _HINT_M_H_
#define _HINT_M_H_

#include "../def_global.h"
#include "../containers/relation.h"
#include "../containers/offsets.h"
#include "../containers/offsets_templates.cpp"
#include "../indices/hierarchicalindex.h"
#include <boost/dynamic_bitset.hpp>


// Base HINT^m, no optimizations activated
class HINT_M : public HierarchicalIndex
{
private:
    Relation **pOrgs, **pReps;
    RecordId **pOrgs_sizes;
    size_t   **pReps_sizes;
    
    // Construction
    inline void updateCounters(const Record &r);
    inline void updatePartitions(const Record &r);
    
public:
    // Construction
    HINT_M(const Relation &R, const unsigned int numBits, const unsigned int maxBits);
    void print(const char c);
    void getStats();
    ~HINT_M();
    
    // Querying
    vector<Record> executeTopDown_gOverlaps(RangeQuery Q, int k);
    vector<Record> executeBottomUp_gOverlaps(RangeQuery Q, int k);
};


// HINT^m with subs+sort optimization activated
class HINT_M_SubsSort : public HierarchicalIndex
{
private:
    Relation **pOrgsIn;
    Relation **pOrgsAft;
    Relation **pRepsIn;
    Relation **pRepsAft;
    RecordId **pOrgsIn_sizes, **pOrgsAft_sizes;
    size_t   **pRepsIn_sizes, **pRepsAft_sizes;
    
    // Construction
    inline void updateCounters(const Record &r);
    inline void updatePartitions(const Record &r);
    
public:
    // Construction
    HINT_M_SubsSort(const Relation &R, const unsigned int numBits, const unsigned int maxBits);
    void getStats();
    ~HINT_M_SubsSort();
    
    // Querying
    vector<Record> executeBottomUp_gOverlaps(RangeQuery Q, int k);
};


// Comparators
inline bool CompareTimestampPairsByStart(const pair<Timestamp, Timestamp> &lhs, const pair<Timestamp, Timestamp> &rhs)
{
    return (lhs.first < rhs.first);
}

inline bool CompareTimestampPairsByEnd(const pair<Timestamp, Timestamp> &lhs, const pair<Timestamp, Timestamp> &rhs)
{
    return (lhs.second < rhs.second);
}
#endif // _HINT_M_H_
