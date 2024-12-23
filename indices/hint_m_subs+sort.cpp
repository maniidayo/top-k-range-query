/******************************************************************************
 * Project:  hint
 * Purpose:  Indexing interval data
 * Author:   Panagiotis Bouros, pbour@github.io
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

#include "hint_m.h"
#include <queue>

struct CompareWeight {
    bool operator()(const Record& a, const Record& b) {
        return a.weight > b.weight;  // Lightweight intervals should have higher priority
    }
};

inline void HINT_M_SubsSort::updateCounters(const Record &r)
{
    int level = 0;
    Timestamp a = r.start >> (this->maxBits-this->numBits);
    Timestamp b = r.end   >> (this->maxBits-this->numBits);
    Timestamp prevb;
    int firstfound = 0, lastfound = 0;
    
    
    while (level < this->height && a <= b)
    {
        if (a%2)
        { //last bit of a is 1
            if (firstfound)
            {
                if ((a == b) && (!lastfound))
                {
                    this->pRepsIn_sizes[level][a]++;
                    lastfound = 1;
                }
                else
                    this->pRepsAft_sizes[level][a]++;
            }
            else
            {
                if ((a == b) && (!lastfound))
                    this->pOrgsIn_sizes[level][a]++;
                else
                    this->pOrgsAft_sizes[level][a]++;
                firstfound = 1;
            }
            a++;
        }
        if (!(b%2))
        { //last bit of b is 0
            prevb = b;
            b--;
            if ((!firstfound) && b < a)
            {
                if (!lastfound)
                    this->pOrgsIn_sizes[level][prevb]++;
                else
                    this->pOrgsAft_sizes[level][prevb]++;
            }
            else
            {
                if (!lastfound)
                {
                    this->pRepsIn_sizes[level][prevb]++;
                    lastfound = 1;
                }
                else
                {
                    this->pRepsAft_sizes[level][prevb]++;
                }
            }
        }
        a >>= 1; // a = a div 2
        b >>= 1; // b = b div 2
        level++;
    }
}


inline void HINT_M_SubsSort::updatePartitions(const Record &r)
{
    int level = 0;
    Timestamp a = r.start >> (this->maxBits-this->numBits);
    Timestamp b = r.end   >> (this->maxBits-this->numBits);
    Timestamp prevb;
    int firstfound = 0, lastfound = 0;
    
    
    while (level < this->height && a <= b)
    {
        if (a%2)
        { //last bit of a is 1
            if (firstfound)
            {
                if ((a == b) && (!lastfound))
                {
                    this->pRepsIn[level][a][this->pRepsIn_sizes[level][a]] = r;
                    this->pRepsIn_sizes[level][a]++;
                    lastfound = 1;
                }
                else
                {
                    this->pRepsAft[level][a][this->pRepsAft_sizes[level][a]] = r;
                    this->pRepsAft_sizes[level][a]++;
                }
            }
            else
            {
                if ((a == b) && (!lastfound))
                {
                    this->pOrgsIn[level][a][this->pOrgsIn_sizes[level][a]] = r;
                    this->pOrgsIn_sizes[level][a]++;
                }
                else
                {
                    this->pOrgsAft[level][a][this->pOrgsAft_sizes[level][a]] = r;
                    this->pOrgsAft_sizes[level][a]++;
                }
                firstfound = 1;
            }
            a++;
        }
        if (!(b%2))
        { //last bit of b is 0
            prevb = b;
            b--;
            if ((!firstfound) && b < a)
            {
                if (!lastfound)
                {
                    this->pOrgsIn[level][prevb][this->pOrgsIn_sizes[level][prevb]] = r;
                    this->pOrgsIn_sizes[level][prevb]++;
                }
                else
                {
                    this->pOrgsAft[level][prevb][this->pOrgsAft_sizes[level][prevb]] = r;
                    this->pOrgsAft_sizes[level][prevb]++;
                }
            }
            else
            {
                if (!lastfound)
                {
                    this->pRepsIn[level][prevb][this->pRepsIn_sizes[level][prevb]] = r;
                    this->pRepsIn_sizes[level][prevb]++;
                    lastfound = 1;
                }
                else
                {
                    this->pRepsAft[level][prevb][this->pRepsAft_sizes[level][prevb]] = r;
                    this->pRepsAft_sizes[level][prevb]++;
                }
            }
        }
        a >>= 1; // a = a div 2
        b >>= 1; // b = b div 2
        level++;
    }
}


HINT_M_SubsSort::HINT_M_SubsSort(const Relation &R, const unsigned int numBits, const unsigned int maxBits) : HierarchicalIndex(R, numBits, maxBits)
{
    // Initialize statistics
    this->numOriginalsIn = this->numOriginalsAft = this->numReplicasIn = this->numReplicasAft = 0;

    // Step 1: one pass to count the contents inside each partition.
    this->pOrgsIn_sizes  = (RecordId **)malloc(this->height*sizeof(RecordId *));
    this->pOrgsAft_sizes = (RecordId **)malloc(this->height*sizeof(RecordId *));
    this->pRepsIn_sizes  = (size_t **)malloc(this->height*sizeof(size_t *));
    this->pRepsAft_sizes = (size_t **)malloc(this->height*sizeof(size_t *));

    for (auto l = 0; l < this->height; l++)
    {
        auto cnt = (int)(pow(2, this->numBits-l));
        
        //calloc allocates memory and sets each counter to 0
        this->pOrgsIn_sizes[l]  = (RecordId *)calloc(cnt, sizeof(RecordId));
        this->pOrgsAft_sizes[l] = (RecordId *)calloc(cnt, sizeof(RecordId));
        this->pRepsIn_sizes[l]  = (size_t *)calloc(cnt, sizeof(size_t));
        this->pRepsAft_sizes[l] = (size_t *)calloc(cnt, sizeof(size_t));
    }
    
    if (!R.empty()) {
        for (const Record &r : R)
            this->updateCounters(r);    
    }

    // Step 2: allocate necessary memory.
    this->pOrgsIn  = new Relation*[this->height];
    this->pOrgsAft = new Relation*[this->height];
    this->pRepsIn  = new Relation*[this->height];
    this->pRepsAft = new Relation*[this->height];
    for (auto l = 0; l < this->height; l++)
    {
        auto cnt = (int)(pow(2, this->numBits-l));
        
        this->pOrgsIn[l]  = new Relation[cnt];
        this->pOrgsAft[l] = new Relation[cnt];
        this->pRepsIn[l]  = new Relation[cnt];
        this->pRepsAft[l] = new Relation[cnt];
        
        for (auto j = 0; j < cnt; j++)
        {
            this->pOrgsIn[l][j].resize(this->pOrgsIn_sizes[l][j]);
            this->pOrgsAft[l][j].resize(this->pOrgsAft_sizes[l][j]);
            this->pRepsIn[l][j].resize(this->pRepsIn_sizes[l][j]);
            this->pRepsAft[l][j].resize(this->pRepsAft_sizes[l][j]);
        }
    }
    for (auto l = 0; l < this->height; l++)
    {
        auto cnt = (int)(pow(2, this->numBits-l));
        
        memset(this->pOrgsIn_sizes[l], 0, cnt*sizeof(RecordId));
        memset(this->pOrgsAft_sizes[l], 0, cnt*sizeof(RecordId));
        memset(this->pRepsIn_sizes[l], 0, cnt*sizeof(size_t));
        memset(this->pRepsAft_sizes[l], 0, cnt*sizeof(size_t));
    }
    

    // Step 3: fill partitions.
    for (const Record &r : R)
        this->updatePartitions(r);
    

    // Step 4: sort partition contents.
    for (auto l = 0; l < this->height; l++)
    {
        auto cnt = (int)(pow(2, this->numBits-l));
        for (auto j = 0; j < cnt; j++)
        {
            this->pOrgsIn[l][j].sortByStart();
            this->pOrgsAft[l][j].sortByStart();
            this->pRepsIn[l][j].sortByEnd();
        }
    }
    
    // Free auxiliary memory.
    for (auto l = 0; l < this->height; l++)
    {
        free(pOrgsIn_sizes[l]);
        free(pOrgsAft_sizes[l]);
        free(pRepsIn_sizes[l]);
        free(pRepsAft_sizes[l]);
    }
    free(pOrgsIn_sizes);
    free(pOrgsAft_sizes);
    free(pRepsIn_sizes);
    free(pRepsAft_sizes);
}


HINT_M_SubsSort::~HINT_M_SubsSort()
{
    for (auto l = 0; l < this->height; l++)
    {
        delete[] this->pOrgsIn[l];
        delete[] this->pOrgsAft[l];
        delete[] this->pRepsIn[l];
        delete[] this->pRepsAft[l];
    }
    delete[] this->pOrgsIn;
    delete[] this->pOrgsAft;
    delete[] this->pRepsIn;
    delete[] this->pRepsAft;
}


void HINT_M_SubsSort::getStats()
{
    size_t sum = 0;
    for (auto l = 0; l < this->height; l++)
    {
        auto cnt = pow(2, this->numBits-l);
        
        this->numPartitions += cnt;
        for (int j = 0; j < cnt; j++)
        {
            this->numOriginalsIn  += this->pOrgsIn[l][j].size();
            this->numOriginalsAft += this->pOrgsAft[l][j].size();
            this->numReplicasIn   += this->pRepsIn[l][j].size();
            this->numReplicasAft  += this->pRepsAft[l][j].size();
            if ((this->pOrgsIn[l][j].empty()) && (this->pOrgsAft[l][j].empty()) && (this->pRepsIn[l][j].empty()) && (this->pRepsAft[l][j].empty()))
                this->numEmptyPartitions++;
        }
    }
    
    this->avgPartitionSize = (float)(this->numIndexedRecords+this->numReplicasIn+this->numReplicasAft)/(this->numPartitions-numEmptyPartitions);
}


// Generalized predicates, ACM SIGMOD'22 gOverlaps
vector<Record> HINT_M_SubsSort::executeBottomUp_gOverlaps(RangeQuery Q, int k)
{
    priority_queue<Record, vector<Record>, CompareWeight> minHeap;

    Relation::iterator iterBegin;
    RelationIterator iter, iterEnd;
    Timestamp a = Q.start >> (this->maxBits-this->numBits); // prefix
    Timestamp b = Q.end   >> (this->maxBits-this->numBits); // prefix
    Record qdummyE(0, Q.end+1, Q.end+1);
    Record qdummyS(0, Q.start, Q.start);
    bool foundzero = false;
    bool foundone = false;
    
    for (auto l = 0; l < this->numBits; l++)
    {
        if (foundone && foundzero)
        {
            // Partition totally covers lowest-level partition range that includes query range
            // all contents are guaranteed to be results
            
            // Handle the partition that contains a: consider both originals and replicas
            iterBegin = this->pRepsIn[l][a].begin();
            iterEnd = this->pRepsIn[l][a].end();

            for (iter = iterBegin; iter != iterEnd; iter++)
            {
                auto interval = Record(iter->id, iter->start, iter->end, iter->weight);
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                }
            }
            iterBegin =this->pRepsAft[l][a].begin();
            iterEnd = this->pRepsAft[l][a].end();
            for (iter = iterBegin; iter != iterEnd; iter++)
            {
                
                auto interval = Record(iter->id, iter->start, iter->end, iter->weight);
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                }
            }
            
            // Handle rest: consider only originals
            for (auto j = a; j <= b; j++)
            {
                iterBegin = this->pOrgsIn[l][j].begin();
                iterEnd = this->pOrgsIn[l][j].end();
                for (iter = iterBegin; iter != iterEnd; iter++)
                {
                    
                auto interval = Record(iter->id, iter->start, iter->end, iter->weight);
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                }
                }
                iterBegin = this->pOrgsAft[l][j].begin();
                iterEnd = this->pOrgsAft[l][j].end();
                for (iter = iterBegin; iter != iterEnd; iter++)
                {
                    
                auto interval = Record(iter->id, iter->start, iter->end, iter->weight);
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                }
                }
            }
        }
        else
        {
            // Comparisons needed
            
            // Handle the partition that contains a: consider both originals and replicas, comparisons needed
            if (a == b)
            {
                // Special case when query overlaps only one partition, Lemma 3
                if (!foundzero && !foundone)
                {
                    iterBegin = this->pOrgsIn[l][a].begin();
                    iterEnd = lower_bound(iterBegin, this->pOrgsIn[l][a].end(), qdummyE);
                    for (iter = iterBegin; iter != iterEnd; iter++)
                    {
                        if (Q.start <= iter->end)
                        {
                            
                auto interval = Record(iter->id, iter->start, iter->end, iter->weight);
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                }
                        }
                    }
                    iterBegin = this->pOrgsAft[l][a].begin();
                    iterEnd = lower_bound(iterBegin, this->pOrgsAft[l][a].end(), qdummyE);
                    for (iter = iterBegin; iter != iterEnd; iter++)
                    {
                        
                auto interval = Record(iter->id, iter->start, iter->end, iter->weight);
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                }
                    }
                }
                else if (foundzero)
                {
                    iterBegin = this->pOrgsIn[l][a].begin();
                    iterEnd = lower_bound(iterBegin, this->pOrgsIn[l][a].end(), qdummyE);
                    for (iter = iterBegin; iter != iterEnd; iter++)
                    {
                        
                auto interval = Record(iter->id, iter->start, iter->end, iter->weight);
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                }
                    }
                    iterBegin = this->pOrgsAft[l][a].begin();
                    iterEnd = lower_bound(iterBegin, this->pOrgsAft[l][a].end(), qdummyE);
                    for (iter = iterBegin; iter != iterEnd; iter++)
                    {
                        
                auto interval = Record(iter->id, iter->start, iter->end, iter->weight);
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                }
                    }
                }
                else if (foundone)
                {
                    iterBegin = this->pOrgsIn[l][a].begin();
                    iterEnd = this->pOrgsIn[l][a].end();
                    for (iter = iterBegin; iter != iterEnd; iter++)
                    {
                        if (Q.start <= iter->end)
                        {
                            
                auto interval = Record(iter->id, iter->start, iter->end, iter->weight);
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                }
                        }
                    }
                    iterBegin = this->pOrgsAft[l][a].begin();
                    iterEnd = this->pOrgsAft[l][a].end();
                    for (iter = iterBegin; iter != iterEnd; iter++)
                    {
                        
                auto interval = Record(iter->id, iter->start, iter->end, iter->weight);
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                }
                    }
                }
            }
            else
            {
                // Lemma 1
                if (!foundzero)
                {
                    iterBegin = this->pOrgsIn[l][a].begin();
                    iterEnd = this->pOrgsIn[l][a].end();
                    for (iter = iterBegin; iter != iterEnd; iter++)
                    {
                        if (Q.start <= iter->end)
                        {
                            
                auto interval = Record(iter->id, iter->start, iter->end, iter->weight);
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                }
                        }
                    }
                }
                else
                {
                    iterBegin = this->pOrgsIn[l][a].begin();
                    iterEnd = this->pOrgsIn[l][a].end();
                    for (iter = iterBegin; iter != iterEnd; iter++)
                    {
                        
                auto interval = Record(iter->id, iter->start, iter->end, iter->weight);
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                }
                    }
                }
                iterBegin = this->pOrgsAft[l][a].begin();
                iterEnd = this->pOrgsAft[l][a].end();
                for (iter = iterBegin; iter != iterEnd; iter++)
                {
                    
                auto interval = Record(iter->id, iter->start, iter->end, iter->weight);
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                }
                }
            }
            
            // Lemma 1, 3
            if (!foundzero)
            {
                iterEnd = this->pRepsIn[l][a].end();
                iterBegin = lower_bound(this->pRepsIn[l][a].begin(), this->pRepsIn[l][a].end(), qdummyS, CompareRecordsByEnd);
                for (iter = iterBegin; iter != iterEnd; iter++)
                {
                    
                auto interval = Record(iter->id, iter->start, iter->end, iter->weight);
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                }
                }
            }
            else
            {
                iterBegin = this->pRepsIn[l][a].begin();
                iterEnd = this->pRepsIn[l][a].end();
                for (iter = iterBegin; iter != iterEnd; iter++)
                {
                    
                auto interval = Record(iter->id, iter->start, iter->end, iter->weight);
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                }
                }
            }
            iterBegin = this->pRepsAft[l][a].begin();
            iterEnd = this->pRepsAft[l][a].end();
            for (iter = iterBegin; iter != iterEnd; iter++)
            {
                
                auto interval = Record(iter->id, iter->start, iter->end, iter->weight);
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                }
            }
            
            if (a < b)
            {
                if (!foundone)
                {
                    // Handle the rest before the partition that contains b: consider only originals, no comparisons needed
                    for (auto j = a+1; j < b; j++)
                    {
                        iterBegin = this->pOrgsIn[l][j].begin();
                        iterEnd = this->pOrgsIn[l][j].end();
                        for (iter = iterBegin; iter != iterEnd; iter++)
                        {
                            
                auto interval = Record(iter->id, iter->start, iter->end, iter->weight);
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                } iter->id;
                        }
                        iterBegin = this->pOrgsAft[l][j].begin();
                        iterEnd = this->pOrgsAft[l][j].end();
                        for (iter = iterBegin; iter != iterEnd; iter++)
                        {
                            
                auto interval = Record(iter->id, iter->start, iter->end, iter->weight);
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                }
                        }
                    }
                    
                    // Handle the partition that contains b: consider only originals, comparisons needed
                    iterBegin = this->pOrgsIn[l][b].begin();
                    iterEnd = lower_bound(iterBegin, this->pOrgsIn[l][b].end(), qdummyE);
                    for (iter = iterBegin; iter != iterEnd; iter++)
                    {
                        
                auto interval = Record(iter->id, iter->start, iter->end, iter->weight);
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                }
                    }
                    iterBegin = this->pOrgsAft[l][b].begin();
                    iterEnd = lower_bound(iterBegin, this->pOrgsAft[l][b].end(), qdummyE);
                    for (iter = iterBegin; iter != iterEnd; iter++)
                    {
                        
                auto interval = Record(iter->id, iter->start, iter->end, iter->weight);
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                }
                    }
                }
                else
                {
                    for (auto j = a+1; j <= b; j++)
                    {
                        iterBegin = this->pOrgsIn[l][j].begin();
                        iterEnd = this->pOrgsIn[l][j].end();
                        for (iter = iterBegin; iter != iterEnd; iter++)
                        {
                            
                auto interval = Record(iter->id, iter->start, iter->end, iter->weight);
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                }
                        }
                        iterBegin = this->pOrgsAft[l][j].begin();
                        iterEnd = this->pOrgsAft[l][j].end();
                        for (iter = iterBegin; iter != iterEnd; iter++)
                        {
                            
                auto interval = Record(iter->id, iter->start, iter->end, iter->weight);
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                }
                        }
                    }
                }
            }
            
            if ((!foundone) && (b%2)) //last bit of b is 1
                foundone = 1;
            if ((!foundzero) && (!(a%2))) //last bit of a is 0
                foundzero = 1;
        }
        a >>= 1; // a = a div 2
        b >>= 1; // b = b div 2
    }
    
    // Handle root.
    if (foundone && foundzero)
    {
        // All contents are guaranteed to be results
        iterBegin = this->pOrgsIn[this->numBits][0].begin();
        iterEnd = this->pOrgsIn[this->numBits][0].end();
        for (iter = iterBegin; iter != iterEnd; iter++)
        {            
            auto interval = Record(iter->id, iter->start, iter->end, iter->weight);
            if ((int)minHeap.size() < k) {
                minHeap.push(interval);
            } else if (interval.weight > minHeap.top().weight) {                
                minHeap.pop();
                minHeap.push(interval);
            }
        }
    }
    else
    {
        // Comparisons needed
        iterBegin = this->pOrgsIn[this->numBits][0].begin();
        iterEnd = lower_bound(iterBegin, this->pOrgsIn[this->numBits][0].end(), qdummyE);
        for (iter = iterBegin; iter != iterEnd; iter++)
        {
            if (Q.start <= iter->end)
            {
                
                auto interval = Record(iter->id, iter->start, iter->end, iter->weight);
                if ((int)minHeap.size() < k) {
                    minHeap.push(interval);
                } else if (interval.weight > minHeap.top().weight) {                
                    minHeap.pop();
                    minHeap.push(interval);
                }
            }
        }
    }

    vector<Record> result;

    while (!minHeap.empty()) {
        result.emplace_back(minHeap.top());
        minHeap.pop();
    }

    return result;
}