/*
 * 
 *  Copyright (c) 2008, 2011, 2012, 2013
 *  David Hernandez, Patrice Francois, Jacques Schrenzel
 * 
 *  This file is part of EDENA.
 *
 *  EDENA is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EDENA is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with EDENA.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef READSLAYOUTS_H
#define READSLAYOUTS_H

#include <cstddef>
#include <iostream>
#include <string>
#include <vector>
#include <deque>

#include "readsStorage.h"
using namespace std;

class ReadsStorage;
class Pairing;
class CountMatrix;

class ReadsLayout {
    
public:
    ReadsLayout();
    virtual ~ReadsLayout();

    void init(ReadsStorage *R, unsigned int n);
    void allocateOverHangingLinks();
    void cleanMemory();

    void save(ostream&);
    bool load(istream&, ReadsStorage*);

    void initLayout(size_t layout, bool dir, unsigned int nodeId);
    void setLayoutNodeId(size_t layout, unsigned int nodeId);

    string getDirectRead(size_t i) const;
    string getReverseRead(size_t i) const;
    void getPCharRead(char *dest, size_t i, bool dir) const;
    unsigned int getMultiplicity(size_t i);

    //list related accessor
    size_t getBegin(size_t index) const;
    size_t getEnd(size_t index) const;
    size_t reverseComplement(size_t index);
    size_t merge(size_t l1, size_t l2, bool direct, int size);
    size_t link(size_t l1, size_t l2, bool direct, int ovsize);

    bool isUniDirectional(size_t i) const;
    unsigned int getSequenceLength(size_t listIndex) const;
    unsigned int getNReads(size_t listIndex) const;
    
    inline unsigned int getNext(size_t i) const {return next[i];}
    inline unsigned int getPrevious(size_t i) const {return previous[i];}
    inline unsigned int getNextOverHanging(size_t i) const{return nextOverHanging[i];}
    inline unsigned int getPreviousOverHanging(size_t i) const {return previousOverHanging[i];}
    
    int getPosition(size_t i) const;
    bool getDirection(size_t i) const;
    unsigned int getNodeId(size_t i) const;
    
private:

    inline void setNext(size_t i, unsigned int id);
    inline void setPrevious(size_t i, unsigned int id);
    inline void setPosition(size_t i, int pos);
    inline void setDirection(size_t i, bool dir);
    inline void setPosition(size_t i, int pos, bool dir);
    void setNodeId(size_t i, unsigned int id);
    void unChain(size_t i);

public:

    unsigned int getLastIdentical(size_t i) const;
    void setLastIdentical(size_t i, unsigned int lastId);

    string getSequence(size_t listIndex);
    string getReverseSequence(size_t listIndex);
    void getVcoverage(size_t listIndex, bool direction, vector<unsigned int> &cov);

    void print(size_t index, ostream &out, bool dir, unsigned int start, unsigned int distance, Pairing *P);

    unsigned int sampleOH(size_t listIndex, bool dir, unsigned int maxD, unsigned int *distr1);

    double getMeanOverlapSize(size_t i);
    bool checkLayout(size_t index); //dev, destroy the lastIdentical tab!!
    void writeReadsAsFasta(size_t index, ostream &out);
    
    void flagReads(size_t index, char state);
    
    void buildOverHangingLinks(size_t index);
  
private:

    ReadsStorage *R;
    unsigned int tabSize; // > number of reads+1
    unsigned int *next;
    unsigned int *previous;
    
    //used for directly accessing next overhanging reads, and this bypassing identical aligned reads
    unsigned int *nextOverHanging;
    unsigned int *previousOverHanging;
    
    unsigned int *nodeId;
    int *position; //encode relative position and direction
    unsigned char *flags;
    
    unsigned int *lastIdentical; //last identical read (zero if first one)
};

#endif // READSLAYOUTS_H
