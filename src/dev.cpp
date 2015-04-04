#include "node.h"
#include "overlapGraph.h"

extern unsigned int global_count;

void OverlapsGraph::computeIrreductibeEdges()
{
    minOverlap = R->minOvSize;
    nNodes = R->getN_nrReads();
    Node::nNodes = nNodes;
    Node::setEdgeSortedFlag(true);
    allocateNodeTab(nNodes);
    for (unsigned int i = 1; i <= nNodes; i++)
    {
        nodesTab[i].layout = L->getLastIdentical(i);
    }
    g_count=0;
    nodeQueue.clear();
    unsigned int node;

    for (unsigned int i = 1; i <= nNodes; i++)
    {
        if (nodesTab[i].ovComputed()==true)
            continue; //

        nodesTab[i].computeOverlaps();
        nodesTab[i].anotherMarkTransitiveEdges();
        //  nodeQueue.push_front(i);

        while (!nodeQueue.empty())
        {
            node = nodeQueue.back();//*(nodeQueue.end() - 1);
            nodeQueue.pop_back();
            nodesTab[node].anotherMarkTransitiveEdges();
        }
    }
}

void Node::anotherMarkTransitiveEdges()
{

    //used by the "on the fly" reduction procedure
    unsigned int id,idB, size=0,sizeB;
    bool dir,dirB;
    unsigned int longestOH;
    unsigned int overhang;
    unsigned int nEdges;
    int absoluteIndexShift;
    int targetSize;
    unsigned int edgeIndex;

    stack.clear();
    
    //Linear expected algorithm adapted from:
    //E.W. Myers, The fragment assembly graph, Bioinformatics vol 21 Suppl.2 2005 pages ii79-ii85
    NodeMem nM, nMB;
    
    bool currentDir=false;
    
    do
    {
        currentDir=!currentDir;
        
        absoluteIndexShift=0;
        if (!currentDir)
            absoluteIndexShift=nOvRight;
       
        nEdges=getNEdges(currentDir);
        
        for (unsigned int i=0; i<nEdges; i++)
        {
            getNeighbor(currentDir,i,id,size,dir);
            
            if (N[id].isInplay(dir))
                N[id].setMultipleEdges(dir);
            else
                N[id].setInplay(dir);

            nM.id=id;
            nM.incomIndex=i+absoluteIndexShift;
            nM.dir=dir;
            nM.overHang=R->getReadsLength() - size;
            nM.isDelegated=false;
            stack.push_front(nM);
            
//            if (N[id].ovComputed()==false)
//            {
//                N[id].computeOverlaps();
//                //enqueue for transitive reduction
//                OverlapsGraph::nodeQueue.push_front(id);
//            }
        }

        longestOH = R->getReadsLength() - size; //longest overhanging

        while (stack.empty() == false)
        {
            nM=stack.back();//*(stack.end()-1);
            stack.pop_back();
            id=nM.id;
            dir=nM.dir;
            overhang=nM.overHang;
            
            if (edgeBitsFlag[nM.incomIndex] == true &&
                    nM.isDelegated==false) //edges already marked ->multiple edges
            {
                continue;
            }

            if (N[id].isEliminated(dir) == false || nM.isDelegated)
            {

                if (N[id].ovComputed() == false)
                {
                    N[id].computeOverlaps();
                    //enqueue for transitive reduction
                    OverlapsGraph::nodeQueue.push_front(id);
                }
                
                if (N[id].isReduced())
                { //delegate to neighbor

                    nEdges=N[id].getNEdges(dir);

                    for (unsigned int i = 0; i < nEdges; i++)
                    {
                        global_count++;

                        N[id].getNeighbor(dir, i, idB, sizeB, dirB);

                        if ((R->getReadsLength() - sizeB) + overhang > longestOH)
                            break;

                        if (N[idB].hasMultipleEdges(dirB))//means "traversed AND multiple edges
                        {
                            targetSize = sizeB - overhang;

                            if (currentDir)
                                edgeIndex = getEdgeIndex(true, idB, targetSize, dirB);
                            else
                                edgeIndex = getEdgeIndex(false, idB, targetSize, !dirB);

                            edgeBitsFlag[edgeIndex] = true;
                        }
                        else if (N[idB].isInplay(dirB))
                        {
                            N[idB].setEliminated(dirB);
                        }

                        nMB.id = idB;
                        nMB.dir = dirB;
                        nMB.incomIndex = i + absoluteIndexShift;
                        nMB.isDelegated = true;
                        nMB.overHang = overhang + R->getReadsLength() - sizeB;
                        stack.push_back(nMB);
                    }
                }
                else
                { //Myers algorithm 

                    nEdges=N[id].getNEdges(dir);

                    for (unsigned int i = 0; i < nEdges; i++)
                    {
                        global_count++;

                        N[id].getNeighbor(dir, i, idB, sizeB, dirB);
                        
                        if ((R->getReadsLength() - sizeB) + overhang > longestOH)
                            break;

                        if (N[idB].hasMultipleEdges(dirB))
                        {
                            //this point as been added from Myers's algorithm which do not take account for
                            //multiple edges between two nodes

                            targetSize = sizeB - overhang;
                            if (currentDir)
                                edgeIndex = getEdgeIndex(true, idB, targetSize, dirB);
                            else
                                edgeIndex = getEdgeIndex(false, idB, targetSize, !dirB);

                            edgeBitsFlag[edgeIndex] = true;
                        }
                        else if (N[idB].isInplay(dirB))
                        {
                            N[idB].setEliminated(dirB);
                        }
                    }
                } //end Myers's algorithm
            }
        } //end stack pop_back

        //init flags and mark edges
       
        nEdges=getNEdges(currentDir);
        
        for (unsigned int i=0; i<nEdges; i++)
        {
            getNeighbor(currentDir,i,id,size,dir);
            if(N[id].isEliminated(dir))
                edgeBitsFlag[i+absoluteIndexShift]=true;
        }
        for (unsigned int i=0; i<nEdges; i++)
        {
            getNeighbor(currentDir,i,id,size,dir);
            //N[id].flagsBit &= ~126;//do not change flag 1 and 128
            initTransitiveReductionFlags();
        }
        
    } while (currentDir == true);
    
    removeMarkedEdgeUnrec(255);
    
    setReduced();
    reallocMemory();
}