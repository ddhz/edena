/* 
 * File:   DevShell.cpp
 * Author: david
 * 
 * Created on November 26, 2012, 2:16 PM
 */

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

#include "DevShell.h"
#include "stat.h"
#include <limits>
using namespace std;

DevShell::DevShell() {
}

DevShell::DevShell(const DevShell& orig) {
}

DevShell::~DevShell() {
}

void DevShell::init(OverlapsGraph *g, Pairing *p, ReadsStorage *r, ReadsLayout *l) {
    G = g;
    P = p;
    R = r;
    L = l;
    checkDistances = true;
    minNPair = 5;
    minRatio = 0.95;
    if (!BKW.isAllocated())
        BKW.init(G);
}

void DevShell::prompt() {
    string str, buffer;

    bool finished = false;

    do {
        cout << "edenaDev> ";
        args.clear();
        getline(cin, buffer);
        iss.clear();
        iss.str(buffer);
        do {
            iss >> str;
            if (!iss.fail())
                args.push_back(str);
        } while (!iss.fail());

        if (args.size() == 0)
            continue;

        if (args[0] == "g") {
            drawGraph();
        } else if (args[0] == "t") {
            growSearchTree();
        } else if (args[0] == "e") {
            elongateSearchTree();
        } else if (args[0] == "t2s") {
            Tree2seq();
        } else if (args[0] == "p2s") {
            path2Seq();
        } else if (args[0] == "rp") {
            reversePath();
        } else if (args[0] == "l") {
            printReadsLayout();
        } else if (args[0] == "h" || args[0] == "?" || args[0] == "help") {
            printHelp(cerr);
        } else if (args[0] == "findRead") {
            locateRead();
        } else if (args[0] == "checkD") {
            setCheckDistances();
        } else if (args[0] == "p3") {
            p3();
        } else if (args[0] == "q" || args[0] == "Q" || args[0] == "quit") {
            cerr << "bye" << endl;
            finished = true;
        } else if (args[0] == "dev") {
            dev();
        } else if (args[0] == "cc") {
            sampleEdgeContext();
        } else {
            whatDoYouMean(cerr);
        }

    } while (finished == false);
}

void DevShell::printHelp(ostream &out) {
    cout << "DEV SHELL COMMANDS:\n"
            << "   g <PATH> <maxDepth>  draw graph around <PATH> (could be a single node)\n"
            << "                        (require graphviz). 'maxDepth' is optional.\n"
            << "                        The graph is written to the file \"OSG.ps\".\n"
            << "   t <PATH>             Initialize the search Tree with a path.\n"
            << "                        PATH can be prefixed with \'R\' to specify\n"
            << "                        the path reverse-complemented. The tree is\n"
            << "                        written to the file \"searchTree.ps\".\n"
            << "   t <no argument>      Grow the search tree by one more stage.\n"
            << "   e <treenode_ID>      Elongate the determined path using treenode_ID.\n"
            << "                        treenode_ID must be one of the candidate.\n"
            << "                        This command may updates the file \"searchTree.ps\".\n"
            << "   e <no argument>      Same as above, but using the candidate determined\n"
            << "                        by the extension algorithm (if even).\n"
            << "   t2s                  Extract the sequence corresponding to the backward\n"
            << "                        walker determined path.\n"
            << "   checkD yes/no        Set the 'check PE distance' flag:\n"
            << "                        yes: pairs must be connected within allowed range\n"
            << "                             to be matched.\n"
            << "                        no: pairs are matched as long as a path between the\n"
            << "                        the two reads exists.\n"
            << "   p3 <P1[P2]P3>        Get primers to PCR that flanks the sub-path P2\n"
            << "                        This command provides a sequence including the\n"
            << "                        '[...]' marker, as suited to be pasted in the\n"
            << "                        primer3 Web site (http://primer3.wi.mit.edu/)\n"
            << "                        See below for the syntax of <P1[P2]P3>\n"
            << "   rp <PATH>            Display the PATH reverse-complemented.\n"
            << "   p2s <PATH>           Display the sequence corresponding\n"
            << "                        to a path in the graph. The path is required to\n"
            << "                        be valid and unambiguous to be displayed.\n"
            << "   l <ID>               Display the reads layout in node <ID>.\n"
            << "                        !!can be long to complete for long nodes...!!\n"
            << "   l <ID> <maxD>        Display the layout up to a distance of maxD.\n"
            << "   l <ID> <begin> <end> Display the layout for coordinates begin..end.\n"
            << "   findRead <readID>    Display the Node ID and position of readID.\n"
            << "   q                    bye!\n\n"
            << "   <ID> can be specified as:\n"
            << "               1234        (node 1234)\n"
            << "              !1234        (same but reverse comp)\n"
            << "   <readID> is the read number, as displayed by the 'l' (layout) command\n\n"
            << "   <PATH> is a '-' separated list of nodeIDs as displayed in the '.lay'\n"
            << "                   file. A single node <ID> is also considered as a path\n"
            << "   <P1[P2]P3> is a path made up of three sub-paths P1, P2 and P3 with P2\n"
            << "              possibly empty.\n"
            << "              Examples of valid syntaxes:\n"
            << "              12-!23[3-!2]!4\n"
            << "              12-!23[]3-!2-!4\n"
            << "              12[!23]3-!2-!4\n"
            << "   To allow pasting line longer that 4096 char: type \"stty -icanon\" before\n"
            << "   launching Edena\n";
}

void DevShell::growSearchTree() {
    if (args.size() > 1) {
        G->nodeListToEPath(args.at(1), ePath);
        BKW.initWalker(ePath, minNPair, minRatio, checkDistances, 10000);
        BKW.stepElongate(0, chosen);
        BKW.dot("searchTree");
        if (chosen != 0)
            cerr << "Chosen candidate: " << chosen << endl;
        else
            cerr << "no determined elongation" << endl;

    } else if (args.size() == 1) {
        if (BKW.isInitialized() == false) {
            cerr << "you must first initialize the backward walker with a path" << endl;
            return;
        }
        BKW.stepElongate(0, chosen);
        BKW.dot("searchTree");

        if (chosen != 0)
            cerr << "Chosen candidate: " << chosen << endl;
        else
            cerr << "no determined elongation" << endl;
    }


};

void DevShell::elongateSearchTree() {

    unsigned int candidate;
    if (BKW.isInitialized() == false) {
        cerr << "You must first initialize the search tree with a path" << endl;
        return;
    }

    if (args.size() == 1) {
        if (chosen != 0) {
            cerr << "Chosen candidate: " << chosen << endl;
            BKW.acceptChoice(chosen);
            chosen = 0;
            BKW.dot("searchTree");
        } else {//try another tree stage
            BKW.stepElongate(0, chosen);
            if (chosen != 0) {
                cerr << "Chosen candidate: " << chosen << endl;
                BKW.acceptChoice(chosen);
                chosen = 0;
            } else {
                cerr << "no chosen candidate" << endl;
            }

            BKW.dot("searchTree");
        }
    } else if (args.size() > 1) {
        iss.clear();
        iss.str(args.at(1));
        iss >> candidate;

        if (BKW.acceptChoice(candidate) == 0) {
            cout << "tree_node: " << candidate << " is not a candidate" << endl;
        } else {
            cerr << "Chosen candidate: " << candidate << endl;
            BKW.dot("searchTree");
        }
    }

    args.clear(); //workaround to avoid tree initializationb
    args.push_back("fake");
    growSearchTree();
};

void DevShell::Tree2seq() {

    vector<unsigned int> epath;
    vector<unsigned int> path;
    vector<bool> dir;
    string s;
    ostringstream oss;
    ofstream out;

    if (BKW.isInitialized() == false) {
        cerr << "You must first initialize the search tree with a path" << endl;
        return;
    }

    BKW.getPath(epath);
    G->ePathToNodePath(epath, pathNode, pathDir);
    s = G->ePathToSeq(epath);

    for (size_t i = 1; i < pathNode.size(); i++) {
        if (!pathDir[i])
            oss << '!';
        oss << pathNode[i];
        if (i < pathNode.size() - 1)
            oss << '-';
    }

    if (s != "") {
        out.open("out.dev");
        cerr << "SearchTreeToSeq (l=" << pathNode[0] << ") " << oss.str();
        out << ">SearchTreeToSeq (l=" << pathNode[0] << ") " << oss.str();
        cerr << '\n';
        out << '\n';
        lineWrap(cerr, s, 70);
        lineWrap(out, s, 70);
        out.close();
    }

}

void DevShell::drawGraph() {

    if (args.size() < 2)
        return;

    vector<unsigned int> ePath;
    vector<unsigned int> nn;
    vector<bool> dd;
    G->nodeListToEPath(args.at(1), ePath);
    G->ePathToNodePath(ePath, nn, dd);

    if (ePath.empty()) {
        whatDoYouMean(cerr);
        return;
    }

    unsigned int gDepth = 8;

    if (args.size() == 3) {
        iss.clear();
        iss.str(args.at(2));
        iss >> gDepth;
    }

    ofstream out("OSG.dot");
    Node::initDot(out);

    for (size_t i = 1; i < nn.size(); i++) {
        G->nodesTab[nn.at(i)].setTargeted();
    }
    for (size_t i = 1; i < nn.size(); i++) {
        G->nodesTab[nn.at(i)].dotAround(out, gDepth);
    }
    Node::closeDot(out);

    out.close();
    out.clear();

    //executing dot command
    string command = "dot -Tps -o ";
    command += "OSG.ps OSG.dot";
    cout << "executing " << command << " ... " << flush;
    int unused __attribute__((unused)); //get rid of gcc warning
    unused = system(command.c_str());
    cout << "done" << endl;
    // G->nodesTab[nodeID].dotLocalGraph(gDepth, "OSG");
}

void DevShell::path2Seq() {

    if (args.size() < 2)
        return;

    G->nodeListToEPath(args.at(1), ePath);
    string str = G->ePathToSeq(ePath);
    ofstream out;
    if (str != "") {
        out.open("out.dev");
        cerr << "pathToSeq (l=" << str.size() << ") " << args.at(1);
        out << ">pathToSeq (l=" << str.size() << ") " << args.at(1);
        cerr << '\n';
        out << '\n';
        lineWrap(cerr, str, 70);
        lineWrap(out, str, 70);
        out.close();
    } else
        cout << "invalid path" << endl;

}

void DevShell::reversePath() {

    string pp = "R";
    vector<unsigned int> ePath;
    vector<unsigned int> nn;
    vector<bool> dd;
    pp += args.at(1);
    G->nodeListToEPath(pp, ePath);
    G->ePathToNodePath(ePath, nn, dd);
    for (unsigned int i = 1; i < nn.size(); i++) {
        if (!dd[i])
            cout << "!";
        cout << nn[i];
        if (i < nn.size() - 1)
            cout << "-";
    }
    cout << endl;

};

void DevShell::printReadsLayout() {

    parseNodeId(args.at(1));

    if (nodeID == 0) {
        whatDoYouMean(cerr);
        return;
    }
    unsigned int maxD = numeric_limits<unsigned int>::max();
    unsigned int start = 0;

    if (args.size() == 3) {
        iss.clear();
        iss.str(args.at(2));
        iss >> maxD;
        start = 0;
    } else if (args.size() == 4) {
        iss.clear();
        iss.str(args.at(2));
        iss >> start;
        iss.clear();
        iss.str(args.at(3));
        iss >> maxD;
    }

    L->print(G->nodesTab[nodeID].getLayout(), cerr, nodeDir, start, maxD, P);
    cout << endl;
};

void DevShell::locateRead() {


    unsigned int readId = 0;
    bool dir = true;

    if (args.size() < 2) {
        whatDoYouMean(cerr);
        return;
    }
    iss.clear();
    iss.str(args.at(1));
    iss >> readId; //read ID
    if (iss.fail()) {
        whatDoYouMean(cerr);
        return;
    }

    if (readId < 1 || readId >= R->getNumberOfReads() + 1) {
        cout << "<readID> out of range (1.." << R->getNumberOfReads() << ")\n";
        return;
    }
    if (L->getDirection(readId) == false)
        dir = false;
    cerr << L->getNodeId(readId);
    if (!dir)
        cerr << " < ";
    else
        cerr << " > ";
    cerr << "pos: " << L->getPosition(readId) << endl;
}

void DevShell::setCheckDistances() {
    if (args.size() != 2) {
        cerr << "possible paramters are 'yes' or 'no'" << endl;
        return;
    }

    if (args[1] == "yes") {
        checkDistances = true;
        cerr << "check PE distance flag: YES" << endl;
    } else if (args[1] == "no") {
        checkDistances = false;
        cerr << "check PE distance flag: NO" << endl;
    } else
        cerr << "possible paramters are 'yes' or 'no'" << endl;

    BKW.setCheckDistances(checkDistances);
}

void DevShell::p3() {
    
    //check annotated path string (search for '[' and ']'
    //and set p1 and p2 member according to the position of '[' and ']'
    //in the path
    if (!set_p1p2_ePath()) {
        cerr << "invalid syntax" << endl;
        return;
    }

    //   annotatedPath = "n1[n2]n3
    //   this will produce annotated DNA string, excluding the
    //   overlap region, which may be bad if neighbor node is 
    //   not unique.
    
    //   .................... n1
    //                 .............. n2   
    //                         ..................... n2
    //                       [   ]
    //                       <---> = may be a single nucleotide!             
    unsigned int id, size;
    bool dir, currentDir;
    unsigned int currentNode;
    bool checkDir;
    string s, tmp;

    currentNode = ePath.at(1);
    if (ePath.at(2) < G->nodesTab[currentNode].getNRightOverlap())
        currentDir = true;
    else
        currentDir = false;

    s = G->nodesTab[currentNode].getSequence(currentDir);

    for (unsigned int i = 2; i < ePath.size(); i++) {
        //check path validity
        if (ePath.at(i) < G->nodesTab[currentNode].getNRightOverlap())
            checkDir = true;
        else
            checkDir = false;

        if (checkDir != currentDir) {
            cout << "problem" << endl;
            return;
        }

        G->nodesTab[currentNode].getOverlap(ePath.at(i), id, size, dir);
        
        if (p1 == i - 1) {
            //insert '[' before the overlaps
            s.insert(s.length() - size, "[");
        }
        if (p2 == i - 1) {
            s += ']';
        }

        if (!dir)
            currentDir = !currentDir;
        tmp = G->nodesTab[id].getSequence(currentDir);
        tmp = tmp.substr(size);
        s += tmp;
        currentNode = id;
    }

    ofstream out("out.dev");
    lineWrap(cerr, s, 70);
    out << "edenaDev> p3 " << args.at(1) << '\n';
    lineWrap(out, s, 70);
    out.close();
}

bool DevShell::set_p1p2_ePath() {
    p1 = p2 = 0;

    if (args.size() < 2)
        return false;

    nodeList = args.at(1);
    unsigned int count = 0, mem_i = 0;
    bool open = false;
    bool closed = false;

    for (size_t i = 0; i < nodeList.size(); i++) {
        if (nodeList[i] == '-')
            count++;
        else if (nodeList[i] == '[') {
            if (open == true || closed == true || i == 0)
                return false;

            open = true;
            //convert '[' to '-' for the path to be a valid nodepath
            nodeList[i] = '-';
            count++;
            p1 = count;
            mem_i = i;
        } else if (nodeList[i] == ']') {
            if (open == false || closed == true || i == nodeList.size() - 1)
                return false;

            closed = true;


            if (i > mem_i + 1) // x-x-x[x-x-x]x-x-x
            {
                count++;
                nodeList[i] = '-';
            } else //x-x-x[]x-x-x
            {
                nodeList.erase(nodeList.begin() + i);
            }

            p2 = count;
        }
    }

    if (open == false || closed == false)
        return false;

    G->nodeListToEPath(nodeList, ePath);

    if (ePath.size() > 0)
        return true;
    return false;

}

void DevShell::dev() {
    
    unsigned int readId;
    iss.clear();
    iss.str(args.at(1));
    iss >> readId;
    cout << "multiplicity: " << L->getMultiplicity(readId) << endl;
    cout << "paired ID: " << P->getPairing(readId) << endl;
//    double minCoverage;
//    string s;
//    iss.clear();
//    iss.str(args.at(1));
//    iss >> minCoverage; //read ID
//    if (iss.fail()) {
//        whatDoYouMean(cerr);
//        return;
//    }
//
//    ofstream out("out.dev");
//
//    for (unsigned int i = 1; i <= G->getNNodes(); i++) {
//        if (G->nodesTab[i].getCoverage() >= minCoverage) {
//            out << ">node_" << i
//                    << " length=" << G->nodesTab[i].getSequenceLength()
//                    << " cov=" << G->nodesTab[i].getCoverage()
//                    << '\n';
//            s = G->nodesTab[i].getSequence(true);
//            lineWrap(out, s, 70);
//
//        }
//    }
}

void DevShell::parseNodeId(string buffer) {
    if (buffer[0] == '!')//reverse
    {
        nodeDir = false;
        buffer = buffer.substr(1);
    } else
        nodeDir = true;

    iss.clear();
    iss.str(buffer);
    iss >> nodeID;
    if (iss.fail())
        nodeID = 0;

    if (nodeID < 1 || nodeID > G->getNNodes()) {
        cout << "<nodeID> out of range (1.." << G->getNNodes() << ")\n";
        nodeID = 0;
    }
}

void DevShell::sampleEdgeContext() {

    unsigned int maxD = 200;
    G->nodeListToEPath(args.at(1), ePath);
    if (ePath.size() != 3) {
        cout << "Path must be made of two nodes (=one edge)\n";
        return;
    }


    unsigned int node1 = ePath[1];
    unsigned int node2 = G->nodesTab[node1].getOverlapId(ePath[2]);

    bool node1Dir, node2Dir;

    if (ePath[0] == 0)
        node1Dir = false;
    else
        node1Dir = true;

    node2Dir = G->nodesTab[node1].getOverlapDirection(ePath[2]);



    //reliableCutoff: used to flag edges that are, or have at least one reliable brother

    //sample requirements. 
    unsigned int minSampleSize = 25;
    unsigned int minSampledLength = 0;

    unsigned int distr1[512];
    unsigned int distr2[512];

    unsigned int nOhSampled1 = 0;
    unsigned int nOhSampled2 = 0;
    unsigned int oh;
    double prob,prob2, mean;
    unsigned int rl = R->getReadsLength();

    oh = G->nodesTab[node1].getOverlapSize(ePath[2]);
    oh = G->getReadLength() - oh;
    cout << "overhanging length OH: " << oh << endl;

    memset(distr1, 0, 512 * sizeof (unsigned int));

    nOhSampled1 = G->nodesTab[node1].sampleOverHangs(!node1Dir, maxD, distr1);
    memset(distr2, 0, 512 * sizeof (unsigned int));
    nOhSampled2 = G->nodesTab[node2].sampleOverHangs(node2Dir, maxD, distr2);

    unsigned int nOhSampled = nOhSampled1 + nOhSampled2;

    unsigned int sampleSize = 0, sum = 0;

    for (unsigned int i = 0; i < rl; i++) {
        distr1[i] += distr2[i];
        sampleSize += distr1[i];
        sum += i * distr1[i];
    }
    ofstream out("out.dev");
    out << "oh\tnumber\trawNumber\n";
    out << "0\t" << distr1[0] << "\t" << distr1[511] + distr2[511] << '\n';
    for (unsigned int i = 1; i < rl; i++) {
        out << i << '\t' << distr1[i] << '\t' << distr1[i] << '\n';
    }
    out.close();


    cout << "nReads: " << nOhSampled << '\n';
    cout << "sampleSize: " << sampleSize << '\n';
    cout << "sum: " << sum << '\n';
    //note: sum is equivalent to sampled length




    if (sampleSize < minSampleSize || sum < minSampledLength) {
        prob = prob2 = 2.0; //canceled
    } else {

        mean = (double) sum / sampleSize;
      
        prob2 = cdfShiftedGeom(1.0 / (mean + 1), oh);
        //  prob = cdfOH(1.0 / (mean + 1), oh, 1.0);
    }
   // cout << std::setprecision(10) << std::scientific <<  "prob: " << prob << endl;
     cout << std::setprecision(10) << std::scientific << "prob: " << prob2 << endl;


}