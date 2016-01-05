/**
trie_cuckoo_noise_07_21.cpp
create by: Yunhong
create time: 07/21/2015
*/


#include "Trie_cuckoo_noise_07_21.h"

ofstream outfileR;

float maxHaoOver = 0.0;
int indexHaoOver = 0.0;
float minHaoOver = 1.0;
int indexMinHaoOver = 0.0;

int main(int argc, char * argv[])
{
    cout << "\n";
    cout << "PRIME_OPENMP\n";
    cout << "  C++/OpenMP version\n";

    cout << "\n";
    cout << "  Number of processors available = " << omp_get_num_procs ( ) << "\n";
    cout << "  Number of threads =              " << omp_get_max_threads ( ) << "\n";

    // ---------------------------
    // Global variables
    CUCKOO_BLACK_SIZE = strtof(argv[6],NULL); // black table size
    FLOW_EST_SIZE = 0;       // flow estimator size
    float overTarget = 0.01;    // target overselection rate

    // -----------------------------
    /* initialize random seed: */
    srand (time(NULL));

    // time v.r.
    struct timeval gen_start, gen_end; /* gettimeofday stuff */
    struct timezone tzp;
    long timeInv = 0;

    // ---------------------------------------------------
    /*define mask*/
    vector<size_t> mask;
    size_t maskIP;
    mask.clear();
    cout<<"* Compute prefixes main ... ..."<<endl<<endl;;
    for(int i = 8; i <= 32; i++)
    {
        maskIP = (size_t(pow(2,i))-1)<<(32-i);
        mask.push_back(maskIP);
    }

    std::string inFileName = argv[1];
    intss& vuniquePrefix;
    intss& vuniqueAggPrefix;

    vuniquePrefix.clear();
    vuniqueAggPrefix.clear();

    loadKeys2Filter(inFileName, mask, vuniquePrefix, vuniqueAggPrefix);



    // ------------------------------------------------
    // File name for caida trace
    const char * fp[] = {"mapped-caida-1",
                         "mapped-caida-6",
                         "mapped-caida-11",
                         "mapped-caida-16",
                         "mapped-caida-21",
                         "mapped-caida-26",
                         "mapped-caida-31",
                         "mapped-caida-36",
                         "mapped-caida-41",
                         "mapped-caida-46",
                         "mapped-caida-51",
                         "mapped-caida-56"
                        };

    // -------------------------------------------------
    // Open out file, and write result into it
    std::ofstream outfile0;
    string outfileName = ("outfile0_simple_"+ string(argv[2]) + '_' +
                          string(argv[4]))+ "_tstNum_"+ string(argv[5])+"_b"+string(argv[6])+".csv";
    outfile0.open(outfileName.c_str());

    // wirte resource allocation
    //std::ofstream outfileR;
    outfileName = ("resouce_assign_"+ string(argv[2]) + '_' +
                          string(argv[4]))+ "_tstNum_"+ string(argv[5])+"_b"+string(argv[6])+".csv";
    outfileR.open(outfileName.c_str());

    // Open out file, and write result into it
    std::ofstream keyCountfile;
    outfileName = ("keyCountfile_simple_"+ string(argv[2]) + '_' +
                          string(argv[4]))+ "_tstNum_"+ string(argv[5])+"_b"+string(argv[6])+".csv";
    keyCountfile.open(outfileName.c_str());

    // Open out file, and write result into it
    std::ofstream keyCountfileEst;
    outfileName = ("keyCountfileEst_simple_"+ string(argv[2]) + '_' +
                          string(argv[4]))+ "_tstNum_"+ string(argv[5])+"_b"+string(argv[6])+".csv";
    keyCountfileEst.open(outfileName.c_str());

    // Open out file, and write result into it
    std::ofstream keyCountfileDiff;
    outfileName = ("keyCountfileDiff_simple_"+ string(argv[2]) + '_' +
                          string(argv[4]))+ "_tstNum_"+ string(argv[5])+"_b"+string(argv[6])+".csv";
    keyCountfileDiff.open(outfileName.c_str());

        // Open out file, and write result into it
    //std::ofstream keyCountfileAll;
    //outfileName = ("keyCountfileAll_simple_"+ string(argv[2]) + '_' +\
                          string(argv[4]))+ "_tstNum_"+ string(argv[5])+"_b"+string(argv[6])+".csv";
   // keyCountfileAll.open(outfileName.c_str());



    strings keycs;
    vector<int> keyActioncs;

    strings keycsM;
    vector<int> keyActioncsM;
    vector<size_t> keyCountcsM;  // est
    vector<size_t> keyCountcs0M; // accurate
    vector<size_t> keyCountDiffsM;
    cuckooFilter.returnKey(keycsM,keyActioncsM,mL0, keyCountcs,keyCountcs0,keyCountDiffs,keyCountcsM,
                                       keyCountcs0M,keyCountDiffsM);

    // create RLearn for recvs
    RLearn* rLearn[switchSize][actionSize];
    for(int si = 0; si < switchSize; si++)
    {
        for(int i = 0; i < actionSize; i++)
        {
            rLearn[si][i] = new RLearn();
            initRLearn(rLearn[si][i]);

        }
    }



    // -----------------------------------------------
    float falsePos = 0;
    float falsePos0 = 0;
    float haoFalsePos = 0;
    float haoFalsePos0 = 0;
    float haoFalsePosTotal;
    float haoFalsePos0Total;
    float overAggr = 0.0f;

    vector<string> overBigKeys;
    vector<size_t> overBigKeyNos;
    vector<string> blackKeys;
    vector<size_t> blackKeyNos;
    vector<int> blackActions;
    vector<float> haoOvers;
    haoOvers.assign(actionSize,0);

    vector<float> haoOversInv;
    haoOversInv.assign(actionSize,0);

    size_t countNum = 0;
    size_t countNum0 = 0;
    double countIP = 0;
    double countIP0 = 0;
    size_t countBlack = 0;

    double countIPTotal = 0.0f;
    double countIP0Total = 0.0f;
    double keySum = 0.0f;
    double pktSum = 0.0f;
    double keySumTotal = 0.0f;
    double pktSumTotal = 0.0f;
    double countIPTotalOff = 0.0f;
    double countIP0TotalOff = 0.0f;
    double keySumTotalOff = 0.0f;
    double pktSumTotalOff = 0.0f;
    double aggrSum= 0.0f;


    floats keySums;
    floats countIPs;
    keySums.assign(actionSize,0);
    countIPs.assign(actionSize,0);

    floats keySumsInv;
    floats countIPsInv;
    keySumsInv.assign(actionSize,0);
    countIPsInv.assign(actionSize,0);

    uint64_t line = 0;
    size_ts slotNums;
    slotNums.assign(actionSize,0);

    // ---------------------------
    // Define a trie
    Trie *trie;            // define tree
    trie = new Trie();            // define tree

    // -------------------------
    //std::ifstream infile;

    // ----------------------------
    // load packets from file
    for (int fi = 0; fi < 6; fi++)
    {

        // --------------------------
        // Open trace file
        string pathname = fp[fi];
        std::string inFileName = argv[3];
        inFileName += pathname;
        infile.open(inFileName.c_str());
        cout<<inFileName.c_str()<<endl;
        cout<<inFileName<<endl;
        if(!infile)
            std::cout << "* TestIP File Error " << std::endl;

        // ------------------------------
        countIPTotalOff += countIP;
        countIP0TotalOff += countIP0;
        keySumTotalOff += keySum;
        pktSumTotalOff += pktSum;

        //line = 0;
        countBlack = 0;
        countIP = 0.0f;
        countIP0 = 0.0f;
        keySum = 0.0f;
        pktSum = 0.0f;
        //aggrSum = 0.0f;
        //keySums.clear();
        //countIPs.clear();
        //keySums.assign(actionSize,0);
        //countIPs.assign(actionSize,0);

        keySumsInv.clear();
        countIPsInv.clear();
        keySumsInv.assign(actionSize,0);
        countIPsInv.assign(actionSize,0);


        int nthreads, tid;

        size_t  ei;
        bool isEndFlag = 0;
        size_t updateInvDis = 10000; // interval for display
        size_t readNum = strtof(argv[5],NULL);    // interval for feedback

        while(!isEndFlag )
        {
            // init overselection rate each cycle
            if(line < 100000)
            {
                keySumsInv.clear();
                countIPsInv.clear();
                keySumsInv.assign(actionSize,0);
                countIPsInv.assign(actionSize,0);

                for(int ai = 0; ai < actionSize; ai++)
                rLearn[ai]->clearList();
            }

            // -------------------------------
            // read file
            strings flows;
            size_ts flowNos;
            readFile0(infile, flows, flowNos, readNum, isEndFlag);
            size_t flowNum = flows.size();

            # pragma omp parallel for \
            shared ( infile,isEndFlag, updateInvDis, line, actionSize, mask, keySum, pktSum, aggrSum, keySums, countIPs, countIP, countNum, countIP0, countNum0, \
                     countBlack, overBigKeys,overBigKeyNos,  trie, bigNonCounts, timeCounts, uniquePrefix,uniqueAggPrefix, cuckooFilter, \
                     cuckooFilterInit0,cuckooBlackKeyTable,cuckooTableKey,cuckooAggrKeyTable,cuckooFilterFlowEst) \
            private ( ei )

            for(ei = 0; ei < flowNum; ei++)
            {
                // read 1 flow
                bool readFlag1 = 0;
                size_t flowNoInt = 1;
                uint32_t flowInt;

                #pragma omp critical
                {
                    /*if(!readFlag1)
                    {
                        readFlag1 = (infile>>flowInt>>flowSizeInt);
                    }*/
                    flowInt = parseIPV4string(flows[ei].c_str());
                    flowNoInt = flowNos[ei];
                    readFlag1 = 1;

                }

                if(readFlag1)
                {
                    // return position
                    long hBuckhit = -1;
                    int slot_ihit = -1;
                    long hBuck = -1;
                    int slot_i = -1;
                    // ---------------------------
                    // ---------------------------
                    // lookup  key
                    {
                        int prefix;
                        int flag_look,flag_look0, flag_lookkey;
                        bool flag_lookblack = 0;
                        uint32_t subIP;
                        string flowstr;
                        string flowIpv4 = parsedec2IPV4(flowInt);

                        uint32_t ip = flowInt; //convert IP to int
                        string keyType_cur = "0";
                        string flow_action_str;

                        int keySumLocal = 0;
                        int aggrSumLocal = 0;
                        int countIPLocal = 0;
                        int countNumLocal = 0;
                        int countIP0Local = 0;
                        int countNum0Local = 0;
                        vector<int> keySumsLocal;
                        keySumsLocal.clear();
                        keySumsLocal.assign(actionSize,0);
                        vector<int> countIPsLocal;
                        countIPsLocal.clear();
                        countIPsLocal.assign(actionSize,0);

                        //pkt_sum += flowNoInt;

                        // -------------------------------
                        // look up blackkey
                        int iflowaction;
                        if(CUCKOO_BLACK_SIZE > 0)
                        {
                            int mi = 24;
                            //for(int mi = 0; mi <= 32-8; mi++)
                            {

                                subIP = ip & mask[mi]; // mask
                                flowstr = parsedec2IPV4(ip);//s.str();
                                prefix = mi+8;

                                flag_lookblack = cuckooBlackKeyTable.LookUpKeyAction(flowstr,prefix,iflowaction);
                                if (flag_lookblack)
                                {
                                    countBlack += flowNoInt;
                                    //break;
                                }
                            }


                        }

                        // ------------------------------
                        flag_look = 0;
                        if (flag_lookblack == 0) // not a blackkey
                        {
                            // look up key
                            for(int mi = uniquePrefix.size()-1; mi >=0; mi--)
                            {
                                subIP = ip & mask[uniquePrefix[mi]-8];
                                flowstr = parsedec2IPV4(subIP);
                                prefix = uniquePrefix[mi];

                                flag_lookkey = cuckooTableKey.LookUpKeyActionCount(flowstr,
                                                                                   prefix,iflowaction,flowNoInt);
                                if (flag_lookkey)
                                {
                                    //cout<<hBuck<<" "<<slot_i<<endl;
                                    //keyCountcs0[hBuck][slot_i] += flowNoInt;
                                    //keyCountDiffs[hBuck][slot_i] = keyCountcs[hBuck][slot_i] - keyCountcs0[hBuck][slot_i];
                                    keyType_cur = "1";
                                    flow_action_str = num2str(iflowaction);
                                    //#pragma omp atomic
                                    keySumLocal += (flowNoInt);
                                    for(int ai = 0; ai <actionSize; ai++)
                                    {
                                        if(iflowaction == ai)
                                        {
                                            //#pragma omp atomic
                                            keySumsLocal[ai] += (flowNoInt);
                                        }

                                    }
                                    break;

                                }
                            }



                            // -------------------------------
                            // lookup key from cuckoo filter
                            vector<int> iactions;
                            bool isAKey = 0;

                            for(int mi = 0; mi < uniqueAggPrefix.size(); mi++)
                            {
                                bool overbigFlag = 0;
                                subIP = ip & mask[uniqueAggPrefix[mi]-8];
                                flowstr = parsedec2IPV4(subIP)+"/"+num2str(uniqueAggPrefix[mi]);
                                iactions.clear();
                                hBuckhit = -1;
                                slot_ihit = -1;
                                flag_look = cuckooFilter.LookUpKeyActionsCount(flowstr,iactions,flowNoInt, mL0,
                                                                                keyCountcs,  keyCountcs0, keyCountDiffs,hBuckhit,slot_ihit);
                                if(hBuckhit!= -1 && slot_ihit!= -1)
                                {
                                    hBuck = hBuckhit;
                                    slot_i = slot_ihit;
                                }
                                if (flag_look != 0)
                                {
                                    isAKey = 1;
                                    countNumLocal += flag_look;
                                    countIPLocal += flag_look*(flowNoInt);

                                    // each action lookups
                                    for(int ai = 0; ai <actionSize; ai++)
                                    {
                                        for(int li = 0; li < iactions.size(); li++)
                                        {
                                            if(iactions[li] == ai)
                                            {
                                                countIPsLocal[ai] += (flowNoInt);
                                            }

                                        }
                                    }

                                    // overselection count
                                    if(!flag_lookkey && !overbigFlag)
                                    {
                                        overbigFlag = 1;
                                        #pragma omp critical
                                        {
                                            trie->addWordCountNum(DecToBin(flowInt),32, isnonkey, iactions[0], countIPLocal);
                                        }

                                    }
                                }
                            }

                                                        // -------------------------------
                            // look up key
                            for(int mi = uniquePrefix.size()-1; mi >=0; mi--)
                            {
                                subIP = ip & mask[uniquePrefix[mi]-8];
                                flowstr = parsedec2IPV4(subIP);
                                prefix = uniquePrefix[mi];

                                flag_lookkey = cuckooTableKey.LookUpKeyAction(flowstr,
                                                                                   prefix,iflowaction);
                                if (flag_lookkey)
                                {
                                    //cout<<hBuck<<" "<<slot_i<<endl;
                                    keyCountcs0[hBuck][slot_i] += flowNoInt;
                                    keyCountDiffs[hBuck][slot_i] = keyCountcs[hBuck][slot_i] - keyCountcs0[hBuck][slot_i];
                                    //keyType_cur = "1";
                                    //flow_action_str = num2str(iflowaction);
                                    //#pragma omp atomic
                                    //keySumLocal += (flowNoInt);
                                    /*for(int ai = 0; ai <actionSize; ai++)
                                    {
                                        if(iflowaction == ai)
                                        {
                                            //#pragma omp atomic
                                            keySumsLocal[ai] += (flowNoInt);
                                        }

                                    }*/
                                    break;

                                }
                            }

                            // -------------------------------
                            // look up aggregate key
                            bool isAggregatekey = 0;
                            //if(cuckooAggrKeyTable.mm > 1)
                            {
                                for(int mi = uniqueAggPrefix.size()-1; mi >=0; mi--)
                                {
                                    subIP = ip & mask[uniqueAggPrefix[mi]-8];
                                    flowstr = parsedec2IPV4(subIP);
                                    prefix = uniqueAggPrefix[mi];
                                    isAggregatekey = cuckooAggrKeyTable.LookUpKey(flowstr,prefix);
                                    if (isAggregatekey && !flag_lookkey)
                                    {
                                        //#pragma omp atomic
                                        aggrSumLocal += (flowNoInt);
                                        break;

                                    }
                                }
                            }

                            // ---------------------------------------------------
                            // add to flow est filter
                            if(isAKey == 0)  // not inside cuckooFilter
                            {
                                // add to filter
                                //if(FLOW_EST_SIZE>0)
                               // addFlowEst(flowIpv4, bigNonCounts, timeCounts, overbigKey );

                            }

                        }

                        // -------------------------------
                        // lookup key from cuckoo filterInit0
                        vector<int> iactions;
                        for(int mi = 0; mi < uniquePrefix.size(); mi++)
                        {
                            subIP = ip & mask[uniquePrefix[mi]-8];
                            flowstr = parsedec2IPV4(subIP)+"/"+num2str(uniquePrefix[mi]);
                            iactions.clear();
                            flag_look0 =cuckooFilterInit0.LookUpKeyActions(flowstr,iactions);
                            if (flag_look0)
                            {
                                //#pragma omp atomic
                                countNum0Local += flag_look0;
                                //#pragma omp atomic
                                countIP0Local += flag_look0*(flowNoInt);

                            }

                        }
                        // lookup key end
                        // ---------------------------------------------
                        // ---------------------------------------------

                        // ---------------------------------------------
                        // update global variable

                        {
                            line ++;
                            #pragma omp atomic
                            pktSum += flowNoInt;
                            #pragma omp atomic
                            keySum += keySumLocal;
                            #pragma omp atomic
                            aggrSum += aggrSumLocal;
                            #pragma omp atomic
                            countIP += countIPLocal;
                            #pragma omp atomic
                            countNum += countNumLocal;
                            #pragma omp atomic
                            countIP0 += countIP0Local;
                            #pragma omp atomic
                            countNum0 += countNum0Local;
                            for(int ai = 0; ai <actionSize; ai++)
                            {
                                #pragma omp atomic
                                keySums[ai] += keySumsLocal[ai];
                                #pragma omp atomic
                                countIPs[ai] += countIPsLocal[ai];

                                #pragma omp atomic
                                keySumsInv[ai] += keySumsLocal[ai];
                                #pragma omp atomic
                                countIPsInv[ai] += countIPsLocal[ai];
                            }
                        }

                    }
                }
                else
                {
                    isEndFlag = 1; // end read flag
                }
                //tid = omp_get_thread_num();

                // ---------------------------------
                // update rates and write to file

                #pragma omp critical
                {
                    if(size_t(line)%200 == 0)
                    {
                        // update total overselects value
                        countIPTotal = countIPTotalOff + countIP;
                        countIP0Total = countIP0TotalOff + countIP0;
                        keySumTotal = keySumTotalOff + keySum;
                        pktSumTotal = pktSumTotalOff + pktSum;

                        // --------------------------------------------
                        // overselection rate
                        // ----------------------------------------------
                        if ((pktSum - keySum) != 0)
                        {
                            falsePos = float(countIP-keySum)/float(pktSum-keySum);
                            falsePos0 = float(countIP0-keySum)/float(pktSum-keySum);
                        }


                        if(keySum != 0)
                        {
                            haoFalsePos = float(countIP-keySum)/float(keySum);
                            haoFalsePos0 = float(countIP0-keySum)/float(keySum);
                        }

                        for(int ai = 0; ai < actionSize; ai++)
                        {
                            if(keySums[ai] != 0)
                            {
                                haoOvers[ai] = float(countIPs[ai]-keySums[ai])/float(keySums[ai]);
                            }
                            if(keySumsInv[ai] != 0)
                            {
                                haoOversInv[ai] = float(countIPsInv[ai]-keySumsInv[ai])/float(keySumsInv[ai]);
                            }
                        }

                        if(keySumTotal != 0)
                        {
                            haoFalsePosTotal = (countIPTotal-keySumTotal)/(keySumTotal);
                            haoFalsePos0Total = (countIP0Total-keySumTotal)/(keySumTotal);
                        }

                        overAggr = aggrSum/keySumTotal;
                    }

                    if(size_t(line)%updateInvDis == 0 || line%flowNum == 0)
                    {

                        // ------------------------------------------------
                        // Display and write to file
                        outfile0<<"line,"<<line<<",total,"<<haoFalsePosTotal<<",total0,"<<haoFalsePos0Total<<",false,"<<
                                falsePos<<",over,"<<haoFalsePos<<",false0,"<<falsePos0<<",over0,"<<
                                haoFalsePos0<<",overaggr,"<<overAggr<<",overcuckoo,"<<(haoFalsePosTotal-overAggr)<<",";
                        for(int ai = 0; ai < actionSize; ai++)
                             outfile0<<"over_ai,"<<haoOversInv[ai]<<",";
                        for(int ai = 0; ai < actionSize; ai++)
                             outfile0<<"over_ai,"<<haoOvers[ai]<<",";

                        outfile0<<"countIP,"<<countIPTotal<<",countIP0,"<<countIP0Total<<",keysum,"<<keySumTotal<<
                                ",pktSum,"<<pktSumTotal<<",aggrSum,"<<aggrSum<<",blackkey_num,"<<countBlack<<",feedback,"<<blackBackSize<<",feedsumportion,"
                                <<feedSumPortion<<",finger0,"<<finger0<<",finger,"<<finger<<",time,"<<timeInv<<endl;


                        cout<<endl<<"line "<<line<<" total "<<haoFalsePosTotal<<" total0 "<<haoFalsePos0Total<<" false "<<
                                falsePos<<" over "<<haoFalsePos<<" false0 "<<falsePos0<<" over0 "<<
                                haoFalsePos0<<" overaggr "<<overAggr<<" overcuckoo "<<(haoFalsePosTotal-overAggr)<<endl;

                        for(int ai = 0; ai < actionSize; ai++)
                             cout<<"over_ai "<<haoOversInv[ai]<<" ";
                        for(int ai = 0; ai < actionSize; ai++)
                             cout<<"over_ai,"<<haoOvers[ai]<<",";

                        cout<<"countIP "<<countIPTotal<<" countIP0 "<<countIP0Total<<" keysum "<<keySumTotal<<
                                " pktSum "<<pktSumTotal<<" aggrSum "<<aggrSum<<" blackkey_num "<<countBlack<<" feedback "<<blackBackSize<<" feedsumportion "
                                <<feedSumPortion<<" finger0 "<<finger0<<" finger "<<finger<<" time "<<timeInv<<endl<<endl;
                    }
                }

            }

            #pragma omp barrier

            // ---------------------------------------------
            // feed back portionFeedBack% overselections and overselection rates for actions
            //if(size_t(line)%2000000 == 0 )
            {

                // ---------------------------
                // the overselects for feedback
                vector<char> word;
                trie->getLeaf(trie->root,word,blackKeys,blackKeyNos, blackActions);
                trie->deleteChild(trie->root);
                delete trie;
                trie = new Trie();            // define tree

                //trie->deleteChild(trie->root);

                int init = 0;
                float overTotalSum = accumulate(blackKeyNos.begin(), blackKeyNos.end(), init);

                // sort balckkeys
                keySort3(blackKeys,blackKeyNos, blackActions);

                float portionFeedBack = strtof(argv[4],NULL);
                blackBackSize = portionFeedBack*blackKeys.size();
                /*if(blackKeys.size()<100)
                {
                    blackBackSize = blackKeys.size();
                }
                else if(blackKeys.size()<200 && portionFeedBack<=0.5)
                {
                    blackBackSize = 0.5*blackKeys.size();
                }
                else if(blackKeys.size()<1000 && portionFeedBack <=0.1)
                {
                    blackBackSize = 0.1*blackKeys.size();
                }*/
                blackKeys.erase(blackKeys.begin()+blackBackSize, blackKeys.end());
                blackKeyNos.erase(blackKeyNos.begin()+blackBackSize, blackKeyNos.end());
                blackActions.erase(blackActions.begin()+blackBackSize, blackActions.end());

                float feedSum = accumulate(blackKeyNos.begin(), blackKeyNos.end(), init);
                feedSumPortion = feedSum/overTotalSum;

                // -----------------------------
                // Get instant reward and Q learn
                //cout<<"* Q learn!"<<endl;
                //if(pktSumTotal > 9.59646e+08)
                //if(line>60000)
                {
                    // cout<<"* compare: "<<countIPsInv[0]<<" "<<countIPs[0]<<endl;
                    maxHaoOver = haoOvers[0];
                    indexHaoOver = 0;
                    minHaoOver = haoOvers[0];
                    indexMinHaoOver = 0;

                    for(int ri = 0; ri < actionSize; ri++)
                    {
                        //cout<<rLearn[ri]<<endl;
                        rLearn[ri]->update(slotNums[ri],haoOvers[ri]);
                        //cout<<"* qleran!"<<endl;

                        // find maximum haoover
                        if(haoOvers[ri] > maxHaoOver)
                        {
                            maxHaoOver = haoOvers[ri];
                            indexHaoOver = ri;
                        }

                        // find maximum haoover
                        if(haoOvers[ri] < minHaoOver)
                        {
                            minHaoOver = haoOvers[ri];
                            indexMinHaoOver = ri;
                        }

                        // compute ebuse0
                        if(line<1000000)
                        {
                            EPSILON0 = 1.0;
                        }
                        else
                        {
                            EPSILON0 = 0.20;
                        }
                        rLearn[ri]->qLearn();
                        printQList(rLearn[ri]);
                    }

                    // ------------------------------------
                    // feedback
                    // -------------------------
                    //time before calling function
                    gettimeofday(&gen_start, &tzp);
                    // ---------------------------------
                    // call function
                    //feedbackBlackkey(blackKeys);
                    cout<<"* feedback bks!"<<endl;
                    feedbackBlackkeyRL(blackKeys, blackActions, rLearn, actionSize,slotNums, line);
                    cout<<"* feedback bks end!"<<endl;
                    /*if(10 < blackKeys.size())
                    {
                        for(size_t i = 0; i < 10; i++)
                        {
                             cout<<"* blackkey: "<<blackKeys[i]<<" num: "<<blackKeyNos[i]<<" actions: "<<blackActions[i]<<endl;
                        }
                    }*/


                    // ------------------------------
                    // time after calling function
                    gettimeofday(&gen_end, &tzp);
                    // time interval
                    timeInv = print_elapsed("Aggr: ", &gen_start,&gen_end, 1);
                }

                // write to file
                strings().swap(keycsM);
                ints().swap(keyActioncsM);
                size_ts().swap(keyCountcsM);
                size_ts().swap(keyCountcs0M);
                size_ts().swap(keyCountDiffsM);

                cuckooFilter.returnKey(keycsM,keyActioncsM,mL0, keyCountcs,keyCountcs0,keyCountDiffs,keyCountcsM,
                                       keyCountcs0M,keyCountDiffsM);
                /*for(int i = 0; i < 10; i++)
                cout<<keyCountcsM[i]<<" "<<keyCountcs0M[i]<<" "<<keyCountDiffsM[i]<<"     ";
                cout<<endl;*/
                /*for(size_t i = 0; i < keycsM.size(); i++)
                {
                    keyCountfileEst<<keyCountcsM[i]<<" ";
                    //keyCountfileAll<<keyCountcsM[i]<<" ";
                }
                //keyCountfileAll<<endl;
                for(size_t i = 0; i < keycsM.size(); i++)
                {
                    keyCountfile<<keyCountcs0M[i]<<" ";
                    //keyCountfileAll<<keyCountcs0M[i]<<" ";
                }
                //keyCountfileAll<<endl;
                for(size_t i = 0; i < keycsM.size(); i++)
                {
                    keyCountfileDiff<<keyCountDiffsM[i]<<" ";
                    //keyCountfileAll<<keyCountDiffsM[i]<<" ";
                    / *if(keyCountDiffsM[i]!=0)
                    {
                        cout<<"Diff: "<<keyCountDiffsM[i]<<endl;
                    }* /
                }
                //keyCountfileAll<<endl;
                keyCountfile<<endl;
                keyCountfileEst<<endl;
                keyCountfileDiff<<endl;*/

                // ------------------------------------
                // print haoOvers
                /*for(int i = 0; i< actionSize; i++)
                    cout<<"Action: "<<i <<" haoOver: "<<" "<<haoOvers[i]<<" ";
                cout<<endl;*/

                // ------------------------------------
                // aggregation
                //isInit = 1;
                /*if(line >= 20000000 && line%20000000 == 0 && haoFalsePosTotal > overTarget)
                {

                //trie->deleteChild(trie->root);
                //delete trie;
                //trie = new Trie();            // define tree

                    aggregation(keys,keyprefixlengths,keyActions, \
                            mask, actionSize, storage, isInit, finger,uniqueAggPrefix, blackKeys, blackKeyNos, haoOvers, overTarget);
                }*/


                strings().swap(blackKeys);
                vector<size_t>().swap(blackKeyNos);
                vector<int>().swap(blackActions);
            }

            strings().swap(flows);
            vector<size_t>().swap(flowNos);

        }
        //ifstream().swap(infile);
        infile.clear();
        infile.close();
    } //for fi

    //-------------------------------------------------------------
    /*clear the data structure*/
    mask.clear();

    cuckooFilter.ClearTable();
    cuckooBlackKeyTable.ClearTable();
    cuckooTableKey.ClearTable();
    cuckooAggrKeyTable.ClearTable();

    (keyCountfile.close());
    keyCountfileDiff.close();
    keyCountfileEst.close();
    (outfile0.close());
}

void feedbackBlackkey(vector<string>& overBigKeys)
{
    // Add blackkeys to cuckooTable
    GLOBAL_BIGNONKEYNUM = 100.0;
    vector<string> blackKeys;
    blackKeys.clear();

    //cout<<"* overBigKeys No: "<<overBigKeys.size()<<endl;

    for(int i = 0; i < overBigKeys.size(); i++)
    {
        //if(overBigKeyNos[i]>GLOBAL_BIGNONKEYNUM)
        {
            blackKeys.push_back(overBigKeys[i]);
        }
    }

    // Add blackkey to cuckooTable
    //cout<<"* Add balckkey to cuckooTable!"<<endl;
    float loadFactor = 0.9f;
    int slotNo = 4;
    size_t blackKeySize = CUCKOO_BLACK_SIZE;
    size_t bucketSize = int(blackKeySize/(loadFactor*slotNo))+1;
    int fingerprintNew = 12;
    long MaxKickoutNum = 1000;
    cuckooBlackKeyTable.ClearTable();
    cuckooBlackKeyTable.CuckooTableInit(bucketSize,fingerprintNew,slotNo, \
                                        MaxKickoutNum);

    for(size_t i = 0; i < blackKeys.size(); i++)
    {
        if(i<CUCKOO_BLACK_SIZE)
        {
            bool addFlag = cuckooBlackKeyTable.AddKeyPrefix(blackKeys[i],32, 4);
            if(!addFlag)
            {
                cout<<"1 add fail..."<<endl;
            }
        }
    }
    size_t newBlackSize = blackKeys.size();

    // Import previous blackkeys
    ifstream blackKeyFileIn;
    blackKeyFileIn.open(BLACKFILENAME.c_str());
    string blackKeyStr;
    int cuckooBlackSize = CUCKOO_BLACK_SIZE;
    int prefix = 32;
    while(blackKeyFileIn>>blackKeyStr>>prefix && blackKeys.size()<cuckooBlackSize)
    {
        bool blackFlag = 0;
        int iflowaction;
        blackFlag = cuckooBlackKeyTable.LookUpKeyAction(blackKeyStr,prefix,iflowaction);

        if(!blackFlag)
        blackKeys.push_back(blackKeyStr);
    }
    blackKeyFileIn.clear();
    blackKeyFileIn.close();

    // -----------------------------------
    // Write blackkey to file
    //cout<<"* Write blackkey to file!"<<endl;
    ofstream blackKeyFileOut;
    blackKeyFileOut.open(BLACKFILENAME.c_str());
    for(int i = 0; i < blackKeys.size(); i++)
        blackKeyFileOut<<blackKeys[i]<<" "<<32<<endl;
    blackKeyFileOut.clear();
    blackKeyFileOut.close();

    // --------------------------------------
    //cout<<"balckkey  size: "<<blackKeys.size()<<endl;
    for(size_t i = newBlackSize; i < blackKeys.size(); i++)
    {
        bool addFlag = cuckooBlackKeyTable.AddKeyPrefix(blackKeys[i],32, 4);
        if(!addFlag)
        {
            cout<<"2 add fail..."<<endl;
        }
    }



}


void feedbackBlackkey1(strings& overBigKeys)
{

    float loadFactor = 0.9f;
    int slotNo = 4;
    int cuckooBlackSize = CUCKOO_BLACK_SIZE;
    int blackKeySize = cuckooBlackSize;
    int bucketSize = int(blackKeySize/(loadFactor*slotNo))+1;
    int fingerprintNew = 12;
    int MaxKickoutNum = 1000;
    cuckooBlackKeyTable.ClearTable();
    cuckooBlackKeyTable.CuckooTableInit(bucketSize,fingerprintNew,slotNo,
                                        MaxKickoutNum);
    int prefix = 32;
    // write to file
    ofstream blackkeyFile;
    blackkeyFile.open("blackkeyFile");

    sort(overBigKeys.begin(),overBigKeys.end());
    vector<string>::iterator it;
    it = unique (overBigKeys.begin(), overBigKeys.end());
    overBigKeys.resize( distance(overBigKeys.begin(),it) );
    cout<<"* overBigKeys No: "<<overBigKeys.size()<<endl;
    int mini = 0;
    if(overBigKeys.size()>cuckooBlackSize)
    {
        mini = overBigKeys.size()-cuckooBlackSize;
        overBigKeys.erase(overBigKeys.begin(),overBigKeys.begin()+mini);
    }
    for(int i = overBigKeys.size()-1; i > 0 ; i--)
    {

        cuckooBlackKeyTable.AddKeyPrefix(overBigKeys[i],prefix, 4);
        blackkeyFile<<overBigKeys[i]<<endl;
    }
    blackkeyFile.close();

}

bool readFile0(ifstream& infile, vector<string> &flow, vector<size_t> &flow_cnt, size_t readNum, bool& isEndFlag)
{
    Trie *bTrie = new Trie();
    uint32_t flowInt;
    size_t flowNo;
    size_t in_num = 0;

    while((infile >> flowInt >>flowNo) && (in_num < readNum))
    {

        bTrie->addWordCount(DecToBin(flowInt),32, isnonkey, flowNo);
        in_num ++;
        //if(in_num%1000000 == 0)
            //cout<<"loading file ...."<<in_num<<"  ";
    }
    vector<char> word;
    vector<int> flowActions;
    bTrie->getLeaf(bTrie->root,word,flow,flow_cnt,flowActions);


    if(in_num < readNum)
    {
        isEndFlag = 1;
    }
    //cout<<"* Flow size: "<<flow.size()<<endl;

    delete bTrie;
    return true;
}


void initRLearn(RLearn* rLearn)
{
    // build a RL for rev1

    // create a qlist and initialize it
    float minState = 20;
    float maxState = 820;   // ovs
    float minAction = 20;
    float maxAction = 820; // #slots

    float state = minState;
    float action = minAction;
    float ovs = 0;
    float ovsTarget = 0.01;

    size_t numS = 8;
    size_t numA = 8;

    //float initState = 500;
    //float initAction = 500;


    //rLearn->rLearnInit();

    // initialize the table
    cout<<"* init qlist..."<<endl;
    rLearn->initQtable(minState,maxState,minAction,maxAction, numS,  numA,  ovsTarget);


    // initialized state
    cout<<"* init rlearn state... "<<endl;
    rLearn->update(state, ovs);

    printQList(rLearn);

    // get suggected action

    // get immedite reward


    // Q leran
    //rLearn.qLearn();
}

void updateBlacklist(vector<string>& overBigKeys, vector<int>& overActions, RLearn* rLearn, int actionSeq
, vector<string>& blackkeyPres, vector<int>& blackActionPres, ofstream& blackKeyFileOut, size_t& slotNum)
{
    // for each recv
    // get the suggested action
    //slotNum = rLearn->selectActionSuggest();
    cout<<"* suggested action++++++++++++++++++++++++++++++++++: "<<slotNum<<endl;

    // write to file
    outfileR<<slotNum<<",";

    // -----------------------------------------------
    // identify the blackkeys for each action
    vector<string> blackKeysRecv;
    vector<int> blackActRecv;
    blackKeysRecv.clear();
    blackActRecv.clear();

    size_t blackSize = overBigKeys.size();
    for(size_t i = 0; i < blackSize; i++)
    {
        if(i < slotNum)
        {
            if(overActions[i] == actionSeq)
            {
                blackKeysRecv.push_back(overBigKeys[i]);
                blackActRecv.push_back(overActions[i]);
                //overBigKeys.erase();
            }
        }
        else // blacklist is full
        {
            break;
        }

    }

    size_t blackKeysRecvSize = blackKeysRecv.size();
    //cout<<"* action: "<<actionSeq<<" #cur blackkeys: "<<blackKeysRecvSize<<endl;

    // --------------------------------------------
    // put current blackkeys to cuckoo table
    for(size_t i = 0; i < blackKeysRecvSize; i++)
    {
        if(i<slotNum)
        {
            bool addFlag = cuckooBlackKeyTable.AddKeyPrefix(blackKeysRecv[i],32, 4);
            if(!addFlag)
            {
                cout<<"1 add fail..."<<endl;
            }
        }
    }

    // -------------------------------------------
    // load the previous blackkeys for each action

    size_t keysNum = blackKeysRecvSize;
    size_t keyPresNum = blackkeyPres.size();


    // record old keys
    vector<string> blackkeyPresOld;
    vector<int> actionPresOld;
    size_t keysPreSeq = 0;

    for(size_t i = 0; i < keyPresNum; i++)
    {
        if(keysNum < slotNum && blackActionPres[i] == actionSeq) // IF with the same action
        {
            bool blackFlag = 0;
            int iflowaction;
            int prefix = 32;
            blackFlag = cuckooBlackKeyTable.LookUpKeyAction(blackkeyPres[i],prefix,iflowaction);

            if(!blackFlag)
            {
                blackKeysRecv.push_back(blackkeyPres[i]);
                blackActRecv.push_back(blackActionPres[i]);
                keysNum ++;
                keysPreSeq ++;

                // put it into the blacklist
                bool addFlag = cuckooBlackKeyTable.AddKeyPrefix(blackkeyPres[i],32, 4);
                if(!addFlag)
                {
                    cout<<"2 add fail..."<<endl;
                }
            }

        }

    }
    for(size_t i = keysPreSeq-1; i < keyPresNum; i++)
    {
        if(blackActionPres[i] == actionSeq && keysNum<CUCKOO_BLACK_SIZE)
        {
            blackkeyPresOld.push_back(blackkeyPres[i]);
            actionPresOld.push_back(blackActionPres[i]);
        }
    }


    // -----------------------------------------------
    // wirte blackkeys for each action to blacklist
    size_t bkSize = blackKeysRecv.size();

    for(int i = 0; i < bkSize; i++)
        blackKeyFileOut<<blackKeysRecv[i]<<" "<<32<<" "<<blackActRecv[i]<<endl;

    // write other keys
    size_t bkOldSize = blackkeyPresOld.size();

    for(int i = 0; i < bkOldSize; i++)
        blackKeyFileOut<<blackkeyPresOld[i]<<" "<<32<<" "<<actionPresOld[i]<<endl;



}
void feedbackBlackkeyRL(vector<string>& overBigKeys, vector<int>& overActions, RLearn* rLearn[],
int actionSize, size_ts& slotNums, size_t line)
{

    // -------------------------------------------
    // load the previous blackkeys for each action
    ifstream blackKeyFileIn;
    blackKeyFileIn.open(BLACKFILENAME.c_str());

    vector<string> blackkeyPres;
    vector<int> actionPres;
    string blackKeyStr;
    int prefix = 32;
    int actionInt = 0;

    while(blackKeyFileIn>>blackKeyStr>>prefix>>actionInt)
    {
        blackkeyPres.push_back(blackKeyStr);
        actionPres.push_back(actionInt);

    }
    blackKeyFileIn.clear();
    blackKeyFileIn.close();


    // -------------------------------------------
    // Add blackkey to cuckooTable
    //cout<<"* Add balckkey to cuckooTable!"<<endl;
    float loadFactor = 0.9f;
    int slotNo = 4;
    size_t blackKeySize = CUCKOO_BLACK_SIZE;
    size_t bucketSize = int(blackKeySize/(loadFactor*slotNo))+1;
    int fingerprintNew = 12;
    long MaxKickoutNum = 1000;
    cuckooBlackKeyTable.ClearTable();
    cuckooBlackKeyTable.CuckooTableInit(bucketSize,fingerprintNew,slotNo, \
                                        MaxKickoutNum);

    // write to history files
    ofstream blackKeyFileOut;
    blackKeyFileOut.open(BLACKFILENAME.c_str());

    /*if(line < 80000)
    {
        for(int i = 0; i < actionSize; i++)
        slotNums[i] = CUCKOO_BLACK_SIZE/4;
    }
    else*/
        selectAction(rLearn, actionSize, slotNums);

    printVec(slotNums);

    for(int i = 0; i < actionSize; i++)
    {

        updateBlacklist(overBigKeys, overActions, rLearn[i], i, blackkeyPres, actionPres,
        blackKeyFileOut, slotNums[i]);
    }

    outfileR<<endl;

    blackKeyFileOut.clear();
    blackKeyFileOut.close();


}

void printQList(RLearn* rLearn)
{
    cout<<"Print Qlist!"<<endl;
    for(int i = 0; i < rLearn->_numS; i++)
    {
        cout<<rLearn->qList[i][0].state<<" ";
        for(int j = 0; j < COLNUM; j++)
        {
            // print aValue
            cout<<rLearn->qList[i][j].qValue<<" ";
        }
        cout<<endl;
    }
}

void selectAction(RLearn* rLearn[], int actionSize, size_ts& slotNums)
{
    //if()
    vector<vector<QEntry> > qVecs;

    //range
    int iMin[actionSize], iMax[actionSize];

    for(int i = 0; i < actionSize; i++)
    {
        int stateIndex = rLearn[i]->findIndex(rLearn[i]->_states,rLearn[i]->_state);
        //cout<<"* stateIndex: "<<stateIndex<<endl;

        qVecs.push_back(rLearn[i]->qList[stateIndex]);

        int randNum = (int(rand())%100);
        float randx = float(randNum)/100.0;
        cout<<"* randNum: "<<randNum<<endl;

        // random action
        if(randx < EPSILON0)
        {
            int randa = -1;
            while(rLearn[i]->qList[stateIndex][randa].qValue==-100 || randa == -1 || randa == 0)
            {
                //if(stateQStart > 0 && stateQStart < _numS-1)
                randa =rand()%COLNUM;
                iMin[i] = randa;
                iMax[i] = randa + 1;

            }

        }
        else  // max action
        {
            iMin[i] = 0;
            iMax[i] = COLNUM;
        }
    }

    if(maxHaoOver>0.01)
    {
        iMin[indexHaoOver] = COLNUM-1;
        iMax[indexHaoOver] = COLNUM;
    }

    if(minHaoOver < 0.01 && qVecs[indexMinHaoOver][1].qValue != -100)
    {
        iMin[indexMinHaoOver] = 1;
        iMax[indexMinHaoOver] = 2;
    }
    else if (minHaoOver < 0.01 && qVecs[indexMinHaoOver][2].qValue != -100)
    {
        iMin[indexMinHaoOver] = 2;
        iMax[indexMinHaoOver] = 3;
    }

    vector<QSum> qSums;
    for (int i1 = iMin[0]; i1 < iMax[0]; i1++)
    {

        for(int i2 = iMin[1]; i2 < iMax[1]; i2++)
        {
            for(int i3 = iMin[2]; i3 < iMax[2]; i3++)
            {
                for(int i4 = iMin[3]; i4 < iMax[3]; i4++)
                {
                    bool isValid = (qVecs[0][i1].qValue != -100) && (qVecs[1][i2].qValue != -100) &&
                    (qVecs[2][i3].qValue != -100) && (qVecs[3][i4].qValue != -100);

                    // no qvalue is -100
                    if(isValid)
                    {

                        QSum qSum;

                        // state
                        qSum.states[0] = qVecs[0][i1].state;
                        qSum.states[1] = qVecs[1][i2].state;
                        qSum.states[2] = qVecs[2][i3].state;
                        qSum.states[3] = qVecs[3][i4].state;

                        // action
                        qSum.actions[0] = qVecs[0][i1].action;
                        qSum.actions[1] = qVecs[1][i2].action;
                        qSum.actions[2] = qVecs[2][i3].action;
                        qSum.actions[3] = qVecs[3][i4].action;

                        // actionsum
                        qSum.actionSum = qVecs[0][i1].action + qVecs[1][i2].action
                        + qVecs[2][i3].action + qVecs[3][i4].action;

                        //sum
                        qSum.sum = qVecs[0][i1].qValue + qVecs[1][i2].qValue + qVecs[2][i3].qValue + qVecs[3][i4].qValue;

                        // all actions should be constrainted by the blacklist volume
                        if(qSum.actionSum <= CUCKOO_BLACK_SIZE)
                            qSums.push_back(qSum);
                    }
                }

            }
        }
    }

    // find the peak assignment
    if(qSums.size() == 0) // the old srategy
    {
        for (int i = 0; i < actionSize; i++)
            slotNums[i] = qVecs[i][0].action;
    }
    else
    findMax(qSums, slotNums, actionSize);

    vector<QSum>().swap(qSums);
    vector<vector<QEntry> >().swap(qVecs);
}

void findMax(vector<QSum>& qSums, size_ts& slotNums, int actionSize )
{
    size_t qSumSize = qSums.size();

    if(qSumSize == 0)
    {
        cout<<"No feasible solutions+++++++++++++++++++++++++++++++++"<<endl;
    }

    // the 1st biggest haoover
    /*if(maxHaoOver > 0.01)
    {

    }*/

    float qMax = qSums[qSumSize-1].sum;

    size_t index = qSumSize-1;

    for(size_t i = 0; i < qSumSize; i++)
    {
        if(qSums[i].sum > qMax)
        {
            qMax = qSums[i].sum;
            index = i;
        }
    }

    for(int i = 0; i < actionSize; i++)
        slotNums[i] = qSums[index].actions[i];
}

void printVec(vector<size_t>& vec)
{
    size_t vecSize = vec.size();
    for(size_t i = 0; i < vecSize; i++)
    {
        cout<<vec[i]<<" ";
    }
    cout<<endl;
}

void loadKeys2Filter( string& inFileName, vector<size_t>& mask, intss& vuniquePrefix, intss& vuniqueAggPrefix)
{
    // ------------------------------
    int switchNum = 3;
    int actionSize = 2;

    // ------------------------------
    /*load key file*/
    for(int si = 0; si < switchNum; si++)
    {
        inFileName = inFileName + num2str(si+1);
        std::ifstream infile(inFileName.c_str());
        if(!infile)
            std::cout << "Train File Error " << std::endl;

        int flowPrefixInt;
        string flowStr;
        vector<int> keyprefixlengths;
        vector<string> keys;
        while(infile >> flowPrefixInt >> flowStr)
        {
            keyprefixlengths.push_back(flowPrefixInt);
            keys.push_back(flowStr);
        }
        cout<<"* Key size: "<<keys.size()<<endl;
        infile.clear();
        infile.close();

        // --------------------------------
        // get unique prefix length
        vector<int> uniquePrefix;
        vector<int> uniqueAggPrefix;
        prefixNum(keyprefixlengths, uniquePrefix);

        // ------------------------------
        /*assign actions*/

        vector<int> keyActions;
        assignAction(keys,keyActions,actionSize);

        // --------------------------------
        /*cuckoo table*/
        // load factor
        // m: key number, f: fingerprint width, bc: slots num per bucket,
        // MaxNumKicks: max kickout number in a cuckoo filter
        cout<<"* Init cuckoo table ... ..."<<endl<<endl;
        float a;
        int f,bc;
        a = 0.9;
        bc = 4;
        long m = long(keys.size()/(a*bc))+1;
        f = 12;
        long MaxNumKicks = 1000;
        cuckooTableKey[si].ClearTable();
        cuckooTableKey[si].CuckooTableInit(m,f,bc,MaxNumKicks);

        // --------------------------------
        // Add original key to cuckoo table
        cout<<"* Add key to cuckoo table ... ..."<<endl;
        bool isAddTable;
        for (int i = 0; i < keys.size(); i++)
        {
            isAddTable = cuckooTableKey[si].AddKeyPrefix(keys[i],(keyprefixlengths[i]),keyActions[i]);

            if(isAddTable == 0)
            {
                cout<<"* Flag_add fail"<<endl;
                cout<<"* Order: "<<i<<"  ";
            }

        }

        // -----------------------------------------------
        // init cuckooFilter for flow estimation
        long flowEstSize = FLOW_EST_SIZE;
        m = flowEstSize/(a*bc)+1;
        f = 14;
        cuckooFilterFlowEst[si].ClearTable();
        cuckooFilterFlowEst[si].cuckooFilterInit(m,f,bc,MaxNumKicks);

        // init black table
        m = CUCKOO_BLACK_SIZE/(a*bc)+1;
        cuckooBlackKeyTable[si].ClearTable();
        cuckooBlackKeyTable[si].CuckooTableInit(m,f,bc,
                                            MaxNumKicks);
        // -----------------------------------
        // Init blackkey file
        cout<<"* Write blackkey to file!"<<endl;
        ofstream blackKeyFileOut;
        BLACKFILENAME = "blackkeyfile_" + string(argv[2]) + '_' +
                              string(argv[4])+ "_tstNum_"+ string(argv[5])+"_b"+string(argv[6]);
        blackKeyFileOut.open(BLACKFILENAME.c_str());
        blackKeyFileOut.clear();
        blackKeyFileOut.close();
        uint16_t blackBackSize = 0;
        float feedSumPortion = 0.0;

        // -----------------------------------
        // Init aggr file
        cout<<"* Write aggr to file!"<<endl;
        ofstream aggrFileOut;
        AGGRFILENAME = "aggrfile" + string(argv[2]) + '_' +
                              string(argv[4])+ "_tstNum_"+ string(argv[5])+"_b"+string(argv[6]);
        aggrFileOut.open(AGGRFILENAME.c_str());
        aggrFileOut.clear();
        aggrFileOut.close();

        // -----------------------------------------------
        // Two counters for estimator
        int2s bigNonCounts;
        long2s timeCounts;
        bigNonCounts= vector<vector<int> > (m, vector<int>(bc, 0));
        timeCounts = vector<vector<long> > (m, vector<long>(bc, 0));

        // ---------------------------------------------
        /* init cuckoo filter without aggregation*/
        cout<<"* Init cuckoo filter ... ..."<<endl;

        float storage = strtof(argv[2],NULL); // storage size
        int finger = 0;
        long but = 200000/3;
        char mL0[but][4][20];//
        bzero(&mL0,sizeof(mL0));
        //cout<<sizeof(mL0)<<endl;
        //return 0;
        vector<vector<size_t> > keyCountcs = vector<vector<size_t> > (but, vector<size_t>(4, 0));;
        vector<vector<size_t> > keyCountcs0= vector<vector<size_t> > (but, vector<size_t>(4, 0));;
        vector<vector<size_t> > keyCountDiffs= vector<vector<size_t> > (but, vector<size_t>(4, 0));;
        initCuckoo(keys,keyprefixlengths,keyActions,storage, finger, mL0, cuckooFilter[si],cuckooFilterInit0[si]);
        int finger0 = finger;

        /*for(size_t i = 0; i < 100; i++)
        {
            for(int j=0; j <4; j++)
                cout<<mL0[i][j]<<" ";
        }
        cout<<endl;*/


        // ------------------------------------------------
        // Init aggregation
        bool isInit = 1;
        initAggregation(keys,keyprefixlengths,keyActions,
                        mask, actionSize, storage, isInit, finger,uniqueAggPrefix,mL0,cuckooFilter[si], cuckooAggrKeyTable[si]);

        // output
        vuniquePrefix.push_back(uniquePrefix);
        vuniqueAggPrefix.push_back(uniqueAggPrefix);
    }// end si
}




