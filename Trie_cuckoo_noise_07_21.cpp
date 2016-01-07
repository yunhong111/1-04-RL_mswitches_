/**
trie_cuckoo_noise_07_21.cpp
create by: Yunhong
create time: 07/21/2015
*/


#include "Trie_cuckoo_noise_07_21.h"

ofstream outfileR[3];

int main(int argc, char * argv[])
{
    cout << "\n";
    cout << "PRIME_OPENMP\n";
    cout << "  C++/OpenMP version\n";

    cout << "\n";
    cout << "  Number of processors available = " << omp_get_num_procs ( ) << "\n";
    cout << "  Number of threads =              " << omp_get_max_threads ( ) << "\n";

    // -----------------------------------------------------------------------80
    // Global variables
    CUCKOO_BLACK_SIZE = strtof(argv[6],NULL); // black table size
    FLOW_EST_SIZE = 0;       // flow estimator size
    float overTarget = 0.01;    // target overselection rate

    // -----------------------------------------------------------------------80
    /* initialize random seed: */
    srand (time(NULL));

    // time v.r.
    struct timeval gen_start, gen_end; /* gettimeofday stuff */
    struct timezone tzp;
    long timeInv = 0;

    // -----------------------------------------------------------------------80
    int switchNum = 2;
    const int actionSize = ACTIONSIZE;

    // -----------------------------------------------------------------------80
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

    // -----------------------------------------------------------------------80
    // load key files and add keys to cuckoo filter
    std::string inFileName = argv[1];
    VUPrefix vuniquePrefix;
    VUPrefix vuniqueAggPrefix;

    vuniquePrefix.clear();
    vuniqueAggPrefix.clear();

    long but = 180000/3;
    char mL0[but][actionSize][20];//
    bzero(&mL0,sizeof(mL0));

    int finger = 0;
    int finger0 = 0;

    loadKeys2Filter(inFileName, mask, vuniquePrefix, vuniqueAggPrefix, mL0, argv, finger, finger0, switchNum, actionSize);

    // ----------------------------------------------------------------------80
    vector<vector<size_t> > keyCountcs = vector<vector<size_t> > (but, vector<size_t>(actionSize, 0));
    vector<vector<size_t> > keyCountcs0= vector<vector<size_t> > (but, vector<size_t>(actionSize, 0));
    vector<vector<size_t> > keyCountDiffs= vector<vector<size_t> > (but, vector<size_t>(actionSize, 0));

    uint16_t blackBackSize = 0;
    float feedSumPortion = 0.0;

    // -----------------------------------------------------------------------80
    // File name for caida trace
    const char * fp[] = {"mapped-caida-1","mapped-caida-6", "mapped-caida-11",
                         "mapped-caida-16","mapped-caida-21","mapped-caida-31",
                         "mapped-caida-36","mapped-caida-41","mapped-caida-46",
                         "mapped-caida-51","mapped-caida-56"
                        };

    // -----------------------------------------------------------------------80
    // Open out file, and write result into it
    std::ofstream outfile0[switchNum];
    for(int si = 0; si < switchNum; si++)
    {
        string outfileName = ("outfile0_simple_"+ string(argv[2]) + '_' +
                              string(argv[4]))+ "_tstNum_"+ string(argv[5])+"_b"+string(argv[6])+"_s"+num2str(si)+".csv";
        outfile0[si].open(outfileName.c_str());

        // wirte resource allocation
        //std::ofstream outfileR;
        outfileName = ("resouce_assign_"+ string(argv[2]) + '_' +
                       string(argv[4]))+ "_tstNum_"+ string(argv[5])+"_b"+string(argv[6])+"_s"+num2str(si)+".csv";
        outfileR[si].open(outfileName.c_str());
    }

    // -----------------------------------------------------------------------80
    int ***sec = new int**[switchNum];
    // create RLearn for recvs
    RLearn*** rLearn = new RLearn**[switchNum];
    for (int i = 0; i < switchNum; ++i)
    {
        rLearn[i] = new RLearn*[actionSize];
    }
    //RLearn* rLearn[switchNum][actionSize];
    for(int si = 0; si < switchNum; si++)
    {
        for(int i = 0; i < actionSize; i++)
        {
            rLearn[si][i] = new RLearn();
            initRLearn(rLearn[si][i]);

        }
    }

    // -----------------------------------------------------------------------80
    floats falsePos;
    floats falsePos0;
    floats haoFalsePos;
    floats haoFalsePos0;
    floats haoFalsePosTotal;
    floats haoFalsePos0Total;
    floats overAggr;

    falsePos.assign(switchNum,0);
    falsePos0.assign(switchNum,0);
    haoFalsePos.assign(switchNum,0);
    haoFalsePos0.assign(switchNum,0);
    haoFalsePosTotal.assign(switchNum,0);
    haoFalsePos0Total.assign(switchNum,0);
    overAggr.assign(switchNum,0);


    vector<string> overBigKeys;
    vector<size_t> overBigKeyNos;
    vector<string> blackKeys;
    vector<size_t> blackKeyNos;
    vector<int> blackActions;
    VBkInfo vBkInfo;

    floatss haoOvers;
    haoOvers.resize(switchNum);
    for(int si = 0; si < switchNum; si++)
    {
        haoOvers[si].assign(actionSize, 0);
    }

    floatss haoOversInv;
    haoOversInv.resize(switchNum);
    for(int si = 0; si < switchNum; si++)
    {
        haoOversInv[si].assign(actionSize, 0);
    }

    size_ts countNum;
    size_ts countNum0;
    doubles countIP;
    doubles countIP0;
    size_ts countBlack;

    doubles countIPTotal;
    doubles countIP0Total;
    doubles keySum;
    doubles pktSum;
    doubles keySumTotal;
    doubles pktSumTotal;
    doubles countIPTotalOff;
    doubles countIP0TotalOff;
    doubles keySumTotalOff;
    doubles pktSumTotalOff;
    doubles aggrSum;

    countNum.assign(switchNum, 0);
    countNum0.assign(switchNum, 0);
    countIP.assign(switchNum, 0);
    countIP0.assign(switchNum, 0);
    countBlack.assign(switchNum, 0);

    countIPTotal.assign(switchNum, 0);
    countIP0Total.assign(switchNum, 0);
    keySum.assign(switchNum, 0);
    pktSum.assign(switchNum, 0);
    keySumTotal.assign(switchNum, 0);
    pktSumTotal.assign(switchNum, 0);
    countIPTotalOff.assign(switchNum, 0);
    countIP0TotalOff.assign(switchNum, 0);
    keySumTotalOff.assign(switchNum, 0);
    pktSumTotalOff.assign(switchNum, 0);
    aggrSum.assign(switchNum, 0);


    floatss keySums;
    floatss countIPs;
    keySums.resize(switchNum);
    countIPs.resize(switchNum);

    for(int si = 0; si < switchNum; si++)
    {
        keySums[si].assign(actionSize,0);
        countIPs[si].assign(actionSize,0);
    }

    floatss keySumsInv;
    floatss countIPsInv;
    keySumsInv.resize(switchNum);
    countIPsInv.resize(switchNum);

    for(int si = 0; si < switchNum; si++)
    {
        keySumsInv[si].assign(actionSize,0);
        countIPsInv[si].assign(actionSize,0);
    }

    uint64_t line = 0;
    size_tss slotNums;
    slotNums.resize(switchNum);

    for (int si = 0; si < switchNum; si++)
    {
        slotNums[si].assign(actionSize,0);  // resize each of the contained vectors
    }

    // -----------------------------------------------------------------------80
    // Define a trie
    Trie *trie[3];            // define tree
    for(int si = 0; si < switchNum; si++)
        trie[si] = new Trie();            // define tree

    // -----------------------------------------------------------------------80
    std::ifstream infile;

    // -----------------------------------------------------------------------80
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

        // -------------------------------------------------------------------80
        for(int si = 0; si < switchNum; si++)
        {
            countIPTotalOff[si] += countIP[si];
            countIP0TotalOff[si] += countIP0[si];
            keySumTotalOff[si] += keySum[si];
            pktSumTotalOff[si] += pktSum[si];

            //line = 0;
            countBlack[si] = 0;
            countIP[si] = 0.0f;
            countIP0[si] = 0.0f;
            keySum[si] = 0.0f;
            pktSum[si] = 0.0f;

            keySumsInv[si].clear();
            countIPsInv[si].clear();
            keySumsInv[si].assign(actionSize,0);
            countIPsInv[si].assign(actionSize,0);
        }

        size_t  ei;
        bool isEndFlag = 0;
        size_t updateInvDis = 10000; // interval for display
        size_t readNum = strtof(argv[5],NULL);    // interval for feedback

        while(!isEndFlag )
        {
            // init overselection rate each cycle
            if(line < 100000)
            {
                for (int si = 0; si < switchNum; si++)
                {
                    keySumsInv[si].clear();
                    countIPsInv[si].clear();
                    keySumsInv[si].assign(actionSize,0);
                    countIPsInv[si].assign(actionSize,0);

                    for(int ai = 0; ai < actionSize; ai++)
                        rLearn[si][ai]->clearList();
                }
            }

            // ---------------------------------------------------------------80
            // read file
            strings flows;
            size_ts flowNos;
            readFile0(infile, flows, flowNos, readNum, isEndFlag);
            size_t flowNum = flows.size();

            # pragma omp parallel for \
            shared ( infile,isEndFlag, updateInvDis, line, mask, keySum, pktSum, aggrSum, keySums, countIPs, countIP, countNum, countIP0, countNum0, \
                     countBlack,  trie, vuniquePrefix,vuniqueAggPrefix, cuckooFilter, \
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

                    // -------------------------------------------------------80
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

                        size_ts keySumLocal;
                        size_ts aggrSumLocal;
                        size_ts countIPLocal;
                        size_ts countNumLocal;
                        size_ts countIP0Local;
                        size_ts countNum0Local;

                        keySumLocal.assign(switchNum,0);
                        aggrSumLocal.assign(switchNum,0);
                        countIPLocal.assign(switchNum,0);
                        countNumLocal.assign(switchNum,0);
                        countIP0Local.assign(switchNum,0);
                        countNum0Local.assign(switchNum,0);

                        size_tss keySumsLocal;
                        keySumsLocal.clear();
                        keySumsLocal.resize(switchNum);
                        size_tss countIPsLocal;
                        countIPsLocal.clear();
                        countIPsLocal.resize(switchNum);
                        for(int si = 0; si <switchNum; si++)
                        {
                            keySumsLocal[si].assign(actionSize,0);
                            countIPsLocal[si].assign(actionSize,0);
                        }

                        // ---------------------------------------------------80
                        // look up blackkey
                        int iflowaction;
                        for(int si = 0; si < switchNum; si++)
                        {
                            if(CUCKOO_BLACK_SIZE > 0)
                            {
                                int mi = 24;
                                //for(int mi = 0; mi <= 32-8; mi++)
                                {

                                    subIP = ip & mask[mi]; // mask
                                    flowstr = parsedec2IPV4(ip);//s.str();
                                    prefix = mi+8;

                                    flag_lookblack = cuckooBlackKeyTable[si].LookUpKeyAction(flowstr,prefix,iflowaction);
                                    if (flag_lookblack)
                                    {
                                        countBlack[si] += flowNoInt;
                                        //break;
                                    }
                                }


                            }
                        }

                        // ------------------------------
                        flag_look = 0;
                        if (flag_lookblack == 0) // not a blackkey
                        {
                            // look up key
                            for(int si = 0; si < switchNum; si++)
                            {
                                for(int mi = vuniquePrefix[si].uPres.size()-1; mi >=0; mi--)
                                {

                                    subIP = ip & mask[vuniquePrefix[si].uPres[mi]-8];
                                    flowstr = parsedec2IPV4(subIP);
                                    prefix = vuniquePrefix[si].uPres[mi];


                                    flag_lookkey = cuckooTableKey[si].LookUpKeyActionCount(flowstr,
                                                   prefix,iflowaction,flowNoInt);
                                    if (flag_lookkey)
                                    {
                                        keyType_cur = "1";
                                        flow_action_str = num2str(iflowaction);
                                        //#pragma omp atomic
                                        keySumLocal[si] += (flowNoInt);
                                        for(int ai = 0; ai <actionSize; ai++)
                                        {
                                            if(iflowaction == ai)
                                            {
                                                //#pragma omp atomic
                                                keySumsLocal[si][ai] += (flowNoInt);
                                            }

                                        }
                                        break;
                                    }
                                }
                            }

                            // -------------------------------
                            // lookup key from cuckoo filter
                            vector<int> iactions;
                            bool isAKey = 0;

                            for(int si = 0; si < switchNum; si++)
                            {
                                for(int mi = 0; mi < vuniqueAggPrefix[si].uPres.size(); mi++)
                                {
                                    bool overbigFlag = 0;
                                    subIP = ip & mask[vuniqueAggPrefix[si].uPres[mi]-8];
                                    flowstr = parsedec2IPV4(subIP)+"/"+num2str(vuniqueAggPrefix[si].uPres[mi]);
                                    iactions.clear();
                                    hBuckhit = -1;
                                    slot_ihit = -1;
                                    flag_look = cuckooFilter[si].LookUpKeyActionsCount(flowstr,iactions,flowNoInt, mL0,
                                                keyCountcs,  keyCountcs0, keyCountDiffs,hBuckhit,slot_ihit);
                                    if(hBuckhit!= -1 && slot_ihit!= -1)
                                    {
                                        hBuck = hBuckhit;
                                        slot_i = slot_ihit;
                                    }
                                    if (flag_look != 0)
                                    {
                                        isAKey = 1;
                                        countNumLocal[si] += flag_look;
                                        countIPLocal[si] += flag_look*(flowNoInt);

                                        // each action lookups
                                        for(int ai = 0; ai <actionSize; ai++)
                                        {
                                            for(int li = 0; li < iactions.size(); li++)
                                            {
                                                if(iactions[li] == ai)
                                                {
                                                    countIPsLocal[si][ai] += (flowNoInt);
                                                }

                                            }
                                        }

                                        // overselection count
                                        if(!flag_lookkey && !overbigFlag)
                                        {
                                            overbigFlag = 1;
                                            #pragma omp critical
                                            {
                                                trie[si]->addWordCountNum(DecToBin(flowInt),32, isnonkey, iactions[0], countIPLocal[si]); // ?
                                            }

                                        }
                                    }
                                }
                            }

                            // -------------------------------
                            // look up aggregate key
                            bool isAggregatekey = 0;
                            //if(cuckooAggrKeyTable.mm > 1)
                            for(int si = 0; si < switchNum; si++)

                            {
                                for(int mi = vuniqueAggPrefix[si].uPres.size()-1; mi >=0; mi--)
                                {
                                    subIP = ip & mask[vuniqueAggPrefix[si].uPres[mi]-8];
                                    flowstr = parsedec2IPV4(subIP);
                                    prefix = vuniqueAggPrefix[si].uPres[mi];
                                    isAggregatekey = cuckooAggrKeyTable[si].LookUpKey(flowstr,prefix);
                                    if (isAggregatekey && !flag_lookkey)
                                    {
                                        //#pragma omp atomic
                                        aggrSumLocal[si] += (flowNoInt);
                                        break;

                                    }
                                }
                            }

                        }

                        // -------------------------------
                        // lookup key from cuckoo filterInit0
                        vector<int> iactions;
                        for(int si = 0; si < switchNum; si++)
                        {
                            for(int mi = 0; mi < vuniquePrefix[si].uPres.size(); mi++)
                            {
                                subIP = ip & mask[vuniquePrefix[si].uPres[mi]-8];
                                flowstr = parsedec2IPV4(subIP)+"/"+num2str(vuniquePrefix[si].uPres[mi]);
                                iactions.clear();
                                flag_look0 =cuckooFilterInit0[si].LookUpKeyActions(flowstr,iactions);
                                if (flag_look0)
                                {
                                    //#pragma omp atomic
                                    countNum0Local[si] += flag_look0;
                                    //#pragma omp atomic
                                    countIP0Local[si] += flag_look0*(flowNoInt);

                                }

                            }
                        }
                        // lookup key end

                        // ---------------------------------------------
                        // update global variable
                        line ++;

                        #pragma omp critical
                        {
                            for(int si = 0; si < switchNum; si++)
                            {
                                //#pragma omp atomic
                                pktSum[si] += flowNoInt;
                                //#pragma omp atomic
                                keySum[si] += keySumLocal[si];
                                //#pragma omp atomic
                                aggrSum[si] += aggrSumLocal[si];
                                //#pragma omp atomic
                                countIP[si] += countIPLocal[si];
                                //#pragma omp atomic
                                countNum[si] += countNumLocal[si];
                                //#pragma omp atomic
                                countIP0[si] += countIP0Local[si];
                                //#pragma omp atomic
                                countNum0[si] += countNum0Local[si];
                                for(int ai = 0; ai <actionSize; ai++)
                                {
                                    //#pragma omp atomic
                                    keySums[si][ai] += keySumsLocal[si][ai];
                                    //#pragma omp atomic
                                    countIPs[si][ai] += countIPsLocal[si][ai];

                                    //#pragma omp atomic
                                    keySumsInv[si][ai] += keySumsLocal[si][ai];
                                    //#pragma omp atomic
                                    countIPsInv[si][ai] += countIPsLocal[si][ai];
                                }
                            }
                        }

                    }
                }
                else
                {
                    isEndFlag = 1; // end read flag
                }

                // ---------------------------------
                // update rates and write to file
                #pragma omp critical
                {
                    for(int si = 0; si < switchNum; si ++)
                    {
                        if(size_t(line)%200 == 0)
                        {
                            // update total overselects value
                            countIPTotal[si] = countIPTotalOff[si] + countIP[si];
                            countIP0Total[si] = countIP0TotalOff[si] + countIP0[si];
                            keySumTotal[si] = keySumTotalOff[si] + keySum[si];
                            pktSumTotal[si] = pktSumTotalOff[si] + pktSum[si];

                            // --------------------------------------------
                            // overselection rate
                            // ----------------------------------------------
                            if ((pktSum[si] - keySum[si]) != 0)
                            {
                                falsePos[si] = float(countIP[si]-keySum[si])/float(pktSum[si]-keySum[si]);
                                falsePos0[si] = float(countIP0[si]-keySum[si])/float(pktSum[si]-keySum[si]);
                            }


                            if(keySum[si] != 0)
                            {
                                haoFalsePos[si] = float(countIP[si]-keySum[si])/float(keySum[si]);
                                haoFalsePos0[si] = float(countIP0[si]-keySum[si])/float(keySum[si]);
                            }

                            for(int ai = 0; ai < actionSize; ai++)
                            {
                                if(keySums[si][ai] != 0)
                                {
                                    haoOvers[si][ai] = float(countIPs[si][ai]-keySums[si][ai])/float(keySums[si][ai]);
                                }
                                if(keySumsInv[si][ai] != 0)
                                {
                                    haoOversInv[si][ai] = float(countIPsInv[si][ai]-keySumsInv[si][ai])/float(keySumsInv[si][ai]);
                                }
                            }

                            if(keySumTotal[si] != 0)
                            {
                                haoFalsePosTotal[si] = (countIPTotal[si]-keySumTotal[si])/(keySumTotal[si]);
                                haoFalsePos0Total[si] = (countIP0Total[si]-keySumTotal[si])/(keySumTotal[si]);
                            }

                            overAggr[si] = aggrSum[si]/keySumTotal[si];
                        }

                        if(size_t(line)%updateInvDis == 0 || line%flowNum == 0)
                        {

                            // ------------------------------------------------
                            // Display and write to file
                            outfile0[si]<<"line,"<<line<<",total,"<<haoFalsePosTotal[si]<<",total0,"<<haoFalsePos0Total[si]<<",false,"<<
                                        falsePos[si]<<",over,"<<haoFalsePos[si]<<",false0,"<<falsePos0[si]<<",over0,"<<
                                        haoFalsePos0[si]<<",overaggr,"<<overAggr[si]<<",overcuckoo,"<<(haoFalsePosTotal[si]-overAggr[si])<<",";
                            for(int ai = 0; ai < actionSize; ai++)
                                outfile0[si]<<"over_ai,"<<haoOversInv[si][ai]<<",";
                            for(int ai = 0; ai < actionSize; ai++)
                                outfile0[si]<<"over_ai,"<<haoOvers[si][ai]<<",";

                            outfile0[si]<<"countIP,"<<countIPTotal[si]<<",countIP0,"<<countIP0Total[si]<<",keysum,"<<keySumTotal[si]<<
                                        ",pktSum,"<<pktSumTotal[si]<<",aggrSum,"<<aggrSum[si]<<",blackkey_num,"<<countBlack[si]<<",feedback,"<<blackBackSize<<",feedsumportion,"
                                        <<feedSumPortion<<",finger0,"<<finger0<<",finger,"<<finger<<",time,"<<timeInv<<endl;


                            cout<<endl<<"line "<<line<<" total "<<haoFalsePosTotal[si]<<" total0 "<<haoFalsePos0Total[si]<<" false "<<
                                falsePos[si]<<" over "<<haoFalsePos[si]<<" false0 "<<falsePos0[si]<<" over0 "<<
                                haoFalsePos0[si]<<" overaggr "<<overAggr[si]<<" overcuckoo "<<(haoFalsePosTotal[si]-overAggr[si])<<endl;

                            for(int ai = 0; ai < actionSize; ai++)
                                cout<<"over_ai "<<haoOversInv[si][ai]<<" ";
                            for(int ai = 0; ai < actionSize; ai++)
                                cout<<"over_ai,"<<haoOvers[si][ai]<<",";

                            cout<<"countIP "<<countIPTotal[si]<<" countIP0 "<<countIP0Total[si]<<" keysum "<<keySumTotal[si]<<
                                " pktSum "<<pktSumTotal[si]<<" aggrSum "<<aggrSum[si]<<" blackkey_num "<<countBlack[si]<<" feedback "<<blackBackSize<<" feedsumportion "
                                <<feedSumPortion<<" finger0 "<<finger0<<" finger "<<finger<<" time "<<timeInv<<endl<<endl;
                        }
                    }
                }

            }

            // ---------------------------------------------
            // feed back portionFeedBack% overselections and overselection rates for actions
            //if(size_t(line)%2000000 == 0 )
            {

                // ---------------------------
                // the overselects for feedback

                for(int si = 0; si < switchNum; si++)
                {
                    blackKeys.clear();
                    blackKeyNos.clear();
                    blackActions.clear();

                    vector<char> word;
                    trie[si]->getLeaf(trie[si]->root,word,blackKeys,blackKeyNos, blackActions);
                    trie[si]->deleteChild(trie[si]->root);
                    delete trie[si];
                    trie[si] = new Trie();            // define tree

                    int init = 0;
                    float overTotalSum = accumulate(blackKeyNos.begin(), blackKeyNos.end(), init);

                    // sort balckkeys
                    keySort3(blackKeys,blackKeyNos, blackActions);

                    float portionFeedBack = strtof(argv[4],NULL);
                    blackBackSize = portionFeedBack*blackKeys.size();

                    blackKeys.erase(blackKeys.begin()+blackBackSize, blackKeys.end());
                    blackKeyNos.erase(blackKeyNos.begin()+blackBackSize, blackKeyNos.end());
                    blackActions.erase(blackActions.begin()+blackBackSize, blackActions.end());

                    float feedSum = accumulate(blackKeyNos.begin(), blackKeyNos.end(), init);
                    feedSumPortion = feedSum/overTotalSum;

                    // output
                    BkInfo bkInfo;
                    bkInfo.bkActions = blackActions;
                    bkInfo.bkPres.assign(blackKeys.size(),32);
                    bkInfo.bks = blackKeys;
                    bkInfo.bkSizes = blackKeyNos;

                    vBkInfo.push_back(bkInfo);
                }


                // -----------------------------
                // Get instant reward and Q learn
                //cout<<"* Q learn!"<<endl;
                for(int si = 0; si < switchNum; si++)
                {
                    for(int ri = 0; ri < actionSize; ri++)
                    {
                        //cout<<rLearn[ri]<<endl;
                        rLearn[si][ri]->update(slotNums[si][ri],haoOvers[si][ri]);
                        //cout<<"* qleran!"<<endl;

                        // compute ebuse0
                        if(line<1000000)
                        {
                            EPSILON0 = 1.0;
                        }
                        else
                        {
                            EPSILON0 = 0.20;
                        }
                        rLearn[si][ri]->qLearn();
                        printQList(rLearn[si][ri]);
                    }
                }


                // ------------------------------------
                // feedback
                // ------------------------------------

                //time before calling function
                gettimeofday(&gen_start, &tzp);

                // ---------------------------------
                // call function
                cout<<"* feedback bks!"<<endl;
                feedbackBlackkeyRL(vBkInfo, rLearn, actionSize,switchNum,slotNums, line);
                cout<<"* feedback bks end!"<<endl;

                // ------------------------------
                // time after calling function
                gettimeofday(&gen_end, &tzp);

                // time interval
                timeInv = print_elapsed("Aggr: ", &gen_start,&gen_end, 1);

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

    //cuckooFilter.ClearTable();
    //cuckooBlackKeyTable.ClearTable();
    //cuckooTableKey.ClearTable();
    //cuckooAggrKeyTable.ClearTable();

    // (outfile0.close());
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
                     , vector<string>& blackkeyPres, vector<int>& blackActionPres, ofstream& blackKeyFileOut, size_t& slotNum, int si)
{
    // for each recv
    // get the suggested action
    //slotNum = rLearn->selectActionSuggest();
    cout<<"* suggested action++++++++++++++++++++++++++++++++++: "<<slotNum<<endl;

    // write to file
    outfileR[si]<<slotNum<<",";

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
            bool addFlag = cuckooBlackKeyTable[si].AddKeyPrefix(blackKeysRecv[i],32, 4);
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
            blackFlag = cuckooBlackKeyTable[si].LookUpKeyAction(blackkeyPres[i],prefix,iflowaction);

            if(!blackFlag)
            {
                blackKeysRecv.push_back(blackkeyPres[i]);
                blackActRecv.push_back(blackActionPres[i]);
                keysNum ++;
                keysPreSeq ++;

                // put it into the blacklist
                bool addFlag = cuckooBlackKeyTable[si].AddKeyPrefix(blackkeyPres[i],32, 4);
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
void feedbackBlackkeyRL(VBkInfo& vBkInfo, RLearn** rLearn[],
                        int actionSize, int switchNum, size_tss& slotNums, size_t line)
{

    selectAction(rLearn, actionSize, switchNum,slotNums);

    // -------------------------------------------
    // load the previous blackkeys for each action
    for(int si = 0; si < switchNum; si++)
    {
        cout<<"* feedbackBlackkeyRL:si: "<<si<<endl;
        string blackFileName = BLACKFILENAME + num2str(si);
        ifstream blackKeyFileIn;
        blackKeyFileIn.open(blackFileName.c_str());

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
        //cuckooBlackKeyTable[si].ClearTable();
        cuckooBlackKeyTable[si].CuckooTableInit(bucketSize,fingerprintNew,slotNo, \
                                                MaxKickoutNum);

        // write to history files
        ofstream blackKeyFileOut;
        blackKeyFileOut.open(blackFileName.c_str());


        /*if(line < 80000)
        {
            for(int i = 0; i < actionSize; i++)
            slotNums[i] = CUCKOO_BLACK_SIZE/4;
        }
        else*/


        //printVec(slotNums);

        for(int i = 0; i < actionSize; i++)
        {
            strings overBigKeys = vBkInfo[si].bks;
            ints overActions = vBkInfo[si].bkActions;

            updateBlacklist(overBigKeys, overActions, rLearn[si][i], i, blackkeyPres, actionPres,
                            blackKeyFileOut, slotNums[si][i], si);
        }

        outfileR[si]<<endl;

        blackKeyFileOut.clear();
        blackKeyFileOut.close();
    }


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

void selectAction(RLearn** rLearn[], int actionSize, int switchNum, size_tss& slotNums)
{
    floatss actionVec;
    actionVec.resize(switchNum*actionSize);
    for(int si = 0; si < actionVec.size(); si ++)
    {
        actionVec[si].assign(COLNUM,0);
    }

    //if()
    vector<vector<vector<QEntry> > > qVecs;
    qVecs.resize(switchNum);

    for(int si = 0; si < switchNum; si++)
    {
        qVecs[si].resize(actionSize);
    }

    //range
    int iMin[switchNum][actionSize], iMax[switchNum][actionSize];
    for(int si = 0; si < switchNum; si++)
    {
        for(int i = 0; i < actionSize; i++)
        {
            cout<<"* si: "<<si<<" ai: "<<i<<endl;
            int stateIndex = rLearn[si][i]->findIndex(rLearn[si][i]->_states,rLearn[si][i]->_state);
            cout<<"* stateIndex: "<<stateIndex<<endl;

            qVecs[si][i] = (rLearn[si][i]->qList[stateIndex]);

            cout<<"* actionVec: "<<endl;
            for(int ci = 0; ci < COLNUM; ci++)
            {
                actionVec[si*actionSize+i][ci]  = rLearn[si][i]->qList[stateIndex][ci].action;
                //cout<<"* actionVec[si*actionSize+i][ci]: "<<actionVec[si*actionSize+i][ci]<<" ";

            }
            //cout<<endl;

            int randNum = (int(rand())%100);
            float randx = float(randNum)/100.0;
            cout<<"* randNum: "<<randNum<<endl;

            // random action
            if(randx < EPSILON0)
            {
                int randa = -1;
                while(rLearn[si][i]->qList[stateIndex][randa].qValue==-100 || randa == -1/* || randa == 0*/)
                {
                    //if(stateQStart > 0 && stateQStart < _numS-1)
                    randa =rand()%COLNUM;
                    cout<<"* randa: "<<randa<<endl;
                    iMin[si][i] = randa;
                    iMax[si][i] = randa + 1;

                }

            }
            else  // max action
            {
                iMin[si][i] = 0;
                iMax[si][i] = COLNUM;
            }
        }
    }

    /*if(maxHaoOver>0.01)
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
    }*/
    cout<<"* find max actions!"<<endl;

    vector<vector<QSum> > qSums;
    qSums.resize(switchNum);

    for(int si = 0; si < switchNum; si ++)
    {
        for (int i1 = iMin[si][0]; i1 < iMax[si][0]; i1++)
        {

            for(int i2 = iMin[si][1]; i2 < iMax[si][1]; i2++)
            {

                bool isValid = (qVecs[si][0][i1].qValue != -100) && (qVecs[si][1][i2].qValue != -100) ;
                // no qvalue is -100
                if(isValid)
                {

                    QSum qSum;

                    // state
                    qSum.states[0] = qVecs[si][0][i1].state;
                    qSum.states[1] = qVecs[si][1][i2].state;

                    // action
                    qSum.actions[0] = qVecs[si][0][i1].action;
                    qSum.actions[1] = qVecs[si][1][i2].action;

                    // actionsum
                    qSum.actionSum = qVecs[si][0][i1].action + qVecs[si][1][i2].action;

                    //sum
                    qSum.sum = qVecs[si][0][i1].qValue + qVecs[si][1][i2].qValue;

                    // all actions should be constrainted by the blacklist volume
                    if(qSum.actionSum <= CUCKOO_BLACK_SIZE)
                        qSums[si].push_back(qSum);
                }
            }
        }
    }

    // find the peak assignment
    cout<<"* find slotNums!"<<endl;
    for(int si = 0; si < switchNum; si++)
    {
        if(qSums[si].size() == 0) // the old srategy
        {
            for (int i = 0; i < actionSize; i++)
                slotNums[si][i] = qVecs[si][i][0].action;
        }
        else
        {
            findMax(qSums[si], slotNums[si], actionSize, switchNum);
        }
    }

    // find the index of the selected action
    cout<<"* find the index of selected action!"<<endl;

    for(int si = 0; si < switchNum; si ++)
    {
        printVec(slotNums[si]);
        //printVec(actionVec[si]);
        for(int ai = 0; ai < actionSize; ai++)
        {
            rLearn[si][ai]->_actionSuggestIndex = rLearn[si][ai]->findIndex(actionVec[si*actionSize+ai], (float)slotNums[si][ai]);
            cout<<"* rLearn[si][ai]->_actionSuggestIndex: "<<rLearn[si][ai]->_actionSuggestIndex<<endl;
        }
    }

    //vector<QSum>().swap(qSums);
    //vector<vector<QEntry> >().swap(qVecs);
}

void findMax(vector<QSum>& qSums, size_ts& slotNums, int actionSize, int switchNum )
{
    size_t qSumSize = qSums.size();

    if(qSumSize == 0)
    {
        cout<<"No feasible solutions+++++++++++++++++++++++++++++++++"<<endl;
    }

    float qMax = qSums[0].sum;

    size_t index = 0;

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

void loadKeys2Filter(string& inFileName, vector<size_t>& mask, VUPrefix& vuniquePrefix,
                     VUPrefix& vuniqueAggPrefix, char mL0[][ACTIONSIZE][20], char * argv[], int& finger, int& finger0,int switchNum, int actionSize)
{

    // ------------------------------
    /*load key file*/
    for(int si = 0; si < switchNum; si++)
    {
        cout<<"* switch: "<<si<<"++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;

        string inFileName0 = inFileName + num2str(si+1);
        std::ifstream infile(inFileName0.c_str());
        if(!infile)
            std::cout << "Train File Error " << std::endl;

        int flowPrefixInt;
        string flowStr;
        int actionInt;

        vector<int> keyprefixlengths;
        vector<string> keys;
        vector<int> keyActions;

        while(infile >> flowPrefixInt >> flowStr>>actionInt)
        {
            keyprefixlengths.push_back(flowPrefixInt);
            keys.push_back(flowStr);
            keyActions.push_back(actionInt);
        }
        cout<<"* Key size: "<<keys.size()<<endl;
        infile.clear();
        infile.close();

        // --------------------------------
        // get unique prefix length
        vector<int> uniquePrefix;
        vector<int> uniqueAggPrefix;
        prefixNum(keyprefixlengths, uniquePrefix);

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
        //int finger = 0;
        long but = 200000/3;
        //cout<<sizeof(mL0)<<endl;
        //return 0;

        initCuckoo(keys,keyprefixlengths,keyActions,storage, finger, mL0, cuckooFilter[si],cuckooFilterInit0[si]);
        finger0 = finger;

        cout<<"* finger: "<<finger<<endl;

        // ------------------------------------------------
        // Init aggregation
        bool isInit = 1;

        initAggregation(keys,keyprefixlengths,keyActions,
                        mask, actionSize, storage, isInit, finger,uniqueAggPrefix,mL0,cuckooFilter[si], cuckooAggrKeyTable[si]);

        // output
        UPrefix uPrefix;
        uPrefix.uPres  = uniquePrefix;
        vuniquePrefix.push_back(uPrefix);

        uPrefix.uPres.clear();
        uPrefix.uPres  = uniqueAggPrefix;

        vuniqueAggPrefix.push_back(uPrefix);
    }// end si
}




