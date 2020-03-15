/* <===<===<===<===<===<===<===<===<===~===>===>===>===>===>===>===>===>===>
 * File Name:    AutoTTTSAlg.cpp
 * Author:       Xin-Xin MA, Hao-Kai SUN
 * Created:      2019-11-09 Sat 15:49:02 CST
 * <<=====================================>>
 * Last Updated: 2020-01-07 Tue 14:05:12 CST
 *           By: Hao-Kai SUN
 *     Update #: 209
 * <<======== COPYRIGHT && LICENSE =======>>
 *
 * Copyright © 2019 SUN Hao-Kai <spin.hk@outlook.com>. All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs.  If not, see <https://www.gnu.org/licenses/>.
 *
 * ============================== CODES ==============================>>> */
#include "AutoTTTSAlg/AutoTTTSAlg.h"

#include "EmcRecEventModel/RecEmcShower.h"
#include "EvTimeEvent/RecEsTime.h"

#include "EventModel/Event.h"
#include "EventModel/EventHeader.h"
#include "EventModel/EventModel.h"
#include "EventNavigator/EventNavigator.h"
#include "EvtRecEvent/EvtRecEvent.h"

#include "GaudiKernel/IJobOptionsSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "EvtRecEvent/EvtRecVeeVertex.h"
#include "HadronInfo/ParticleInf.h"
#include "TupleSvc/DecayChainSvc.h"
#include "AutoTTTSAlg/AutoTTTSAlg.h"
#include "AutoTTTSAlg/selector/SignalCandidate.h"
#include "HadronInfo/LambdaInfo.h"

#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
using namespace std;

bool GetBestCandidate(CDDecayList& signalList, CDDecayList::iterator& best) {
    best = signalList.particle_end();
    if (signalList.empty()) {
        return false;
    }
    double minQ = 1E6, Q = 1E6;
    for (CDDecayList::iterator itr = signalList.particle_begin();
         itr != signalList.particle_end(); ++itr) {
        const CDCandidate& signal = (*itr).particle();
        LambdaInfo lambdaInfo(signal);
        Q = lambdaInfo.chisq();
        ///  cout << "Q = " << Q << endl;
        if (Q < minQ) {
            minQ = Q;
            best = itr;
        }
    }
    return true;
}


AutoTTTSAlg::AutoTTTSAlg(const std::string& name, ISvcLocator* pSvcLocator)
    : Algorithm(name, pSvcLocator) {
    declareProperty("UseMatch", m_useMatch = false);
    declareProperty("Mother", m_Mother = 443);
    declareProperty("Ecm", m_Ecm = 3.097);
    declareProperty("ReadBeamE", m_readBeamE = false);
    declareProperty("FillMCInfo", m_fillMCInfo = false);
    declareProperty("FillMCParAll", m_fillMCParA = true);
    declareProperty("MinChargedTracks", m_minTrks = 0);
    declareProperty("MaxChargedTracks", m_maxTrks = 100);
    declareProperty("MinShowers", m_minShowers = 0);
    declareProperty("MaxShowers", m_maxShowers = 100);
    declareProperty("InfoLevel", m_InfoLvl = 0);
    declareProperty("SigFid", m_sigFid);
    declareProperty("TagFid", m_tagFid);
    m_tagFid.clear();
    m_sigFid.clear();
}

AutoTTTSAlg::~AutoTTTSAlg() { 
}

StatusCode AutoTTTSAlg::initialize() {
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "in initialize()" << endmsg;
    time_t now = time(0);
    char* dt = ctime(&now);
    std::cout << "The date : " << dt << std::endl;
    this->InitTuple();
    this->LoadSvc(log);
    
    string tail;
    for (int i=0; i< m_decayTreeList.size(); ++i) {
        stringstream ss;
        ss << i;
        ss >> tail;
        m_sigTuple[i] = this->BookTuple("sig" + tail);
        if (m_sigTuple[i]) {
            cout << "Book tuple 'sig" << tail 
                << "' successful!" << endl;
            cout << "Tag mode: [";
            for(vector<int>::iterator itr = m_tagFid.begin();
                    itr != m_tagFid.end(); 
                    ++itr){
                cout << *itr << ", ";
            }
            cout << "]" << endl;
            cout << "Sig mode: [";
            vector<int> sig_fid = m_decayTreeList[i].GetFID();
            for(vector<int>::iterator itr = sig_fid.begin();
                    itr != sig_fid.end(); 
                    ++itr){
                cout << *itr << ", ";
            }
            cout << "]" << endl;
        }
        m_SigTupleSvc[i].BindTuple(m_sigTuple[i]);
        m_SigTupleSvc[i].Register(*gBeamInfoSvc);

        gDecayChainSvc.SetDecayTree(DecayTree(m_tagFid));
        gDecayChainSvc.SetTitle("_tag");
        gDecayChainSvc.BulkRegister(m_SigTupleSvc[i]);
        // m_sigTuple[i]->addItem("findSignal", m_findSignal[i]);

        m_sigTuple[i]->addItem("tagmrec", m_recMass[i]);
        gDecayChainSvc.SetTitle("_sig");
        gDecayChainSvc.SetDecayTree(m_decayTreeList[i]);
        gDecayChainSvc.BulkRegister(m_SigTupleSvc[i]);
        // m_SigTupleSvc[i].Register(*gPhotonConverSvc);
        gMCTruthInfoSvc->SetName("");
        m_SigTupleSvc[i].Register(*gMCTruthInfoSvc);
    }

    if (m_fillMCInfo) {
        mctp = this->BookTuple("mc");
        if (mctp) {
            cout << "Book Tuple 'mc' successful!" << endl;
        }
        m_mcTupleSvc.BindTuple(mctp);
        m_mcTupleSvc.Register(*gMCTruthInfoSvc);
    }
    log << MSG::INFO << "successfully return from initialize()" << endmsg;
    return StatusCode::SUCCESS;
}

void AutoTTTSAlg::InitTuple(){
    m_decayTreeList.clear();
    // get the decaytree
    vector<int> tmp;
    for(vector<int>::iterator itr= m_sigFid.begin();
        itr != m_sigFid.end();
        ++itr){
        if ((*itr) == 0) {
            if (tmp.empty()) continue;
            m_decayTreeList.push_back(DecayTree(tmp)); 
            tmp.clear();
            continue;
        }
        tmp.push_back(*itr);
    }
    if (!tmp.empty()) {
        m_decayTreeList.push_back(DecayTree(tmp));
    }
    int nTuple = m_decayTreeList.size();
    cout << "Info in <AutoTTTSAlg::InitTuple()>: " 
        << "#signal channels " << nTuple << endl;
    // m_SigTupleSvc = new TupleSvc[nTuple];
    // m_mcTupleSvc = TupleSvc(10, 10, 10, 10);
    // m_findSignal = new NTuple::Item<int> [nTuple];
}

vector<int> AutoTTTSAlg::getChannelCC(const vector<int>& fid) {
    vector<int> fidcc;
    fidcc.clear();
    for (vector<int>::const_iterator itr = fid.begin(); itr != fid.end(); ++itr) {
    }
    return fidcc;
}

StatusCode AutoTTTSAlg::finalize() {
    MsgStream logmess(msgSvc(), name());
    cout << __func__ << " " << __LINE__ << endl;
    // delete[] m_SigTupleSvc;
    cout << __func__ << " " << __LINE__ << endl;
    // delete[] m_findSignal;
    cout << __func__ << " " << __LINE__ << endl;
    return StatusCode::SUCCESS;
}

NTuple::Tuple* AutoTTTSAlg::BookTuple(const string& name,
                                     const string& comment) {
    MsgStream log(msgSvc(), Algorithm::name());
    NTuplePtr tmptuple(ntupleSvc(), string("FILE/") + name);
    if (tmptuple) {
        return tmptuple;
    }

    return ntupleSvc()->book(string("FILE/") + name, CLID_ColumnWiseTuple,
                             comment);
}

StatusCode AutoTTTSAlg::LoadSvc(MsgStream& log) {
    // load MCTruthInfo service
    IMCTruthInfo* i_matruthInfo;
    StatusCode sc_MCInfo = service("MCTruthInfo", i_matruthInfo);
    if (sc_MCInfo.isFailure()) {
        std::cout << "could not load MCTruthInfo!" << std::endl;
        return 0;
    }
    gMCTruthInfoSvc = dynamic_cast<MCTruthInfo*>(i_matruthInfo);
    gMCTruthInfoSvc->SetDecayTree(DecayTree(m_tagFid));
    
    // load BeamInfoSvc
    IBeamInfoSvc* i_BeamInfoSvc;
    StatusCode sc_BeamInfoSvc = service("BeamInfoSvc", i_BeamInfoSvc);
    if (sc_BeamInfoSvc.isFailure()) {
        std::cout << "could not load BeamInfoSvc!" << std::endl;
        return 0;
    }
    gBeamInfoSvc = dynamic_cast<BeamInfoSvc*>(i_BeamInfoSvc);
    if (!m_readBeamE) {
        gBeamInfoSvc->SetEcm(m_Ecm);
    }
    // load PhotonConverSvc
    IPhotonConverSvc* i_PhotonConverSvc;
    StatusCode sc_PhotonConverSvc = service("PhotonConverSvc", i_PhotonConverSvc);
    if (sc_PhotonConverSvc.isFailure()) {
        std::cout << "could not load PhotonConverSvc!" << std::endl;
        return 0;
    }
    gPhotonConverSvc = dynamic_cast<PhotonConverSvc*>(i_PhotonConverSvc);
    gPhotonConverSvc->SetDecayTree(DecayTree(m_tagFid));

    // load MCTruthMatchSvc
    IMCTruthMatchSvc* i_matchSvc;
    StatusCode sc_MC = service("MCTruthMatchSvc", i_matchSvc);
    if (sc_MC.isFailure()) {
        std::cout << "could not load MCTruthMatchSvc!" << std::endl;
        return 0;
    }
    m_matchSvc = dynamic_cast<MCTruthMatchSvc*>(i_matchSvc);

    // load McDecayModeSvc
    IMcDecayModeSvc* i_McDecayModeSvc;
    StatusCode sc_McDecayModeSvc = service("McDecayModeSvc", i_McDecayModeSvc);
    if (sc_McDecayModeSvc.isFailure()) {
        log << MSG::FATAL << "could not load McDecayModeSvc" << endmsg;
        return sc_McDecayModeSvc;
    }
    m_McDecayModeSvc = dynamic_cast<McDecayModeSvc*>(i_McDecayModeSvc);
    
    IDCListSvc* i_DCListSvc;
    StatusCode sc_DCListSvc = service("DCListSvc", i_DCListSvc);
    if (sc_DCListSvc.isFailure()) {
        log << MSG::FATAL << "could not load DCListSvc" << endmsg;
        return sc_DCListSvc;
    }
    m_DCListSvc = dynamic_cast<DCListSvc*>(i_DCListSvc);
    return StatusCode::SUCCESS;
}

void AutoTTTSAlg::FillInfo(const CDCandidate& tag, 
                           const DecayTree& decayTree,
                           TupleSvc& tupleSvc) {
    // fill beam status info
    if (m_InfoLvl > 0) {
        cout << "AutoTTTSAlg::FillInfo() " << endl;
        cout << "\t fillInfo begin " << endl;
    }
    tupleSvc << *gBeamInfoSvc;
    // fill the kininatic infomation
    ///   cout << "gDecayChainSvc.Express" << endl;
    if (m_InfoLvl > 0) {
        cout << "\t gDecayChainSvc.Express() " << endl;
    }
    gDecayChainSvc.SetDecayTree(decayTree);
    gDecayChainSvc.Express(tag, tupleSvc);
    if (m_InfoLvl > 0) {
        cout << "\t Fill MCtrthInfo " << endl;
    }
    tupleSvc << *gMCTruthInfoSvc;
}

void AutoTTTSAlg::SearchSignal(int index, const DecayTree& signalTree,
        TupleSvc& tupleSvc) {
    if (m_readBeamE) {
        gBeamInfoSvc->GetEcm(m_Ecm);
        localSignalCandidate.SetEcm(m_Ecm);
        //   gPhotonConverSvc->SetEcm(m_Ecm);
    } else {
        localSignalCandidate.SetEcm(m_Ecm);
        // gPhotonConverSvc->SetEcm(m_Ecm);
    }
    CDDecayList tagList(localSignalCandidate);
    tagList = m_DCListSvc->DecayList(m_tagFid);
    if (tagList.empty()) return;

    if (m_InfoLvl >= 2) {
        if (tagList.size() > 1) {
            std::cout << "*** Multi Signal: " << tagList.size() << " ***"
                      << std::endl;
        }
        // print(signal);
    }
    if (m_InfoLvl >= 2) {
        cout << "Search for best candidate!" << endl;
    }
    CDDecayList::iterator best = tagList.particle_end();
    bool status = GetBestCandidate(tagList, best);
    if (best == tagList.particle_end()) {
        std::cout << "Error: didnot find a best candidate" << std::endl;
        return;
    }
    const CDCandidate& tag = (*best).particle();
    if (m_InfoLvl >= 1) {
        cout << "\tFill tag Info!" << endl;
        cout << "\tp4 = " << tag.p4() << endl;
    }
    HepLorentzVector p4Beam(0.011*m_Ecm, 0, 0, m_Ecm);
    m_recMass[index]= (p4Beam - tag.p4()).m();
    gDecayChainSvc.SetTitle("_tag");
    FillInfo(tag, DecayTree(m_tagFid), tupleSvc);

    CDDecayList sigList = m_DCListSvc->DecayList(signalTree.GetFID());
    if (sigList.empty()) {
        tupleSvc.Write();
        return;
    }

    CDDecayList::iterator sigbest = sigList.particle_end();
    for (CDDecayList::iterator itr = sigList.particle_begin();
            itr != sigList.particle_end(); ++itr) {
        const CDCandidate& sig = (*itr).particle();
        if (sig.overlap(tag)) continue;
        sigbest = itr;
        // keep the first
        break;
    }
    if (sigbest == sigList.particle_end()){
        tupleSvc.Write();
        return;
    }
    
    const CDCandidate& sig = (*sigbest).particle();
    if (m_InfoLvl >= 1) {
        cout << "\tFill sig Info!" << endl;
        cout << "\tp4 = " << sig.p4() << endl;
    }
    gDecayChainSvc.SetTitle("_sig");
    FillInfo(sig, signalTree, tupleSvc);
    tupleSvc.Write();
    if (m_InfoLvl >= 1) {
        cout << "Write" << endl;
    }
}

StatusCode AutoTTTSAlg::execute() {
    SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(),
                                          "/Event/EvtRec/EvtRecEvent");
    int chrTrks = evtRecEvent->totalCharged();
    if (chrTrks < m_minTrks || chrTrks > m_maxTrks) {
        return StatusCode::SUCCESS;
    }

    int showers = evtRecEvent->totalNeutral();
    if (showers < m_minShowers || showers > m_maxShowers) {
        return StatusCode::SUCCESS;
    }
    if (m_fillMCInfo) {
        m_mcTupleSvc << *gMCTruthInfoSvc;
        m_mcTupleSvc.Write();
    }
    m_DCListSvc->init();
    for (int i =0; i< m_decayTreeList.size(); ++i){
        SearchSignal(i, m_decayTreeList[i], m_SigTupleSvc[i]);
    }
    return StatusCode::SUCCESS;
}
