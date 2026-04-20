//
// Created by Nadrino on 22/07/2021.
//


#include "SampleSet.h"
#include "GundamGlobals.h"

#include "Logger.h"

#include <memory>


void SampleSet::configureImpl(){

  auto sampleListConfig = GenericToolbox::Json::fetchValue(_config_, {{"sampleList"}, {"fitSampleList"}}, JsonType());
  LogDebugIf(GundamGlobals::isDebug()) << sampleListConfig.size() << " samples defined in the config." << std::endl;

  if( _sampleList_.empty() ){
    // from scratch
    _sampleList_.reserve( sampleListConfig.size() );
    int iSample{0};
    for( auto& sampleConfig : sampleListConfig ){
      _sampleList_.emplace_back();
      _sampleList_.back().setIndex( iSample++ );
      _sampleList_.back().configure( sampleConfig );

      LogDebugIf(GundamGlobals::isDebug()) << "Defined sample: " << _sampleList_.back().getName() << std::endl;

      // remove from the list if not enabled
      if( not _sampleList_.back().isEnabled() ){
        LogDebugIf(GundamGlobals::isDebug()) << "-> removing this sample as it is disabled." << std::endl;
        _sampleList_.pop_back(); iSample--;
      }
    }
  }
  else{
    // for temporary config overrides of propagators,
    // we want to read the config without removing the content of the samples

    // need to check how many samples are enabled. It should match the list.
    size_t nSamples{0};
    for(const auto & sampleConfig : sampleListConfig){
      if( not GenericToolbox::Json::fetchValue(sampleConfig, "isEnabled", true) ) continue;
      nSamples++;
    }
    LogThrowIf(nSamples != _sampleList_.size(), "Can't reload config with different number of samples");

    for( size_t iSample = 0 ; iSample < _sampleList_.size() ; iSample++ ){
      if( not GenericToolbox::Json::fetchValue(sampleListConfig[iSample], "isEnabled", true) ) continue;
      _sampleList_[ iSample ].configure( sampleListConfig[iSample] ); // read the config again
    }
  }

  LogDebugIf(GundamGlobals::isDebug()) << sampleListConfig.size() << " samples were defined." << std::endl;
}
void SampleSet::initializeImpl() {
  for( auto& sample : _sampleList_ ){ sample.initialize(); }
}

void SampleSet::clearEventLists(){
  for( auto& sample : _sampleList_ ){ sample.getEventList().clear(); }
}

std::vector<std::string> SampleSet::fetchRequestedVariablesForIndexing() const{
  std::vector<std::string> out;
  for (auto &sample: _sampleList_) {
    for (auto &binContext: sample.getHistogram().getBinContextList()) {
      for (auto &edges: binContext.bin.getEdgesList()) { GenericToolbox::addIfNotInVector(edges.varName, out); }
    }
  }
  return out;
}
void SampleSet::copyEventsFrom(const SampleSet& src_){
  LogThrowIf(
      src_.getSampleList().size() != this->getSampleList().size(),
      "Can't copy events from mismatching sample lists. src(" << src_.getSampleList().size() << ")"
      << "dst(" << this->getSampleList().size() << ")."
  );

  for( size_t iSample = 0 ; iSample < src_.getSampleList().size() ; iSample++ ){
    this->getSampleList()[iSample].getEventList() = src_.getSampleList()[iSample].getEventList();
  }
}
size_t SampleSet::getNbOfEvents() const {
  return std::accumulate(
      _sampleList_.begin(), _sampleList_.end(), size_t(0),
      [](size_t sum_, const Sample& s_){ return sum_ + s_.getEventList().size(); });
}

void SampleSet::printConfiguration() const {

  LogInfo << _sampleList_.size() << " samples defined." << std::endl;
  for( auto& sample : _sampleList_ ){ sample.printConfiguration(); }

}
std::string SampleSet::getSampleBreakdown() const{
  GenericToolbox::TablePrinter t;

  t << "Sample" << GenericToolbox::TablePrinter::NextColumn;
  t << "# of binned event" << GenericToolbox::TablePrinter::NextColumn;
  t << "total rate (weighted)" << GenericToolbox::TablePrinter::NextLine;

  for( auto& sample : _sampleList_ ){
    t << sample.getName() << GenericToolbox::TablePrinter::NextColumn;
    t << sample.getNbBinnedEvents() << GenericToolbox::TablePrinter::NextColumn;
    t << sample.getSumWeights() << GenericToolbox::TablePrinter::NextLine;
  }

  return t.generateTableString();
}

void SampleSet::throwStatErrors_SimFitToy(bool ThrowMCEvent, bool ThrowBinContent){

  LogInfo << "Throwing statistical error for SimFitToy..." << std::endl;
  LogInfo << "- ThrowMCEvent: " << (ThrowMCEvent?"true":"false") << std::endl;
  LogInfo << "- ThrowBinContent: " << (ThrowBinContent?"true":"false") << std::endl;

  std::vector<Sample>& vec_samples = getSampleList();

  int _nTotalBins{0};
  std::vector<int> _sampleIndicesForEachBin;
  std::vector<int> _localBinIndicesForEachBin;

  int NBinsForEachSample[vec_samples.size()];
  for(unsigned int i_sample=0; i_sample<vec_samples.size(); i_sample++){
    NBinsForEachSample[i_sample] = vec_samples[i_sample].getHistogram().getBinContentList().size();
    _nTotalBins += NBinsForEachSample[i_sample];
  }
  for(unsigned int i_sample=0; i_sample<vec_samples.size(); i_sample++){
    for(int i = 0; i < NBinsForEachSample[i_sample]; ++i) {
      _sampleIndicesForEachBin.push_back( i_sample );
      _localBinIndicesForEachBin.push_back( i );
    }
  }

  TMatrixTSym<double> StatCovValues;
  StatCovValues.ResizeTo(_nTotalBins, _nTotalBins);
  StatCovValues.Zero();
  for(unsigned int idx_global_i=0; idx_global_i<_nTotalBins; idx_global_i++){

    int idx_sample_i = _sampleIndicesForEachBin[idx_global_i];
    const auto& sample_i = vec_samples[idx_sample_i];
    int idx_local_i = _localBinIndicesForEachBin[idx_global_i];

    std::vector<Event*> vec_EvtList_i = sample_i.getHistogram().getBinContextList()[idx_local_i].eventPtrList;

    for(unsigned int idx_global_j=idx_global_i; idx_global_j<_nTotalBins; idx_global_j++){

      int idx_sample_j = _sampleIndicesForEachBin[idx_global_j];
      const auto& sample_j = vec_samples[ idx_sample_j ];
      int idx_local_j = _localBinIndicesForEachBin[idx_global_j];

      std::vector<Event*> vec_EvtList_j = sample_j.getHistogram().getBinContextList()[idx_local_j].eventPtrList;

      double this_BinContent = 0.;
      for(Event* EvtList_i: vec_EvtList_i){
        EventUtils::Indices& EvtIndices_i = EvtList_i->getIndices();
        for(Event* EvtList_j: vec_EvtList_j){
          EventUtils::Indices& EvtIndices_j = EvtList_j->getIndices();
          if(EvtIndices_i.entry==EvtIndices_j.entry){

            bool Is_i_Thrown = EvtList_i->StatThrown;
            bool Is_j_Thrown = EvtList_j->StatThrown;

            if(!Is_i_Thrown && !Is_j_Thrown){

              double this_Poisson_Weight = double(gRandom->Poisson(1)) * EvtList_i->getEventWeight();

              if(ThrowMCEvent){
                EvtList_i->getWeights().current = this_Poisson_Weight;
                EvtList_i->StatThrown = true;

                EvtList_j->getWeights().current = this_Poisson_Weight;
                EvtList_j->StatThrown = true;
              }
            }

            this_BinContent += EvtList_i->getWeights().current;

          } // Found common event between two bins
        } // END Loop over events j
      } // END Loop over events i

      StatCovValues(idx_global_i, idx_global_j) = this_BinContent;
      StatCovValues(idx_global_j, idx_global_i) = this_BinContent;

    }

  }

  // Reset flag
  for(unsigned int idx_global_i=0; idx_global_i<_nTotalBins; idx_global_i++){
    int idx_sample_i = _sampleIndicesForEachBin[idx_global_i];
    const auto& sample_i = vec_samples[idx_sample_i];
    int idx_local_i = _localBinIndicesForEachBin[idx_global_i];
    std::vector<Event*> vec_EvtList_i = sample_i.getHistogram().getBinContextList()[idx_local_i].eventPtrList;
    for(Event* EvtList_i: vec_EvtList_i){
      EvtList_i->StatThrown = false;
    }
  }

  // Now based on BinContent
  // Throw based on BinContent
  // We do off-diagonal
  if(!ThrowBinContent) return;

  for(unsigned int idx_global_i=0; idx_global_i<_nTotalBins; idx_global_i++){

    int idx_sample_i = _sampleIndicesForEachBin[idx_global_i];
    const auto& sample_i = vec_samples[idx_sample_i];
    int idx_local_i = _localBinIndicesForEachBin[idx_global_i];

    std::vector<Event*> vec_EvtList_i = sample_i.getHistogram().getBinContextList()[idx_local_i].eventPtrList;

    for(unsigned int idx_global_j=0; idx_global_j<_nTotalBins; idx_global_j++){

      if(idx_global_i==idx_global_j) continue;

      int idx_sample_j = _sampleIndicesForEachBin[idx_global_j];
      const auto& sample_j = vec_samples[ idx_sample_j ];
      int idx_local_j = _localBinIndicesForEachBin[idx_global_j];

      if(idx_sample_i==idx_sample_j) continue;

      std::vector<Event*> vec_EvtList_j = sample_j.getHistogram().getBinContextList()[idx_local_j].eventPtrList;

      // binContent

      double BinContent_BeforeThrow = StatCovValues(idx_global_i, idx_global_j);
      if(BinContent_BeforeThrow==0.) continue;
      double BinContent_AfterThrow = double( gRandom->Poisson( BinContent_BeforeThrow ) );

      double ThisStatThrowSF = BinContent_AfterThrow/BinContent_BeforeThrow;
      //printf("[JSKIMDEBUG] (i, j) = (%d, %d), Before: %f, After: %f -> SF = %f\n", idx_global_i, idx_global_j, BinContent_BeforeThrow, BinContent_AfterThrow, ThisStatThrowSF);

      // Update data event list weight
      for(Event* EvtList_i: vec_EvtList_i){
        EventUtils::Indices& EvtIndices_i = EvtList_i->getIndices();
        for(Event* EvtList_j: vec_EvtList_j){
          EventUtils::Indices& EvtIndices_j = EvtList_j->getIndices();
          if(EvtIndices_i.entry==EvtIndices_j.entry){

            if( !(EvtList_i->StatThrown) ){
              EvtList_i->getWeights().current *= ThisStatThrowSF;
              EvtList_i->StatThrown = true;
            }
            if( !(EvtList_j->StatThrown) ){
              EvtList_j->getWeights().current *= ThisStatThrowSF;
              EvtList_j->StatThrown = true;
            }

          } // Found common event between two bins
        } // END Loop over events j
      } // END Loop over events i

    }
  }

  for(unsigned int i_sample=0; i_sample<vec_samples.size(); i_sample++){

    Histogram& h = vec_samples[i_sample].getHistogram();
    std::vector<Histogram::BinContent>& vec_binContent = h.getBinContentList();

    for(unsigned int i_bin=0; i_bin<vec_binContent.size(); i_bin++){

      Histogram::BinContent& binContent = vec_binContent[i_bin];

      binContent.sumWeights = 0;
      binContent.sqrtSumSqWeights = 0;

      for (auto *eventPtr: h.getBinContextList()[i_bin].eventPtrList) {
        double weight{eventPtr->getEventWeight()};
        binContent.sumWeights += weight;
        binContent.sqrtSumSqWeights += weight * weight;
      }

      binContent.sqrtSumSqWeights = sqrt(binContent.sqrtSumSqWeights);
    }

  }

}
