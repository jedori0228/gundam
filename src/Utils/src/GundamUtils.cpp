//
// Created by Adrien Blanchet on 25/05/2023.
//

#include "GundamUtils.h"
#include "VersionConfig.h" // the only place it is included
#include "SourceConfig.h" // the only place it is included

#include "GenericToolbox.Fs.h"
#include "Logger.h"

#include <sstream>

#ifndef DISABLE_USER_HEADER
LoggerInit([]{ Logger::getUserHeader() << "[" << FILENAME << "]"; });
#endif


namespace GundamUtils {

  // forward version config that was auto generated by CMake
  std::string getVersionStr(){ return GundamVersionConfig::getVersionStr(); }
  std::string getVersionFullStr(){
    std::stringstream ss;
    ss << GundamVersionConfig::getVersionStr();
    if( not GundamVersionConfig::getVersionPostCommitNb().empty() ){
      ss << "+" << GundamVersionConfig::getVersionPostCommitNb() << "-" << GundamVersionConfig::getVersionPostCommitHash();
      ss << "/" << GundamVersionConfig::getVersionBranch();
    }
    return ss.str();
  }
  std::string getSourceCodePath(){
    return {SOURCE_CONFIG_H};
  }
  bool isNewerOrEqualVersion( const std::string& minVersion_ ){
    if( GundamUtils::getVersionStr() == "X.X.X" ){
      LogAlert << "Can't check version requirement. Assuming OK." << std::endl;
      return true;
    }
    auto minVersionSplit = GenericToolbox::splitString(minVersion_, ".");
    LogThrowIf(minVersionSplit.size() != 3, "Invalid version format: " << minVersion_);

    // stripping "f" tag
    if( minVersionSplit[2].back() == 'f' ){ minVersionSplit[2].pop_back(); }

    if( std::stoi(GundamVersionConfig::getVersionMajor()) > std::stoi(minVersionSplit[0]) ) return true; // major is GREATER -> OK
    if( std::stoi(GundamVersionConfig::getVersionMajor()) < std::stoi(minVersionSplit[0]) ) return false; // major is LOWER -> NOT OK
    // major is equal -> next

    if( std::stoi(GundamVersionConfig::getVersionMinor()) > std::stoi(minVersionSplit[1]) ) return true; // minor is GREATER -> OK
    if( std::stoi(GundamVersionConfig::getVersionMinor()) < std::stoi(minVersionSplit[1]) ) return false; // minor is LOWER -> NOT OK
    // minor is equal -> next

    if( std::stoi(GundamVersionConfig::getVersionMicro()) > std::stoi(minVersionSplit[2]) ) return true; // revision is GREATER -> OK
    if( std::stoi(GundamVersionConfig::getVersionMicro()) < std::stoi(minVersionSplit[2]) ) return false; // revision is LOWER -> NOT OK
    // minor is equal -> OK
    return true;
  }

  std::string generateFileName(const CmdLineParser& clp_, const std::vector<AppendixEntry>& appendixDict_){
    std::vector<std::string> appendixList{};

    int maxArgLength{64};

    for( const auto& appendixDictEntry : appendixDict_ ){
      if( clp_.isOptionTriggered( appendixDictEntry.optionName ) ){
        appendixList.emplace_back( appendixDictEntry.appendix );
        if( clp_.getNbValueSet(appendixDictEntry.optionName) > 0 ){

          auto args = clp_.getOptionValList<std::string>(appendixDictEntry.optionName);
          for( auto& arg : args ){
            // strip potential slashes and extensions
            arg = GenericToolbox::getFileName(arg, false);
            if( arg.size() > maxArgLength ){
              // print dotdot if too long
              arg = arg.substr(0, maxArgLength);
              arg[arg.size()-1] = '.';
              arg[arg.size()-2] = '.';
              arg[arg.size()-3] = '.';
            }

            // cleanup from special chars
            arg = GenericToolbox::generateCleanFileName(arg);

            if( not arg.empty() ){
              if( not appendixList.back().empty() ){ appendixList.back() += "_"; }
              appendixList.back() += arg;
            }
          }
        }

      } // option triggered?
    }

    return GenericToolbox::joinVectorString(appendixList, "_");
  }

  bool ObjectReader::quiet{false};
  bool ObjectReader::throwIfNotFound{false};
  bool ObjectReader::readObject( TDirectory* f_, const std::string& objPath_){ return readObject<TObject>(f_, objPath_); }

}
