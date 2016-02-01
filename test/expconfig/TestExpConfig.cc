#include "catch.hpp"
#include "catch_config.h"

#include "unpacker/Unpacker.h"
#include "expconfig/ExpConfig.h"

#include "expconfig/detectors/CB.h"

#include <iostream>

using namespace ant;
using namespace std;

void getdetector();
void getlastfound();
void getall();

TEST_CASE("ExpConfig Get (all)", "[expconfig]") {
    ExpConfig::Setup::Cleanup();
    getall();
}

TEST_CASE("ExpConfig GetLastFound", "[expconfig]") {
    ExpConfig::Setup::Cleanup();
    getlastfound();
}

TEST_CASE("ExpConfig GetDetector", "[expconfig]") {
    ExpConfig::Setup::Cleanup();
    getdetector();
}

void getall() {
    auto setupnames = ExpConfig::Setup::GetNames();
    for(auto setupname : setupnames) {
        REQUIRE(ExpConfig::Setup::Get(setupname) != nullptr);
    }
}

void getdetector() {
    ExpConfig::Setup::ManualName = "Setup_Test";
    REQUIRE_NOTHROW(Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-big.dat.xz"));

    auto tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    REQUIRE(tagger.get() != nullptr);
    auto tagger_fromtype = ExpConfig::Setup::GetDetector(Detector_t::Type_t::EPT);
    REQUIRE(tagger == tagger_fromtype);

    auto cb = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();
    REQUIRE(cb.get() != nullptr);
    REQUIRE(cb->GetNChannels() == 720);

    REQUIRE_THROWS_AS(ExpConfig::Setup::GetDetector(Detector_t::Type_t::Tagger), ExpConfig::Exception);
}

void getlastfound() {
    REQUIRE(ExpConfig::Setup::GetLastFound().get() == nullptr);
    ExpConfig::Setup::ManualName = "Setup_Test";
    REQUIRE(ExpConfig::Setup::GetLastFound().get() != nullptr);
}
