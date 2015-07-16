#ifndef ANT_TEVENT_H
#define ANT_TEVENT_H

#include "TDataRecord.h"
#include "TTrack.h"
#include "TCluster.h"
#include "TTaggerHit.h"

#ifndef __CINT__
#include <iomanip>
#include <sstream>
#endif

namespace ant {

struct TEvent : TDataRecord
{
  TEvent() : TDataRecord() {}
  TEvent(const TID& id) :
    TDataRecord(id)
  {}
  virtual ~TEvent() {}

  std::vector<ant::TTrack> Tracks;
  std::vector<ant::TTaggerHit> TaggerHits;

#ifndef __CINT__
  virtual std::ostream& Print( std::ostream& s) const override {
    s << "TEvent:";
    s << " " << TaggerHits.size() << " Taggerhits:\n";
    for(auto& th : TaggerHits) {
        s << "  " << th << "\n";
    }
    s << " " << Tracks.size() << " Tracks\n";
    for(auto& t : Tracks) {
        s << "  " << t << "\n";
        for(auto& c : t.Clusters) {
            s << "   " << c << "\n";
            for(auto& h : c.Hits) {
                s << "    " << h << "\n";
            }
        }
    }
    return s;
  }
#endif

  ClassDef(TEvent, ANT_UNPACKER_ROOT_VERSION)

};

}

#endif // ANT_TEVENT_H
