#include "Tagger.h"

#include "detail/Tagger_2007_1508_elements.h"
#include "detail/Tagger_2010_03_450_elements.h"
#include "detail/Tagger_2015_450_elements.h"

#include "tree/TID.h"

#include <limits>
#include <cassert>


using namespace std;
using namespace ant;
using namespace ant::expconfig::detector;

const std::string Tagger::ScalerName = "Tagger_Scalers";

void Tagger::BuildMappings(
        vector<UnpackerAcquConfig::hit_mapping_t> &hit_mappings,
        vector<UnpackerAcquConfig::scaler_mapping_t>& scaler_mappings) const
{

    vector<UnpackerAcquConfig::scaler_mapping_t::entry_t> scaler_entries;

    for(const Element_t& element : elements) {

        // TDC/scaler information is most important
        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::Timing,
                                  element.Channel,
                                  element.TDC
                                  );

        // build the scaler entries
        scaler_entries.emplace_back(element.Channel,
                                    element.Scaler);

        // ADC information are rarely present for the EPT
        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::Integral,
                                  element.Channel,
                                  element.ADC
                                  );
    }

    // map the scalers
    scaler_mappings.emplace_back(Tagger::ScalerName, scaler_entries);
}



