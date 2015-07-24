#include "Timing.h"
#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"

#include "tree/TDetectorRead.h"

#include <cstdint>

using namespace std;
using namespace ant;
using namespace ant::calibration;

void Timing::ProcessEvent(const Event &)
{

}

void Timing::Finish()
{

}

void Timing::ShowResult()
{

}

Timing::Timing(Detector_t::Type_t detectorType,
        Calibration::Converter::ptr_t converter,
        double defaultOffset,
        const interval<double>& timeWindow, // default {-inf, inf}
        const double defaultGain, // default gain is 1.0
        const std::vector< TKeyValue<double> >& gains
        ) :
    Calibration::Module(
        std_ext::formatter()
        << Detector_t::ToString(detectorType)
        << "_Timing"
           ),
    DetectorType(detectorType),
    Converter(move(converter)),
    TimeWindow(timeWindow),
    DefaultOffset(defaultOffset),
    Offsets(),
    DefaultGain(defaultGain),
    Gains()
{
    if(Converter==nullptr)
        throw std::runtime_error("Given converter should not be nullptr");

    // fill a gain vector from given key-value pairs
    // for faster access (if some are given at all)
    if(gains.empty())
        return;
    unsigned maxkey = 0;
    for(const auto& gain : gains)
        maxkey = gain.Key>maxkey ? gain.Key : maxkey;
    Gains.resize(maxkey+1, DefaultGain);
    for(const auto& gain : gains)
        Gains[gain.Key] = gain.Value;
}

void Timing::ApplyTo(const map< Detector_t::Type_t, list< TDetectorReadHit* > >& hits)
{
    // search for to be calibrated timings
    const auto it_dethits = hits.find(DetectorType);
    if(it_dethits == hits.end())
        return;

    const auto& dethits = it_dethits->second;

    // now calibrate the timings (ignore any other kind of hits)
    for(TDetectorReadHit* dethit : dethits) {
        if(dethit->GetChannelType() != Channel_t::Type_t::Timing)
            continue;

        // the Converter is smart enough to account for reference timings!
        const auto& values = Converter->Convert(dethit->RawData);
        dethit->Values.reserve(values.size());

        // apply gain/offset to each of the values (might be multihit)
        for(double value : values) {
            if(Gains.empty())
                value *= DefaultGain;
            else
                value *= Gains[dethit->Channel];

            if(Offsets.empty())
                value -= DefaultOffset;
            else
                value -= Offsets[dethit->Channel];

            if(!TimeWindow.Contains(value))
                continue;

            dethit->Values.push_back(value);
        }
    }
}

