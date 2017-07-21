#include "Setup.h"

#include "detectors/Trigger.h"
#include "detectors/CB.h"
#include "detectors/PID.h"
#include "detectors/TAPS.h"
#include "detectors/TAPSVeto.h"
#include "detectors/EPT.h"

namespace ant {
namespace expconfig {
namespace setup {

/**
 * @brief Common base class for all Setups of the 2014 Ent Point Tagger beam times
 */
class Setup_2014_EPT : public Setup
{
protected:
    const bool MCTaggerHits;
    const bool cherenkovInstalled;
    const std::shared_ptr<detector::Trigger_2014> Trigger;
    const std::shared_ptr<detector::EPT_2014> EPT;
    const std::shared_ptr<detector::CB> CB;
    const std::shared_ptr<detector::PID_2014> PID;
    const std::shared_ptr<detector::TAPS_2013_11> TAPS;
    const std::shared_ptr<detector::TAPSVeto_2014> TAPSVeto;

public:

    Setup_2014_EPT(const std::string& name, OptionsPtr opt);

    virtual double GetElectronBeamEnergy() const override;

    void BuildMappings(std::vector<hit_mapping_t>& hit_mappings,
                       std::vector<scaler_mapping_t>& scaler_mappings) const override;

    virtual candidatebuilder_config_t GetCandidateBuilderConfig() const override;

    virtual triggersimu_config_t GetTriggerSimuConfig() const override;

    virtual UnpackerA2GeantConfig::promptrandom_config_t GetPromptRandomConfig() const override;
};

}}} // namespace ant::expconfig::setup
