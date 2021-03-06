#pragma once

#include "base/OptionsList.h"

#include <map>
#include <list>
#include <memory>
#include <functional>

namespace ant {
namespace expconfig {

class Setup;

/**
 * @brief The SetupRegistry class semi-automatically registers Setups
 *
 * \note don't forget to use AUTO_REGISTER_SETUP
 * \note linking order is important to get this registry properly working, see CMakeLists.txt
 */
class SetupRegistry
{
friend class SetupRegistration;

private:
    using Creator = std::function<std::shared_ptr<Setup>(const std::string& name, OptionsPtr)>;
    using setup_creators_t = std::map<std::string, Creator>;
    using setups_t = std::map<std::string, std::shared_ptr<Setup> >;
    setup_creators_t setup_creators;
    setups_t setups;
    OptionsPtr options;

    void RegisterSetup(Creator, std::string);
    static SetupRegistry& get_instance();

    SetupRegistry();
    ~SetupRegistry();
public:
    static std::shared_ptr<Setup> GetSetup(const std::string& name);
    static std::list<std::string> GetNames();
    static void AddSetup(const std::string& name, std::shared_ptr<Setup> setup);
    static void SetSetupOptions(OptionsPtr opt);
    static void Cleanup();
};

/**
 * @brief The SetupRegistration class is instantiated for each setup
 */
class SetupRegistration
{
public:
    SetupRegistration(SetupRegistry::Creator, std::string);
};

template<class T>
std::shared_ptr<Setup> setup_factory(const std::string& name, OptionsPtr opt)
{
    return std::make_shared<T>(name, opt);
}

#define AUTO_REGISTER_SETUP(setup) \
    SetupRegistration _setup_registration_ ## setup(ant::expconfig::setup_factory<setup>, #setup);

}} // namespace ant::expconfig
