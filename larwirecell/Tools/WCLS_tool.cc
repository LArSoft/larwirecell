#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"

#include "larwirecell/Interfaces/IArtEventVisitor.h"
#include "larwirecell/Interfaces/MainTool.h"

#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetWalker.h"
#include "fhiclcpp/coding.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/OptionalDelegatedParameter.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

#include "WireCellApps/Main.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/String.h"

#include "WireCellUtil/NamedFactory.h"

#include <cstdio>
#include <map>
#include <string>
#include <vector>

namespace {

  // Render the non-atom (sequence or table) top-level values of a
  // ParameterSet as JSON text, which is valid Jsonnet code.  Atoms
  // inside them become JSON strings, numbers, booleans or null.
  class JsonnetCoder : public fhicl::ParameterSetWalker {
  public:
    // Top-level key to JSON text.
    std::map<std::string, std::string> codes;

  private:
    struct level_t {
      bool is_table;
      bool first;
    };
    std::vector<level_t> m_stack;
    std::string m_key, m_code;

    void begin_value(key_t const& key)
    {
      if (m_stack.empty()) {
        m_key = key;
        m_code.clear();
        return;
      }
      auto& lvl = m_stack.back();
      if (!lvl.first) { m_code += ","; }
      lvl.first = false;
      if (lvl.is_table) { m_code += quote(key) + ":"; }
    }
    void end_value()
    {
      if (m_stack.empty()) { codes[m_key] = m_code; }
    }

    void enter_table(key_t const& key, any_t const&) override
    {
      begin_value(key);
      m_code += "{";
      m_stack.push_back({true, true});
    }
    void exit_table(key_t const&, any_t const&) override
    {
      m_stack.pop_back();
      m_code += "}";
      end_value();
    }
    void enter_sequence(key_t const& key, any_t const&) override
    {
      begin_value(key);
      m_code += "[";
      m_stack.push_back({false, true});
    }
    void exit_sequence(key_t const&, any_t const&) override
    {
      m_stack.pop_back();
      m_code += "]";
      end_value();
    }
    void atom(key_t const& key, any_t const& value) override
    {
      if (m_stack.empty()) { return; } // top-level atoms are taken as Jsonnet code
      begin_value(key);
      m_code += atom_code(key, value);
    }

    // FHiCL holds atoms in canonical form: strings double-quoted,
    // numbers normalized (no leading '+', hex/octal decoded).
    static std::string atom_code(key_t const& key, any_t const& value)
    {
      if (fhicl::detail::is_nil(value)) { return "null"; }
      auto const& text = std::any_cast<std::string const&>(value);
      if (!text.empty() && text[0] == '"') {
        std::string str;
        fhicl::detail::decode(value, str);
        return quote(str);
      }
      if (text.find("infinity") != std::string::npos || (!text.empty() && text[0] == '(')) {
        throw cet::exception("WireCellLArSoft")
          << "structs value \"" << key << "\" = " << text << " has no Jsonnet form\n";
      }
      return text; // number, true or false
    }

    static std::string quote(std::string const& str)
    {
      std::string ret = "\"";
      for (unsigned char c : str) {
        switch (c) {
        case '"': ret += "\\\""; break;
        case '\\': ret += "\\\\"; break;
        case '\n': ret += "\\n"; break;
        case '\t': ret += "\\t"; break;
        case '\r': ret += "\\r"; break;
        default:
          if (c < 0x20) {
            char buf[8];
            snprintf(buf, sizeof buf, "\\u%04x", c);
            ret += buf;
          }
          else {
            ret += c;
          }
        }
      }
      return ret + "\"";
    }
  };
}

namespace wcls {

  // https://cdcvs.fnal.gov/redmine/projects/fhicl-cpp/wiki/Fhiclcpp_types_in_detail#TableltT-KeysToIgnoregt
  struct WCLSKeysToIgnore {
    std::set<std::string> operator()()
    {
      // Ignore these for validation.
      return {"params"};
    }
  };

  // https://cdcvs.fnal.gov/redmine/projects/art/wiki/Configuration_validation_and_description
  struct WCLSConfig {
    typedef fhicl::Sequence<std::string> string_list_t;
    typedef fhicl::OptionalSequence<std::string> optional_string_list_t;
    //typedef fhicl::OptionalTable<fhicl::ParameterSet> generic_pset_t;
    typedef fhicl::OptionalDelegatedParameter generic_pset_t;

    // These are WCT config items to pass through
    string_list_t configs{fhicl::Name("configs"),
                          fhicl::Comment("List of one or more WCT configuration files."
                                         "\nThey are located w.r.t. the WCT load path.")};
    string_list_t apps{fhicl::Name("apps"),
                       fhicl::Comment("List of one or more WCT application objects to execute.")};
    string_list_t plugins{fhicl::Name("plugins"),
                          fhicl::Comment("List of WCT component plugin libraries to load.\n"
                                         "They are located w.r.t. the OS library load path")};

    optional_string_list_t paths{
      fhicl::Name("paths"),
      fhicl::Comment("Optional list of file system paths to add to the WCT "
                     "configuration file load path."
                     "\nThis augments the WIRECELL_PATH environment variable.")};
    generic_pset_t params{
      fhicl::Name("params"),
      fhicl::Comment("Optional table giving external variables to inject into WCT configuration.")};
    generic_pset_t structs{
      fhicl::Name("structs"),
      fhicl::Comment(
        "Optional table giving external Jsonnet code to inject into WCT configuration.\n"
        "A string value is Jsonnet code; a sequence or table value is converted to JSON.")};

    // These are items needed by the tool
    optional_string_list_t inputers{
      fhicl::Name("inputers"),
      fhicl::Comment("List of WCT components which act as WCT sources.\n"
                     "They are called before WCT executes on each Art Event object")};
    optional_string_list_t outputers{
      fhicl::Name("outputers"),
      fhicl::Comment("List of WCT components which act as WCT sinks.\n"
                     "They are called after WCT executes on each Art Event object.")};

    optional_string_list_t logsinks{
      fhicl::Name("logsinks"),
      fhicl::Comment("List of WCT log sinks.\n"
                     "Eg the strings 'stdout', 'stderr' or a file name.\n"
                     "An optional log level may be appended with ':<level>'.")};
    optional_string_list_t loglevels{
      fhicl::Name("loglevels"),
      fhicl::Comment("List of minimum WCT logger levels.\n"
                     "Specify as '<logger>:<level>' or as just '<level>' for default.")};
  };

  class WCLS : public MainTool {
  public:
    using Parameters = art::ToolConfigTable<WCLSConfig, WCLSKeysToIgnore>;

    explicit WCLS(Parameters const& ps);
    virtual ~WCLS() {}

    void produces(art::ProducesCollector& collector)
    {
      for (auto iaev : m_outputers) {
        iaev->produces(collector);
      }
    }
    void process(art::Event& event);

  private:
    WireCell::Main m_wcmain;
    wcls::IArtEventVisitor::vector m_inputers, m_outputers;
    // for c2: m_prod is not used
    // art::EDProducer* m_prod;
  };
}

wcls::WCLS::WCLS(wcls::WCLS::Parameters const& params) : m_wcmain()
{
  const auto& wclscfg = params();
  WCLSConfig::optional_string_list_t::value_type slist;

  if (wclscfg.logsinks(slist)) {
    for (auto logsink : slist) {
      //std::cerr << "Log sink: \"" << logsink << "\"\n";
      auto ls = WireCell::String::split(logsink, ":");
      if (ls.size() == 2) { m_wcmain.add_logsink(ls[0], ls[1]); }
      else {
        m_wcmain.add_logsink(ls[0]);
      }
    }
  }
  slist.clear();
  if (wclscfg.loglevels(slist)) {
    for (auto loglevel : slist) {
      //std::cerr << "Log level: \"" << loglevel << "\"\n";
      auto ll = WireCell::String::split(loglevel, ":");
      if (ll.size() == 2) { m_wcmain.set_loglevel(ll[0], ll[1]); }
      else {
        m_wcmain.set_loglevel("", ll[0]);
      }
    }
  }
  slist.clear();
  WireCell::Log::set_pattern("[%H:%M:%S.%03e] %L [%^%=8n%$] %v");

  // transfer configuration

  // required

  for (auto cfg : wclscfg.configs()) {
    m_wcmain.add_config(cfg);
  }

  for (auto app : wclscfg.apps()) {
    m_wcmain.add_app(app);
  }

  for (auto plugin : wclscfg.plugins()) {
    m_wcmain.add_plugin(plugin);
  }

  // optional

  if (wclscfg.paths(slist)) {
    for (auto path : slist) {
      m_wcmain.add_path(path);
    }
  }
  slist.clear();

  {
    fhicl::ParameterSet wcps;
    if (wclscfg.params.get_if_present(wcps)) {
      for (auto key : wcps.get_names()) {
        auto value = wcps.get<std::string>(key);
        m_wcmain.add_var(key, value);
      }
    }
  }
  {
    fhicl::ParameterSet wcps;
    if (wclscfg.structs.get_if_present(wcps)) {
      // An atom value is taken as Jsonnet code text.  A sequence or
      // table value is converted to the equivalent JSON.
      JsonnetCoder coder;
      wcps.walk(coder);
      for (auto key : wcps.get_names()) {
        auto it = coder.codes.find(key);
        auto value = it == coder.codes.end() ? wcps.get<std::string>(key) : it->second;
        m_wcmain.add_code(key, value);
      }
    }
  }

  //std::cerr << "Initialize Wire Cell\n";
  try {
    m_wcmain.initialize();
  }
  catch (WireCell::Exception& e) {
    std::cerr << "Wire Cell Toolkit threw an exception\n";
    auto msg = errstr(e);
    std::cerr << msg << std::endl;
    throw cet::exception("WireCellLArSoft") << msg;
  }

  if (wclscfg.inputers(slist)) {
    for (auto inputer : slist) {
      auto iaev = WireCell::Factory::find_tn<IArtEventVisitor>(inputer);
      m_inputers.push_back(iaev);
      std::cerr << "Inputer: \"" << inputer << "\"\n";
    }
  }
  slist.clear();
  if (wclscfg.outputers(slist)) {
    for (auto outputer : slist) {
      auto iaev = WireCell::Factory::find_tn<IArtEventVisitor>(outputer);
      m_outputers.push_back(iaev);
      std::cerr << "Outputer: \"" << outputer << "\"\n";
    }
  }
  slist.clear();
}

void wcls::WCLS::process(art::Event& event)
{
  for (auto iaev : m_inputers) {
    //std::cerr << "pre visit\n";
    iaev->visit(event);
  }

  //std::cerr << "Running Wire Cell Toolkit...\n";
  m_wcmain();
  //std::cerr << "... Wire Cell Toolkit done\n";

  for (auto iaev : m_outputers) {
    //std::cerr << "post visit\n";
    iaev->visit(event);
  }
}

DEFINE_ART_CLASS_TOOL(wcls::WCLS)

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
