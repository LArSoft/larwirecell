#include "MultiChannelNoiseDB.h"

#include "art/Framework/Principal/Event.h"

#include "WireCellUtil/NamedFactory.h"

#include <unordered_set>

WIRECELL_FACTORY(wclsMultiChannelNoiseDB,
                 wcls::MultiChannelNoiseDB,
                 wcls::IArtEventVisitor,
                 WireCell::IChannelNoiseDatabase,
                 WireCell::IConfigurable)

using namespace WireCell;

wcls::MultiChannelNoiseDB::MultiChannelNoiseDB() {}

wcls::MultiChannelNoiseDB::~MultiChannelNoiseDB() {}

void wcls::MultiChannelNoiseDB::MultiChannelNoiseDB::visit(art::Event& event)
{
  for (auto one : m_rules) {
    if (one.check(event)) {
      m_pimpl = one.chndb;
      // std::cerr << "wclsMultiChannelNoiseDB: using: " << type(*one.chndb)
      //           << " @" << (void*)one.chndb.get()
      //           << " from MultiChannelNoiseDB @" << (void*)this
      //           <<"\n";
      m_pimpl_visitor = one.visitor;
      if (m_pimpl_visitor) { // okay to be nullptr if not a wclsChannelNoiseDB
        m_pimpl_visitor->visit(event);
      }
      return;
    }
  }
  THROW(KeyError() << errmsg{
          "MultiChannelNoiseDB: no matching rule for event, consider 'bool' catch all in config"});
}

struct ReturnBool {
  bool ok;
  ReturnBool(Json::Value jargs) // simple value, true or false
    : ok(jargs.asBool())
  {}
  bool operator()(art::Event& /*event*/) { return ok; }
};

struct RunRange {
  int first, last;
  RunRange(Json::Value jargs)
  { // {first: NNN, last: MMM}
    first = jargs["first"].asInt();
    last = jargs["last"].asInt();
  }
  bool operator()(art::Event& event)
  {
    int run = event.run();
    return first <= run and run <= last;
  }
};

struct RunList {
  std::unordered_set<int> runs;
  RunList(Json::Value jargs)
  { // [1001, 1002, 1003]
    for (auto j : jargs) {
      runs.insert(j.asInt());
    }
  }
  bool operator()(art::Event& event)
  {
    int run = event.run();
    return runs.find(run) != runs.end();
  }
};

struct RunStarting {
  int run;
  RunStarting(Json::Value jargs)
  { // 1001
    run = jargs.asInt();
  }
  bool operator()(art::Event& event)
  {
    int thisrun = event.run();
    return thisrun >= run;
  }
};
struct RunBefore {
  int run;
  RunBefore(Json::Value jargs)
  { // 1001
    run = jargs.asInt();
  }
  bool operator()(art::Event& event)
  {
    int thisrun = event.run();
    return thisrun < run;
  }
};

void wcls::MultiChannelNoiseDB::MultiChannelNoiseDB::configure(const WireCell::Configuration& cfg)
{
  for (auto jrule : cfg["rules"]) {
    auto rule = jrule["rule"].asString();
    auto tn = jrule["chndb"].asString();
    std::cerr << "\tMultiChannelNoiseDB: " << tn << " using rule: " << rule << std::endl;
    auto chndb = Factory::find_tn<IChannelNoiseDatabase>(tn);
    if (!chndb) {
      THROW(KeyError() << errmsg{"Failed to find (sub) channel noise DB object: " + tn});
    }
    IArtEventVisitor::pointer visitor = Factory::find_maybe_tn<IArtEventVisitor>(tn);
    auto jargs = jrule["args"];
    if (rule == "runrange") {
      m_rules.push_back(SubDB(RunRange(jargs), chndb, visitor));
      continue;
    }
    if (rule == "runlist") {
      m_rules.push_back(SubDB(RunList(jargs), chndb, visitor));
      continue;
    }
    if (rule == "runstarting") {
      m_rules.push_back(SubDB(RunStarting(jargs), chndb, visitor));
      continue;
    }
    if (rule == "runbefore") {
      m_rules.push_back(SubDB(RunBefore(jargs), chndb, visitor));
      continue;
    }

    if (rule == "bool") {
      m_rules.push_back(SubDB(ReturnBool(jargs), chndb, visitor));
      continue;
    }
    THROW(KeyError() << errmsg{"Unknown multi channel noise DB rule: " + rule});
  }
}

/* example:
 rules: [
   {
     rule: "runlist",
     chndb: "wclsChannelNoiseDB:prehwfix",
     args: [1001, 1003, 1005],
   },
   {
     rule: "runrange",
     chndb: "wclsChannelNoiseDB:",
     args: { first: 1001, second: 1003},
   },
   {  // provide a catch all or risk an runtime exception
     rule: "bool",
     chndb: "OmniChannelNoiseDB:",
     args: true
   },
   // also "runstarting" and "runbefore"
 ]
*/

WireCell::Configuration wcls::MultiChannelNoiseDB::MultiChannelNoiseDB::default_configuration()
  const
{
  Configuration cfg;
  cfg["rules"] = Json::arrayValue;
  return cfg;
}
