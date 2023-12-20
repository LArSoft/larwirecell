#ifndef PRIVATE_SEDDUMPER_H
#define PRIVATE_SEDDUMPER_H

#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include <algorithm>
#include <fstream>
#include <string>

// Get SimEnergyDeposit at label from event and dump to filename.
inline void sed_dumper(const art::Event& event,
                       const std::string& label,
                       const std::string& filename,
                       const std::string& prefix = "")
{
  if (label.empty() or filename.empty()) { return; }

  art::Handle<std::vector<sim::SimEnergyDeposit>> sedvh;
  const bool okay = event.getByLabel(label, sedvh);
  if (not okay) { return; }

  std::fstream fstr(filename, std::fstream::out | std::fstream::app);
  const size_t ndepos = sedvh->size();
  for (size_t index = 0; index < ndepos; ++index) {
    auto const& sed = sedvh->at(index);
    fstr << prefix << label << " sed"
         << " " << ndepos << " " << index << " " << event.run() << " " << event.subRun() << " "
         << event.event() << " " << sed.TrackID() << " " << sed.OrigTrackID() << " "
         << sed.PdgCode() << " " << sed.T() << " " << sed.X() << " " << sed.Y() << " " << sed.Z()
         << " " << sed.E() << " " << sed.NumElectrons() << " " << sed.NumPhotons() << "\n";
  }
}

inline int depoid(const WireCell::IDepo::pointer& depo)
{
  auto p = depo->prior();
  if (p) { return p->id(); }
  return depo->id();
}

inline bool idorder(const WireCell::IDepo::pointer& a, const WireCell::IDepo::pointer& b)
{
  return depoid(a) < depoid(b);
}

inline void depo_dumper(WireCell::IDepo::vector depos,
                        const std::string& label,
                        const std::string& filename,
                        const std::string& prefix = "")
{
  if (label.empty() or filename.empty()) { return; }

  std::sort(depos.begin(), depos.end(), idorder);

  auto ndepos = depos.size();
  std::fstream fstr(filename, std::fstream::out | std::fstream::app);
  for (size_t index = 0; index < ndepos; ++index) {
    auto const& child = depos.at(index);
    auto parent = child->prior();

    std::vector<WireCell::IDepo::pointer> two = {child, parent};

    for (size_t ind = 0; ind < 2; ++ind) {
      const auto& depo = two[ind];
      if (!depo) { continue; }

      auto const& pos = depo->pos();
      fstr << prefix << label << " depo"
           << " " << ind << " " << ndepos << " " << index << " " << depoid(depo) << " "
           << depo->id() << " " << depo->pdg() << " " << depo->charge() << " " << depo->energy()
           << " " << depo->time() << " " << pos.x() << " " << pos.y() << " " << pos.z() << " "
           << depo->extent_long() << " " << depo->extent_tran() << "\n";
    }
  }
}

#endif
