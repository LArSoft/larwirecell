#ifndef WIRECELL_QLMATCH_OPFLASH
#define WIRECELL_QLMATCH_OPFLASH

#include "WireCellIface/ITensorSet.h"
#include <memory>
#include <set>
#include <vector>

namespace WireCell::QLMatch {
  class Opflash {
  public:
    typedef std::shared_ptr<Opflash> pointer;
    /// @brief Construct a Opflash object from a 2D ITensor with doubles.
    /// @param ten: ncol need to larger than nchan+1. 0: time, 1-nchan: PE
    /// each row is a flash
    Opflash(const ITensor::pointer ten,
            const int idx,
            const double threshold,
            const int nchan = 32);
    ~Opflash();

    void set_flash_id(int value) { flash_id = value; };
    void set_flash_type(int value) { type = value; };

    int get_flash_id() { return flash_id; };
    double get_time() { return time; };
    double get_total_PE() { return total_PE; };
    const std::vector<double>& get_PEs() const { return PE; };
    double get_PE(int ch) { return PE[ch]; };
    double get_PE_err(int ch) { return PE_err[ch]; };
    bool get_fired(int ch);
    int get_num_fired() { return fired_channels.size(); };
    int get_type() { return type; }
    double get_low_time() { return low_time; };
    double get_high_time() { return high_time; };
    int get_num_channels() const { return m_nchan; };
    double get_threshold() const { return m_threshold; };

    // void swap_channels();

  protected:
    int m_nchan;
    double m_threshold;

    int type;
    int flash_id;
    double low_time;
    double high_time;
    double time;
    double total_PE;

    std::vector<int> fired_channels;
    std::vector<double> PE;
    std::vector<double> PE_err;
  };

  struct OpFlashCompare {
    bool operator()(Opflash* a, Opflash* b) const
    {
      if (a->get_time() < b->get_time()) { return true; }
      else {
        return false;
      }
    }
  };

  typedef std::vector<Opflash*> OpflashSelection;
  typedef std::set<Opflash*, OpFlashCompare> OpFlashSet;

} // namespace WireCell::QLMatch

#endif
