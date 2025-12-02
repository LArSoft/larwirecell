#include "Opflash.h"

#include "WireCellAux/SimpleTensor.h"
#include "WireCellUtil/Exceptions.h"

WireCell::QLMatch::Opflash::Opflash(const ITensor::pointer ten,
                                    const int idx,
                                    const double threshold,
                                    const int nchan)
{
  if (ten->shape().size() != 2) {
    raise<ValueError>("input tensor dim %d != 2", ten->shape().size());
  }
  const int nrow = ten->shape()[0];
  const int ncol = ten->shape()[1];
  if (nrow < idx + 1) { raise<ValueError>("input tensor nrow %d < idx+1", nrow); }
  if (ncol < nchan + 1) { raise<ValueError>("input tensor ncol %d < nchan+1", ncol); }

  flash_id = idx;
  m_nchan = nchan;
  m_threshold = threshold;

  /// FIXME: understand and implement the following
  type = 0;
  low_time = 0;
  high_time = 0;

  typedef boost::multi_array<double, 2> MultiArray;
  boost::array<MultiArray::index, 2> shape = {nrow, ncol};
  boost::multi_array_ref<double, 2> ten_data((double*)ten->data(), shape);

  time = ten_data[idx][0];
  total_PE = 0;
  PE.resize(nchan, 0);
  PE_err.resize(nchan, 1);

  for (int i = 0; i < nchan; ++i) {
    PE[i] = ten_data[idx][i + 1];
    if (PE[i]<1){
      PE_err[i] = 0.3;
    }
    else{
      PE_err[i] = 0.3*PE[i]; // pow(0.1*PE[i],2);
    }
    total_PE += PE[i];
    if (PE[i] > threshold) { fired_channels.push_back(i); }
  }
}

WireCell::QLMatch::Opflash::~Opflash() {}

bool WireCell::QLMatch::Opflash::get_fired(int ch)
{
  if (std::find(fired_channels.begin(), fired_channels.end(), ch) == fired_channels.end()) {
    return false;
  }
  else {
    return true;
  }
}
