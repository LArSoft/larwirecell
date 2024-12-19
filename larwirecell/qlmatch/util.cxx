#include "util.h"

#include "WireCellUtil/Exceptions.h"
#include "WireCellImg/PointCloudFacade.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"


#include <fstream>

using namespace WireCell;
using namespace WireCell::PointCloud::Facade;
using namespace WireCell::PointCloud::Tree;

void WireCell::QLMatch::dump_bee_3d(const Points::node_t& root, const std::string& fn)
{
    using spdlog::debug;
    using WireCell::PointCloud::Facade::float_t;
    using WireCell::PointCloud::Facade::int_t;

    Json::Value bee;
    bee["runNo"] = 0;
    bee["subRunNo"] = 0;
    bee["eventNo"] = 0;
    bee["geom"] = "uboone";
    bee["type"] = "cluster";

    std::vector<float_t> x;
    std::vector<float_t> y;
    std::vector<float_t> z;
    std::vector<float_t> q;
    std::vector<int_t> cluster_id;
    int_t cid = 0;
    for (const auto cnode : root.children()) {  // this is a loop through all clusters ...
        Scope scope = {"3d", {"x", "y", "z"}};
        const auto& sv = cnode->value.scoped_view(scope);

        const auto& spcs = sv.pcs();  // spcs 'contains' all blobs in this cluster ...
        // int npoints = 0;

        for (const auto& spc : spcs) {  // each little 3D pc --> (blobs)   spc represents x,y,z in a blob
            if (spc.get().get("x") == nullptr) {
                debug("No x in point cloud, skip");
                continue;
            }

            // assume others exist
            const auto& x_ = spc.get().get("x")->elements<float_t>();
            const auto& y_ = spc.get().get("y")->elements<float_t>();
            const auto& z_ = spc.get().get("z")->elements<float_t>();
            const size_t n = x_.size();
            x.insert(x.end(), x_.begin(), x_.end());  // Append x_ to x
            y.insert(y.end(), y_.begin(), y_.end());
            z.insert(z.end(), z_.begin(), z_.end());
            q.insert(q.end(), n, 1.0);
            cluster_id.insert(cluster_id.end(), n, cid);
            // npoints += n;
        }

        // spc.kd() // kdtree ...
        // const auto& skd = sv.kd();
        // std::cout << "xin6: " << sv.npoints() << " " << npoints << " " << spcs.size() << " " <<
        // skd.points().size() << std::endl;

        ++cid;
    }

    Json::Value json_x(Json::arrayValue);
    for (const auto& val : x) {
        json_x.append(val / units::cm);
    }
    bee["x"] = json_x;

    Json::Value json_y(Json::arrayValue);
    for (const auto& val : y) {
        json_y.append(val / units::cm);
    }
    bee["y"] = json_y;

    Json::Value json_z(Json::arrayValue);
    for (const auto& val : z) {
        json_z.append(val / units::cm);
    }
    bee["z"] = json_z;

    Json::Value json_q(Json::arrayValue);
    for (const auto& val : q) {
        json_q.append(val);
    }
    bee["q"] = json_q;

    Json::Value json_cluster_id(Json::arrayValue);
    for (const auto& val : cluster_id) {
        json_cluster_id.append(val);
    }
    bee["cluster_id"] = json_cluster_id;

    // Write cfg to file
    std::ofstream file(fn);
    if (file.is_open()) {
        Json::StreamWriterBuilder writer;
        writer["indentation"] = "    ";
        writer["precision"] = 6; // significant digits
        std::unique_ptr<Json::StreamWriter> jsonWriter(writer.newStreamWriter());
        jsonWriter->write(bee, &file);
        file.close();
    }
    else {
        raise<ValueError>("Failed to open file: " + fn);
    }
}

#include <json/json.h>

void WireCell::QLMatch::dump_bee_flash(const ITensorSet::pointer& ts, const std::string& fn)
{
    using spdlog::debug;
    const auto& tens = ts->tensors();
    if (tens->size() != 1) {
        raise<ValueError>("Expected 1 tensor, got %d", tens->size());
    }
    const auto& ten = tens->at(0);
    typedef boost::multi_array<double, 2> MultiArray;
    boost::array<MultiArray::index, 2> shape = {(int)ten->shape()[0], (int)ten->shape()[1]};
    boost::multi_array_ref<double, 2> mar((double*)ten->data(), shape);

    Json::Value data;
    data["runNo"] = 0;
    data["subRunNo"] = 0;
    data["eventNo"] = 0;
    data["geom"] = "sbnd";
    data["op_t"] = Json::Value(Json::arrayValue);
    data["op_pes"] = Json::Value(Json::arrayValue);
    data["op_pes_pred"] = Json::Value(Json::arrayValue);
    data["op_peTotal"] = Json::Value(Json::arrayValue);
    data["cluster_id"] = Json::Value(Json::arrayValue);
    data["op_nomatching_cluster_ids"] = Json::Value(Json::arrayValue);

    debug("shape: {} {}", shape[0], shape[1]);
    for (size_t i = 0; i < (size_t)shape[0]; ++i) {
        data["op_t"].append(mar[i][0]);
        data["cluster_id"].append(i);
        double op_peTotal = 0;
        auto op_pes = Json::Value(Json::arrayValue);
        for (size_t j = 1; j < (size_t)shape[1]; ++j) {
            op_peTotal += mar[i][j];
            op_pes.append(mar[i][j]);
            // std::cout << mar[i][j] << " ";
        }
        data["op_peTotal"].append(op_peTotal);
        data["op_pes"].append(op_pes);
        data["op_pes_pred"].append(op_pes);
        // std::cout << std::endl;
    }

    // Write cfg to file
    std::ofstream file(fn);
    if (file.is_open()) {
        Json::StreamWriterBuilder writer;
        writer["indentation"] = "    ";
        writer["precision"] = 6; // significant digits
        std::unique_ptr<Json::StreamWriter> jsonWriter(writer.newStreamWriter());
        jsonWriter->write(data, &file);
        file.close();
    }
    else {
        raise<ValueError>("Failed to open file: " + fn);
    }
}