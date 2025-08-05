// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

// In extraction/LocalPropertyHdf5Output.cc
// suitable data types for HDF5

#ifdef USE_HDF5

#include "extraction/LocalPropertyHdf5Output.h"
#include "net/IOCommunicator.h"
#include <vector>
#include <cassert>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <filesystem>
#include <variant>

// Hdf5 check
#define H5_CHECK(err) (h5_check((err), __FILE__, __LINE__))
inline herr_t h5_check(herr_t err, const char* file, int line) {
    if (err < 0) {
        fprintf(stderr, "HDF5 Error at %s:%d\n", file, line);
        H5Eprint(H5E_DEFAULT, stderr);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    return err;
}

namespace {
    template<typename T>
    void WriteDataset(hid_t group_id, const std::string& name, const std::vector<T>& buffer, 
                      hid_t hdf_type, hsize_t global_rows, hsize_t local_rows, unsigned field_len, 
                      const hemelb::net::IOCommunicator& comms, MPI_Comm mpi_comm) {
        
        hsize_t global_dims[2] = {global_rows, field_len};
        hid_t filespace_id = H5Screate_simple(2, global_dims, NULL);
        H5_CHECK(filespace_id);
        
        hid_t dataset_id = H5Dcreate2(group_id, name.c_str(), hdf_type, filespace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5_CHECK(dataset_id);

        hsize_t local_dims[2] = {local_rows, field_len};
        hid_t memspace_id = H5Screate_simple(2, local_dims, NULL);
        H5_CHECK(memspace_id);

        hsize_t start[2] = {(hsize_t)comms.Scan(local_rows, MPI_SUM) - local_rows, 0};
        hsize_t count[2] = {local_rows, field_len};
        H5_CHECK(H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, start, NULL, count, NULL));
        
        hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
        H5_CHECK(dxpl);
        H5_CHECK(H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE));
        H5_CHECK(H5Dwrite(dataset_id, hdf_type, memspace_id, filespace_id, dxpl, buffer.data()));
        
        H5_CHECK(H5Pclose(dxpl));
        H5_CHECK(H5Sclose(memspace_id));
        H5_CHECK(H5Dclose(dataset_id));
        H5_CHECK(H5Sclose(filespace_id));
    }
    // type mapping for HDF5
    template<typename T>
    hid_t hdf5_type_from_variant();

    template<>
    hid_t hdf5_type_from_variant<float>() { return H5T_NATIVE_FLOAT; }

    template<>
    hid_t hdf5_type_from_variant<double>() { return H5T_NATIVE_DOUBLE; }

    template<>
    hid_t hdf5_type_from_variant<std::int32_t>() { return H5T_NATIVE_INT32; }

    template<>
    hid_t hdf5_type_from_variant<std::uint32_t>() { return H5T_NATIVE_UINT32; }

    template<>
    hid_t hdf5_type_from_variant<std::int64_t>() { return H5T_NATIVE_INT64; }

    template<>
    hid_t hdf5_type_from_variant<std::uint64_t>() { return H5T_NATIVE_UINT64; }
}

namespace hemelb::extraction
{
  LocalPropertyHdf5Output::LocalPropertyHdf5Output(IterableDataSource& dataSource,
                                                     const PropertyOutputFile& outputSpec,
                                                     const net::IOCommunicator& comms)
    : LocalPropertyOutput(dataSource, outputSpec, comms), mpi_comm(comms) {}

  LocalPropertyHdf5Output::~LocalPropertyHdf5Output() {
    if (file_id >= 0) {
      WriteXDMFFile();
      H5Fclose(file_id);
    }
  }

  //  WriteXDMFFile for proper data type and precision
    void LocalPropertyHdf5Output::WriteXDMFFile()
    {
        if (comms.Rank() != 0) return;

        using namespace hemelb::extraction;
        using hemelb::extraction::code::type_to_enum;
        using hemelb::extraction::code::type_to_size;

        std::string h5_filename = outputSpec.filename.string();
        std::string xmf_filename = h5_filename;
        size_t pos = xmf_filename.rfind(".h5");
        if (pos != std::string::npos)
            xmf_filename.replace(pos, 3, ".xmf");
        else
            xmf_filename += ".xmf";

        std::filesystem::path h5_path(h5_filename);
        std::string h5_basename = h5_path.filename().string();

        std::ofstream xmf_file(xmf_filename);
        xmf_file << "<?xml version=\"1.0\" ?>\n";
        xmf_file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        xmf_file << "<Xdmf Version=\"3.0\">\n";
        xmf_file << "  <Domain>\n";
        xmf_file << "    <Grid Name=\"HemeLB Collection\" GridType=\"Collection\" CollectionType=\"Temporal\">\n";

        for (unsigned long ts : written_timesteps)
        {
            xmf_file << "      <Grid Name=\"step_" << ts << "\" GridType=\"Uniform\">\n";
            xmf_file << "        <Time Value=\"" << ts << "\"/>\n";
            xmf_file << "        <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << global_site_count << "\"/>\n";

            // Geometry
            xmf_file << "        <Geometry GeometryType=\"XYZ\">\n";
            xmf_file << "          <DataItem Format=\"HDF\" Dimensions=\"" << global_site_count << " 3\" DataType=\"UInt\" Precision=\"4\">"
                        << h5_basename << ":/step_" << ts << "/geometry</DataItem>\n";
            xmf_file << "        </Geometry>\n";

            // Attributes
            for (const auto& field : outputSpec.fields)
            {
                unsigned field_len = GetFieldLength(field.src);
                std::string attr_type = (field_len > 1) ? "Vector" : "Scalar";

                std::string data_type_str = "Float";
                std::string precision_str = "4";

                // Determine the data type and precision based on the field typecode
                auto tc = code::type_to_enum(field.typecode);
                auto prec = code::type_to_size(field.typecode);

                // Set the data type and precision based on the type code
                switch (tc) {
                case io::formats::extraction::TypeCode::FLOAT:
                case io::formats::extraction::TypeCode::DOUBLE:
                    data_type_str = "Float";
                    break;
                case io::formats::extraction::TypeCode::INT32:
                case io::formats::extraction::TypeCode::INT64:
                    data_type_str = "Int";
                    break;
                case io::formats::extraction::TypeCode::UINT32:
                case io::formats::extraction::TypeCode::UINT64:
                    data_type_str = "UInt";
                    break;
                default:
                    data_type_str = "Unknown";
                }

                precision_str = std::to_string(prec);

                xmf_file << "        <Attribute Name=\"" << field.name << "\" AttributeType=\"" << attr_type << "\" Center=\"Node\">\n";
                xmf_file << "          <DataItem Format=\"HDF\" Dimensions=\"" << global_site_count << " " << field_len
                        << "\" DataType=\"" << data_type_str << "\" Precision=\"" << precision_str << "\">"
                        << h5_basename << ":/step_" << ts << "/" << field.name << "</DataItem>\n";
                xmf_file << "        </Attribute>\n";
            }

            xmf_file << "      </Grid>\n";
        }

        xmf_file << "    </Grid>\n";
        xmf_file << "  </Domain>\n";
        xmf_file << "</Xdmf>\n";
        xmf_file.close();
    }



    // Write
    void LocalPropertyHdf5Output::Write(unsigned long timestepNumber, unsigned long totalSteps)
    {
    if (!ShouldWrite(timestepNumber)) return;

    if (file_id < 0) {
        hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
        H5_CHECK(H5Pset_fapl_mpio(fapl, mpi_comm, MPI_INFO_NULL));
        file_id = H5Fcreate(outputSpec.filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
        H5_CHECK(H5Pclose(fapl));
    }

    written_timesteps.push_back(timestepNumber);

    std::string group_name = "step_" + std::to_string(timestepNumber);
    hid_t group_id = H5Gcreate(file_id, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5_CHECK(group_id);

    std::vector<std::vector<uint32_t>> geometry_data;
    std::map<std::string, std::vector<char>> field_buffers;
    std::map<std::string, hid_t> field_types;
    std::map<std::string, unsigned> field_lengths;

    dataSource.Reset();
    while (dataSource.ReadNext()) {
        if (!outputSpec.geometry->Include(dataSource, dataSource.GetPosition()))
        continue;

        const auto& pos = dataSource.GetPosition();
        geometry_data.push_back({(uint32_t)pos.x(), (uint32_t)pos.y(), (uint32_t)pos.z()});

        for (const auto& field_spec : outputSpec.fields) {
        auto& buffer = field_buffers[field_spec.name];
        auto len = GetFieldLength(field_spec.src);
        field_lengths[field_spec.name] = len;

        std::visit([&](auto tag) {
            using T = decltype(tag);
            T values[9] = {};
            if constexpr (std::is_same_v<T, float>) {
            if (std::holds_alternative<source::Pressure>(field_spec.src)) values[0] = static_cast<T>(dataSource.GetPressure());
            if (std::holds_alternative<source::ShearStress>(field_spec.src)) values[0] = static_cast<T>(dataSource.GetShearStress());
            } else if constexpr (std::is_same_v<T, double>) {
            if (std::holds_alternative<source::Velocity>(field_spec.src)) {
                auto v = dataSource.GetVelocity();
                values[0] = static_cast<T>(v.x());
                values[1] = static_cast<T>(v.y());
                values[2] = static_cast<T>(v.z());
            }
            }
            buffer.insert(buffer.end(), reinterpret_cast<char*>(values), reinterpret_cast<char*>(values + len));
            field_types[field_spec.name] = hdf5_type_from_variant<T>();
        }, field_spec.typecode);
        }
    }

    std::vector<uint32_t> flat_coords;
    for (auto& p : geometry_data)
        flat_coords.insert(flat_coords.end(), p.begin(), p.end());

    hsize_t geo_dims[2] = {global_site_count, 3};
    hid_t geo_space = H5Screate_simple(2, geo_dims, NULL);
    hid_t geo_dset = H5Dcreate2(group_id, "geometry", H5T_NATIVE_UINT32, geo_space,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5_CHECK(geo_dset);
    hsize_t start[2] = {comms.Scan(local_site_count, MPI_SUM) - local_site_count, 0};
    hsize_t count[2] = {local_site_count, 3};
    H5Sselect_hyperslab(geo_space, H5S_SELECT_SET, start, NULL, count, NULL);
    hid_t memspace = H5Screate_simple(2, count, NULL);
    hid_t xfer = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(xfer, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(geo_dset, H5T_NATIVE_UINT32, memspace, geo_space, xfer, flat_coords.data());
    H5Pclose(xfer); H5Sclose(memspace); H5Sclose(geo_space); H5Dclose(geo_dset);

    for (const auto& field_spec : outputSpec.fields) {
        const std::string& name = field_spec.name;
        const auto& raw_buffer = field_buffers[name];
        auto h5type = field_types[name];
        unsigned len = field_lengths[name];
        hsize_t dims[2] = {global_site_count, len};
        hid_t space = H5Screate_simple(2, dims, NULL);
        hid_t dset = H5Dcreate2(group_id, name.c_str(), h5type, space,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hsize_t start[2] = {comms.Scan(local_site_count, MPI_SUM) - local_site_count, 0};
        hsize_t count[2] = {local_site_count, len};
        H5Sselect_hyperslab(space, H5S_SELECT_SET, start, NULL, count, NULL);
        hid_t memspace = H5Screate_simple(2, count, NULL);
        hid_t xfer = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(xfer, H5FD_MPIO_COLLECTIVE);
        H5Dwrite(dset, h5type, memspace, space, xfer, raw_buffer.data());
        H5Pclose(xfer); H5Sclose(memspace); H5Sclose(space); H5Dclose(dset);
    }

    H5_CHECK(H5Gclose(group_id));
    }
}
#endif