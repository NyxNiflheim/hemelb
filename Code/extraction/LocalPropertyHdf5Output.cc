// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

// LocalPropertyHdf5Output.cc
// In extraction/LocalPropertyHdf5Output.cc
#ifdef USE_HDF5

#include "extraction/LocalPropertyHdf5Output.h"
#include "net/IOCommunicator.h"
#include <vector>
#include <cassert>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <filesystem>

// check HDF5 error handling
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
    // formatting helper function
    template <typename... Ts>
    std::string safe_fmt(std::string const& pattern, Ts... args) {
        int sz = std::snprintf(nullptr, 0, pattern.c_str(), args...);
        if (sz < 0) throw std::runtime_error("Formatting error");
        std::string ans(sz + 1, '\0');
        std::snprintf(ans.data(), ans.size(), pattern.c_str(), args...);
        ans.pop_back();
        return ans;
    }
    
    // writes a dataset to the HDF5 file
    void WriteDataset(hid_t group_id, const std::string& name, const std::vector<double>& buffer, 
                      hsize_t global_rows, hsize_t local_rows, unsigned field_len, 
                      const hemelb::net::IOCommunicator& comms, MPI_Comm mpi_comm) {
        
        hsize_t global_dims[2] = {global_rows, field_len};
        hid_t filespace_id = H5Screate_simple(2, global_dims, NULL);
        H5_CHECK(filespace_id);
        
        hid_t dataset_id = H5Dcreate2(group_id, name.c_str(), H5T_NATIVE_DOUBLE, filespace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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
        H5_CHECK(H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id, dxpl, buffer.data()));
        
        H5_CHECK(H5Pclose(dxpl));
        H5_CHECK(H5Sclose(memspace_id));
        H5_CHECK(H5Dclose(dataset_id));
        H5_CHECK(H5Sclose(filespace_id));
    }
}

namespace hemelb::extraction
{
  LocalPropertyHdf5Output::LocalPropertyHdf5Output(IterableDataSource& dataSource,
                                                     const PropertyOutputFile& outputSpec,
                                                     const net::IOCommunicator& comms)
    : LocalPropertyOutput(dataSource, outputSpec, comms), mpi_comm(comms)
  {
      std::string_view p = this->outputSpec.filename.native();
      auto i_pcd = p.rfind(".h5");
      auto beginning = (i_pcd != std::string_view::npos) ? p.substr(0, i_pcd) : p;
      
      output_file_pattern += beginning;
      output_file_pattern += "_%06ld.h5";
  }

  void LocalPropertyHdf5Output::WriteXDMFFile(const std::string& h5_filename, unsigned long timestep)
  {
    if (comms.Rank() != 0) return;
    
    std::string xmf_filename = h5_filename;
    size_t pos = xmf_filename.rfind(".h5");
    if (pos != std::string::npos) {
        xmf_filename.replace(pos, 3, ".xmf");
    } else {
        xmf_filename += ".xmf";
    }

    std::filesystem::path h5_path(h5_filename);
    std::string h5_basename = h5_path.filename().string();
    
    std::ofstream xmf_file(xmf_filename);
    xmf_file << "<?xml version=\"1.0\" ?>\n";
    xmf_file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    xmf_file << "<Xdmf Version=\"3.0\">\n";
    xmf_file << "  <Domain>\n";
    xmf_file << "    <Grid Name=\"HemeLB Grid\" GridType=\"Uniform\">\n";
    xmf_file << "      <Time Value=\"" << timestep << "\"/>\n";
    xmf_file << "      <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << global_site_count << "\"/>\n";
    xmf_file << "      <Geometry GeometryType=\"XYZ\">\n";
    xmf_file << "        <DataItem Format=\"HDF\" Dimensions=\"" << global_site_count << " 3\" DataType=\"Float\" Precision=\"8\">" << h5_basename << ":/geometry</DataItem>\n";
    xmf_file << "      </Geometry>\n";

    for (const auto& field : outputSpec.fields) {
        unsigned field_len = GetFieldLength(field.src);
        std::string type = (field_len > 1) ? "Vector" : "Scalar";
        xmf_file << "      <Attribute Name=\"" << field.name << "\" AttributeType=\"" << type << "\" Center=\"Node\">\n";
        xmf_file << "        <DataItem Format=\"HDF\" Dimensions=\"" << global_site_count << " " << field_len << "\" DataType=\"Float\" Precision=\"8\">" << h5_basename << ":/" << field.name << "</DataItem>\n";
        xmf_file << "      </Attribute>\n";
    }
    xmf_file << "    </Grid>\n";
    xmf_file << "  </Domain>\n";
    xmf_file << "</Xdmf>\n";
    xmf_file.close();
  }

  void LocalPropertyHdf5Output::Write(unsigned long timestepNumber, unsigned long totalSteps)
  {
    if (!ShouldWrite(timestepNumber)) return;

    std::string current_h5_filename = safe_fmt(output_file_pattern, 6, timestepNumber);
    
    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5_CHECK(fapl);
    H5_CHECK(H5Pset_fapl_mpio(fapl, mpi_comm, MPI_INFO_NULL));
    hid_t file_id = H5Fcreate(current_h5_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
    H5_CHECK(H5Pclose(fapl));
    
    std::vector<double> coords_buffer;
    coords_buffer.reserve(local_site_count * 3);
    
    std::vector<std::vector<double>> field_buffers(outputSpec.fields.size());
    for(size_t i = 0; i < outputSpec.fields.size(); ++i) {
        field_buffers[i].reserve(local_site_count * GetFieldLength(outputSpec.fields[i].src));
    }

    // traverse the data only once
    dataSource.Reset();
    while (dataSource.ReadNext()) {
      if (outputSpec.geometry->Include(dataSource, dataSource.GetPosition())) {
        // Get position
        const auto& pos = dataSource.GetPosition();
        coords_buffer.push_back(static_cast<double>(pos.x()));
        coords_buffer.push_back(static_cast<double>(pos.y()));
        coords_buffer.push_back(static_cast<double>(pos.z()));

        // Get and check vadidity of all fields
        for (size_t i = 0; i < outputSpec.fields.size(); ++i) {
          std::visit([&](auto&& src_type) {
             if constexpr (std::is_same_v<std::decay_t<decltype(src_type)>, source::Pressure>) {
                double val = dataSource.GetPressure();
                field_buffers[i].push_back(std::isinf(val) ? -1.0 : val);       // store in buffer
             } else if constexpr (std::is_same_v<std::decay_t<decltype(src_type)>, source::Velocity>) {
                const auto& vel = dataSource.GetVelocity();
                field_buffers[i].push_back(std::isinf(vel.x()) ? 0.0 : vel.x());
                field_buffers[i].push_back(std::isinf(vel.y()) ? 0.0 : vel.y());
                field_buffers[i].push_back(std::isinf(vel.z()) ? 0.0 : vel.z());
             } else if constexpr (std::is_same_v<std::decay_t<decltype(src_type)>, source::ShearStress>) {
                double val = dataSource.GetShearStress();
                field_buffers[i].push_back(std::isinf(val) ? -1.0 : val);
             }
          }, outputSpec.fields[i].src);
        }
      }
    }
    // Write dataset once all collected
    WriteDataset(file_id, "geometry", coords_buffer, global_site_count, local_site_count, 3, comms, mpi_comm);

    for (size_t i = 0; i < outputSpec.fields.size(); ++i) {
        const auto& field_spec = outputSpec.fields[i];
        WriteDataset(file_id, field_spec.name, field_buffers[i], global_site_count, local_site_count, GetFieldLength(field_spec.src), comms, mpi_comm);
    }

    H5_CHECK(H5Fclose(file_id));
    
    WriteXDMFFile(current_h5_filename, timestepNumber);
  }
}
#endif


