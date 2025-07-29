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

  void LocalPropertyHdf5Output::WriteXDMFFile()
  {
    if (comms.Rank() != 0) return;
    
    std::string h5_filename = outputSpec.filename.string();
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
    xmf_file << "    <Grid Name=\"HemeLB Collection\" GridType=\"Collection\" CollectionType=\"Temporal\">\n";

    for (unsigned long ts : written_timesteps) {
        xmf_file << "      <Grid Name=\"step_" << ts << "\" GridType=\"Uniform\">\n";
        xmf_file << "        <Time Value=\"" << ts << "\"/>\n";
        xmf_file << "        <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << global_site_count << "\"/>\n";
        xmf_file << "        <Geometry GeometryType=\"XYZ\">\n";
        xmf_file << "          <DataItem Format=\"HDF\" Dimensions=\"" << global_site_count << " 3\" DataType=\"Float\" Precision=\"8\">" << h5_basename << ":/step_" << ts << "/geometry</DataItem>\n";
        xmf_file << "        </Geometry>\n";

        for (const auto& field : outputSpec.fields) {
            unsigned field_len = GetFieldLength(field.src);
            std::string type = (field_len > 1) ? "Vector" : "Scalar";
            
            std::string precision_str = "4";
            std::string datatype_str = std::visit([](auto t) -> std::string {
                using T = decltype(t);
                if constexpr (std::is_same_v<T, int> || std::is_same_v<T, unsigned int> || std::is_same_v<T, long int> || std::is_same_v<T, long unsigned int>) return "Int";
                if constexpr (std::is_same_v<T, float>) return "Float";
                if constexpr (std::is_same_v<T, double>) return "Float";
                return "Unknown";
            }, field.typecode);

            if (datatype_str == "Float") {
                 precision_str = std::visit([](auto t) { using T = decltype(t); return std::to_string(sizeof(T)); }, field.typecode);
            }

            xmf_file << "        <Attribute Name=\"" << field.name << "\" AttributeType=\"" << type << "\" Center=\"Node\">\n";
            xmf_file << "          <DataItem Format=\"HDF\" Dimensions=\"" << global_site_count << " " << field_len << "\" DataType=\"" << datatype_str << "\" Precision=\"" << precision_str << "\">" << h5_basename << ":/step_" << ts << "/" << field.name << "</DataItem>\n";
            xmf_file << "        </Attribute>\n";
        }
        xmf_file << "      </Grid>\n";
    }

    xmf_file << "    </Grid>\n";
    xmf_file << "  </Domain>\n";
    xmf_file << "</Xdmf>\n";
    xmf_file.close();
  }

  void LocalPropertyHdf5Output::Write(unsigned long timestepNumber, unsigned long totalSteps) {
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
    
    // get buffer sizes
    std::vector<uint32_t> coords_buffer;
    std::vector<float> pressure_buffer;
    std::vector<double> velocity_buffer;
    std::vector<float> shearstress_buffer;
    
    // traverse once and collect data
    dataSource.Reset();
    while (dataSource.ReadNext()) {
      if (outputSpec.geometry->Include(dataSource, dataSource.GetPosition())) {
        const auto& pos = dataSource.GetPosition();
        coords_buffer.push_back(pos.x());
        coords_buffer.push_back(pos.y());
        coords_buffer.push_back(pos.z());

        for (const auto& field_spec : outputSpec.fields) {
            if (std::holds_alternative<source::Pressure>(field_spec.src)) {
                double val = dataSource.GetPressure();
                pressure_buffer.push_back(static_cast<float>(val));
                // pressure_buffer.push_back(static_cast<float>((std::isinf(val) || std::isnan(val)) ? -1.0 : val));
            } else if (std::holds_alternative<source::Velocity>(field_spec.src)) {
                const auto& vel = dataSource.GetVelocity();
                velocity_buffer.push_back(static_cast<double>(vel.x()));
                velocity_buffer.push_back(static_cast<double>(vel.y()));
                velocity_buffer.push_back(static_cast<double>(vel.z()));
                // velocity_buffer.push_back(static_cast<double>((std::isinf(vel.x()) || std::isnan(vel.x())) ? 0.0 : vel.x()));
                // velocity_buffer.push_back(static_cast<double>((std::isinf(vel.y()) || std::isnan(vel.y())) ? 0.0 : vel.y()));
                // velocity_buffer.push_back(static_cast<double>((std::isinf(vel.z()) || std::isnan(vel.z())) ? 0.0 : vel.z()));
            } else if (std::holds_alternative<source::ShearStress>(field_spec.src)) {
                double val = dataSource.GetShearStress();
                shearstress_buffer.push_back(static_cast<float>(val));
                // shearstress_buffer.push_back(static_cast<float>((std::isinf(val) || std::isnan(val)) ? -1.0 : val));
            }
        }
      }
    }
  
    // 3. 并行写入所有数据集
    WriteDataset(group_id, "geometry", coords_buffer, H5T_NATIVE_UINT32, global_site_count, local_site_count, 3, comms, mpi_comm);
    
    for (const auto& field_spec : outputSpec.fields) {
        if (std::holds_alternative<source::Pressure>(field_spec.src)) {
            WriteDataset(group_id, field_spec.name, pressure_buffer, H5T_NATIVE_FLOAT, global_site_count, local_site_count, 1, comms, mpi_comm);
        } else if (std::holds_alternative<source::Velocity>(field_spec.src)) {
            WriteDataset(group_id, field_spec.name, velocity_buffer, H5T_NATIVE_DOUBLE, global_site_count, local_site_count, 3, comms, mpi_comm);
        } else if (std::holds_alternative<source::ShearStress>(field_spec.src)) {
            WriteDataset(group_id, field_spec.name, shearstress_buffer, H5T_NATIVE_FLOAT, global_site_count, local_site_count, 1, comms, mpi_comm);
        }
    }

    H5_CHECK(H5Gclose(group_id));
    // H5_CHECK(H5Fclose(file_id));
    // file_id = -1; // Reset file_id to indicate that the file is closed
  }
}
#endif